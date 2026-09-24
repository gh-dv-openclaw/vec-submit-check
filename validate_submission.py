#!/usr/bin/env python3
"""提交前的本地格式校验：把一个 .h5ad 对着某块板的契约逐条检查。

赛方的 format check 不限次也不消耗配额，但每次都要上传几百 MB 再等回包。
本地先过一遍，把「基因顺序错了」这类问题在上传前挡掉。

    python validate_submission.py --board T1:val --input pred.h5ad

契约来源 `panels/`（`bash fetch_panels.sh` 从官网拉取，无需登录）。
检查项对应官网 /challenge/data 的 "Gene panel, per board" 与各任务的 File format 表。

注意 panel 不等于「目标文件的基因集」，而是「该板加载的所有 stage 的有序交集」——
embryo 板因此只有 498 个基因（E6.75/E7.25 没测 Casp4 和 Pnliprp1），
所以基因表必须从 panels 文件取，不能从任何一个数据文件取。

另外查 .X 的尺度：原始计数、没取 log、归一化到别的总量，这些文件都能过格式检查、照常出分，
但打分器按 log1p(CP10k) 去比，分数没有意义，而且不会有任何提示。

退出码：0 全过；1 有 FAIL，上传会被拒；3 收得下，但 .X 不像 log1p(CP10k)。
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

PANELS = Path(__file__).parent / "panels"

# 赛方发布的每个文件，.X 都是 log1p(CP10k)：每个细胞的计数先缩放到总和 10000，再取 log1p。
# 单个值因此不会超过 log1p(10000)——那是一个细胞的全部计数都落在同一个基因上的情形。
TARGET_SUM = 1e4
LOG1P_CP10K_MAX = float(np.log1p(TARGET_SUM))   # 9.2103
# 「各细胞行和相等」的容差。真正「归一化 → 变换」出来的文件，行和只差浮点误差（发布文件 ≤ 1e-6）；
# 模型在 log 空间乘一个常数也会让某个逆变换下的行和「差不多」相等 —— ×1.44 ≈ 1/ln2 时只差 0.2–0.8%，
# 看起来就像 log2。所以卡在 1e-4，两边各留二十倍以上。
CONST_TOL = 1e-4

# 把 .X 还原成「每个细胞归一化后的计数」的几种逆变换。发布的文件在 expm1 下每个细胞加起来恒为
# 10000，其余几种对应常见的错法：哪一种让各细胞的行和几乎相等，就说明 .X 是怎么来的。
# 每个函数都满足 f(0) = 0，所以稀疏矩阵只需对非零元做。
INVERSES = (
    ("log1p", np.expm1),
    ("没取 log", lambda v: v),
    ("log2(1+x)", lambda v: np.exp2(v) - 1),
    ("log10(1+x)", lambda v: np.power(10.0, v) - 1),
    ("log1p 做了两次", lambda v: np.expm1(np.expm1(v))),
)


class Report:
    """收集检查结果。FAIL 让退出码为 1；.X 尺度可疑让退出码为 3；WARN 只提示。"""

    def __init__(self) -> None:
        self.failed = 0
        self.suspect = 0

    def _emit(self, tag: str, name: str, detail: str) -> None:
        print(f"  [{tag:^4}] {name:<22} {detail}")

    def ok(self, name: str, detail: str = "") -> None:
        self._emit("PASS", name, detail)

    def warn(self, name: str, detail: str) -> None:
        self._emit("WARN", name, detail)

    def suspicious(self, name: str, detail: str) -> None:
        """收得下、会被打分，但分数没有意义。比 FAIL 更贵：它会用掉一次评分。"""
        self.suspect += 1
        self._emit("WARN", name, detail)

    def skip(self, name: str, detail: str) -> None:
        self._emit("SKIP", name, detail)

    def fail(self, name: str, detail: str) -> None:
        self.failed += 1
        self._emit("FAIL", name, detail)


def check_genes(r: Report, adata: ad.AnnData, panel: list[str]) -> None:
    """var_names 必须与 panel 逐元素同名同序。顺序错也拒。"""
    got = list(adata.var_names)
    if got == panel:
        r.ok("gene panel", f"{len(panel)} 个基因，顺序一致")
        return

    if len(got) != len(panel):
        r.fail("gene panel", f"基因数 {len(got)}，应为 {len(panel)}")
    if set(got) == set(panel):
        first = next(i for i, (a, b) in enumerate(zip(got, panel)) if a != b)
        r.fail("gene panel", f"基因集合正确但顺序不同，首个差异在第 {first} 位："
                             f"{got[first]!r} ≠ {panel[first]!r} —— 按 panel 重排即可")
        return

    missing = [g for g in panel if g not in set(got)]
    extra = [g for g in got if g not in set(panel)]
    if missing:
        r.fail("gene panel", f"缺 {len(missing)} 个：{missing[:5]}{' …' if len(missing) > 5 else ''}")
    if extra:
        r.fail("gene panel", f"多 {len(extra)} 个：{extra[:5]}{' …' if len(extra) > 5 else ''}")


def check_cells(r: Report, adata: ad.AnnData, board: dict) -> None:
    """细胞数要落在区间内。上界是评分主机的内存约束，下界是统计有效性。"""
    n, lo, hi = adata.n_obs, board["min_cells"], board["max_cells"]
    if lo <= n <= hi:
        note = ""
        if n > 2500:
            note = "（提示：MMD 只取 2000、energy/variogram 各取 1500，多出的细胞会被丢弃）"
        r.ok("cell count", f"{n} 在 [{lo}, {hi}] 内{note}")
    else:
        r.fail("cell count", f"{n} 不在 [{lo}, {hi}] 内")


def check_X(r: Report, adata: ad.AnnData) -> bool:
    """.X 必须有限且非负。出现负值通常意味着上游做过中心化。返回数值本身是否可用。"""
    X = adata.X
    dense = X.data if sparse.issparse(X) else np.asarray(X).ravel()

    n_bad = int((~np.isfinite(dense)).sum())
    if n_bad:
        r.fail(".X finite", f"{n_bad} 个 NaN/Inf —— 缺失指标按 skill=0 计，不豁免")
    else:
        r.ok(".X finite", "无 NaN/Inf")

    n_neg = int((dense < 0).sum())
    if n_neg:
        r.fail(".X non-negative", f"{n_neg} 个负值（最小 {dense.min():.4g}）—— 通常是被中心化过")
    else:
        # 稀疏矩阵只存非零元，dense 拿到的是 X.data；有隐式 0 时真正的最小值就是 0
        has_implicit_zeros = sparse.issparse(X) and X.nnz < X.shape[0] * X.shape[1]
        lo = 0.0 if has_implicit_zeros else (float(dense.min()) if dense.size else 0.0)
        r.ok(".X non-negative", f"最小 {lo:.4g}")

    dt = X.dtype
    if dt != np.float32:
        r.warn(".X dtype", f"{dt}，官网标的是 float32（不致拒收，但会放大文件）")
    else:
        r.ok(".X dtype", "float32")

    kind = "sparse" if sparse.issparse(X) else "dense"
    r.ok(".X layout", f"{kind}（稀疏/稠密都收，评分前会转稠密）")
    return n_bad == 0 and n_neg == 0


def scale_stats(X, n_rows: int = 2000, seed: int = 0) -> dict:
    """抽最多 n_rows 个细胞，量出能把 log1p(CP10k) 和常见错法分开的几个数。

    vmax           全矩阵最大值
    nnz            抽到的细胞里非零值的个数
    frac_int       非零值里恰为整数的比例              原始计数、二值化 → 1
    frac_logint    非零值里 expm1 之后恰为整数的比例    log1p(原始计数) → 1
    rowsum_*       每个细胞 Σ expm1(x) 的 10/50/90 分位   log1p(CP10k) → 10000
    const_inverse  (名字, 行和)：INVERSES 里第一个让各细胞行和相差不到 CONST_TOL 的逆变换

    「恰为整数」按存储精度判：回写后与原值相差不超过 2 个 ulp。
    """
    n = X.shape[0]
    if n == 0:
        return {"vmax": 0.0, "nnz": 0}
    rows = np.sort(np.random.default_rng(seed).choice(n, size=min(n, n_rows), replace=False))
    if sparse.issparse(X):
        sub = sparse.csr_matrix(X[rows], copy=True)
    else:                                       # 稠密的 T1 文件整块转 float64 很大，分块转成稀疏
        sub = sparse.vstack([sparse.csr_matrix(np.asarray(X[rows[i:i + 256]]))
                             for i in range(0, len(rows), 256)], format="csr")
    sub.eliminate_zeros()                       # 就地操作；sub 是副本，碰不到调用方的矩阵
    vals = sub.data

    out = {"vmax": float(X.max()), "nnz": int(vals.size)}
    if not vals.size:
        return out
    prec = vals.dtype if vals.dtype.kind == "f" else np.dtype(np.float64)
    v = vals.astype(np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        k = np.round(np.expm1(v))
        back = np.log1p(k).astype(prec).astype(np.float64)
    ulp = np.spacing(np.abs(vals).astype(prec)).astype(np.float64)
    out["frac_int"] = float(np.mean(v == np.round(v)))
    out["frac_logint"] = float(np.mean((k >= 1) & (np.abs(back - v) <= 2 * ulp)))

    nonempty = np.diff(sub.indptr) > 0
    for name, inverse in INVERSES:
        e = sub.astype(np.float64)
        with np.errstate(over="ignore", invalid="ignore"):
            e.data = inverse(e.data)
        sums = np.asarray(e.sum(axis=1)).ravel()[nonempty]
        if name == "log1p":
            p10, p50, p90 = np.quantile(sums, [0.1, 0.5, 0.9], method="lower")   # 不插值，inf 也能排
            out.update(rowsum_p10=float(p10), rowsum_med=float(p50), rowsum_p90=float(p90))
        lo, hi = sums.min(), sums.max()
        if "const_inverse" not in out and lo > 0 and np.isfinite(hi) and hi <= lo * (1 + CONST_TOL):
            out["const_inverse"] = (name, float(np.median(sums)))
    return out


def _num(x: float) -> str:
    """一百万以内写全（CPM 的 999,999.99 与 1,000,000.01 显示一样），再大用科学计数。"""
    return f"{x:,.0f}" if abs(x) < 1e7 else f"{x:.3g}"


def scale_verdict(s: dict) -> tuple[bool, str]:
    """(是否可疑, 说明)。可疑 = 服务器照收、照打分，但分数没有意义。

    顺序有讲究：没取 log 的 MERFISH 文件 expm1 之后也「恰为整数」（数值大到 float64 只剩整数），
    所以「行和恒定」与「数值大到不可能」要排在 log1p(原始计数) 那条前面。
    """
    if s["nnz"] == 0:
        return True, "抽到的细胞全是 0 —— 没有可比的表达"
    if s["frac_int"] >= 0.9:
        if s["vmax"] <= 1:
            return True, "非零值全是 1 —— 像二值化的「有/无」，不是表达量"
        return True, (f"非零值 {s['frac_int']:.0%} 是整数、最大 {s['vmax']:.0f} —— 像原始计数。"
                      f"先把每个细胞缩放到总和 10000 再 log1p")
    name, total = s.get("const_inverse") or (None, None)
    if name == "log1p":
        if abs(total / TARGET_SUM - 1) <= 0.02:
            return False, f"log1p(CP10k)：每个细胞 expm1 行和 = {_num(total)}，与发布文件同一尺度"
        return True, (f"每个细胞 expm1 行和都是 {_num(total)}，不是 10,000 —— log1p 之前归一化到了"
                      f"别的总量（scanpy 的 normalize_total 不传 target_sum 时取中位文库大小）")
    if name == "没取 log":
        return True, f"每个细胞的值加起来都是 {_num(total)} —— 做了归一化，但没取 log1p"
    if name in ("log2(1+x)", "log10(1+x)"):
        return True, (f"按 {name} 还原后每个细胞加起来都是 {_num(total)} —— "
                      f"发布文件用的是自然对数 log1p")
    if name == "log1p 做了两次":
        return True, f"连做两次 expm1 后每个细胞加起来都是 {_num(total)} —— log1p 做了两次"
    med = s.get("rowsum_med", 0.0)
    if not np.isfinite(med) or med > 1e20:
        return True, f"最大 {s['vmax']:.4g}，expm1 行和中位 {med:.3g} —— 看起来没取 log"
    if s["frac_logint"] >= 0.9:
        return True, (f"expm1 之后 {s['frac_logint']:.0%} 的值是整数 —— 像 log1p(原始计数)，"
                      f"log1p 之前没按文库大小归一化")
    detail = (f"没有计数 / 未取 log / 未归一化的迹象；每个细胞 expm1 行和不恒定"
              f"（p10–p90 = {_num(s['rowsum_p10'])}–{_num(s['rowsum_p90'])}，"
              f"预测改过数值或截掉过基因时都会这样）")
    if s["vmax"] > LOG1P_CP10K_MAX + 1e-3:
        detail += (f"；最大 {s['vmax']:.4g} 超过 log1p(10000) = {LOG1P_CP10K_MAX:.2f}，"
                   f"模型在 log 空间放大过数值时正常")
    return False, detail


def check_scale(r: Report, adata: ad.AnnData) -> None:
    """原始计数、没取 log、归一化到别的总量 —— 都收得下、照常出分，而分数没有意义。"""
    suspicious, detail = scale_verdict(scale_stats(adata.X))
    (r.suspicious if suspicious else r.ok)(".X scale", detail)


def check_obsm(r: Report, adata: ad.AnnData, board: dict) -> None:
    """T2/T3 必须带 3D 坐标，缺了在评分前就被拒。"""
    for key in board["obsm_required"]:
        if key not in adata.obsm:
            r.fail(f"obsm[{key!r}]", "缺失 —— 有表达无坐标的文件在评分前直接拒")
            continue
        arr = np.asarray(adata.obsm[key])
        if arr.ndim != 2 or arr.shape[1] < 3:
            r.fail(f"obsm[{key!r}]", f"形状 {arr.shape}，需要 [n_obs, ≥3]")
        elif not np.isfinite(arr[:, :3]).all():
            r.fail(f"obsm[{key!r}]", "前三列含 NaN/Inf")
        else:
            extra = f"，第 3 列之后被忽略" if arr.shape[1] > 3 else ""
            span = arr[:, :3].max(0) - arr[:, :3].min(0)
            r.ok(f"obsm[{key!r}]", f"{arr.shape}{extra}，包围盒 {np.round(span, 2).tolist()}")

    if not board["needs_coords"] and adata.obsm:
        r.warn("obsm", f"该板不读坐标，{list(adata.obsm)} 是无用负重")


def check_obs(r: Report, adata: ad.AnnData) -> None:
    """细胞类型标签从不提交，赛方用冻结分类器自己标。"""
    if "celltype" in adata.obs:
        r.warn("obs['celltype']", "提交不带标签，赛方用自己的冻结分类器标 —— 删掉可减小文件")
    else:
        r.ok("obs", "未携带 celltype（正确）")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--board", required=True, help="板名，如 T1:val / T2:heart:val_extrap / T3:gata4")
    p.add_argument("--input", required=True, help="待校验的 .h5ad")
    p.add_argument("--panels", default=str(PANELS / "panels_index.json"))
    a = p.parse_args(argv)

    panels = json.loads(Path(a.panels).read_text())
    if a.board not in panels:
        p.error(f"未知的板 {a.board!r}；可选：{', '.join(panels)}")
    board = panels[a.board]

    panel_path = Path(a.panels).parent / board["genes_file"]
    if not panel_path.exists():
        p.error(f"缺基因表 {panel_path}，先跑 bash fetch_panels.sh")
    panel = panel_path.read_text().split()

    path = Path(a.input)
    size_mb = path.stat().st_size / 2**20
    adata = ad.read_h5ad(path)

    print(f"\n{a.board}  ({board['label']})")
    print(f"  文件 {path.name}  {size_mb:.1f} MB  →  {adata.n_obs} cells × {adata.n_vars} genes\n")

    r = Report()
    check_genes(r, adata, panel)
    check_cells(r, adata, board)
    if check_X(r, adata):
        check_scale(r, adata)
    else:
        r.skip(".X scale", "先修上面的 NaN/Inf 或负值")
    check_obsm(r, adata, board)
    check_obs(r, adata)

    if size_mb > 1200:
        r.fail("file size", f"{size_mb:.0f} MB 超过单文件 1200 MB 上限")

    print()
    if r.failed:
        print(f"  {r.failed} 项不合格 —— 上传会被拒\n")
        return 1
    if r.suspect:
        print("  格式检查全过，但 .X 不像 log1p(CP10k)：这份文件会被接收、会用掉一次评分，")
        print("  打分器却按 log1p(CP10k) 去比，分数没有意义。确认是有意为之再上传（退出码 3）\n")
        return 3
    print("  格式检查全过（这不预示分数，只说明文件收得下）\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
