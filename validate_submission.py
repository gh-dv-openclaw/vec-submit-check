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


class Report:
    """收集检查结果。FAIL 会让退出码非零；WARN 只提示。"""

    def __init__(self) -> None:
        self.failed = 0

    def _emit(self, tag: str, name: str, detail: str) -> None:
        print(f"  [{tag:^4}] {name:<22} {detail}")

    def ok(self, name: str, detail: str = "") -> None:
        self._emit("PASS", name, detail)

    def warn(self, name: str, detail: str) -> None:
        self._emit("WARN", name, detail)

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


def check_X(r: Report, adata: ad.AnnData) -> None:
    """.X 必须有限且非负。出现负值通常意味着上游做过中心化。"""
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


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--board", required=True, help="板名，如 T1:val / T2:heart:val_extrap / T3:gata4")
    p.add_argument("--input", required=True, help="待校验的 .h5ad")
    p.add_argument("--panels", default=str(PANELS / "panels_index.json"))
    a = p.parse_args()

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
    check_X(r, adata)
    check_obsm(r, adata, board)
    check_obs(r, adata)

    if size_mb > 1200:
        r.fail("file size", f"{size_mb:.0f} MB 超过单文件 1200 MB 上限")

    print()
    if r.failed:
        print(f"  {r.failed} 项不合格 —— 上传会被拒\n")
        return 1
    print("  格式检查全过（这不预示分数，只说明文件收得下）\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
