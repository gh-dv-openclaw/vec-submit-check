#!/usr/bin/env python3
"""validate_submission.py 的测试。全部用合成数据，不含任何比赛数据。

    python -m pytest -q
    python test_validate_submission.py        # 不装 pytest 也能跑

合成数据按发布文件的做法构造：整数计数 → 每个细胞缩放到总和 10000 → log1p。
每种错法各造一份，检查退出码和提示里点名的原因。尺度相关的检查在两种文库大小下各测一遍：
scRNA 那样的几千个计数，和 MERFISH 那样的一两百个（小文库下 expm1(.X) 恰为整数的比例更高，
正确的文件也会有一部分，阈值要放得过）。
"""
from __future__ import annotations

import contextlib
import io
import json
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate_submission as vs  # noqa: E402

N_GENES = 120
GENES = [f"G{i:04d}" for i in range(N_GENES)]
BOARDS = {
    "TEST:scrna": {"label": "synthetic, no coordinates", "genes_file": "test.genes.txt",
                   "min_cells": 100, "max_cells": 5000, "obsm_required": [], "needs_coords": False},
    "TEST:merfish": {"label": "synthetic, needs coordinates", "genes_file": "test.genes.txt",
                     "min_cells": 100, "max_cells": 5000, "obsm_required": ["spatial_3D"],
                     "needs_coords": True},
}
REGIMES = {"scrna": 5000, "merfish": 150}          # 每个细胞的平均计数


# ---------- 造数据 ----------

def counts(n_cells: int = 400, mean_lib: int = 5000, seed: int = 0) -> sparse.csr_matrix:
    """整数计数：四种「细胞类型」各有自己的表达谱（少数基因高、多数稀疏），文库大小逐细胞不同。

    类型要有几种：所有细胞共用一个表达谱时，各细胞的行和天然更接近，测不出对编辑的误报。
    """
    rng = np.random.default_rng(seed)
    profiles = rng.dirichlet(np.full(N_GENES, 0.3), size=4)
    kind = rng.integers(0, 4, n_cells)
    lib = rng.lognormal(np.log(mean_lib), 0.5, n_cells).astype(int) + 20
    return sparse.csr_matrix(np.stack([rng.multinomial(n, profiles[k])
                                       for n, k in zip(lib, kind)]).astype(np.float64))


def normalise(C: sparse.csr_matrix, target) -> sparse.csr_matrix:
    lib = np.asarray(C.sum(axis=1)).ravel()
    return sparse.csr_matrix(C.multiply(np.asarray(target, dtype=float) / lib[:, None]))


def log1p_cp10k(C: sparse.csr_matrix) -> sparse.csr_matrix:
    return normalise(C, 1e4).log1p()


def f32(M) -> sparse.csr_matrix:
    return sparse.csr_matrix(M, dtype=np.float32)


def adata(X, coords: bool = False, genes=GENES) -> ad.AnnData:
    a = ad.AnnData(X=X)
    a.var_names = list(genes)
    a.obs_names = [f"c{i}" for i in range(X.shape[0])]
    if coords:
        a.obsm["spatial_3D"] = np.random.default_rng(1).normal(size=(X.shape[0], 3)).astype(np.float32)
    return a


def run(a: ad.AnnData, board: str = "TEST:scrna") -> tuple[int, str]:
    """写盘、跑一遍 CLI 入口，返回 (退出码, 输出)。"""
    with tempfile.TemporaryDirectory(prefix="vsc-test-") as d:
        d = Path(d)
        (d / "test.genes.txt").write_text("\n".join(GENES) + "\n")
        (d / "index.json").write_text(json.dumps(BOARDS))
        a.write_h5ad(d / "pred.h5ad")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            code = vs.main(["--board", board, "--input", str(d / "pred.h5ad"),
                            "--panels", str(d / "index.json")])
    return code, buf.getvalue()


def expect(a: ad.AnnData, code: int, *phrases: str, board: str = "TEST:scrna") -> str:
    got, out = run(a, board)
    assert got == code, f"退出码 {got}，应为 {code}\n{out}"
    for ph in phrases:
        assert ph in out, f"输出里没有 {ph!r}\n{out}"
    return out


# ---------- 正确的文件必须放行 ----------

def test_released_format_passes():
    for regime, lib in REGIMES.items():
        out = expect(adata(f32(log1p_cp10k(counts(mean_lib=lib)))), 0, "[PASS] .X scale", "log1p(CP10k)")
        assert "[WARN]" not in out, f"{regime}: 不该有任何警告\n{out}"


def test_dense_passes():
    expect(adata(f32(log1p_cp10k(counts())).toarray()), 0, "[PASS] .X scale")


def test_resampled_real_cells_pass():
    """挑真细胞、有放回重采样：每行仍是一个真实细胞，行和仍是 10000。"""
    X = f32(log1p_cp10k(counts()))
    idx = np.random.default_rng(2).choice(X.shape[0], size=600, replace=True)
    expect(adata(X[idx]), 0, "[PASS] .X scale", "log1p(CP10k)")


def test_model_edits_in_log_space_pass():
    """模型在 log 空间改数值是正常的预测，不该被当成错法——哪怕最大值超过 9.21。"""
    rng = np.random.default_rng(3)
    for regime, lib in REGIMES.items():
        X = f32(log1p_cp10k(counts(mean_lib=lib)))
        scaled = f32(X * 1.5)
        out = expect(adata(scaled), 0, "[PASS] .X scale", "模型在 log 空间放大过数值时正常")
        assert "[WARN] .X scale" not in out, regime
        jitter = X.copy()
        jitter.data *= rng.uniform(0.8, 1.25, size=jitter.data.size).astype(np.float32)
        expect(adata(jitter), 0, "[PASS] .X scale")


def test_float64_only_warns_dtype():
    out = expect(adata(sparse.csr_matrix(log1p_cp10k(counts()), dtype=np.float64)), 0, "[WARN] .X dtype")
    assert "[PASS] .X scale" in out


# ---------- 收得下、但分数没有意义的：退出码 3，并点名原因 ----------

def test_raw_counts_flagged():
    for lib in REGIMES.values():
        expect(adata(f32(counts(mean_lib=lib))), 3, "[WARN] .X scale", "原始计数")


def test_binarised_flagged():
    expect(adata(f32(counts() > 0)), 3, "二值化")


def test_normalised_but_not_logged_flagged():
    for lib in REGIMES.values():
        expect(adata(f32(normalise(counts(mean_lib=lib), 1e4))), 3, "没取 log1p")


def test_log1p_of_raw_counts_flagged():
    for lib in REGIMES.values():
        expect(adata(f32(counts(mean_lib=lib).log1p())), 3, "log1p(原始计数)")


def test_cpm_flagged():
    expect(adata(f32(normalise(counts(), 1e6).log1p())), 3, "1,000,000", "不是 10,000")


def test_scanpy_default_median_target_flagged():
    """normalize_total 不传 target_sum 时用中位文库大小。MERFISH 上那是 1e4 的百分之一。"""
    for lib in REGIMES.values():
        C = counts(mean_lib=lib)
        med = np.median(np.asarray(C.sum(axis=1)).ravel())
        expect(adata(f32(normalise(C, med).log1p())), 3, "normalize_total")


def test_wrong_log_base_flagged():
    base2 = f32(log1p_cp10k(counts()) / np.log(2))
    expect(adata(base2), 3, "log2(1+x)", "自然对数")
    base10 = f32(log1p_cp10k(counts()) / np.log(10))
    expect(adata(base10), 3, "log10(1+x)")


def test_double_log_flagged():
    expect(adata(f32(log1p_cp10k(counts()).log1p())), 3, "log1p 做了两次")


def test_all_zero_flagged():
    expect(adata(sparse.csr_matrix((400, N_GENES), dtype=np.float32)), 3, "全是 0")


# ---------- 原有的格式检查：退出码 1 ----------

def test_gene_order_names_first_difference():
    swapped = GENES.copy()
    swapped[7], swapped[8] = swapped[8], swapped[7]
    expect(adata(f32(log1p_cp10k(counts())), genes=swapped), 1, "首个差异在第 7 位", "'G0008'")


def test_cell_count_out_of_range():
    expect(adata(f32(log1p_cp10k(counts(n_cells=50)))), 1, "[FAIL] cell count")


def test_nan_fails_and_scale_is_skipped():
    X = f32(log1p_cp10k(counts()))
    X.data[:5] = np.nan
    expect(adata(X), 1, "[FAIL] .X finite", "[SKIP] .X scale")


def test_centred_matrix_fails():
    X = f32(log1p_cp10k(counts())).toarray()
    expect(adata(X - X.mean(axis=0)), 1, "[FAIL] .X non-negative")


def test_missing_coordinates_fail_on_spatial_board():
    X = f32(log1p_cp10k(counts()))
    expect(adata(X), 1, "[FAIL] obsm['spatial_3D']", board="TEST:merfish")
    expect(adata(X, coords=True), 0, "[PASS] obsm['spatial_3D']", board="TEST:merfish")


# ---------- 实现细节 ----------

def test_scale_stats_does_not_touch_the_input():
    """csr_matrix(M) 与 M 共用索引数组，就地 eliminate_zeros 会把调用方的矩阵改错位。"""
    X = f32(log1p_cp10k(counts()))
    X.data[::7] = 0                                     # 显式存储的零
    before = (X.data.copy(), X.indices.copy(), X.indptr.copy())
    vs.scale_stats(X)
    assert all(np.array_equal(a, b) for a, b in zip(before, (X.data, X.indices, X.indptr)))


def test_scale_stats_on_released_format_is_exact():
    s = vs.scale_stats(f32(log1p_cp10k(counts())))
    name, total = s["const_inverse"]
    assert name == "log1p" and abs(total - 1e4) < 1, s
    assert s["frac_int"] == 0.0 and s["vmax"] <= vs.LOG1P_CP10K_MAX + 1e-3, s


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  ok    {name}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {name}\n{e}\n")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
