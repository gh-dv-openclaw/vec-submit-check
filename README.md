# vec-submit-check · Virtual Embryo Challenge 提交自查工具与中文入门

给 [Virtual Embryo Challenge](https://virtualembryo.ai/challenge)（NeurIPS 2026）参赛者的
提交格式校验器，外加一份中文入门说明。

赛方提供**不限次、不消耗配额**的在线格式检查，但每次都要上传几百 MB 再等回包。
这个工具把同样的检查放在本地做，几秒钟给出结果，还会指出具体是哪个基因、哪一位错了。

```bash
git clone https://github.com/<你的用户名>/vec-submit-check && cd vec-submit-check
bash fetch_panels.sh                 # 拉基因 panel，公开的，不需要登录
pip install anndata scipy numpy      # 或 uv pip install
python validate_submission.py --board T1:val --input pred.h5ad
```

输出长这样：

```
T1:val  (T1 · validation (E10.5))
  文件 pred.h5ad  84.2 MB  →  3000 cells × 32285 genes

  [PASS] gene panel             32285 个基因，顺序一致
  [PASS] cell count             3000 在 [1000, 5118] 内
  [PASS] .X finite              无 NaN/Inf
  [PASS] .X non-negative        最小 0
  [PASS] .X dtype               float32
  [PASS] .X layout              sparse（稀疏/稠密都收，评分前会转稠密）
  [PASS] obs                    未携带 celltype（正确）

  格式检查全过（这不预示分数，只说明文件收得下）
```

板名有五个：`T1:val`、`T2:heart:val_interp`、`T2:heart:val_extrap`、
`T2:embryo:val_interp`、`T3:gata4`。

---

## 提交契约里真正会绊人的几处

下面每一条都能在官网 [/challenge/data](https://virtualembryo.ai/challenge/data) 找到出处，
但分散在很长的页面里，这里按「容易踩」的顺序排。

### 1. 基因 panel 不等于「目标文件里的基因」

这是最容易错的一条。panel 是**该块板加载的所有 stage 的有序交集**，不是任何单个文件的基因集。

embryo 那块板因此只有 **498** 个基因而不是 500：`Casp4` 和 `Pnliprp1` 在 E6.75 和 E7.25
没有测，所以它们从每一块 embryo 板里消失了——**包括目标 stage E7.5 明明测了它们的那块**。
理由是「拿一个参照 stage 从未观测过的基因，没法给预测打分」。

所以基因表要从 `panels/` 里取，**不要从任何一个数据文件里取**。

### 2. 顺序错也会被拒

`var_names` 是逐元素比对的：同名、同序。集合对了但顺序不同一样拒收。
本工具会直接告诉你第一个错位在第几位：

```
[FAIL] gene panel   基因集合正确但顺序不同，首个差异在第 1247 位：'Gata4' ≠ 'Gata6' —— 按 panel 重排即可
```

修法就是 `adata = adata[:, panel]`。

### 3. 细胞数：有区间，但多了没用

每块板有 `min_cells` / `max_cells`，在 `panels_index.json` 里。

- **上限是评分主机的内存约束**，不是生物学。T1 的答案会在评分机上展开成
  `n_cells × 32285` 的 float32。
- **下限是统计有效性**：每个主指标都在估一个分布，太少了估计量被自身采样噪声主导。
- **往高了报没有任何好处**：三个最重的指标各自再做一次下采样——unbiased MMD 取 2000，
  energy distance 和 variogram 各取 1500。**超过约 2000 个细胞会在被测量之前就被丢弃。**
- 细胞数本身不计分，也不与任何东西比较。它不是「你觉得胚胎有多大」的预测。

几千个细胞在每块板上都是好的提交。

### 4. T2 / T3 必须带 3D 坐标

`obsm["spatial_3D"]`，形状 `[n_obs, 3]`，第三列之后被忽略。
**有表达没坐标的文件在评分开始之前就被拒**，不会得到任何分数。

坐标系随你定：所有空间指标对平移和旋转不变。
**已知盲点是手性**——镜像的胚胎和正确的得分相同，官网明说了这一点。

### 5. 不要提交细胞类型标签

赛方对每一份提交都用同一个冻结分类器自己标类型，所以 `obs["celltype"]` 既不需要也不会被读。
带着它只是让文件更大。

### 6. 缺失的指标按「零技能」计，不豁免

一个你的输出产生不了的指标、或者算出来是 NaN 的指标，在它所属的问题里**记 0 分**，
不会被悄悄从平均里剔除。否则等于让一个模型在更容易的子集上被评判。

所以 `.X` 里出现 NaN/Inf 不是小事。负值通常意味着矩阵在上游某处被中心化过了——
提交要求是 log 归一化后的非负值。

### 7. 单文件不超过 1200 MB

---

## 数据怎么下

注册（GitHub 或 Google 登录）之后，训练 stage 就可以下了。
浏览器里点着下也行；要在服务器上下，接口是这两个：

```bash
# 登录态在 localStorage 的 ve-challenge-session 里，取 .token 字段
curl -H "Authorization: Bearer $TOKEN" \
  https://kg.virtualembryo.ai/challenge/data/manifest

# 逐个文件换预签名链接（TTL 900 秒，所以要下的时候再换，别批量换完慢慢下）
curl -H "Authorization: Bearer $TOKEN" \
  "https://kg.virtualembryo.ai/challenge/data/link?key=T2/E6.75.h5ad"
```

manifest 会给出每个文件的 `bytes`，下完对一下尺寸。训练数据总共约 2.26 GiB。

> 把 token 写进文件（`chmod 600`）再用 `curl -K` 读，别内联到命令行里——
> 内联会让它出现在 `ps` 的输出中。

**验证集和测试集的真值不会发放**，提交回来的是分数不是答案。
到最终阶段（2026-10-20）验证集答案会全部释放，测试集答案永不释放。

---

## 本地跑分

赛方开源了打分器 [`aristoteleo/veckit`](https://github.com/aristoteleo/veckit)，
可以离线对着**你自己提供的**参照文件打出指标面板：

```bash
pip install veckit
veckit --task T1 --input pred.h5ad --target T1/9.5.h5ad --reference T1/8.5.h5ad
```

两点提醒：

1. `--reference`（T3 是 `--wt`）不传会默认成 `--input` 自己，此时 `de_score` 和
   `de_direction` 恒为 0。这**只在测 copy_last / wt_identity 这类「不做改变」的基线时才正确**，
   给真模型打分一定要显式传真的参照 stage。
2. veckit 对 `--target` 是照单全收的。赛方建板时会把每个 stage **下采样到 10%** 再比，
   目标还要对半切（一半当真值、一半用来估 ceiling）。本地要跟官方尺度可比的话，
   这一步得自己做。

---

## 出处

本文每一条都来自官网，主要是
[/challenge/data](https://virtualembryo.ai/challenge/data)、
[/challenge/evaluation](https://virtualembryo.ai/challenge/evaluation)、
[/challenge/rules](https://virtualembryo.ai/challenge/rules)。
有出入以官网为准，也欢迎开 issue 指出。

本仓库**不包含任何比赛数据**。`panels/` 下的基因表是官网公开发布、无需登录即可下载的文件。
比赛数据的使用条款见 [Terms of Use 第 14 条](https://virtualembryo.ai/challenge/terms)：
数据尚未发表，任何对外使用需事先获得 UCSD Neil Chi 组的批准。

MIT License。
