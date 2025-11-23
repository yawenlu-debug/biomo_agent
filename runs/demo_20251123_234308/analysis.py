import os
import scanpy as sc
import matplotlib.pyplot as plt

# Step 1: 读取数据
print("Step 1: 读取数据完成")
adata = sc.read_h5ad("/hpc2hdd/home/ylu752/jhupload/runs/singlecellexpress.h5ad")

# Step 2: 预处理和 QC
print("Step 2: 预处理和 QC 完成")

# 自动识别线粒体基因
mt_mask = adata.var_names.str.upper().str.startswith("MT-")
if mt_mask.sum() > 0:
    adata.var["mt"] = mt_mask
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
else:
    sc.pp.calculate_qc_metrics(adata, inplace=True)

# 准备要绘制的 QC 指标列表
qc_metrics = []
for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
    if col in adata.obs.columns:
        qc_metrics.append(col)

# 绘制小提琴图
os.makedirs("runs/demo_20251123_234308/output/figures", exist_ok=True)
plt.figure(figsize=(8, 6))
sc.pl.violin(adata, qc_metrics, groupby=None, show=False)
plt.savefig("runs/demo_20251123_234308/output/figures/qc_violin.png", dpi=300, bbox_inches="tight")
plt.close()

# 绘制散点图
plt.figure(figsize=(6, 6))
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
plt.savefig("runs/demo_20251123_234308/output/figures/qc_scatter.png", dpi=300, bbox_inches="tight")
plt.close()

# 过滤低质量细胞和基因
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 归一化、log1p、高变基因选择
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")

# Step 3: PCA / 邻居图 / 聚类 / UMAP 完成
print("Step 3: PCA / 邻居图 / 聚类 / UMAP 完成")

sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5)
adata.obs["leiden"] = adata.obs["leiden"].astype(str)

sc.tl.umap(adata)

# Step 4: 结果与图像已保存
print("Step 4: 结果与图像已保存")

# 保存 UMAP 聚类图
plt.figure(figsize=(8, 8))
sc.pl.umap(adata, color="leiden", show=False)
plt.savefig("runs/demo_20251123_234308/output/figures/umap_cluster.png", dpi=300, bbox_inches="tight")
plt.close()

# 保存处理后的 adata
os.makedirs("runs/demo_20251123_234308/output/data", exist_ok=True)
adata.write_h5ad("runs/demo_20251123_234308/output/data/processed.h5ad")