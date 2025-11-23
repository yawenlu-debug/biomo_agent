# biomo_agent/nodes/codegen_node.py
from textwrap import dedent
from ..state import AgentState
from ..config import get_llm
from ..datastore.run_manager import save_code
from ..llm_utils import call_llm_with_retry

#代码清洗
def clean_code(raw: str) -> str:
    """
    去掉 ```python / ``` 等 markdown 代码块，以及多余换行。
    """
    lines = raw.splitlines()
    cleaned = []
    for line in lines:
        if line.strip().startswith("```"):
            continue
        cleaned.append(line)
    return "\n".join(cleaned).strip()



#生成代码
def codegen_node(state: AgentState) -> AgentState:
    llm = get_llm()
    data_path = state.data_path
    run_dir = str(state.run_dir)

    # prompt：让模型只输出纯 Python 代码
    system = (
    "You are an expert Python bioinformatics assistant. "
    "You MUST output ONLY pure Python code. "
    "DO NOT include markdown, code fences, triple backticks, or explanations."
    )

    #system = (
    #    "You are an expert Python bioinformatics assistant. "
    #    "Generate a complete Python script using scanpy to analyze scRNA-seq data."
    #)
    user = dedent(f"""
    请帮我生成一段 Python 脚本，要求：

    1. 使用 scanpy 读取这个 h5ad 文件："{data_path}"
    2. 做基础预处理：过滤低基因数细胞和低细胞数基因、归一化、log1p、高变基因选择。
       在预处理和质控阶段：
       - 优先尝试根据基因名自动识别线粒体基因：`mt_mask = adata.var_names.str.upper().str.startswith("MT-")`
       - 如果 `mt_mask.sum() > 0`，则新增一列 `adata.var["mt"] = mt_mask`，并调用
         `sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)`
       - 如果没有检测到任何以 `MT-` 开头的基因，则调用
         `sc.pp.calculate_qc_metrics(adata, inplace=True)`，不要传 `qc_vars`，以避免 KeyError
       - 计算完 QC 指标后，生成两张 QC 图：
         * 准备要绘制的 QC 指标列表：先检查 `adata.obs.columns`，只有存在的列才加入列表（例如 `["n_genes_by_counts", "total_counts"]`，如果 `pct_counts_mt` 存在则加入）。
         * 使用 `sc.pl.violin` 绘制 QC 小提琴图，传入存在的指标列表，参数中设置 `show=False`。
           注意：在这一步绝对不要使用 `groupby='leiden'`，因为此时还没有进行聚类。
         * 使用 `sc.pl.scatter` 绘制 QC 散点图（例如 `total_counts` vs `n_genes_by_counts`），同样设置 `show=False`
         * 在绘图前要 `import matplotlib.pyplot as plt` 和 `import os`，并调用
           `os.makedirs("{run_dir}/output/figures", exist_ok=True)` 确保目录存在
         * 小提琴图保存到 `"{run_dir}/output/figures/qc_violin.png"`，
           散点图保存到 `"{run_dir}/output/figures/qc_scatter.png"`，
           都使用 `plt.savefig(..., dpi=300, bbox_inches="tight")` 并在保存后调用 `plt.close()`
    3. 做 PCA、构建邻居图，使用 leiden 聚类，计算 UMAP
    4. 把聚类结果保存在 `adata.obs["leiden"]` 中（注意转成字符串类型以便绘图）。
    5. 使用 matplotlib 保存 UMAP 聚类图到 "{run_dir}/output/figures/umap_cluster.png"：
       - 复用前面导入的 `matplotlib.pyplot as plt` 和 `os`
       - 调用 `sc.pl.umap(adata, color="leiden", show=False)`（不要使用 `save=` 参数）
       - 再用 `plt.savefig("{run_dir}/output/figures/umap_cluster.png", dpi=300, bbox_inches="tight")` 保存图像
       - 最后 `plt.close()`
    6. 将处理后的 adata 写回到 "{run_dir}/output/data/processed.h5ad"
    7. 在关键步骤（读取数据、预处理、QC 完成、PCA/聚类/UMAP 完成、结果保存完成）适当添加 `print(...)` 日志，方便在命令行中查看脚本执行进度，例如：
       - `print("Step 1: 读取数据完成")`
       - `print("Step 2: 预处理和 QC 完成")`
       - `print("Step 3: PCA / 邻居图 / 聚类 / UMAP 完成")`
       - `print("Step 4: 结果与图像已保存")`

    请只输出可以直接运行的 Python 代码，不要解释，不要 markdown。
    最后请严格遵循：只输出纯 Python 代码，不要 ```python 或任何 markdown。
    """)

    #resp = llm.invoke([("system", system), ("user", user)])
    resp = call_llm_with_retry([("system", system), ("user", user)])
    raw_code = resp.content  # 如果用的是新版 ChatModel，这里可能是 resp.content / resp[0].text
    code = clean_code(raw_code)
    state.code = code
    if state.run_dir:
        save_code(code, state.run_dir)
    return state
