# run_demo.py
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# run_demo.py
from biomo_agent.graph import build_graph
from biomo_agent.state import AgentState
from biomo_agent.datastore.run_manager import create_run_dir

def main():
    # 1. 创建一个 run 目录
    project_id = "demo_project2"
    run_dir = create_run_dir(project_id)
    # 从目录名中推一个 run_id，例如 "demo_project_20251120_162200"
    run_id = run_dir.name

    # 2. 构造初始 state
    state = AgentState(
        user_query="请用 scanpy 对这个数据做基础预处理和聚类",
        data_path="/hpc2hdd/home/ylu752/jhupload/runs/singlecellexpress.h5ad",
        run_id=run_id,
        run_dir=run_dir,
    )

    # 3. 构建并运行图
    app = build_graph()
    final_state = app.invoke(state)

    # 4. 看执行结果
    # langgraph 默认返回的是 dict 状态，这里按 dict 方式读取
    results = final_state.get("results", {})
    error = final_state.get("error")

    print("=== Return code ===")
    print(results.get("returncode"))
    print("=== STDOUT ===")
    print(results.get("stdout", ""))
    print("=== STDERR ===")
    print(results.get("stderr", ""))
    print("=== Error field ===")
    print(error)

    #整体日志流程回放
    print("=== Agent Logs ===")
    for line in final_state.get("logs", []):
        print(line)

    # 5. 你还可以检查：
    # - run_dir / "analysis.py"   （LLM 生成的代码）
    # - run_dir / "figures/umap_cluster.png"
    # - run_dir / "adata_processed.h5ad"

if __name__ == "__main__":
    main()


