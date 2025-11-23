#!/usr/bin/env python
"""
简单的命令行交互脚本，用于调用 biomo_agent。
启动后通过问答方式设置 / 修改数据路径和项目 ID，并流式打印各阶段状态。

示例交互：

    请输入 h5ad 数据文件路径（必填）> /path/to/data.h5ad
    请输入项目 ID（回车默认为 demo）> my_project

进入主对话后，你可以输入自然语言需求，例如：

    我想对这个 h5ad 做预处理并聚类

输入 `exit` / `quit` / `q` 可以退出。

注意：
- 需要提前在环境中设置 DASHSCOPE_API_KEY，用于通义千问 DashScope 接口。
"""

import os
import json
from pathlib import Path
from dataclasses import asdict, is_dataclass

from biomo_agent.graph import build_graph
from biomo_agent.datastore.run_manager import create_run_dir


def _state_to_dict(state):
    """把可能是 dataclass / dict 的状态统一转成 dict。"""
    if state is None:
        return {}
    if isinstance(state, dict):
        return state
    if is_dataclass(state):
        return asdict(state)
    # 兜底：尝试用 __dict__
    return getattr(state, "__dict__", {})


def main():
    # 通过问答方式获取 data_path 与 project_id
    print("=== biomo_agent 命令行交互 ===")
    print("启动前先配置数据文件路径和项目 ID。")

    while True:
        data_path = input("请输入 h5ad 数据文件路径（必填）> ").strip()
        if not data_path:
            print("[提示] 数据路径不能为空。")
            continue
        if not os.path.exists(data_path):
            print(f"[错误] 数据文件不存在：{data_path}")
            continue
        break

    project_id = input("请输入项目 ID（回车默认为 demo）> ").strip() or "demo"

    # 构建 LangGraph 图（只需构建一次）
    app = build_graph()

    print("\n=== 配置完成，进入对话模式 ===")
    print(f"当前数据文件：{os.path.abspath(data_path)}")
    print(f"当前项目 ID：{project_id}")
    print("输入你的自然语言需求（例如：'帮我对这个数据做预处理和聚类'）")
    print("支持的控制指令：")
    print("- '修改数据路径' / '更换数据文件'：重新设置 data_path")
    print("- '修改项目id' / '更换项目id'：重新设置 project_id")
    print("输入 'exit' / 'quit' / 'q' 退出。")

    while True:
        try:
            user_query = input("\n问题> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n已退出。")
            break

        if not user_query:
            continue

        if user_query.lower() in {"exit", "quit", "q"}:
            print("已退出。")
            break

        # 特殊指令：修改数据路径
        if user_query in {"修改数据路径", "修改数据文件", "更换数据", "更换数据文件"}:
            while True:
                new_path = input("请输入新的 h5ad 数据文件路径> ").strip()
                if not new_path:
                    print("[提示] 数据路径不能为空。")
                    continue
                if not os.path.exists(new_path):
                    print(f"[错误] 数据文件不存在：{new_path}")
                    continue
                data_path = new_path
                print(f"[信息] 已更新数据文件：{os.path.abspath(data_path)}")
                break
            continue

        # 特殊指令：修改项目 ID
        if user_query in {"修改项目id", "修改项目ID", "更换项目id", "更换项目"}:
            new_project_id = input("请输入新的项目 ID（回车保持不变）> ").strip()
            if new_project_id:
                project_id = new_project_id
                print(f"[信息] 已更新项目 ID：{project_id}")
            else:
                print(f"[信息] 项目 ID 保持不变：{project_id}")
            continue

        # 为本次对话创建一个 run 目录
        run_dir: Path = create_run_dir(project_id=project_id)
        print(f"[信息] 本次运行目录：{run_dir.resolve()}")

        # 初始状态，交给 LangGraph
        initial_state = {
            "user_query": user_query,
            "data_path": data_path,
            "run_dir": run_dir,
            "run_id": run_dir.name,
        }

        print("[信息] 正在调用 biomo_agent（流式模式），请稍候...")
        final_state = None

        # 使用 LangGraph 的 stream 接口，逐节点打印状态
        for step in app.stream(initial_state):
            for node_name, node_state in step.items():
                s = _state_to_dict(node_state)

                if node_name == "intent":
                    tasks = s.get("tasks") or []
                    print(f"[intent] 解析到任务：{tasks}")
                elif node_name == "plan":
                    plan = s.get("plan") or {}
                    print("[plan] 生成分析计划：")
                    print(json.dumps(plan, ensure_ascii=False, indent=2))
                elif node_name == "codegen":
                    print(f"[codegen] 已生成分析代码：{(run_dir / 'analysis.py').resolve()}")
                elif node_name == "execute":
                    print("[execute] 正在执行分析脚本...")
                elif node_name == "report":
                    print(f"[report] 已生成分析报告：{(run_dir / 'report.md').resolve()}")

                final_state = node_state

        state_dict = _state_to_dict(final_state)

        error = state_dict.get("error")
        results = state_dict.get("results", {}) or {}

        if error:
            print("\n[错误] agent 执行失败：")
            print(error)
        else:
            print("\n[完成] agent 执行结束。")
            returncode = results.get("returncode")
            if returncode is not None:
                print(f"返回码：{returncode}")

            stdout = results.get("stdout")
            if stdout:
                print("\n[stdout]")
                print(stdout)

            stderr = results.get("stderr")
            if stderr:
                print("\n[stderr]")
                print(stderr)

            print("\n[输出文件]")
            print(f"- 计划文件：{(run_dir / 'plan.json').resolve()}")
            print(f"- 分析代码：{(run_dir / 'analysis.py').resolve()}")
            print(f"- 处理后的数据：{(run_dir / 'adata_processed.h5ad').resolve()}")
            print(f"- UMAP 聚类图：{(run_dir / 'figures' / 'umap_cluster.png').resolve()}")
            print(f"- 分析报告：{(run_dir / 'report.md').resolve()}")


if __name__ == "__main__":
    main()


