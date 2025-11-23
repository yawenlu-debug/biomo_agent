# biomo_agent/graph.py
from langgraph.graph import StateGraph, END
from .state import AgentState
from .datastore.run_manager import create_run_dir
from .nodes.intent_node import intent_node
from .nodes.plan_node import plan_node
from .nodes.codegen_node import codegen_node
from .nodes.execute_node import execute_node
from .nodes.fix_code_node import fix_code_node
from .nodes.report_node import report_node

def should_continue(state: AgentState):
    # 如果有错误，且重试次数没用完 -> 转去 fix_code 节点
    # 注意：state.retry_count 初始为 0，最大重试 5 次
    if state.error and state.retry_count < 5:
        return "fix_code"
    # 否则（成功或重试次数耗尽），去生成报告
    return "report"

def build_graph():
    graph = StateGraph(AgentState)

    graph.add_node("intent", intent_node)
    graph.add_node("plan", plan_node)
    graph.add_node("codegen", codegen_node)
    graph.add_node("execute", execute_node)
    graph.add_node("fix_code", fix_code_node)
    graph.add_node("report", report_node)

    graph.set_entry_point("intent")
    graph.add_edge("intent", "plan")
    graph.add_edge("plan", "codegen")
    graph.add_edge("codegen", "execute")
    
    # Conditional edge after execution
    graph.add_conditional_edges(
        "execute",
        should_continue,
        {
            "fix_code": "fix_code",
            "report": "report"
        }
    )
    
    # Fix code 完了之后，重新执行
    graph.add_edge("fix_code", "execute")
    
    graph.add_edge("report", END)

    return graph.compile()
