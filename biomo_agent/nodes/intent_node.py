# biomo_agent/nodes/intent_node.py
from biomo_agent.state import AgentState

def intent_node(state: AgentState) -> AgentState:
    q = state.user_query.lower()
    # 极简规则：看到 cluster 就做 preprocess+cluster
    tasks = []
    if "cluster" in q or "聚类" in q:
        tasks = ["preprocess", "cluster"]
    else:
        # 默认只预处理
        tasks = ["preprocess"]

    state.modality = "scRNA"  # 先写死，之后再区分 spatial
    state.tasks = tasks
    return state
