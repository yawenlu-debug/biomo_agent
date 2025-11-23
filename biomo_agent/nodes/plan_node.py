# biomo_agent/nodes/plan_node.py
from typing import Any, Dict
from ..state import AgentState
from ..datastore.run_manager import save_plan

def plan_node(state: AgentState) -> AgentState:
    steps = []
    if "preprocess" in state.tasks:
        steps.append({"id": "preprocess", "type": "sc_preprocess"})
    if "cluster" in state.tasks:
        steps.append({"id": "cluster", "type": "sc_cluster", "depends_on": ["preprocess"]})

    plan: Dict[str, Any] = {
        "modality": state.modality,
        "steps": steps,
    }
    state.plan = plan

    if state.run_dir:
        save_plan(plan, state.run_dir)
    return state
