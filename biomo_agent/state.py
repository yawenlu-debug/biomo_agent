# biomo_agent/state.py
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

@dataclass
class AgentState:
    user_query: str
    data_path: str

    modality: Optional[str] = None
    tasks: List[str] = field(default_factory=list)
    plan: Dict[str, Any] = field(default_factory=dict)

    code: Optional[str] = None
    run_id: str = ""
    run_dir: Optional[Path] = None

    results: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None
    retry_count: int = 0
    report: Optional[str] = None
