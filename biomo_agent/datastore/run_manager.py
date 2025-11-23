# biomo_agent/datastore/run_manager.py
from pathlib import Path
from datetime import datetime
import json

BASE_RUN_DIR = Path("runs")

def create_run_dir(project_id: str = "demo") -> Path:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = BASE_RUN_DIR / f"{project_id}_{ts}"
    run_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建规范的输出目录结构
    output_dir = run_dir / "output"
    (output_dir / "data").mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(parents=True, exist_ok=True)
    
    return run_dir

def save_plan(plan: dict, run_dir: Path):
    with open(run_dir / "plan.json", "w", encoding="utf-8") as f:
        json.dump(plan, f, indent=2, ensure_ascii=False)

def save_code(code: str, run_dir: Path) -> Path:
    code_path = run_dir / "analysis.py"
    with open(code_path, "w", encoding="utf-8") as f:
        f.write(code)
    return code_path

def save_report(report: str, run_dir: Path) -> Path:
    report_path = run_dir / "report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report)
    return report_path
