from textwrap import dedent
from ..state import AgentState
from ..config import get_llm
from ..datastore.run_manager import save_code
from ..llm_utils import call_llm_with_retry
from .codegen_node import clean_code

def fix_code_node(state: AgentState) -> AgentState:
    state.retry_count += 1
    print(f">> [Auto-Fix] Attempt {state.retry_count}/5: Fixing code based on error...")

    llm = get_llm()
    current_code = state.code
    error_msg = state.error
    
    system = (
        "You are an expert Python bioinformatics assistant. "
        "Your task is to fix the provided Python code based on the execution error. "
        "You MUST output ONLY the fixed pure Python code. "
        "DO NOT include markdown, code fences, or explanations."
    )
    
    user = dedent(f"""
    The following Python code failed to execute:
    
    ```python
    {current_code}
    ```
    
    The error message was:
    
    {error_msg}
    
    Please fix the code to resolve this error. 
    Ensure you keep the original logic (loading data, saving figures to correct paths, etc.) but fix the specific issue causing the crash.
    
    If the error is about missing columns (like 'pct_counts_mt'), please check if the column exists before using it, or skip that part of the plot.
    
    Output ONLY the full fixed Python code.
    """)
    
    resp = call_llm_with_retry([("system", system), ("user", user)])
    fixed_code = clean_code(resp.content)
    
    # Update state
    state.code = fixed_code
    state.error = None  # Clear error for next execution
    
    # Save the fixed code
    if state.run_dir:
        save_code(fixed_code, state.run_dir)
        
    return state

