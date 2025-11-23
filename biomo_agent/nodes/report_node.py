from pathlib import Path
from textwrap import dedent
from ..state import AgentState
from ..llm_utils import call_llm_with_retry
from ..datastore.run_manager import save_report

def report_node(state: AgentState) -> AgentState:
    """
    根据执行结果生成 Markdown 报告。
    """
    if state.error:
        state.report = f"# Execution Failed\n\nError: {state.error}"
        return state

    run_dir = state.run_dir
    results = state.results
    
    # 收集生成了哪些图片
    # 修正路径：从 run_dir/figures 改为 run_dir/output/figures
    figures_dir = run_dir / "output" / "figures"
    figures = []
    if figures_dir.exists():
        figures = [f.name for f in figures_dir.glob("*.png")]

    # 构造 Prompt
    system = (
        "You are a bioinformatics expert. "
        "You will be given the execution logs and a list of generated figures from a single-cell analysis pipeline. "
        "Your goal is to write a professional Markdown report summarizing the analysis."
    )

    user = dedent(f"""
    User Query: "{state.user_query}"
    
    Execution stdout:
    {results.get('stdout', '')}
    
    Execution stderr:
    {results.get('stderr', '')}
    
    Generated Figures (found in 'output/figures/' folder):
    {', '.join(figures) if figures else 'None'}
    
    Please write a report in Markdown format with the following structure:
    
    # Analysis Report
    
    ## Executive Summary
    (Briefly summarize what was done and if it succeeded)
    
    ## Key Results
    (Summarize key findings based on logs. If logs are sparse, just mention the steps completed.)
    
    ## Visualizations
    (For each generated figure, create a subsection. 
     Use standard markdown image syntax: `![Description](output/figures/filename.png)`.
     Add a brief explanation of what this plot typically shows in scRNA-seq analysis.)
    
    ## Conclusion
    (Final wrap-up)
    
    ## Next Steps suggestion
    (Suggest 3 logical next steps as a JSON list inside a ```json block at the very end)
    """)

    resp = call_llm_with_retry([("system", system), ("user", user)])
    report_content = resp.content

    state.report = report_content
    
    if run_dir:
        save_report(report_content, run_dir)
        
    return state
