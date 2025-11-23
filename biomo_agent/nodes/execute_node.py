# biomo_agent/nodes/execute_node.py
import subprocess
import sys
from pathlib import Path
from ..state import AgentState

def execute_node(state: AgentState) -> AgentState:
    if not state.run_dir:
        state.error = "run_dir not set"
        return state

    code_path = Path(state.run_dir) / "analysis.py"
    if not code_path.exists():
        state.error = f"analysis.py not found at {code_path}"
        return state

    print(f"\n>> [Execute] Starting execution of {code_path.name}...\n")

    # 使用 Popen 实时流式输出 stdout/stderr 到终端，同时也收集起来
    process = subprocess.Popen(
        ["python", str(code_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1  # 行缓冲
    )

    stdout_lines = []
    stderr_lines = []

    # 简单的实时打印循环（注：标准做法需要两个线程或 select，简单场景下轮询即可，或者直接让 stdout=None 打印到父进程终端）
    # 这里为了最简单的“看到输出”，我们采用直接打印到终端的方式 (stdout=None)，
    # 但这样会导致 capture 变得困难。
    # 
    # 更好的折中方案：使用 communicate() 之前，无法同时流式打印和捕获。
    # 鉴于用户更看重“实时看到流程”，我们改用简单的“直通模式”：
    # stdout, stderr 直接继承父进程，这样你会直接在屏幕看到所有输出。
    # 缺点：state.results["stdout"] 会是空的，报告里可能就没有日志了。
    
    # 方案优化：使用 tee 的逻辑（一边读一边打）。
    
    # 这里采用最稳妥的实时打印方案：
    # 逐行读取 stdout 并打印
    with process.stdout:
        for line in iter(process.stdout.readline, ''):
            print(line, end='')  # 实时打印到终端
            stdout_lines.append(line)
            
    # 等待结束
    returncode = process.wait()
    
    # 读取剩下的 stderr
    stderr_content = process.stderr.read()
    if stderr_content:
        print(stderr_content, file=sys.stderr) # 实时打印 stderr
        stderr_lines.append(stderr_content)

    state.results["stdout"] = "".join(stdout_lines)
    state.results["stderr"] = "".join(stderr_lines)
    state.results["returncode"] = returncode

    if returncode != 0:
        state.error = state.results["stderr"]
        print(f"\n>> [Execute] Failed with return code {returncode}")
    else:
        print(f"\n>> [Execute] Finished successfully.")

    return state
