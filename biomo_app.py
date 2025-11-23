import streamlit as st
import os
import json
import time
from pathlib import Path
from dataclasses import asdict, is_dataclass
from typing import List, Dict, Any

# 尝试导入 biomo_agent
try:
    from biomo_agent.graph import build_graph
    from biomo_agent.datastore.run_manager import create_run_dir
except ImportError:
    import sys
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from biomo_agent.graph import build_graph
    from biomo_agent.datastore.run_manager import create_run_dir

# 设置页面配置
st.set_page_config(
    page_title="Biomo Agent",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# 初始化 Session State
if "messages" not in st.session_state:
    st.session_state.messages = []
    
@st.cache_resource
def get_agent_graph():
    """缓存编译好的 LangGraph 对象"""
    return build_graph()

def _state_to_dict(state):
    if state is None:
        return {}
    if isinstance(state, dict):
        return state
    if is_dataclass(state):
        return asdict(state)
    return getattr(state, "__dict__", {})

def init_sidebar():
    with st.sidebar:
        st.header("⚙️ 配置")
        
        # 移除了 API Key 输入框
        
        data_path = st.text_input("h5ad 数据路径 (绝对路径)", value="", help="例如: /home/user/data.h5ad")
        project_id = st.text_input("Project ID", value="demo")
        
        st.divider()
        is_ready = False
        if not os.environ.get("DASHSCOPE_API_KEY"):
            st.error("🔴 环境变量 DASHSCOPE_API_KEY 未设置")
            st.info("请在启动 Streamlit 前，在终端运行: export DASHSCOPE_API_KEY=your_key")
        elif not data_path:
            st.warning("🟡 请输入数据路径")
        elif not os.path.exists(data_path):
            st.error(f"🔴 文件不存在: {data_path}")
        else:
            st.success(f"🟢 Ready: {Path(data_path).name}")
            is_ready = True
            
        if st.button("清除对话历史"):
            st.session_state.messages = []
            st.rerun()
            
        return data_path, project_id, is_ready

def scan_artifacts(run_dir: Path) -> List[Dict[str, Any]]:
    """扫描 run_dir 下的生成物"""
    artifacts = []
    if not run_dir or not run_dir.exists():
        return artifacts
    
    # 1. Images (output/figures)
    figures_dir = run_dir / "output" / "figures"
    if figures_dir.exists():
        for img_path in figures_dir.glob("*.png"):
            # image 类型本身不直接支持 content 下载，但在渲染时可以读取
            artifacts.append({
                "type": "image", 
                "path": str(img_path), 
                "name": img_path.name
            })
            
    # 2. Data (output/data) -> 支持 h5ad 下载
    data_dir = run_dir / "output" / "data"
    if data_dir.exists():
        for data_file in data_dir.glob("*.h5ad"):
            artifacts.append({
                "type": "data",
                "path": str(data_file),
                "name": data_file.name
            })

    # 3. Code (根目录)
    code_path = run_dir / "analysis.py"
    if code_path.exists():
        try:
            content = code_path.read_text(encoding="utf-8")
            artifacts.append({"type": "code", "content": content, "path": str(code_path), "name": "analysis.py"})
        except:
            pass

    # 4. Report (根目录)
    report_path = run_dir / "report.md"
    if report_path.exists():
        try:
            content = report_path.read_text(encoding="utf-8")
            artifacts.append({"type": "report", "content": content, "path": str(report_path), "name": "report.md"})
        except:
            pass
            
    return artifacts

def render_message(msg):
    """渲染单条消息"""
    with st.chat_message(msg["role"]):
        # 渲染主文本内容
        if msg.get("content"):
            st.markdown(msg["content"])
        
        # 渲染生成的 Artifacts
        artifacts = msg.get("artifacts", [])
        if artifacts:
            images = [a for a in artifacts if a["type"] == "image"]
            reports = [a for a in artifacts if a["type"] == "report"]
            codes = [a for a in artifacts if a["type"] == "code"]
            data_files = [a for a in artifacts if a["type"] == "data"]
            
            # 1. Images 直接展示 (更直观)
            if images:
                cols = st.columns(min(len(images), 2))
                for idx, img in enumerate(images):
                    with cols[idx % 2]:
                        # use_column_width -> use_container_width
                        st.image(img["path"], caption=img["name"], use_container_width=True)
            
            # 2. 详细文件与代码
            if reports or codes or data_files or images:
                with st.expander("📦 下载与详情 (代码/数据/图片)", expanded=False):
                    
                    # 报告下载
                    for rep in reports:
                        st.markdown("#### 📄 分析报告")
                        # st.markdown(rep["content"]) # 既然外面已经渲染了文本，这里就不重复渲染了，只提供下载
                        st.download_button("⬇️ 下载报告 (MD)", rep["content"], file_name="report.md", key=f"dl_rep_{msg.get('run_dir', '')}_{rep['name']}")
                        st.divider()

                    # 图片下载
                    if images:
                        st.markdown("#### 🖼️ 图片下载")
                        for img in images:
                            # 读取二进制
                            try:
                                with open(img["path"], "rb") as f:
                                    img_bytes = f.read()
                                st.download_button(f"⬇️ 下载图片 ({img['name']})", img_bytes, file_name=img['name'], mime="image/png", key=f"dl_img_{msg.get('run_dir', '')}_{img['name']}")
                            except Exception:
                                st.error(f"无法读取图片: {img['name']}")
                        st.divider()

                    # 数据下载
                    if data_files:
                        st.markdown("#### 💾 数据下载")
                        for df in data_files:
                            # 读取二进制 (h5ad可能很大，但在 streamlit server 端读取通常很快，只要不传给浏览器渲染)
                            # Streamlit download button 会把数据传给前端，大文件可能会慢
                            try:
                                with open(df["path"], "rb") as f:
                                    data_bytes = f.read()
                                st.download_button(f"⬇️ 下载数据 ({df['name']})", data_bytes, file_name=df['name'], key=f"dl_data_{msg.get('run_dir', '')}_{df['name']}")
                            except Exception as e:
                                st.error(f"无法读取数据文件: {e}")
                        st.divider()
                        
                    # 代码下载
                    for code in codes:
                        st.markdown("#### 🐍 分析代码")
                        st.code(code["content"], language="python")
                        st.download_button("⬇️ 下载代码 (PY)", code["content"], file_name="analysis.py", key=f"dl_code_{msg.get('run_dir', '')}_{code['name']}")

def run_agent_stream(user_query, data_path, project_id):
    """运行 Agent 并实时更新 UI"""
    
    # 创建 run 目录
    run_dir = create_run_dir(project_id=project_id)
    
    # 初始化 Agent
    app = get_agent_graph()
    initial_state = {
        "user_query": user_query,
        "data_path": data_path,
        "run_dir": run_dir,
        "run_id": run_dir.name,
    }
    
    status_container = st.status("🚀 Agent 正在启动...", expanded=False)
    
    final_response_text = ""
    artifacts = []
    
    try:
        for step in app.stream(initial_state):
            for node_name, node_state in step.items():
                s = _state_to_dict(node_state)
                
                if node_name == "intent":
                    tasks = s.get("tasks", [])
                    status_container.update(label=f"🎯 Intent: {tasks}", state="running")
                    status_container.write(f"解析到任务：{tasks}")
                    
                elif node_name == "plan":
                    plan = s.get("plan", {})
                    status_container.update(label="📋 Plan: 计划已生成", state="running")
                    status_container.write("生成分析计划：")
                    status_container.json(plan, expanded=False)
                    
                elif node_name == "codegen":
                    status_container.update(label="💻 Codegen: 代码已生成", state="running")
                    status_container.write("代码生成完毕，准备执行...")
                    
                elif node_name == "execute":
                    status_container.update(label="⚙️ Execute: 执行分析中...", state="running")
                    results = s.get("results", {})
                    if results.get("returncode") == 0:
                        status_container.write("✅ 执行成功")
                    else:
                        status_container.error(f"❌ 执行失败: {results.get('stderr')}")
                        
                elif node_name == "report":
                    status_container.update(label="📝 Report: 报告生成完毕", state="running")
                    if s.get("report"):
                        final_response_text = s.get("report")

        status_container.update(label="✅ 分析流程完成", state="complete", expanded=False)
        
        # 扫描生成的文件
        artifacts = scan_artifacts(run_dir)
        
        if not final_response_text:
            final_response_text = "分析已完成。结果如下："

    except Exception as e:
        status_container.update(label="❌ 发生错误", state="error")
        st.error(f"Agent 运行出错: {str(e)}")
        final_response_text = f"执行过程中发生错误: {str(e)}"
    
    return final_response_text, str(run_dir), artifacts

def main():
    data_path, project_id, is_ready = init_sidebar()
    
    st.title("🧬 Biomo Agent Chat")
    
    if data_path:
        st.caption(f"Project: {project_id} | Data: {Path(data_path).name}")
    else:
        st.caption("请先在左侧配置数据路径")

    # 1. 渲染历史消息
    for msg in st.session_state.messages:
        render_message(msg)

    # 2. 处理用户输入
    if prompt := st.chat_input("请输入您的分析需求...", disabled=not is_ready):
        # 添加用户消息
        st.session_state.messages.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)

        # Agent 回复
        with st.chat_message("assistant"):
            response_text, run_dir_str, artifacts = run_agent_stream(prompt, data_path, project_id)
            
            # 显示最终结果
            if response_text:
                st.markdown(response_text)
            
            # 显示 artifacts (本次)
            if artifacts:
                # 图片直接显示
                images = [a for a in artifacts if a["type"] == "image"]
                if images:
                    cols = st.columns(min(len(images), 2))
                    for idx, img in enumerate(images):
                        with cols[idx % 2]:
                            # use_column_width -> use_container_width
                            st.image(img["path"], caption=img["name"], use_container_width=True)
                
                # 其他折叠显示下载
                with st.expander("📦 下载与详情 (代码/数据/图片)", expanded=False):
                    # 这里可以重用 render_message 里的逻辑，但为了简单直接在这写
                    data_files = [a for a in artifacts if a["type"] == "data"]
                    codes = [a for a in artifacts if a["type"] == "code"]
                    reports = [a for a in artifacts if a["type"] == "report"]
                    
                    # Code
                    for code in codes:
                        st.download_button(f"⬇️ 下载代码 {code['name']}", code["content"], file_name=code['name'])
                    
                    # Data
                    for df in data_files:
                        try:
                            with open(df["path"], "rb") as f:
                                st.download_button(f"⬇️ 下载数据 {df['name']}", f.read(), file_name=df['name'])
                        except: pass
                    
                    # Images (download)
                    for img in images:
                        try:
                            with open(img["path"], "rb") as f:
                                st.download_button(f"⬇️ 下载图片 {img['name']}", f.read(), file_name=img['name'], mime="image/png")
                        except: pass

            # 保存到 Session State
            st.session_state.messages.append({
                "role": "assistant",
                "content": response_text,
                "run_dir": run_dir_str,
                "artifacts": artifacts
            })

if __name__ == "__main__":
    main()
