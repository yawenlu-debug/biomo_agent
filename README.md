# Biomo Agent

Biomo Agent 是一个专为生物信息学数据分析设计的智能体系统。它结合了大语言模型（LLM）与自动化工作流，允许用户通过自然语言对话来执行复杂的生物数据分析任务。

Agent 能够自动解析用户意图、制定分析计划、生成并执行 Python 代码、最终输出分析报告及相关图表数据。

## 🚀 核心功能

- **自然语言交互**：通过对话形式提出分析需求（例如："帮我对这个数据做预处理并聚类"）。
- **全自动工作流**：基于 LangGraph 构建，包含意图识别 (Intent)、计划制定 (Plan)、代码生成 (Codegen)、代码执行 (Execute) 和报告生成 (Report) 等环节。
- **多模态产物**：自动生成分析脚本 (`analysis.py`)、可视化图表 (`figures/`)、处理后的数据 (`.h5ad`) 和 Markdown 报告。
- **双重交互界面**：提供直观的 Web UI 和轻量级的命令行 CLI。

## 📂 项目结构

```text
jhupload/
├── biomo_app.py          # Streamlit Web 应用程序入口
├── biomo_cli.py          # 命令行交互脚本入口
├── requirements.txt      # 项目依赖列表
├── biomo_agent/          # Agent 核心逻辑包
│   ├── graph.py          # LangGraph 工作流图构建
│   ├── config.py         # 配置管理
│   ├── state.py          # 状态定义
│   ├── nodes/            # 工作流节点实现 (codegen, execute等)
│   └── datastore/        # 数据存储与运行目录管理
└── README.md             # 项目说明文档
```

## 🛠️ 安装与配置

### 1. 安装依赖
确保你的 Python 环境（建议 Python 3.9+）已安装必要的依赖库：

```bash
pip install -r requirements.txt
```

### 2. 设置环境变量
本项目依赖阿里云通义千问（DashScope）模型服务，运行前**必须**配置 API Key：

**Linux/macOS:**
```bash
export DASHSCOPE_API_KEY="sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
```

**Windows PowerShell:**
```powershell
$env:DASHSCOPE_API_KEY="sk-xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
```

## 💻 运行说明

本项目提供两种运行方式，可根据需求选择。

### 方式一：Web 界面 (Streamlit)
提供完整的可视化交互，支持查看对话历史、实时状态监控以及下载生成的图表和代码。

**启动命令：**
```bash
streamlit run biomo_app.py
```

**使用步骤：**
1. 在左侧边栏输入 `.h5ad` 文件的绝对路径。
2. 输入 Project ID（可选）。
3. 在对话框中输入你的分析需求。
4. 等待 Agent 运行完成，即可在界面上预览和下载结果。

### 方式二：命令行工具 (CLI)
适合在服务器终端快速调试或运行，无图形界面依赖。

**启动命令：**
```bash
python biomo_cli.py
```

**交互指令：**
- 启动后按照提示输入数据路径和 Project ID。
- 输入 `修改数据路径` 或 `修改项目id` 可在运行时变更配置。
- 输入 `q` 或 `exit` 退出程序。

## 📝 输出产物

每次分析运行都会生成一个独立的运行目录（基于 Project ID），包含以下内容：
- `analysis.py`: Agent 编写的具体分析代码。
- `report.md`: 最终的分析报告。
- `output/figures/`: 生成的图片文件（如 UMAP, Heatmap 等）。
- `output/data/`: 处理后的数据文件。
