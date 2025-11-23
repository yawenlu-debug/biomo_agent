# biomo_agent/config.py
import os
from langchain_openai import ChatOpenAI

def get_llm():
    """
    使用通义千问 DashScope 的 OpenAI 兼容接口作为 LangChain 的 LLM。
    依赖环境变量 DASHSCOPE_API_KEY。
    """
    return ChatOpenAI(
        model="qwen-flash",  # 通义千问模型名，和你示例里保持一致
        base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
        openai_api_key=os.environ["DASHSCOPE_API_KEY"],  # 注意这里是环境变量名，不是 key 本身
        temperature=0.0,
    )
