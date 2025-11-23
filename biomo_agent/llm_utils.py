# biomo_agent/llm_utils.py
from typing import List, Tuple
from langchain_core.messages import BaseMessage
from .config import get_llm

def call_llm_with_retry(
    messages: List[Tuple[str, str]],
    max_retries: int = 10,
):
    llm = get_llm()
    last_err = None
    for i in range(max_retries):
        try:
            return llm.invoke(messages)
        except Exception as e:
            err_str = str(e)
            # 硅基流动“系统繁忙”这类报错就重试；其他错误直接抛出也可以
            if "System is too busy" in err_str or "503" in err_str:
                last_err = e
                continue
            raise
    # 重试多次仍失败，抛出最后一次错误
    raise last_err
