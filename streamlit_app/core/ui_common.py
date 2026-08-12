"""
streamlit_app/core/ui_common.py — 跨页面共享 UI 组件

提供项目选择器等公共组件，解决子页面直达时 session_state 为空的问题。
"""

import os
import sys

import streamlit as st

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from streamlit_app.core import config_loader as st_cfg
from streamlit_app.core import state as st_state


def render_project_selector() -> str | None:
    """
    在侧边栏渲染项目选择器，写入 st.session_state['project']。

    Returns:
        选中的 project 名称，无可用项目时返回 None。
    """
    with st.sidebar:
        cfg_projects = [
            p.get("project_name", "")
            for p in st_cfg.get_projects()
            if p.get("project_name")
        ]
        result_projects = st_state.list_projects()
        all_projects = sorted(set(cfg_projects + result_projects))

        if not all_projects:
            st.warning("config.yaml 中无 projects 配置，请先在「⚙️ 项目配置」中创建项目。")
            st.session_state["project"] = ""
            return None

        current = st.session_state.get("project", "")
        idx = all_projects.index(current) if current in all_projects else 0
        project = st.selectbox("项目", all_projects, index=idx, key="_sidebar_project_sel")
        st.session_state["project"] = project
        return project
