"""
streamlit_app/core/config_loader.py — 从 config/config.yaml 读取项目配置

UI 通过此模块获取:
  - 所有项目的 gse_list_csv（SraRunInfo 所在目录）
  - rawdata_dir（SRR.txt 扫描目录）
  - 物种
  - 模块开关
"""

import os

_CONFIG_PATH = "config/config.yaml"
_cache = {}


def _load_yaml(path):
    try:
        import yaml
        with open(path) as f:
            return yaml.safe_load(f) or {}
    except ImportError:
        # 简单 key:value 解析 fallback
        cfg = {}
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and ":" in line:
                    k, _, v = line.partition(":")
                    cfg[k.strip()] = v.strip()
        return cfg


def load_config(config_path=_CONFIG_PATH):
    """加载 config.yaml，返回 dict"""
    if not os.path.exists(config_path):
        return {}
    return _load_yaml(config_path)


def get_projects(config_path=_CONFIG_PATH):
    """返回 config.yaml 中的 projects 列表"""
    cfg = load_config(config_path)
    return cfg.get("projects", [])


def get_project_config(project_name, config_path=_CONFIG_PATH):
    """返回指定 project_name 的配置 dict，未找到返回 None"""
    for proj in get_projects(config_path):
        if proj.get("project_name") == project_name:
            return proj
    return None


def get_sra_info_dir(project_name, config_path=_CONFIG_PATH):
    """
    返回 SRR 信息文件所在目录（workflow/resources/{species}/）。
    从 rawdata_dir 获取，默认回退到 workflow/resources/homo。
    """
    proj = get_project_config(project_name, config_path)
    if not proj:
        return "workflow/resources/homo"
    return proj.get("rawdata_dir", "workflow/resources/homo")


def get_rawdata_dir(project_name, config_path=_CONFIG_PATH):
    proj = get_project_config(project_name, config_path)
    if not proj:
        return "workflow/resources/homo"
    return proj.get("rawdata_dir", "workflow/resources/homo")


def get_species(project_name, config_path=_CONFIG_PATH):
    proj = get_project_config(project_name, config_path)
    if not proj:
        return "homo"
    return proj.get("species", "homo")


def get_planning_cfg(config_path=_CONFIG_PATH):
    """返回 planning 子 dict"""
    cfg = load_config(config_path)
    return cfg.get("planning", {})
