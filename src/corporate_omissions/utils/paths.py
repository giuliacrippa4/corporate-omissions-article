from __future__ import annotations
from pathlib import Path

def repo_root() -> Path:
    # .../src/corporate_omissions/utils/paths.py -> repo root
    return Path(__file__).resolve().parents[3]

def data_dir() -> Path:
    return repo_root() / "data"

def raw_dir() -> Path:
    return data_dir() / "raw"

def processed_dir() -> Path:
    return data_dir() / "processed"

def outputs_dir() -> Path:
    return data_dir() / "outputs"

def figures_dir() -> Path:
    return outputs_dir() / "figures"

def tables_dir() -> Path:
    return outputs_dir() / "tables"

def ensure_dirs() -> None:
    for p in [raw_dir(), processed_dir(), outputs_dir(), figures_dir(), tables_dir()]:
        p.mkdir(parents=True, exist_ok=True)
