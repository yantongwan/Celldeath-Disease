from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.build_database import build_database


@pytest.fixture(scope="session")
def demo_db_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    workdir = tmp_path_factory.mktemp("rcd_atlas_demo")
    db_path = workdir / "rcd_atlas.duckdb"
    build_database(
        project_dir=PROJECT_ROOT,
        raw_dir=PROJECT_ROOT / "data" / "raw" / "demo",
        db_path=db_path,
        processed_dir=workdir / "processed",
        exports_dir=workdir / "exports",
        export_outputs=False,
    )
    return db_path


@pytest.fixture()
def conn(demo_db_path: Path):
    handle = duckdb.connect(str(demo_db_path), read_only=True)
    try:
        yield handle
    finally:
        handle.close()
