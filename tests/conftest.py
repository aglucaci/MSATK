"""Shared pytest fixtures for MSATK."""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def data_dir() -> Path:
    return Path(__file__).parent / "data"
