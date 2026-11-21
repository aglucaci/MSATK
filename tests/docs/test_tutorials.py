"""Tutorials must run as written and the generated pages must match their templates."""

from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
TUT = ROOT / "docs" / "tutorials"
CODE = re.compile(r"```(python|bash)\n(.*?)```", re.S)


def _blocks(text: str) -> list[tuple[str, str]]:
    return CODE.findall(text)


@pytest.mark.parametrize("template", sorted((TUT / "_src").glob("*.md")), ids=lambda p: p.stem)
def test_generated_page_matches_template(template):
    page = TUT / template.name
    assert page.exists(), f"Run: python docs/tutorials/build.py {template.stem}"
    assert _blocks(page.read_text()) == _blocks(template.read_text()), (
        f"{page.name} is stale; run: python docs/tutorials/build.py {template.stem}"
    )


@pytest.mark.skipif(
    any(importlib.util.find_spec(m) is None for m in ("matplotlib", "pandas", "torch")),
    reason="tutorials need matplotlib, pandas and torch",
)
def test_tutorials_execute():
    proc = subprocess.run(
        [sys.executable, str(TUT / "build.py"), "--check"], capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
