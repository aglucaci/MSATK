"""Report generation."""

from .html import render_html_report
from .markdown import render_markdown_report

__all__ = ["render_html_report", "render_markdown_report"]
