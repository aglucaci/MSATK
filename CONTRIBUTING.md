# Contributing to MSATK

Thank you for your interest in contributing to MSATK.

MSATK is a Python and command-line toolkit for profiling multiple sequence alignments. Contributions are welcome as bug reports, documentation improvements, test datasets, feature requests, and code contributions.

## Development Setup

```bash
git clone https://github.com/yourname/msatk.git
cd msatk
python -m venv .venv
pip install -e ".[dev,docs,all]"
pre-commit install
```

## Running Tests

```bash
pytest
pytest tests/unit
pytest tests/integration
pytest -m "not slow"
pytest --cov=msatk --cov-report=term-missing
```

## Code Style

MSATK uses Ruff for linting and formatting, Mypy for optional type checking, and Pytest for tests.

Before opening a pull request, run:

```bash
ruff check src/msatk tests
ruff format --check src/msatk tests
mypy src/msatk
pytest
```

## Adding New Features

Please include unit tests, integration tests for CLI behavior, documentation updates, example usage when user-facing, and output schema notes when outputs change.

## Adding New Test Data

Test data should be small, synthetic, deterministic, documented, and placed under `tests/data/`. Do not add large biological datasets directly to the repository.

## Pull Request Checklist

- [ ] Tests pass locally
- [ ] New functionality has tests
- [ ] Documentation has been updated
- [ ] CLI behavior is documented if changed
- [ ] Output schema changes are described
- [ ] No large files were added accidentally
- [ ] The changelog has been updated if appropriate
