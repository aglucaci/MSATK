# Conda And Bioconda Packaging

MSATK is intended to be installable through PyPI and Conda/Bioconda.

Target user-facing install command:

```bash
mamba install -c bioconda -c conda-forge msatk
```

or:

```bash
conda install -c bioconda -c conda-forge msatk
```

## Local Recipe

The local recipe lives at:

```text
conda-recipe/meta.yaml
```

It builds from the repository checkout so maintainers can test Conda packaging before a PyPI release.

```bash
mamba create -n conda-build-env -c conda-forge conda-build boa anaconda-client
conda activate conda-build-env
conda build conda-recipe
```

Then test installation:

```bash
mamba create -n test-msatk --use-local msatk
conda activate test-msatk
msatk --help
msatk profile --help
```

## Bioconda Submission

After publishing MSATK to PyPI, use:

```text
conda-recipe/meta.bioconda.yaml
```

Copy it to `recipes/msatk/meta.yaml` in a fork of `bioconda/bioconda-recipes`, update the PyPI source SHA256, and open a pull request.

Recommended release path:

```text
PyPI first -> local Conda recipe -> Bioconda PR -> CI-tested Conda builds
```

## Grayskull

After PyPI publication, a draft recipe can also be generated with:

```bash
mamba create -n grayskull -c conda-forge grayskull
conda activate grayskull
grayskull pypi msatk
```

Then manually review dependencies, CLI tests, maintainers, summary, and license metadata.
