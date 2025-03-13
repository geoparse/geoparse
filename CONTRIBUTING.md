# Contributing
Whether you are a novice or experienced software developer, all contributions
and suggestions are welcome!



## Getting Started
GeoParse uses [uv](https://docs.astral.sh/uv/) for Python package and project management.
If you don't have `uv` installed, you need to install it first.
```bash
# On Linux and macOS
curl -LsSf https://astral.sh/uv/install.sh | sh

```

```bash
# On Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

You can also update `uv` using the following command if you already have it installed:
```sh
uv self update
```
After cloning the repository, sync dependencies using:
```sh
uv sync --locked
```
This ensures all dependencies match the locked versions.

## 🎨 Code Style

GeoParse follows [PEP 8](https://peps.python.org/pep-0008/) style
and uses [Ruff](https://github.com/astral-sh/ruff) for linting and formatting.
Set up Ruff as a pre-commit hook to enforce style consistency:.

```sh
uv run pre-commit install
```

## 🧪 Running Tests
Before submitting changes, please ensure all tests pass.
```sh
uv run pytest --cov
```
If adding new functionality, include corresponding tests.

## <img src="https://upload.wikimedia.org/wikipedia/commons/3/38/Jupyter_logo.svg" alt="Jupyter Logo" width="30"> Jupyter Notebooks
If you'd like to contribute by working on the Jupyter Notebook files in the tutorial, please ensure that you clear the output of all cells before submitting your changes. This helps to reduce the file size and keeps the repository clean.

## 🚀 Raising a Pull Requests
Once your changes are ready, push them to your branch and open a pull request.
We will review your contribution and may request updates before merging.

For bug reports and feature requests, please use the GitHub Issues tab.

Thank you for helping improve this project! 🚀  



## Getting Started
If you are looking to contribute to the *GeoParse* codebase, the best place to
start is the GitHub [issues](https://github.com/geoparse/geoparse/issues) tab.
This is also a great place for filing bug reports and making suggestions for
ways in which we can improve the code and documentation.

## Contributing to the Codebase
First create your own fork of GeoParse, then clone it:

```bash
# replace <my-username> with your github username
git clone https://github.com/<my-username>/geoparse.git
```

Once you've obtained a copy of the code, create a development environment that's
separate from your existing Python environment so that you can make and test
changes without compromising your own work environment.

### Environment Setup
This project recommends using [uv](https://docs.astral.sh/uv) to manage the
development environment.

```bash
pip install uv
uv venv
uv sync --all-extras
source .venv/bin/activate
```

#### Build Documentation Locally
```bash
cd docs
make html
```

#### Adding New Dependencies
To add new dependencies to the project, first alter the `pyproject.toml` file. Then to sync the dependencies, run the following command:
```bash
uv sync
```

This will invoke `python scripts/generate_pip_deps_from_conda.py` to convert
`environment.yml` to a `requirements.in` file. This is so that you can install
your development environment using the package manager of your choise like
`pip`, `uv`, `conda`, etc.

Moreover to add new extra dependencies in `pyproject.toml`, it is necessary to
add it to the `dependencies` and `[project.optional-dependencies]` section.


#### Set up `pre-commit`

This project uses [pre-commit](https://pre-commit.com/) to ensure that code
standard checks pass locally before pushing to the remote project repo. Follow
the [installation instructions](https://pre-commit.com/#installation), then
set up hooks with `pre-commit install`. After, `black`, `pylint` and `mypy`
checks should be run with every commit.

Make sure everything is working correctly by running

```bash
pre-commit run --all
```

### Making Changes
Before making changes to the codebase or documentation, create a new branch with:

```bash
git checkout -b <my-branch>
```

We recommend following the branch-naming convention described in [Making Pull Requests](#making-pull-requests).

### DCO-signing Commits

This project enforces the [DCO](https://developercertificate.org/) standard for
contributions, which requires authors to sign off on their commits. This can be
done with the `-s` or `--signoff` flag:

```bash
git commit -s -m 'my commit'
```

Refer to [this guide](https://github.com/src-d/guide/blob/master/developer-community/fix-DCO.md#dco-is-missing)
to add sign-offs retroactivately.

### Run the Full Test Suite Locally

```bash
pytest --cov
```

This will run the full test suite that matches the tests that are in `test

### Run a Specific Test Suite Locally

The above command will run the tests in mamba virtual environments for all of
the supported pandera extras packages, versions of python, pandas, etc. To run
a test for a specific set of versions, first run:

```bash
nox --list
```

You should see an output like this:

```bash
...
* tests(extra='core', pydantic='1.10.11', python='3.9', pandas='2.1.1') -> Run the test suite.
* tests(extra='strategies', pydantic='1.10.11', python='3.9', pandas='2.1.1') -> Run the test suite.
* tests(extra='hypotheses', pydantic='1.10.11', python='3.9', pandas='2.1.1') -> Run the test suite.
...
```

Then run a specific test condition with:

```bash
nox -db uv -s "tests(extra='core', pydantic='1.10.11', python='3.9', pandas='2.1.1')"
```

### Project Releases

Releases are organized under [milestones](https://github.com/pandera-dev/pandera/milestones),
which are be associated with a corresponding branch. This project uses
[semantic versioning](https://semver.org/), and we recommend prioritizing issues
associated with the next release.

### Contributing Documentation

Maybe the easiest, fastest, and most useful way to contribute to this project
(and any other project) is to contribute documentation. If you find an API
within the project that doesn't have an example or description, or could be
clearer in its explanation, contribute yours!

You can also find issues for improving documentation under the
[docs](https://github.com/pandera-dev/pandera/labels/docs) label. If you have
ideas for documentation improvements, you can create a new issue [here](https://github.com/pandera-dev/pandera/issues/new?assignees=&labels=docs&template=documentation-improvement.md&title=)

This project uses Sphinx for auto-documentation and RST syntax for docstrings.
Once you have the code downloaded and you find something that is in need of some
TLD, take a look at the [Sphinx](https://www.sphinx-doc.org/en/1.0/rest.html)
documentation or well-documented
[examples](https://pandera.readthedocs.io/en/stable/_modules/pandera/schemas.html#DataFrameSchema)
within the codebase for guidance on contributing.

You can build the html documentation by running `nox -s docs`. The built
documentation can be found in `docs/_build`.

### Contributing Bugfixes

Bugs are reported under the [bug](https://github.com/pandera-dev/pandera/labels/bug)
label, so if you find a bug create a new issue [here](https://github.com/pandera-dev/pandera/issues/new?assignees=&labels=bug&template=bug_report.md&title=).

### Contributing Enhancements

New feature issues can be found under the
[enhancements](https://github.com/pandera-dev/pandera/labels/enhancement) label.
You can request a feature by creating a new issue [here](https://github.com/pandera-dev/pandera/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=).

### Making Pull Requests

Once your changes are ready to be submitted, make sure to push your changes to
your fork of the GitHub repo before creating a pull request. Depending on the
type of issue the pull request is resolving, your pull request should merge
onto the appropriate branch:

#### Bugfixes

- branch naming convention: `bugfix/<issue number>` or `bugfix/<bugfix-name>`
- pull request to: `main`

#### Documentation

- branch naming convention: `docs/<issue number>` or `docs/<doc-name>`
- pull request to: `release/x.x.x` branch if specified in the issue milestone, otherwise `main`

#### Enhancements

- branch naming convention: `feature/<issue number>` or `feature/<enhancement-name>`
- pull request to: `release/x.x.x` branch if specified in the issue milestone, otherwise `main`

We will review your changes, and might ask you to make additional changes
before it is finally ready to merge. However, once it's ready, we will merge
it, and you will have successfully contributed to the codebase!

### Questions, Ideas, General Discussion

Head on over to the [discussion](https://github.com/pandera-dev/pandera/discussions)
section if you have questions or ideas, want to show off something that you
did with `pandera`, or want to discuss a topic related to the project.

### Dataframe Schema Style Guides

We have guidelines regarding dataframe and schema styles that are encouraged
for each pull request.

If specifying a single column DataFrame, this can be expressed as a one-liner:

```python
DataFrameSchema({"col1": Column(...)})
```

If specifying one column with multiple lines, or multiple columns:

```python
DataFrameSchema(
    {
        "col1": Column(
            int,
            checks=[
                Check(...),
                Check(...),
            ]
        ),
    }
)
```

If specifying columns with additional arguments that fit in one line:

```python
DataFrameSchema(
    {"a": Column(int, nullable=True)},
    strict=True
)
```

If specifying columns with additional arguments that don't fit in one line:

```python
DataFrameSchema(
    {
        "a": Column(
            int,
            nullable=True,
            coerce=True,
            ...
        ),
        "b": Column(
            ...,
        )
    },
    strict=True)
```

## Deprecation policy

This project adopts a rolling policy regarding the minimum supported version of its dependencies, based on [NEP 29](https://numpy.org/neps/nep-0029-deprecation_policy.html):

- **Python**: 42 months
- **NumPy**: 24 months
- **Pandas**: 18 months

This means the latest minor (X.Y) version from N months prior. Patch versions (x.y.Z) are not pinned, and only the latest available at the moment of publishing the release is guaranteed to work.
