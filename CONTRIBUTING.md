# Contributing
---

Whether you are a novice or experienced software developer, all contributions
and suggestions are welcome!

## Getting Started

If you are looking to contribute to the *GeoParse* codebase, the best place to
start is the GitHub [issues](https://github.com/geoparse/geoparse/issues) tab.
This is also a great place for filing bug reports and making suggestions for
ways in which we can improve the code and documentation.


## Contributing to the Codebase
First, create your own fork of GeoParse on [GitHub](https://github.com/geoparse/geoparse/fork). 
Then, clone the forked repository to your local machine:

```bash
git clone https://github.com/<your-username>/geoparse.git
```

Once you have obtained a copy of the code, create a separate development environment 
to make and test changes without affecting your existing Python environment. 
This ensures your work environment remains clean and isolated.

### Environment Setup
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
cd geoparse
uv sync --locked
```
This ensures all dependencies match the locked versions.
Finally you need to activate the Python virtual environment.

```sh
source .venv/bin/activate
```
#### Adding New Dependencies
To add new dependencies to the project, first alter the `pyproject.toml` file. 
It is necessary to add it to the `dependencies` or `dev-dependencies` section.

Then to sync the dependencies, run the following command:
```bash
uv sync
```

#### Set up `pre-commit`
GeoParse uses [pre-commit](https://pre-commit.com/) to ensure that code
standard checks pass locally before pushing to the remote repository. 
This project adheres to [PEP 8](https://peps.python.org/pep-0008/) guidelines
and uses [Ruff](https://github.com/astral-sh/ruff) for linting and formatting.

To set up Ruff as a pre-commit hook and enforce style consistency, run:
```sh
pre-commit install
```

After installation, `Ruff` linting and formatting checks will automatically run with every commit.
To verify that everything is working correctly, run:

```bash
pre-commit run --all
```

#### Build Documentation Locally
```bash
cd docs
make html
```

#### Jupyter Notebooks
If you would like to contribute by working on the Jupyter Notebook tutorials, please ensure that you clear the output of all cells before submitting your changes. This helps to reduce the file size and keeps the repository clean.

### Making Changes
Before modifying the codebase or documentation, create a new branch with:
```sh
git checkout -b <your-branch>
```
Make your changes in this branch before submitting a pull request.

### Running Tests
Before submitting changes, please ensure all tests pass.
```sh
pytest --cov
```
If adding new functionality, include corresponding tests.


### Raising a Pull Requests
Once your changes are ready to be submitted, make sure to push your changes to
your fork of the GitHub repo before creating a pull request.
We will review your contribution and may request updates before merging.

Thank you for helping improve this project!