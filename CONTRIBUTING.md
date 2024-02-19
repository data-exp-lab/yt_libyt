# Contributing

## Report Bugs and Submit Feedback
Report bugs at https://github.com/data-exp-lab/yt_libyt/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Setting Up Development Environment

### Fork and Clone the `yt_libyt` Repo

1. Fork the `yt_libyt` repo on GitHub.
2. Clone your fork locally:
  ```bash
  git clone https://github.com/<your-github-account>/yt_libyt.git
  ```
3. Create a branch for local development:
  ```bash
  git checkout -b name-of-your-bugfix-or-feature
  ```

### Using `tox` to Test, Do Code-Formatting, and Linting

We use [`tox`](https://tox.wiki/en/4.11.3/installation.html) to run:
  - Python unit test ([`pytest`](https://docs.pytest.org/en/7.4.x/))
  - Converting old string to f-string ([`flynt`](https://github.com/ikamensh/flynt#flynt---string-formatting-converter))
  - Code formatting ([`black`](https://black.readthedocs.io/en/stable/))
  - Sort import order ([`isort`](https://pycqa.github.io/isort/index.html))
  - Linting ([`flake8`](https://flake8.pycqa.org/en/latest/))

#### Python Unit test
```bash
tox
```

#### Converting Old String to F-String
```bash
tox -e fstring
```

#### Code Formatting
```bash
tox -e format
```

#### Sort Import Order
```bash
tox -e sort_import
```

#### Linting
```bash
tox -e lint
```

### Pre-Commit

We use [pre-commit](https://pre-commit.com/#install) to check code format and style before committing:
  - Converting old string to f-string ([`flynt`](https://github.com/ikamensh/flynt#flynt---string-formatting-converter))
  - Code formatting ([`black`](https://black.readthedocs.io/en/stable/))
  - Sort import order ([`isort`](https://pycqa.github.io/isort/index.html))
  - Linting ([`flake8`](https://flake8.pycqa.org/en/latest/))

Set up pre-commit for the first time:
```bash
pre-commit install
```

Check every file:
```bash
pre-commit run --all-files
```

Commit your changes and push your branch to GitHub, pre-commit will apply to staged files before committing:
```bash
git add .
git commit -m "Your detailed description of your changes."
git push origin name-of-your-bugfix-or-feature
```
