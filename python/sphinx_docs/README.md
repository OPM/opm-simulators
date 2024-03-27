# Python scripts for building opm-simulators sphinx documentation

## Installation of the python scripts
- Requires python3 >= 3.10

### Using poetry
For development it is recommended to use poetry:

- Install [poetry](https://python-poetry.org/docs/)
- Then run:
```
$ poetry install
$ poetry shell
```

### Installation into virtual environment
If you do not plan to change the code, you can do a regular installation into a VENV:

```
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install .
```
