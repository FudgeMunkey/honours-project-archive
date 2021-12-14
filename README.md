# Electronic Archive

This repository is the electronic archive for Alexander Horvat's 2021 Honours thesis for the Bachelor of Advanced Computing.

## Usage

Clones this repo to your local machine and install the dependencies.

```Bash
git clone <url>
python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
```

To run the test suite,

```Bash
pytest
```

To search for and only run certain files/classes/tests,

```Bash
pytest -k "<name>"
```

## Test Generation

Scripts for test generation were only made for the GSW package. 

1. Copy the GSW documentation into `all_functions.txt` ([link](https://teos-10.github.io/GSW-Python/gsw_flat.html))
2. Run `python3 stub_documentation_functions.py` in `tests/gsw/test_generation` to stub tests for all of the functions in the documentation
3. Edit as needed (i.e. define the strategies)
4. Run using `pytest`