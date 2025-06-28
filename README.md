# Reed-Solomon decoder in Python

## dependencies
- Python 3.12 (high or lower might very well work)
- Numpy and Galois packages

## Installation
clone the git repository
```bash
git clone https://github.com/Raaaainbow/Reed-Solomon-Python
```
Create the virtual enviroment
```
bash scripts/enviroment.sh
```
depending on you python version you can use
```
PYTHON_VERSION=python3.11 bash scripts/environment.sh

```
or use python (if available)
```
PYTHON_VERSION=python bash scripts/environment.sh
```
Activate the virtual enviroment
```
source .venv/bin/activate
```
## Features
This program opperates over a GF(2**m) Galois Field, and has been tested to correct 16 errors in a 255 symbol message. 

## Usage
The code includes 4 examples that you can choose by switching between True and False. You can also create your own messages by modifying the degrees and coefficients. Additionally, you can change how much gets printed by toggling the `verbose` variable.

## Reference
If you wish to know more of the theoretical basis for Reed-Solomon decoding, along with the method used in the repository, you can read the book 'A Course In Error-Correcting Codes' by Jørn Justesen & Tom Høholdt (ISBN: 3-03719-001-9)

## Contributions
This program was largely written by Victor K. Andersen under the name [Wickn](https://github.com/Wickn).

## License
This project is licensed under [LICENSE](LICENSE).
