# Installing the python implementation of isONclust

The Rust port is the recommended way to run isONclust; see the
[README](README.md). This file keeps the original python installation
instructions, which are still the reference the port is checked against.

isONclust is distributed as a python package supported on Linux / OSX with
python v>=3.4 as of version 0.0.2 and above (due to updates in python's
multiprocessing library).

### Using conda
Conda is the preferred way to install isONclust.

1. Create and activate a new environment called isonclust

```
conda create -n isonclust python=3 pip 
source activate isonclust
```

2. Install isONclust 

```
pip install isONclust
```
3. You should now have 'isONclust' installed; try it:
```
isONclust --help
```

Upon start/login to your server/computer you need to activate the conda environment "isonclust" to run isONclust as:
```
source activate isonclust
```

### Using pip 

To install isONclust, run:
```
pip install  isONclust
```
`pip` will install the dependencies automatically for you. `pip` is pythons official package installer and is included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 
### Downloading source from GitHub

#### Dependencies

Make sure the below listed dependencies are installed (installation links below). Versions in parenthesis are suggested as isONclust has not been tested with earlier versions of these libraries. However, isONclust may also work with earliear versions of these libaries.
* [parasail](https://github.com/jeffdaily/parasail-python)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)

In addition, please make sure you use python version >=3.4. isONclust will not work with python 2.

With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/isONclust.git
cd isONclust
./isONclust
```
