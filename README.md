# libconviqt

A library implementing beam convolution on a sphere based on

G. Prézeau and M. Reinecke:  
*Algorithm for the Evaluation of Reduced Wigner Matrices*,  
APJS **190** (2010) 267, [arXiv:1002.1050](https://arxiv.org/abs/1002.1050)  

## Installation

To compile, test and install the C++ library with C-bindings:
```bash
./autogen.sh
./configure
make && make check && make install
```

To install and test the Python wrapper
```bash
cd python
python setup.py install
python setup.py test
```
