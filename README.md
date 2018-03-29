# pyMMDS

Extendible metric MDS in Python (akin to landmark MDS). Classic metric MDS 
transforms a symmetric matrix of metric distances between N objects into a 
set of projections in an orthonormal space. The method is akin to PCA, though 
unlike PCA it doesn't provide a straightforward way to project new objects 
onto the space. This package implements an active metric MDS method
described in "Principal Component and Correspondence Analysis in R" by Dr. 
Herve Abdi (2017) and is thus similar to package bios2mds written in the R 
language. The method allows new objects to be projected onto an existing space 
defined by a set of active objects.

### Requirements

- Python >= 3.5
- NumPy >= 1.14.0
- Pandas >= 0.22.0

### Installation

```
git clone https://github.com/grayfall/pymmds.git
cd pymmds
pip install .
```

or simply 

```
pip install git+https://github.com/grayfall/pymmds.git
```

### Usage

There is only one core object (`mmds.Space`) and two methods to consider:

- `Space.__init__` - takes a symmetric distance matrix of active (landmark)
samples and creates an MDS space
- `Space.project` - takes table of distances between any number of supplementary 
samples and all active samples and projects the former onto the initialised
MDS space.

There is also a utility function `mmds.read_dm` that will help you read the
DMs in case you don't want tweak `pandas.read_csv` yourself.  

All functions and methods are well-documented. You can either use Python's `help` and 
iPython's `?` to view the docs from within your development environment or read the 
[Wiki page](https://github.com/grayfall/pymmds/wiki).
