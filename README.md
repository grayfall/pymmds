# pyMMDS

Extendible metric MDS in Python (akin to landmark MDS). This package implements 
an active metric MDS method described in "Principal Component and Correspondence 
Analysis in R" by Dr. Herve Abdi (2017) and is thus similar to package bios2mds 
written in the R language. The method allows new objects to be projected onto 
an existing eigenbasis defined by a set of active objects and is thus useful 
whenever one needs to have a stable and scalable MDS transformation.

### Requirements

- Python >= 3.5
- NumPy >= 1.14.0
- Pandas >= 0.22.0

### Installation

```
pip install git+https://github.com/grayfall/pymmds.git
```

### Usage

There is only one core object (`mmds.Space`) and two methods to consider:

- `Space.__init__(dm: pandas.DataFrame)` - takes a symmetric distance matrix 
of active (landmark) samples and creates an MDS space
- `Space.project(dm: pandas.DataFrame)` - takes a table of distances between 
any number of supplementary samples and all active samples and projects the 
former onto the initialised MDS space.

There is also a utility function `mmds.read_dm` that will help you read the
DMs in case you don't want to tweak `pandas.read_csv` yourself.  

All functions and methods are well-documented. You can either use Python's `help`,
iPython's `?` or whatever docstring tools your IDE provides to view the docs from 
within your development environment or read the 
[Wiki page](https://github.com/grayfall/pymmds/wiki).
