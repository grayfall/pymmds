# pyMMDS

Extendible metric MDS in Python. Classic metric MDS transforms a 
symmetric matrix of metric distances between N objects into a set of projections 
in an orthonormal space. The method is akin to PCA, though unlike PCA it doesn't 
provide a straightforward way to project new objects onto the space. This 
package implements an active metric MDS method described in "Principal Component 
and Correspondence Analysis in R" by Dr. Herve Abdi (2017) and is thus similar 
to package bios2mds written in the R language. The method allows new objects to 
be projected onto an existing space defined by a set of active objects.

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

### Requirements

- Python >= 3.5
- NumPy >= 1.14.0
- Pandas >= 0.22.0