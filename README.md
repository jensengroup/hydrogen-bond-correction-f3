# Third-Generation Hydrogen-Bonding Correction

Third-Generation Hydrogen-Bonding Corrections
for Semiempirical QM Methods and Force Fields.

## Introduction

Fortran F90 (module) implementation of Martin Korths Hydrogen-bond correction term.


![equation](http://latex.codecogs.com/gif.latex?1%2Bsin%28mc%5E2%29%0D%0A)

% Comment

## Known Limitation


## Please cite as:

Martin Korth, J. Chem. Theory Comput., 2010, 6 (12), pp 3808–3816, DOI: 10.1021/ct100408b

## Installation

Use the make script. Change compiler if you use something else than `gfortran`.

    make

## Parameters

Default parameters are for the PM6-D3H+ model.

    nitrogen = -0.15
    oxygen   = -0.16

change it either via source code, or use a parameter file.


## Usage

### Energy

You will properly want to use the module as part of a QM/SQM package, but
if you just want to try it out, use it as.

    ./f3 -x <structure.xyz>

to calculate the correction energy of the XYZ structure.

### Gradient

If you want to calculate the gradient of a structure you can call

    ./f3 -g -x <structure.xyz>


## Examples

Included is a parameter file and structure.xyz

    ./f3 -p examples/parameter.dat -x examples/structure.xyz


## Implementation

If you want to implement the module in a QM package, you can. It should be simple.
Look in the `program.f90` to see how to properly use the module.



