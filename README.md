# Third-Generation Hydrogen-Bonding Correction


Third-Generation Hydrogen-Bonding Corrections for Semiempirical QM Methods and Force Fields.


Fortran F90 (module) implementation of Martin Korths Hydrogen-bond correction term.

## Please cite as:

Martin Korth, J. Chem. Theory Comput., 2010, 6 (12), pp 3808–3816, DOI: 10.1021/ct100408b

## Installation

Use the make script. Change compiler if you use something else than `gfortran`.

    make

## Usage

### Energy

You will properly want to use the module as part of a QM/SQM package, but
if you just want to try it out, use it as.

    ./f3 -x <structure.xyz>

to calculate the correction energy of the XYZ structure.

### Gradient

If you want to calculate the gradient of a structure you can call

    ./f3 -x <structure.xyz> --gradient


## Implementation

If you want to implement the module in a QM package, you can. It should be simple.
Look in the `program.f90` to see how to properly use the module.

