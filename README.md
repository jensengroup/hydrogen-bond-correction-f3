# Third-Generation Hydrogen-Bonding Correction

Third-Generation Hydrogen-Bonding Corrections for Semiempirical QM Methods and
Force Fields.

## Introduction

Fortran F90 (module) implementation of Martin Korths Hydrogen-bond correction
term.

## Known Limitation

TODO Nitrogen term

## Please cite as:

    Martin Korth,
    J. Chem. Theory Comput., 2010, 6 (12), pp 3808–3816,
    DOI: 10.1021/ct100408b

## Installation

Use the make script. Change compiler if you use something else than `gfortran`.

    make

## Parameters

Default parameters are for the PM6-D3H+ model.

    nitrogen = -0.11
    oxygen   = -0.12

change it either via source code, or use a parameter file.


## Usage

### Energy

You will properly want to use the module as part of a QM/SQM package, but
if you just want to try it out, use it as.

    ./f3_exe [OPTIONS] <structure.xyz>

options:

    -v,        --version        print version information and exit
    -h,        --help           print usage information and exit
    -p <file>, --param <file>   use parameter file <par.dat>
    -d,        --debug          print debug information
    -g,        --gradient       calculate and print gradient


to calculate the correction energy of the XYZ structure.

### Gradient

If you want to calculate the gradient of a structure you can call

    ./f3_exe -g -x <structure.xyz>


## Using a parameter file

Included is a parameter file and structure.xyz

    ./f3_exe -p examples/parameter.dat -x examples/structure.xyz


## Implementation

If you want to implement the module in a QM package, you can. It should be
simple.  Look in the `src/f3.f90` and `src/f3_program.f90` to see how to
properly use the module.



