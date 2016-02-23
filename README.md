<p align="center"><img src="https://github.com/substance33/Spectra/blob/LNS/spectra_logo.png" width="450" alt="Spectra Logo"></p>

## About

**Spectra** is a [MATLAB](http://www.mathworks.com/products/matlab/) library which can be used as a generic experimental data container and is especially suited to neutron scattering datasets. Spectra are stored as spec1d objects which can have mathematical operations applied, where errors are automatically and accuratly calculated.

## Components
This library contains the [SpinW](https://www.github.com/tsdev/spinw) library for calculating spin wave spectra and [ResLibCal](https://www.github.com/McStasMcXtrace/iFit) from the iFit library to calculate instrumental resolution.

Data manipulation such as combining, cutting, re-binning and has been implemented as well as fitting data to arbitary models.

## Disclaimer

Research science contains a few endpoints Software, Publication, Pedagogy, and Reproducibility.  The research science process is very interactive so I personally reserve the ``LNS`` branch (default) for vetted software codes.  I tend to work off of the ``development`` branch while I am implementing new code.

## Getting Started

This software can be obtained by git or as a [zip file](https://github.com/simonward86/Spectra/archive/LNS.zip)

1. ``git clone https://github.com/simonward86/Spectra.git``
2. If you want the additional submodules run ``git clone --recursive https://github.com/simonward86/Spectra.git`` instead.
3. Rename ``startup_example.m`` to ``startup.m`` and place it in the MATLAB home directory or ```fullfile(matlabroot,'toolbox','local')```.
4. If you do the former, you will neet to execute ``startup.m`` when you start MATLAB.

## Information
User guides and examples can be found at the project site https://simonward86.github.io/Spectra look for tagged posts.
