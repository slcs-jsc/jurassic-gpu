# JURASSIC-GPU

The Juelich Rapid Spectral Simulation Code for GPUs (JURASSIC-GPU) is a fast infrared radiative transfer model for the analysis of atmospheric remote sensing measurements.

[![release (latest by date)](https://img.shields.io/github/v/release/slcs-jsc/jurassic-gpu)](https://github.com/slcs-jsc/jurassic-gpu/releases)
[![commits since latest release (by SemVer)](https://img.shields.io/github/commits-since/slcs-jsc/jurassic-gpu/latest)](https://github.com/slcs-jsc/jurassic-gpu/commits/master)
[![last commit](https://img.shields.io/github/last-commit/slcs-jsc/jurassic-gpu.svg)](https://github.com/slcs-jsc/jurassic-gpu/commits/master)
[![code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/jurassic-gpu.svg)](https://github.com/slcs-jsc/jurassic-gpu/tree/master/src)
[![top language](https://img.shields.io/github/languages/top/slcs-jsc/jurassic-gpu.svg)](https://github.com/slcs-jsc/jurassic-gpu/tree/master/src)
[![license](https://img.shields.io/github/license/slcs-jsc/jurassic-gpu.svg)](https://github.com/slcs-jsc/jurassic-gpu/blob/master/COPYING)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.4744174.svg)](https://doi.org/10.5281/zenodo.4744174)

## Features

* The Juelich Rapid Spectral Simulation Code (JURASSIC) is a fast infrared radiative transfer model for the analysis of atmospheric remote sensing measurements.
* JURASSIC uses the emissivity growth approximation (EGA) to conduct the radiative transfer calculations.
* Band transmittances are obtained from pre-calculated look-up tables from line-by-line calculations.
* The model was carefully tested in intercomparisons with the Karlsruhe Optimized and Precise Radiative Transfer Algorithm (KOPRA), the Reference Forward Model (RFM), and the Stand-alone AIRS Radiative Transfer Algorithm (SARTA).
* JURASSIC features an MPI/OpenMP hybrid parallelization for efficient use on supercomputers.
* This version of JURASSIC, referred to as JURASSIC-GPU, is dedicated to host the GPU module for the EGA method within the forward model.
    
## Prerequisites

The following software and tools are required to install JURASSIC-GPU:

    CUDA
    C-compiler (e.g. GCC)
    make

## How to build

Follow these steps to build JURASSIC-GPU:

    cd lib
    ./build.sh
    cd ..
    cd src
    make

A number of test cases are provided along with JURASSIC-GPU. For example, this is showing how to run a test case for the limb geometry:

    cd example/limb
    ./run.sh
    diff rad.tab rad.org

# Further information

The JURASSIC-GPU code is described in the following references:

* Baumeister, P. F., and Hoffmann, L., Fast Infrared Radiative Transfer Calculations Using Graphics Processing Units: JURASSIC-GPU v2.0, Geosci. Model Dev., 2021, submitted.

* Baumeister, P. F., Rombach, B., Hater, T., Griessbach, S., Hoffmann, L,, Bühler, M., and Pleiter, D., Strategies for Forward Modelling of Infrared Radiative Transfer on GPUs, Parallel Computing 2017, Bologna, Italy, 2017

* You can cite the source code of JURASSIC-GPU by using the DOI https://doi.org/10.5281/zenodo.4744174. This DOI represents all versions, and will always resolve to the latest one. Specific DOIs for each release of JURASSIC can be found on the zenodo web site.
    
These are the main references for citing the JURASSIC model in scientific publications:

* Hoffmann, L., and M. J. Alexander, Retrieval of stratospheric temperatures from Atmospheric Infrared Sounder radiance measurements for gravity wave studies, J. Geophys. Res., 114, D07105, https://doi.org/10.1029/2008JD011241, 2009.

* Hoffmann, L., Kaufmann, M., Spang, R., Müller, R., Remedios, J. J., Moore, D. P., Volk, C. M., von Clarmann, T., and Riese, M.: Envisat MIPAS measurements of CFC-11: retrieval, validation, and climatology, Atmos. Chem. Phys., 8, 3671-3688, https://doi.org/10.5194/acp-8-3671-2008, 2008.

* You can cite the source code of JURASSIC by using the DOI https://doi.org/10.5281/zenodo.4572889. This DOI represents all versions, and will always resolve to the latest one. Specific DOIs for each release of JURASSIC can be found on the zenodo web site.

## Contributing

We are interested in sharing JURASSIC-GPU for operational or research applications.

Please do not hesitate to contact us, if you have any further questions or need support.

## License

JURASSIC-GPU is distributed under the [GNU General Public License v3.0](https://github.com/slcs-jsc/jurassic-gpu/blob/jurassic-gpu/COPYING).

## Contact

Paul Baumeister and Lars Hoffmann

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: p.baumeister@fz-juelich.de, l.hoffmann@fz-juelich.de
