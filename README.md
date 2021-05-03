Jurassic-GPU-retrieval

    based on JURASSIC by Lars Hoffmann
    this version is dedicated to host the GPU module for the EGA method within the forward model
    
Prerequisites

    CUDA
    C-compiler (e.g. GCC)
    make

How to build

    cd lib
    ./build.sh
    cd ..
    cd src
    make
    cd ..
    cd example/limb
    ./run.sh
    diff rad.tab rad.org
 
