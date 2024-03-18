OpenBLAS
========

This is a fork of the [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS) and has been improved performance for the FUJITSU A64FX processor.

The following routine is tuned for A64FX.
* SGEMM

# Prerequisites

* GNU Compiler Collection 10.2 or higher

# Installation

## Download the source code
* Clone this repository.

```
$ git clone https://github.com/fujitsu/OpenBLAS.git
```

## Build and install OpenBLAS for A64FX 

### Notice
* To make a library for A64FX, specify `TARGET=A64FX` at make command,
* In the following descriptions, `$INSTALL_PATH` is the name of the directory in which to install the OpenBLAS libraries and header files.

### Build and install the sequential library
```
$ make CC=gcc TARGET=A64FX
$ make CC=gcc TARGET=A64FX install PREFIX=$INSTALL_PATH
```

### Build and install threaded parallel Library parallelized with OpenMP
```
$ make USE_OPENMP=1 CC=gcc TARGET=A64FX
$ make USE_OPENMP=1 CC=gcc TARGET=A64FX install PREFIX=$INSTALL_PATH
```

### Build and install ILP64 sequential library
```
$ make CC=gcc TARGET=A64FX INTERFACE64=1
$ make CC=gcc TARGET=A64FX INTERFACE64=1 install PREFIX=$INSTALL_PATH
```

### Build and install ILP64 threaded parallel Library parallelized with OpenMP
```
$ make USE_OPENMP=1 CC=gcc TARGET=A64FX INTERFACE64=1
$ make USE_OPENMP=1 CC=gcc TARGET=A64FX INTERFACE64=1 install PREFIX=$INSTALL_PATH
```

# Usage

To use the OpenBLAS library, add `-lopenblas` as an option when linking user programs.
Also, use the `-L` and `-I` options to specify the directory where you installed the libraries and header files.

## How to compile

* For Fortran programs

```
$ gfortran a.f -L$INSTALL_PATH/lib -lopenblas
```

* For C programs

```
$ gcc a.c -I$INSTALL_PATH/include -L$INSTALL_PATH/lib  -lopenblas
```

It is recommended to use large pages for performance.
For Technical Computing Suite environments, HPC extension large page library can be used
by adding `-L/opt/FJSVxos/mmm/lib64 -lmpg -Wl,-T/opt/FJSVxos/mmm/util/bss-2mb.lds ` to the options. Specify this option before any other libraries.
```
gcc  a.c -L/opt/FJSVxos/mmm/lib64 -lmpg -Wl,-T/opt/FJSVxos/mmm/util/bss-2mb.lds -I$INSTALL_PATH/include -L$INSTALL_PATH/lib  -lopenblas
```

# Performance

The OpenBLAS library in this product improves the performance of SGEMM as follows:

| Library    | Routine | Parameters                | # of cores | Original OpenBLAS | OpenBLAS tuned for A64FX |
|------------|---------|---------------------------|------------|-------------------|--------------------------|
| Sequential | SGEMM   | No Transpose, M=N=K=5000  | 1 core     | 78 GFlops         | 108 GFlops               | 
| OpenMP     | SGEMM   | No Transpose, M=N=K=10000 | 12 cores   | 827 GFlops        | 1267 GFlops              |

# Restrictions

None.

# License
* See [LICENSE](https://github.com/fujitsu/OpenBLAS/LICENSE) file.

