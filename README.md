# Building Cell Ranger 2.1.0
## Build dependencies
- Python 2.7.13
- rust 1.19.0
- clang 4.0.1
- go 1.9

### Example setup of build dependencies on Ubuntu 16.04.3
```
sudo apt-get install make clang-4.0 golang-1.9 libz-dev

# Add golang to path
export PATH=/usr/lib/go-1.9/bin:$PATH

# Install rustup from https://www.rustup.rs/ . Then:
rustup install 1.19.0
rustup default 1.19.0
```

## Build command
`make`

# Running Cell Ranger
## Runtime dependencies
- Binary dependencies can be found in the Ranger 2.1.0 package (https://support.10xgenomics.com/developers/software/downloads/latest)
  - The Ranger package includes a build of the Martian platform (v2.3.0), which is open source. For more information, go to http://martian-lang.org/ .

## Setting up the environment
```
# Setup Martian and binary dependencies
source /path/to/ranger/sourceme.bash

# Setup Cell Ranger
source /path/to/cellranger/sourceme.bash
```

## Note about Loupe
The binaries required to generate Loupe Cell Browser (.cloupe) and Loupe V(D)J Browser files (.vloupe) are not included in this repository or in the binary dependencies package Ranger. By default, you will get empty .cloupe/.vloupe files when running a version of Cell Ranger built from this repository. The necessary binaries can be obtained from an existing binary version of Cell Ranger by running:
`cp /path/to/cellranger-2.1.0/cellranger-cs/*/lib/bin/{crconverter,vlconverter} /path/to/open-source-cellranger/lib/bin/`

# Support
We do not provide support for building and running this code.

The officially supported release binaries are available at: (https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
