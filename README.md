# Building Cell Ranger 3.0.2
## Build dependencies
- Python 2.7.13
- rust 1.28.0
- clang 6.0
- go 1.11
- node v8.11.4

### Example setup of build dependencies on Ubuntu 16.04.3
```
sudo apt-get install make clang-6.0 golang-1.11-go libz-dev libbz2-dev liblzma-dev

# Add golang to path
export PATH=/usr/lib/go-1.11/bin:$PATH

# Install rustup from https://www.rustup.rs/ . Then:
rustup install 1.28.0
rustup default 1.28.0
```

## Build command
`make`

# Running Cell Ranger
## Runtime dependencies
- Binary dependencies can be found in the Ranger 3.0.2 package (https://support.10xgenomics.com/developers/software/downloads/latest)
  - The Ranger package includes a build of the Martian platform (v3.2.0), which is open source. For more information, go to http://martian-lang.org/ .

## Setting up the environment
```
# Setup Martian and binary dependencies
source /path/to/ranger/sourceme.bash

# Setup Cell Ranger
source /path/to/cellranger/sourceme.bash
```

## Note about Loupe
The binaries required to generate Loupe Cell Browser (.cloupe) and Loupe V(D)J Browser files (.vloupe) are not included in this repository or in the binary dependencies package Ranger. The necessary binaries can be obtained from an existing binary version of Cell Ranger by running:
`cp /path/to/cellranger-3.0.2/cellranger-cs/*/lib/bin/{crconverter,vlconverter} /path/to/open-source-cellranger/lib/bin/`

# Support
We do not provide support for building and running this code.

The officially supported release binaries are available at: (https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest)
