#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# Build cellranger.
#
SHELL := /bin/bash -O extglob

VERSION=$(shell git describe --tags --always --dirty)

### Rust
RUST_SRC_PATH=$(shell pwd)/lib/rust
export CARGO_HOME=$(RUST_SRC_PATH)/.cargo
RUST_BINS=vdj_asm chunk_reads annotate_reads detect_chemistry cr_stage cr_stage_pd

.PHONY: all  clean $(RUST_BINS)   louvain rust-clean  version-files

#
# Targets for development builds.
#
all: $(RUST_BINS) louvain   cython
	make -C tenkit all

clean: rust-clean  louvain-clean cython-clean
	make -C tenkit clean

#
# Targets for cython builds
#

# Find python libs, includes, and default flags.
PYTHON_ROOT:=$(shell python-config --exec-prefix)
PYTHON_LIB_DIR:=$(dir $(lastword $(wildcard $(PYTHON_ROOT)/lib*/libpython*.so*)))
PYTHON_SITE_PKG:=$(lastword $(wildcard $(PYTHON_LIB_DIR)/python*/site-packages))
PYTHON_CFLAGS:=$(shell python-config --cflags)
PYTHON_LDFLAGS:=$(shell python-config --ldflags)

# Find stuff that needs to be cythonized.
CYTHON_SRCS=$(shell find $(PWD)/mro/stages $(PWD)/lib/python -type f -name '*.pyx')
CYTHON_LIBS=$(patsubst %.pyx, %.so, $(CYTHON_SRCS))
CYTHON_BUILDPATH=$(shell pwd)/lib/cython
CYTHON_FLAGS?=--line-directives $(EXTRA_CYTHON_FLAGS)

# Prevent make from automatically deleting intermediate files.
.PRECIOUS: $(CYTHON_BUILDPATH)/%.c $(CYTHON_BUILDPATH)/%.o

.PHONY: cython

$(CYTHON_BUILDPATH)/%.c: $(PWD)/%.pyx
	mkdir -p $(@D) && cython $(CYTHON_FLAGS) -w $(<D) -o $(abspath $@) $(<F)

$(CYTHON_BUILDPATH)/%.o: $(CYTHON_BUILDPATH)/%.c
	$(CC) $(PYTHON_CFLAGS) $(CFLAGS) -g -O3 -c -fPIC -fopenmp \
	    -I$(PYTHON_SITE_PKG)/numpy/core/include \
	    -o $@ \
	    $<

$(PWD)/%.so: $(CYTHON_BUILDPATH)/%.o
	$(CC) -L$(PYTHON_LIB_DIR) $(PYTHON_LDFLAGS) $(LDFLAGS) -shared -fopenmp -fPIC $< -o $@

cython: $(CYTHON_LIBS)

clean-cython:
	rm -rf $(CYTHON_BUILDPATH)
	rm -f $(CYTHON_LIBS)


rust-clean:
	rm -Rf lib/rust/.cargo
	$(foreach dir, $(RUST_BINS), \
	    pushd $(RUST_SRC_PATH)/$(dir) >/dev/null && \
	    cargo clean; \
	    popd >/dev/null; \
	)

louvain:
	make -C lib/louvain
	cp lib/louvain/convert lib/bin
	cp lib/louvain/louvain lib/bin

louvain-clean:
	make -C lib/louvain clean
	rm -f lib/bin/convert
	rm -f lib/bin/louvain

lib/bin:
	mkdir -p lib/bin

$(RUST_BINS): lib/bin
	set -e ; \
	pushd lib/rust >/dev/null; \
	cargo build --release; \
	cp target/release/$@ ../bin/; \
	popd > /dev/null

