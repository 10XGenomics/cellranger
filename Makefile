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
all: $(RUST_BINS) louvain  
	make -C tenkit all

clean: rust-clean  louvain-clean
	make -C tenkit clean

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
	pushd lib/rust/$@ >/dev/null; \
	cargo build --release --bin $@; \
	cp target/release/$@ ../../bin/; \
	popd > /dev/null

