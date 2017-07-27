#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
#

.PHONY: all

all: target/release/tada

target/release/tada: src/*.rs src/*/*.rs
	RUSTFLAGS="${RUSTFLAGS} -A unused_must_use -A unused_variables -A unused_mut" cargo build --bin=tada --release


