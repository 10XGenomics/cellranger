#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#

LIBPY=lib/python
LIBSSW=$(LIBPY)/striped_smith_waterman
#LIBSSW_INIT_PY=$(LIBSSW)/__init__.py
LIBSSW_BUILD_DIR=$(LIBSSW)
LIBSSW_BUILD=$(LIBSSW_BUILD_DIR)/libssw.so
VERSION=$(shell git describe --tags --always --dirty)
GOBINS=barcodeqc godemux q30count

PWD=$(shell pwd)
export GOPATH=$(PWD)/lib/go

.PHONY: test go_generate

all: $(LIBSSW_BUILD) $(GOBINS)

clean: libssw-clean
	rm -rf $(GOPATH)/bin
	rm -rf $(GOPATH)/pkg
	rm -rf lib/bin

test:
	MROPATH=$(PWD)/mro \
		PATH=$(PATH):$(PWD)/bin \
		PYTHONPATH=$(PYTHONPATH):$(PWD)/$(LIBPY) \
		nosetests --with-xunit --with-coverage --cover-erase --cover-html --cover-xml --cover-package=stages,tenkit mro/stages/* $(LIBPY)

$(LIBSSW_BUILD):
	$(MAKE) -C $(LIBSSW_BUILD_DIR) libssw.so

libssw-clean:
	$(MAKE) -C $(LIBSSW_BUILD_DIR) clean

lib/bin:
	mkdir -p lib/bin

go_generate:
	PYTHONPATH=$(PYTHONPATH):$(PWD)/$(LIBPY) \
	go generate tenkit/...

$(GOBINS): lib/bin go_generate
	go install -ldflags "-X $@.__VERSION__='$(VERSION)'" $@
	cp lib/go/bin/$@ lib/bin



