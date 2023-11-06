#!/usr/bin/env python3
#
# Copyright (c) 2017 10x Genomics, Inc. All rights reserved.

"""Parse RunInfo.xml to get flowcell."""

from __future__ import annotations

import sys
import xml.etree.ElementTree


def main(loc):
    tree = xml.etree.ElementTree.parse(loc)

    maxilen = 0
    reads = tree.getroot().findall("./Run/Reads/Read")
    for read in reads:
        cycles = int(read.attrib["NumCycles"])
        isindex = read.attrib["IsIndexedRead"] == "Y"
        if isindex and (cycles > maxilen):
            maxilen = cycles
    print(maxilen)


if __name__ == "__main__":
    main(sys.argv[1])
