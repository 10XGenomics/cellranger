#!/usr/bin/env python3
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.

"""Parse RunInfo.xml to get flowcell."""

from __future__ import annotations

import sys
import xml.etree.ElementTree


def main(loc):
    tree = xml.etree.ElementTree.parse(loc)

    # Get flowcell ID and remove cruft from MiSeq IDs.
    flowcell = tree.getroot().findall("./Run/Flowcell")[0].text
    flowcell = flowcell.replace("000000000-", "")
    print(flowcell)


if __name__ == "__main__":
    main(sys.argv[1])
