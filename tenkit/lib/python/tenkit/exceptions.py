#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Exception classes
#

class NotSupportedException(Exception):
    """ For when a specified method isn't supported yet- e.g. when a distribution is
    specified that hasn't been built yet for the given code.
    """
    pass
