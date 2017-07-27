#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import numpy
import numpy as np
import json
import math

# Yuck yuck yuck yuck yuck!
# The default JSON encoders do bad things if you try to encode
# NaN/Infinity/-Infinity as JSON.  This code takes a nested data structure
# composed of: atoms, dicts, or array-like iteratables and finds all of the
# not-really-a-number floats and converts them to an appropriate string for
# subsequent jsonification.
def json_sanitize(data):
    # This really doesn't make me happy. How many cases we we have to test?
    if (type(data) == float) or (type(data) == numpy.float64):
        # Handle floats specially
        if math.isnan(data):
            return "NaN";
        if (data ==  float("+Inf")):
            return "inf"
        if (data == float("-Inf")):
            return "-inf"
        return data
    elif hasattr(data, 'iterkeys'):
        # Dictionary case
        new_data = {}
        for k in data.keys():
            new_data[k] = json_sanitize(data[k])
        return new_data
    elif hasattr(data, '__iter__'):
        # Anything else that looks like a list. N
        new_data = []
        for d in data:
            new_data.append(json_sanitize(d))
        return new_data
    elif hasattr(data, 'shape') and data.shape == ():
        # Numpy 0-d array
        return np.asscalar(data)
    else:
        return data

def safe_jsonify(data, pretty=False):
    safe_data = json_sanitize(data)
    if pretty:
        return json.dumps(safe_data, indent=4, sort_keys=True)
    else:
        return json.dumps(safe_data)

class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray) and obj.ndim == 1:
                return obj.tolist()
        elif isinstance(obj, np.generic):
            return obj.item()
        return json.JSONEncoder.default(self, obj)

def dump_numpy(data, fp, pretty=False):
    ''' Dump object to json, converting numpy objects to reasonable JSON '''
    if pretty:
        json.dump(data, fp, cls=NumpyAwareJSONEncoder, indent=4, sort_keys=True)
    else:
        json.dump(data, fp, cls=NumpyAwareJSONEncoder)
