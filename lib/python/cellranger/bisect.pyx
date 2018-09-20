
from __future__ import division

import numpy as np
cimport numpy as np

cpdef np.int32_t bisect_left(const np.int32_t[:] a, x, np.ndarray keys):
    cdef int lo = 0, hi = a.shape[0], i
    while hi - lo > 1:
        i = (lo + hi) // 2  # cannot be <1
        y = keys[a[i]]
        if keys[a[i-1]] < x and x <= y:
            return a[i]
        elif y < x:
            lo = i
        elif hi == i + 1:
            hi = i
        else:
            hi = i + 1
    return a[lo]
