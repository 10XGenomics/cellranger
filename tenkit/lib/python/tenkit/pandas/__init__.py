'''
We ran into a problem with a bug in glibc and haswell processors that causes
a lock in some situations. We first found it in bwa and fixed it by linking to
jmalloc. But now it is turning up in pandas in numexpr eval. The following code
makes pandas not use eval when subsetting data frames.
'''

try:
    from pandas.core import computation
except ImportError:  # pandas moved the computation library in 0.20.1
    pass
from pandas import *
computation.expressions.set_use_numexpr(False)
