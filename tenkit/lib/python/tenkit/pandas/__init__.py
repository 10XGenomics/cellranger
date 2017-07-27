'''
We ran into a problem with a bug in glibc and haswell processors that causes
a lock in some situations. We first found it in bwa and fixed it by linking to
jmalloc. But now it is turning up in pandas in numexpr eval. The following code
makes pandas not use eval when subsetting data frames.
'''

from pandas import *
computation.expressions.set_use_numexpr(False)
