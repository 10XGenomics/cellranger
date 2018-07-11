#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#

import copy
import math

import tenkit.stats as tk_stats
import cPickle as pickle
import numpy as np

from abc import ABCMeta, abstractmethod

import tenkit.dict_utils as tk_dict

class SummaryStatistic:
    __metaclass__ = ABCMeta

    @abstractmethod
    def add(self, x):
        pass

    @abstractmethod
    def combine(self, x):
        pass

''' Port of:

 [MP80]  Munro & Paterson, "Selection and Sorting with Limited Storage",
         Theoretical Computer Science, Vol 12, p 315-323, 1980.
 As described in :
 [MRL98] Manku, Rajagopalan & Lindsay, "Approximate Medians and other
         Quantiles in One Pass and with Limited Memory", Proc. 1998 ACM
         SIGMOD, Vol 27, No 2, p 426-435, June 1998.'''
class QuantileEstimator(SummaryStatistic):
    def __init__(self, epsilon, n, b = None, k = None):
        self.epsilon = epsilon
        self.n = n

        if b is None and k is None:
            if epsilon is None or n is None:
                raise Exception(
                        "Must set 'epsilon' and 'n' if 'b' and 'k' are not set"
                        )
            (b, k) = QuantileEstimator.optimalBK(epsilon, n)

        self.b = b
        self.k = k
        self.buffer = [[] for i in xrange(b)]
        self.count = 0

        self.min = None
        self.max = None

        self.indices = []

        self.leaves_sorted = False

    @classmethod
    def _create_b_k(cls, b, k, epsilon = None, n = None):
        return QuantileEstimator(epsilon, n, b, k)

    # get the buffer level that contains the largest index
    def _largest_level(self):
        cur_largest = None
        largest_idx = None
        for i in xrange(len(self.buffer)):
            if self.indices[i] is None or self.indices[i] < 0:
                continue
            val = self.buffer[i][self.indices[i]]
            if val > cur_largest:
                cur_largest = val
                largest_idx = i

        return (largest_idx, cur_largest)


    def _init_indices(self):
        buf0len = len(self.buffer[0])
        buf1len = len(self.buffer[1])

        if not self.leaves_sorted:
            if buf0len > 0:
                self.buffer[0].sort()
            if buf1len > 0:
                self.buffer[1].sort()
            self.leaves_sorted = True

        self.indices = [len(buf) - 1 if len(buf) > 0 else None for buf in self.buffer]

        return


    # Returns a tuple of optimal (b, k) given some epsilonilon and n
    @staticmethod
    def optimalBK(epsilon, n):
        b = 2
        upr = epsilon * n

        # This technically overallocates one buffer, but we use the extra one
        # to make the call in collapse_recursive cleaner
        while (b - 2) * (1L << (b - 2)) + 0.5 <= upr:
            b += 1
        # b -= 1

        # subtract an extra '1' from b to make up for extra b added above
        k = int( math.ceil( float(n) / (1L << (b - 1 - 1)) ) )

        return (b, k)

    def add(self, x):
        self.leaves_sorted = False

        if x < self.min or self.min is None:
            self.min = x

        if x > self.max or self.max is None:
            self.max = x

        if self.count > 0 and self.count % (2 * self.k) == 0:
            self.buffer[0].sort()
            self.buffer[1].sort()
            self.leaves_sorted = True
            level = 1
            self.collapse_recursive(self.buffer[0], level)

        index = 0 if len(self.buffer[0]) < self.k else 1
        self.buffer[index].append(x)
        self.count += 1

    # Slow implementation if a and b are already sorted... will speed up later
    # if an issue - HP
    def collapse(self, a, b):
        tmp = a + b
        tmp.sort()

        del a[:]
        del b[:]

        # what happens if len(a) != len(b)? Shouldn't technically be called...
        # look into - HP

        # start at index 1 to keep consistent w/ paper
        return tmp[1::2]

    def collapse_recursive(self, buf, level):
        # print 'level:\t', level
        merged = self.collapse(self.buffer[level], buf)

        if level + 1 >= len(self.buffer):
            self.buffer.append([])
            self.b += 1

        using_tmp_merge = True
        if len(self.buffer[level + 1]) == 0:
            self.buffer[level + 1] = copy.deepcopy(merged)
            using_tmp_merge = False

        if not using_tmp_merge:
            return

        self.collapse_recursive(merged, level + 1)

    @staticmethod
    def weight(level):
        if level == 0 or level == 1:
            return 1
        return 1L << (level - 1)

    def combine(self, other):
        res = MergedQuantileEstimator([self, other])
        return res

    # should only be called when all buffers are full (except 0 and 1)
    def quantile(self, q):
        if self.count == 0:
            return float('NaN')

        if q < 0.0 or q > 1.0:
            raise Exception('Quantile must be between [0,1]')

        # comparison should be OK as long as user inputs exactly 0.0 or 1.0 (min/max)
        if q == 0.0:
            return self.min
        if q == 1.0:
            return self.max

        self._init_indices()

        target = math.ceil( (1.0 - q) * self.count )
        tot_sum = 0.0

        # Algorithm: start at max and iteratively add the weight of each
        # large element
        while True:
            (i, val)  = self._largest_level()
            self.indices[i] -= 1
            tot_sum += self.weight(i)
            if tot_sum >= target:
                break

        return self.buffer[i][self.indices[i] + 1]


class MergedQuantileEstimator(SummaryStatistic):
    def __init__(self, qes):
        self.qes = qes
        self.count = None

    def add(self, x):
        raise Exception("Can't add to a MergedQuantileEstimator!")

    def _init_indices(self):
        for qe in self.qes:
            qe._init_indices()
        self.count = sum([qe.count for qe in self.qes])
        return

    def min(self):
        return min([qe.min for qe in self.qes])

    def max(self):
        return max([qe.max for qe in self.qes])

    def quantile(self, q):
        if q < 0.0 or q > 1.0:
            raise Exception('Quantile must be between [0,1]')

        # comparison should be OK as long as user inputs exactly 0.0 or 1.0 (min/max)
        if q == 0.0:
            return self.min()
        if q == 1.0:
            return self.max()

        self._init_indices()

        if self.count == 0:
            return float('NaN')

        # print 'count: ', self.count
        target = int(math.ceil( (1.0 - q) * self.count ))
        i_qe = None
        i_max = None
        max_val = None

        tot_sum = 0

        while True:
            max_val = None
            i_max = None
            # get the max amongst all of the quantile estimators
            for i in xrange(len(self.qes)):
                cur_qe = self.qes[i]
                (cur_i_max, cur_max) = cur_qe._largest_level()
                if cur_max is not None and cur_max > max_val:
                    max_val = cur_max
                    i_max = cur_i_max
                    i_qe = i
            self.qes[i_qe].indices[i_max] -= 1
            # print 'adding weight\t', self.qes[i_qe].weight(i_max)
            tot_sum += self.qes[i_qe].weight(i_max)

            if tot_sum >= target:
                break
        qe = self.qes[i_qe]
        i = qe.indices[i_max]

        return qe.buffer[i_max][i + 1]

    def combine(self, qe):
        res = None

        if isinstance(qe, QuantileEstimator):
            res = MergedQuantileEstimator( self.qes + [qe] )
        elif isinstance( qe, MergedQuantileEstimator):
            res = MergedQuantileEstimator( self.qes + qe.qes )
        else:
            raise Exception('Unknown type to combine with MergedQuantileEstimator')

        return res


class DiscreteDistribution(SummaryStatistic):
    def __init__(self, lower, upper):
        self.count = 0
        self.lower = lower
        self.upper = upper
        self.lower_obs = None
        self.upper_obs = None
        self.hist = np.zeros(upper - lower + 1, dtype=int)
        self._prob = None
        self._ecdf = None

    def _x_to_idx(self, x):
        self.check_range(x)
        return x - self.lower

    def _idx_to_x(self, idx):
        if idx < 0 or idx >= len(self.hist):
            raise Exception('Index ' + str(idx) + ' is out of bounds ')
        return idx + self.lower

    def check_range(self, x):
        if x < self.lower or x > self.upper:
            raise Exception('Attempting to add item "' + str(x) +
                    '", which is out of range [' + str(self.lower) + ', ' +
                    str(self.upper) + ']')
        return

    def add(self, x):
        # invalidate these summary variables
        self._prob = None
        self._ecdf = None

        self.check_range(x)
        self.hist[self._x_to_idx(x)] += 1
        self.count += 1

        if x < self.lower_obs or self.lower_obs is None:
            self.lower_obs = x
        if x > self.upper_obs or self.upper_obs is None:
            self.upper_obs = x

        return

    def min(self):
        return self.lower_obs

    def max(self):
        return self.upper_obs

    def init_prob(self):
        if self.count == 0:
            raise Exception('Need to add items before you can call init_prob()')
        self._prob = np.array(self.hist,
                dtype = np.float) / np.float(self.count)
        return

    def init_ecdf(self):
        if self.count == 0:
            raise Exception('Need to add items before you can call init_ecdf()')
        self.init_prob()
        self._ecdf = np.cumsum( self._prob )
        return

    def prob(self, x):
        self.check_range(x)
        if self._prob is None:
            self.init_ecdf()
        return self._prob[self._x_to_idx(x)]

    def ecdf(self, x):
        self.check_range(x)
        if self._ecdf is None:
            self.init_ecdf()
        return self._ecdf[self._x_to_idx(x)]

    def quantile(self, q):
        if self.count == 0:
            return float('NaN')

        if q < 0.0 or q > 1.0:
            raise Exception('Quantile must be between [0,1]')

        if q == 0.0:
            return self.min()
        if q == 1.0:
            return self.max()

        if self._ecdf is None:
            self.init_ecdf()

        lwr = 0
        upr = len(self.hist) - 1
        while True:
            mid = (upr - lwr) / 2 + lwr
            if (q <= self._ecdf[mid] and q > self._ecdf[mid - 1]) or lwr == upr:
                return self._idx_to_x(mid)
            if q > self._ecdf[mid]:
                lwr = min(mid + 1, len(self.hist) - 1)
            else:
                upr = max(mid - 1, 0)

    def combine(self, b):
        lwr = min( self.lower, b.lower )
        upr = max( self.upper, b.upper )
        ret = DiscreteDistribution(lwr, upr)
        ret.count = self.count + b.count
        ret.lower_obs = min(self.min(), b.min())
        ret.upper_obs = max(self.max(), b.max())

        for i in xrange(len(self.hist)):
            ret.hist[ ret._x_to_idx(self._idx_to_x(i)) ] += self.hist[i]

        for i in xrange(len(b.hist)):
            ret.hist[ ret._x_to_idx(b._idx_to_x(i)) ] += b.hist[i]

        return ret

class DictionaryDistribution(SummaryStatistic):
    def __init__(self):
        self.dict = {}

    # x: the key
    # val: the value
    def add(self, x, val = 1):
        if val is not None:
            self.dict[x] = self.dict.get(x, 0) + val
        else:
            self.dict[x] = None

        return self.dict[x]

    def combine(self, dd):
        res = DictionaryDistribution()
        res.dict = tk_dict.add_dicts(self.dict, dd.dict, 1)

        return res

    def get(self, k):
        return self.dict.get(k, 0)

    def pop(self, k):
        if k in self.dict:
            return self.dict.pop(k)
        return None

    def __getitem__(self, k):
        return self.get(k)

    def __setitem__(self, k, v):
        self.pop(k)
        self.dict[k] = v

        return v

''' Port of:
    http://www.johndcook.com/skewness_kurtosis.html
    which is an implementation of the Knuth online variance calculation
    '''
class BasicSummary(SummaryStatistic):
    def __init__(self):
        self.clear()

    def clear(self):
        self.count = 0L
        self.M1 = 0.0
        self.M2 = 0.0
        self.M3 = 0.0
        self.M4 = 0.0
        self.min = None
        self.max = None
        return

    def add(self, x):
        x = float(x)
        n1 = self.count

        self.count += 1

        if x < self.min or self.min is None:
            self.min = x

        if x > self.max or self.max is None:
            self.max = x

        delta = x - self.M1
        delta_n = delta / self.count
        delta_n2 = delta_n * delta_n
        term = delta * delta_n * n1
        self.M1 += delta_n
        self.M4 += term * delta_n2 * \
                ( self.count * self.count - 3*self.count + 3 ) + \
                6 * delta_n2 * self.M2 - 4 * delta_n * self.M3
        self.M3 += term * delta_n * (self.count - 2) - 3 * delta_n * self.M2
        self.M2 += term

        return

    def total(self):
        return self.count * self.mean()

    def min(self):
        return self.min

    def max(self):
        return self.max

    def mean(self):
        return self.M1

    def var(self):
        return self.M2 / (self.count - 1.0)

    # is this unbiased for SD or should multiply by (n-1)/n?
    def sd(self):
        return math.sqrt( self.var() )

    def skewness(self):
        return math.sqrt(self.count) * self.M3 / math.pow( self.M2, 1.5 )

    # NB: the -3 has to do with making the Gaussian have kurtosis 0. Most
    # people compute it this way. R package 'moments' doesn't -3
    def kurtosis(self):
        return float(self.count) * self.M4 / (self.M2 * self.M2) - 3.0

    def combine(self, b):
        res = BasicSummary()

        res.count = self.count + b.count

        # cut out if we have no counts to avoid divide by zero
        if res.count == 0:
            res.M1 = 0.0
            res.M2 = 0.0
            res.M3 = 0.0
            res.M4 = 0.0
            return res


        delta = b.M1 - self.M1
        delta2 = delta * delta
        delta3 = delta * delta2
        delta4 = delta2 * delta2

        res.M1 = tk_stats.robust_divide(self.count * self.M1 + b.count * b.M1, res.count)

        res.M2 = self.M2 + b.M2 + delta2 * self.count * b.count / res.count

        res.M3 = self.M3 + b.M3 + \
                delta3 * self.count * b.count * (self.count - b.count) / ( res.count * res.count ) + \
                3.0 * delta * (self.count * b.M2 - b.count * self.M2) / res.count

        res.M4 = self.M4 + b.M4 + delta4 * self.count * b.count * \
                ( self.count * self.count - self.count * b.count + b.count * b.count ) / \
                ( res.count * res.count * res.count ) + \
                6.0 * delta2 * ( self.count * self.count * b.M2 + b.count * b.count * self.M2) / \
                ( res.count * res.count ) + \
                4.0 * delta * ( self.count * b.M3 - b.count * self.M3 ) / res.count

        return res

class SummaryManager:
    def __init__(self):
        self.summaries = {}
        self.attributes = {}

    def get_summarizer(self, summarizer_name):
        return self.summaries.get(summarizer_name)

    def add_summarizer(self, summarizer, name):
        if not isinstance(summarizer, SummaryStatistic):
            raise Exception('"summarizer" must inheret SummaryStatistic')

        if name in self.summaries:
            raise Exception('Summarizer ' + str(name) + ' already exists in this SummaryManager')

        self.summaries[name] = summarizer

        return self.summaries[name]

    def save(self, fname):
        with open(fname, 'wb') as fhandle:
            #pickle.dump( self.summaries, fhandle )
            pickle.dump( self, fhandle )
        return

    @staticmethod
    def load(fname):
        sm = None

        with open(fname, 'rb') as fhandle:
            sm = pickle.load( fhandle )

        return sm

    def combine(self, sm):
        if set(self.summaries.keys()) != set(sm.summaries.keys()):
            raise Exception("Don't know how to combine summary managers where keys aren't equal!")

        res = SummaryManager()
        for key in self.summaries:
            res.summaries[key] = self.summaries[key].combine( sm.summaries[key] )

        res.attributes.update(sm.attributes)

        return res

    def set_attrib(self, k, v):
        self.attributes[k] = v

        return self.attributes[k]

    def get_attrib(self, k):
        return self.attributes.get(k)

#     class TopKSummary:
#         def __init__(self, k):
#             self.k = k
#             self.items = []
#
#         def accept(self, k, v):
#             for i in range(len(self.items)):
#                 tk,tv = self.items[i]
#                 if x > tk:
#                     self.items.insert(i,(k,v))
#                     if len(self.items) > k:
#                         self.items = self.items[:k]
#                     return
#
#             if len(self.items) < k:
#                 self.items.append((k,v))
#
