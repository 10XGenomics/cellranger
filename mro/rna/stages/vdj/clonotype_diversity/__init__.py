#!/usr/bin/env python
#
# Copyright (c) 2020 10X Genomics, Inc. All rights reserved.
#
# Calculate data required for rarefaction and extrapolation curves.
#

"""A stage for calculating clonotype diversity plot from clonotypes.csv.

Classes:
    Clonotypes: used to load information from clonotypes.csv for the purpose of
    plotting clonotype diversity plots.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import martian
import numpy as np
from scipy.special import loggamma

if TYPE_CHECKING:
    import cellranger.mro_types.filetypes as mro_filetypes

__MRO__ = """
stage CALCULATE_CLONOTYPE_DIVERSITY(
    in  json per_origin_clonotype_hist,
    out json plotly_diversity_chart,
    src py   "stages/vdj/clonotype_diversity",
)
"""


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CALCULATE_CLONOTYPE_DIVERSITY stage inputs."""

    per_origin_clonotype_hist: mro_filetypes.json


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """CALCULATE_CLONOTYPE_DIVERSITY stage outputs."""

    plotly_diversity_chart: mro_filetypes.json


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """CALCULATE_CLONOTYPE_DIVERSITY chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """CALCULATE_CLONOTYPE_DIVERSITY chunk outputs."""


# import plotly.express as px
# print(px.colors.qualitative.D3)
PLOTLY_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
]

# TODO: This is the only element in the web summary created by python code because we
# don't have the rarefaction math implemented in Rust. Would be cleaner to move this
DIVERSITY_CHART_TEMPLATE = {
    "plot": {
        "config": {
            "displayModeBar": True,
            "staticPlot": False,
            "dragmode": "zoom",
            "modeBarButtons": [["toImage"]],
        },
        "data": [],
        "layout": {
            "showlegend": True,
            "margin": {"r": 40, "l": 60, "t": 25},
            "hovermode": "closest",
            "xaxis": {"title": "Number of Cells"},
            "yaxis": {"title": "Number of Unique Clonotypes"},
            "legend": dict(yanchor="top", y=0.99, xanchor="left", x=0.01, bgcolor="rgba(0,0,0,0)"),
        },
    },
    "help": {
        "helpText": "This plot displays rarefaction curves of clonotype diversity. Each curve is from 1 origin.",
        "title": "Clonotype diversity",
    },
}


def main(args: ChunkInputs, outs: ChunkOutputs):
    with open(args.per_origin_clonotype_hist) as f:
        per_origin_clonotype_hist = json.load(f)

    data = []
    for index, (origin, hist) in enumerate(per_origin_clonotype_hist):
        # In the very unlikely case of more than 10 origins, we will cycle through colors
        index = index % len(PLOTLY_COLORS)
        values = hist.values()
        data.extend(
            Diversity(values).calc_rarefaction_curve_plotly(
                origin, PLOTLY_COLORS[index], num_steps=50, stdev=False
            )
        )

    chart = DIVERSITY_CHART_TEMPLATE
    chart["plot"]["data"] = data

    with open(outs.plotly_diversity_chart, "w") as out_file:
        json.dump(chart, out_file, indent=4)


class Curve:
    """Curve information for plotting.

    Attributes:
        x     (list of int): Curve x-axis (number of cells)
        y     (list of float): Curve y-axis (expected value of number of
              unique clonotypes)
        y_std (list of float): Standard deviation of number of unique
              clonotypes in a curve
        y_ciu (list of float): Upper bound of 95% confidence interval
              for number of unique clonotypes in a curve
        y_cil (list of float): Lower bound of 95% confidence interval
              for number of unique clonotypes in a curve
    """

    def __init__(self):
        """Initialize an empty curve."""
        self.x: list[int] = []
        self.y: list[float] = []
        self.y_std: list[float] | list[None] = []
        self.y_ciu: list[float] | list[None] = []
        self.y_cil: list[float] | list[None] = []

    def is_empty(self) -> bool:
        """Check if curve is empty or calculated.

        Returns:
            bool
        """
        if len(self.y) == 0 or len(self.x) == 0:
            return True
        return False

    def is_consistent(self) -> bool:
        """Check if calculated curve is consistent (length match).

        Returns:
            bool
        """
        len_list = [len(i) for i in [self.x, self.y, self.y_cil, self.y_ciu]]
        if len_list.count(len_list[0]) != len(len_list):
            return False
        return True


class Diversity:
    """Represents diversity with rarefaction and extrapolation curves.

    Attributes:
        sorted_hist         (list of tuple): A histogram sorted by abundance
                                            (most abundant has lowest index)
        freq_counts         (dict of int: int): Frequency counts table
        N                   (int): total number of counts in histogram
        rarefaction_curve   (Curve): Rarefaction curve information
        extrapolation_curve (Curve): Extrapolation curve information
    """

    def __init__(self, hist):
        self.sorted_hist = sorted(hist, reverse=True)
        self.N: int = self._calc_n()
        self.freq_counts: dict[int, int] = self._get_freq_counts()
        # Rarefaction curve
        self.rarefaction_curve: Curve = Curve()
        # Extrapolation curve
        self.extrapolation_curve: Curve = Curve()

    def is_diversity_curve_possible(self) -> bool:
        if self.N == 0 or self.f_0_chao1() in [0, -1]:
            return False
        return True

    def _get_freq_counts(self) -> dict[int, int]:
        """Get frequency counts (f_k in Colwell et al 2012).

        Returns:
            dict
        """
        ret_dict: dict[int, int] = {}
        for abund in self.sorted_hist:
            if abund not in ret_dict:
                ret_dict[abund] = 1
            else:
                ret_dict[abund] += 1
        return ret_dict

    def _calc_n(self) -> int:
        """Calculates number of samples.

        (length of input histogram)

        Returns:
            int
        """
        return sum(self.sorted_hist)

    def _alpha(self, n: int, k: int) -> np.float64 | int:
        """Calculates alpha parameter in Colwell et al 2012.

        Args:
            n (int): rarefaction point n
            k (int): Frequency

        Returns:
            int
        """
        if k > self.N - n:
            return 0
        log_alpha: np.float64 = (
            loggamma(self.N - k + 1)
            + loggamma(self.N - n + 1)
            - loggamma(self.N + 1)
            - loggamma(self.N - k - n + 1)
        )
        return np.exp(np.real(log_alpha))

    # Minimum variance unbiased estimator (MVUE) for both hypergeometric
    # and multinomial models
    # eq (4) in Colwell et al 2012
    # eq (5) in Colwell et al 2012
    def _rarefaction(self, n: int) -> tuple[np.float64 | float, np.float64 | float]:
        """Calculates rarefaction and standard deviation of rarefaction for a single point n.

        Args:
            n (int): rarefaction point n

        Returns:
            float, float
        """
        # number of observed clonotypes (S_obs in paper)
        num_clon = float(sum(self.freq_counts.values()))
        sum_helper_exp_val = 0
        for k, count in self.freq_counts.items():
            sum_helper_exp_val += self._alpha(n, k) * count
        exp_val = num_clon - sum_helper_exp_val
        sum_helper_std_dev = 0
        for k, count in self.freq_counts.items():
            sum_helper_std_dev += (1 - self._alpha(n, k)) ** 2 * count - float(
                exp_val**2
            ) / self.assemblage_size_estimate()
        return exp_val, np.sqrt(sum_helper_std_dev)

    def calc_rarefaction_curve(self, num_steps: int = 40):
        """Calculates rarefaction curve for num_steps between 1 and N.

        Args:
            num_steps (int): Number of steps between 1 and N
        """
        step_size = int(self.N // (num_steps - 1))
        if step_size == 0:
            step_size = 1
        rc_x = range(1, self.N, step_size)
        if len(rc_x) < num_steps:
            rc_x = range(1, self.N + step_size, step_size)

        placeholder = 0.0
        rc_y_exp = [placeholder] * len(rc_x)
        rc_y_std = [placeholder] * len(rc_x)
        rc_y_ciu = [placeholder] * len(rc_x)
        rc_y_cil = [placeholder] * len(rc_x)
        for idx, n in enumerate(rc_x):
            rc_y_exp[idx], rc_y_std[idx] = self._rarefaction(n)
            rc_y_cil[idx] = rc_y_exp[idx] - 1.96 * rc_y_std[idx]
            rc_y_ciu[idx] = rc_y_exp[idx] + 1.96 * rc_y_std[idx]
        self.rarefaction_curve.x = list(rc_x)
        self.rarefaction_curve.y = rc_y_exp
        self.rarefaction_curve.y_std = rc_y_std
        self.rarefaction_curve.y_cil = rc_y_cil
        self.rarefaction_curve.y_ciu = rc_y_ciu

    def calc_rarefaction_curve_plotly(
        self, origin: str, color, num_steps: int = 40, stdev: bool = True
    ):
        """Create a plotly curve for rarefaction.

        Args:
            origin (str): Origin string (see vdj inputs)
            color: color of the curve
            num_steps (int): number of extrapolation steps
            stdev (bool): flag for calculation of standard deviation (and confidence intervals)

        Returns:
            [dict]
        """
        self.calc_rarefaction_curve(num_steps)
        lines = []
        lines.append(
            {
                "x": self.rarefaction_curve.x,
                "y": self.rarefaction_curve.y,
                "type": "scatter",
                "name": origin,
                "mode": "lines",
                "line_color": color,
            }
        )
        if stdev:
            lines.append(
                {
                    "x": self.rarefaction_curve.x + self.rarefaction_curve.x[::-1],
                    "y": self.rarefaction_curve.y_ciu + self.rarefaction_curve.y_cil[::-1],
                    "type": "scatter",
                    "name": origin + " (Stdev)",
                    "fill": "toself",
                    "line_color": "rgba(255,255,255,0)",
                    "fillcolor": color,
                    "opacity": 0.2,
                }
            )
        return lines

    # eq (9) in Colwell et al 2012
    # not using the approximation used in eq (9)
    def _extrapolation(self, n_plus_m):
        """Calculates extrapolation for a single point (N + m).

        standard deviation is not implemented.

        Args:
            n_plus_m (int): Distance of extrapolation point to N
                (how many more samples collected)

        Returns:
            float: Expected value
            None: Standard deviation (not implemented)
        """
        # number of observed clonotypes (s_obs in paper)
        s_obs = float(sum(self.freq_counts.values()))
        f_0 = float(self.f_0_chao1())
        f_1 = float(self.freq_counts[1])

        N = float(self.N)  # n in paper
        m = float(n_plus_m - self.N)

        brackets = 1.0 - (1.0 - f_1 / N / f_0) ** m

        exp_val = s_obs + f_0 * brackets
        # Standard deviation for extrapolation is not implemented
        std_dev = None

        return exp_val, std_dev

    def calc_extrapolation_curve(self, num_steps: int = 40, max_n: int = 2000):
        """Calculates extrapolation curve for num_steps between N and max_n.

        Args:
            num_steps (int): Number of steps between N and max_n
            max_n (int): The limit to which the function extrapolates

        Returns:
            None
        """
        step_size = int((max_n - self.N) / (num_steps - 1))
        ec_x = range(self.N, max_n, step_size)
        if len(ec_x) < num_steps:
            ec_x = range(self.N, max_n + step_size, step_size)
        ec_y_exp = [0.0] * len(ec_x)
        ec_y_std = [None] * len(ec_x)
        ec_y_ciu = [None] * len(ec_x)
        ec_y_cil = [None] * len(ec_x)
        for idx, n_plus_m in enumerate(ec_x):
            ec_y_exp[idx], ec_y_std[idx] = self._extrapolation(n_plus_m)
            if ec_y_std[idx] is not None:
                ec_y_cil[idx] = ec_y_exp[idx] - 1.96 * ec_y_std[idx]
                ec_y_ciu[idx] = ec_y_exp[idx] + 1.96 * ec_y_std[idx]
        self.extrapolation_curve.x = list(ec_x)
        self.extrapolation_curve.y = ec_y_exp
        self.extrapolation_curve.y_std = ec_y_std
        self.extrapolation_curve.y_cil = ec_y_cil
        self.extrapolation_curve.y_ciu = ec_y_ciu

    def f_0_chao1(self):
        """Estimate f0 (number of clonotypes that we haven't seen yet) using chao1 estimator.

        Based on: Colwell et al 2012

        Returns:
            float
        """
        if 1 in self.freq_counts and 2 in self.freq_counts:
            return float(self.freq_counts[1] ** 2) / float(2 * self.freq_counts[2])
        elif 1 in self.freq_counts:
            return float(self.freq_counts[1] * (self.freq_counts[1] - 1)) / 2.0
        else:
            # print "Error in f_0_chao1" # TODO log message
            return -1.0

    def assemblage_size_estimate(self):
        """Estimate assemblage size (number of total unique clonotypes) using chao1 estimator.

        Based on: Colwell et al 2012

        Returns:
            float
        """
        return self.N + self.f_0_chao1()

    def create_both_curves(
        self, outfile, rc_steps: int = 40, ec_steps: int = 40, max_n: int = 2000
    ):
        """Calculate rarefaction and extrapolation diversity curves.

        Save them in csv format in outfile.

        Args:
            outfile (str): Path to output file
            rc_steps (int): Rarefaction curve steps between 1 and N
            ec_steps (int): Extrapolation steps between N and max_n
            max_n (int): Maximum number of cells extrapolating to
        """
        # Check if extrapolation_curve is calculated. If not, call the appropriate func to calc
        if self.extrapolation_curve.is_empty():
            self.calc_extrapolation_curve(ec_steps, max_n)
        # Check if rarefaction_curve is calculated. If not, call the appropriate func to calc
        if self.rarefaction_curve.is_empty():
            self.calc_rarefaction_curve(rc_steps)
        if not self.rarefaction_curve.is_consistent():
            return None  # TODO add appropriate logging message
        if not self.extrapolation_curve.is_consistent():
            return None  # TODO add appropriate logging message
        header = "x,y,y_cil,y_ciu,type\n"
        with open(outfile, "w") as outhandle:
            outhandle.write(header)
            for idx, x in enumerate(self.rarefaction_curve.x):
                numbers = [
                    x,
                    self.rarefaction_curve.y[idx],
                    self.rarefaction_curve.y_cil[idx],
                    self.rarefaction_curve.y_ciu[idx],
                ]
                line = ",".join([str(n) for n in numbers] + ["rarefaction"]) + "\n"
                outhandle.write(line)
            for idx, x in enumerate(self.extrapolation_curve.x):
                numbers = [
                    x,
                    self.extrapolation_curve.y[idx],
                    self.extrapolation_curve.y_cil[idx],
                    self.extrapolation_curve.y_ciu[idx],
                ]
                line = ",".join([str(n) for n in numbers] + ["extrapolation"]) + "\n"
                outhandle.write(line)
        return 1
