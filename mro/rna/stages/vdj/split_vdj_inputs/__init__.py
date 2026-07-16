# Copyright (c) 2019 10X Genomics, Inc. All rights reserved.
"""Decide which is T vs B vs T(Gamma/Delta) inputs."""

from __future__ import annotations

from typing import TYPE_CHECKING

import martian

if TYPE_CHECKING:
    import cellranger.mro_types.structs as mro_structs

__MRO__ = """
stage SPLIT_VDJ_INPUTS(
    in  VdjInputs[]    vdj_inputs,
    in  ChemistryDef[] vdj_chemistry_defs,
    in  string[]       vdj_receptors,
    out VdjInputs      vdj_t_input,
    out ChemistryDef   vdj_t_chemistry_def,
    out string         vdj_t_receptor,
    out VdjInputs      vdj_t_gd_input,
    out ChemistryDef   vdj_t_gd_chemistry_def,
    out string         vdj_t_gd_receptor,
    out VdjInputs      vdj_b_input,
    out ChemistryDef   vdj_b_chemistry_def,
    out string         vdj_b_receptor,
    src py             "stages/vdj/split_vdj_inputs",
)
"""

TCR_RECEPTOR = "TR"
TCRGD_RECEPTOR = "TR_GD"
BCR_RECEPTOR = "IG"


class StageInputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SPLIT_VDJ_INPUTS stage inputs."""

    vdj_inputs: list[mro_structs.VdjInputs]
    vdj_chemistry_defs: list[mro_structs.ChemistryDef]
    vdj_receptors: list[str | None]


class StageOutputs(martian.Record):  # pylint: disable=too-few-public-methods
    """SPLIT_VDJ_INPUTS stage outputs."""

    vdj_t_input: None | mro_structs.VdjInputs
    vdj_t_chemistry_def: None | mro_structs.ChemistryDef
    vdj_t_receptor: None | str
    vdj_t_gd_input: None | mro_structs.VdjInputs
    vdj_t_gd_chemistry_def: None | mro_structs.ChemistryDef
    vdj_t_gd_receptor: None | str
    vdj_b_input: None | mro_structs.VdjInputs
    vdj_b_chemistry_def: None | mro_structs.ChemistryDef
    vdj_b_receptor: None | str


class ChunkInputs(StageInputs):  # pylint: disable=too-few-public-methods
    """SPLIT_VDJ_INPUTS chunk inputs."""


class ChunkOutputs(StageOutputs):  # pylint: disable=too-few-public-methods
    """SPLIT_VDJ_INPUTS chunk outputs."""


def main(args: ChunkInputs, outs: ChunkOutputs):
    assert len(args.vdj_inputs) > 0
    assert len(args.vdj_inputs) == len(args.vdj_chemistry_defs)
    assert len(args.vdj_inputs) == len(args.vdj_receptors)

    # this is paranoid, but ensure we've null outputs
    outs.vdj_t_input = None
    outs.vdj_t_gd_input = None
    outs.vdj_b_input = None
    for inp, chem, receptor in zip(args.vdj_inputs, args.vdj_chemistry_defs, args.vdj_receptors):
        if receptor == TCR_RECEPTOR or receptor is None:  # Denovo or TCR
            if outs.vdj_t_input is None:
                outs.vdj_t_input = inp
                outs.vdj_t_chemistry_def = chem
                outs.vdj_t_receptor = receptor
            else:
                assert (
                    outs.vdj_t_chemistry_def == chem
                ), f"mismatched VDJ-T chemistry: {chem} != {outs.vdj_t_chemistry_def}"
                assert (
                    outs.vdj_t_receptor == receptor
                ), f"mismatched VDJ(-T) receptor: {receptor} != {outs.vdj_t_receptor}"
                for k, v in outs.vdj_t_input.items():  # pylint: disable=invalid-name
                    if k == "sample_def":
                        continue
                    assert v == inp[k], f"mismatched VDJ-T {k}: {v} != {inp[k]}"
                outs.vdj_t_input["sample_def"].extend(inp["sample_def"])
        elif receptor == TCRGD_RECEPTOR:  # TCR Gamma/Delta
            if outs.vdj_t_gd_input is None:
                outs.vdj_t_gd_input = inp
                outs.vdj_t_gd_chemistry_def = chem
                outs.vdj_t_gd_receptor = receptor
            else:
                assert (
                    outs.vdj_t_gd_chemistry_def == chem
                ), f"mismatched VDJ-T chemistry: {chem} != {outs.vdj_t_gd_chemistry_def}"
                assert (
                    outs.vdj_t_gd_receptor == receptor
                ), f"mismatched VDJ(-T) receptor: {receptor} != {outs.vdj_t_gd_receptor}"
                for k, v in outs.vdj_t_gd_input.items():  # pylint: disable=invalid-name
                    if k == "sample_def":
                        continue
                    assert v == inp[k], f"mismatched VDJ-T {k}: {v} != {inp[k]}"
                outs.vdj_t_gd_input["sample_def"].extend(inp["sample_def"])
        elif receptor == BCR_RECEPTOR:
            if outs.vdj_b_input is None:
                outs.vdj_b_input = inp
                outs.vdj_b_chemistry_def = chem
                outs.vdj_b_receptor = receptor
            else:
                assert (
                    outs.vdj_b_chemistry_def == chem
                ), f"mismatched VDJ-B chemistry: {chem} != {outs.vdj_b_chemistry_def}"
                assert (
                    outs.vdj_b_receptor == receptor
                ), f"mismatched VDJ(-B) receptor: {receptor} != {outs.vdj_b_receptor}"
                for k, v in outs.vdj_b_input.items():  # pylint: disable=invalid-name
                    if k == "sample_def":
                        continue
                    assert v == inp[k], f"mismatched VDJ-B {k}: {v} != {inp[k]}"
                outs.vdj_b_input["sample_def"].extend(inp["sample_def"])
        else:
            martian.exit(f"Unknown receptor {receptor}")
