#!/usr/bin/env python3
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

from __future__ import annotations

import collections
import glob
import os
import re
import shutil
import subprocess
import xml.etree.ElementTree as etree
from collections.abc import Iterable
from enum import Enum, auto, unique
from typing import TypedDict

import martian


@unique
class Bcl2FastqVersion(Enum):
    V1 = auto()
    V2 = auto()

    def __str__(self):
        # These used to be represented as old/new, preserving for error messages
        if self == Bcl2FastqVersion.V1:
            return "OLD"
        elif self == Bcl2FastqVersion.V2:
            return "NEW"
        else:
            raise RuntimeError()


@unique
class IlluminaSequencer(Enum):
    iSeq = auto()
    MiSeq = auto()
    HiSeq = auto()
    NextSeq = auto()
    NovaSeq = auto()
    NovaSeqX = auto()
    Unknown = auto()


class RTAVersionInformation:
    """Used to capture information about an Illumina run, in particular the type of instrument.

    If the R2 needs to be reverse complemented, and what the version of the RTA software was.
    """

    def __init__(self, input_path: str):
        """Create information about the RTA version used.

        Rather complicated due to all the changes in the Illumina
        output files over the years.  This querie the BCL folder for the RTA version of the run, and optionally also
        finds out whether the I2 read needs to be reverse complemented.

        Args:
            input_path: Path to a RunParameters.xml file
        """
        self._input_path = input_path
        # We'll set this later on if we can
        self.rc_i2 = None

        rp_nextseq_novaseq = os.path.join(input_path, "RunParameters.xml")
        rp_other = os.path.join(input_path, "runParameters.xml")
        if os.path.exists(rp_nextseq_novaseq):
            run_parameters_xml = rp_nextseq_novaseq
        elif os.path.exists(rp_other):
            run_parameters_xml = rp_other
        else:
            martian.exit(
                f"Neither a RunParameters.xml or runParameters.xml file was found in {input_path}"
            )
            raise SystemExit()

        tree = etree.parse(run_parameters_xml)

        # TENKIT-60 NovaSeq / NovaSeqX RunParameters.xml doesn't have the "Setup" node
        root = tree.getroot()
        self._root = root
        setup_node = root.find("Setup")
        if setup_node is not None:
            self.application = setup_node.find("ApplicationName").text
            self.application_version = setup_node.find("ApplicationVersion").text
        else:
            # CSR-477 iSeq has a new runParameters variant!
            application_name = root.find("ApplicationName")
            if application_name is not None:
                self.application = application_name.text
            else:
                self.application = root.find("Application").text

            application_version = root.find("ApplicationVersion")
            if application_version is not None:
                self.application_version = application_version.text
            else:
                # Nova Seq X introduced this
                self.application_version = root.find("SystemSuiteVersion").text
        assert self.application is not None
        assert self.application_version
        # Now we can figure out if we need to RC I2
        self._determine_rc_i2()

        # Sequencers prior to the NovaSeq X have a RTA Version, those after don't.
        if self.sequencer == IlluminaSequencer.NovaSeqX:
            self.rta_version = None
        else:
            # Can we run bcl2fastq 2, or do we need to use the old version?
            if setup_node is not None:
                rta_tag = tree.getroot().find("Setup").find("RTAVersion")
                if rta_tag is None:
                    rta_tag = tree.getroot().find("RTAVersion")
            # TENKIT-60 NovaSeq lowercases the tag
            else:
                rta_tag = tree.getroot().find("RtaVersion")

            self.rta_version = rta_tag.text
            assert self.rta_version is not None
            if self.rta_version.startswith("v"):
                self.rta_version = self.rta_version[1:]

    def log_version_info(self) -> None:
        """Logs information about the RTA version to martian."""
        bcl_params = {
            "ApplicationName": self.application,
            "ApplicationVersion": self.application_version,
        }
        martian.log_info(f"BCL folder RTA Version: {self.rta_version}")
        martian.log_info(f"BCL params: {bcl_params!s}")
        martian.log_info(f"RC'ing i2 read: {self.rc_i2!s}")

    @property
    def sequencer(self) -> IlluminaSequencer:
        """Returns the sequencer we've inferred from the [r|R]unParameters.xml file."""
        if self.application.count("NextSeq") > 0:
            return IlluminaSequencer.NextSeq
        elif self.application.count("MiSeq") > 0:
            return IlluminaSequencer.MiSeq
        # Note that NovaSeqX seems to use <Application>control-software</Application> and so this
        # checking for NovaSeqX is just in case they change this in the future
        elif (
            self.application.count("NovaSeq") > 0
            and self.application.count("NovaSeqX") == 0
            and self.application.count("NovaSeq X") == 0
        ):
            return IlluminaSequencer.NovaSeq
        elif self.application.count("HiSeq") > 0:
            return IlluminaSequencer.HiSeq
        elif self.application.count("iSeq") > 0:
            return IlluminaSequencer.iSeq
        elif (
            self.application.count("control-software") > 0
            or self.application.count("NovaSeqX") > 0
            or self.application.count("NovaSeq X") > 0
        ):
            # Doing a second check to try and future proof this
            # <InstrumentType>NovaSeqXPlus</InstrumentType>
            assert (
                self._root.find("InstrumentType").text.count("NovaSeqX") > 0
            ), "Found `control-software` for application but couldn't confirm this was NovaSeqX"
            return IlluminaSequencer.NovaSeqX
        else:
            return IlluminaSequencer.Unknown

    def _determine_rc_i2(self) -> None:
        """Do we need to Reverse Complement the I2 Read? Figuring this out is very sequencer dependent."""
        # Do we need to RC the I2 read?
        # Our current understanding is that NextSeq and HiSeq X / 4000 require it,
        sequencer = self.sequencer
        if sequencer == IlluminaSequencer.NextSeq:
            self.rc_i2 = True
        elif sequencer == IlluminaSequencer.MiSeq:
            self.rc_i2 = False
        # according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-03.pdf
        # NovaSeq follows MiSeq/HiSeq 2500 workflow for I5 index; this is effectively
        # a noop thanks to rc_i2 being False by default but let's make it explicit
        elif sequencer == IlluminaSequencer.NovaSeq:
            # Note this can set it to None
            self.rc_i2 = _detect_rc_i2_via_recipe_novaseq(
                self._input_path, self.application_version
            )
        elif sequencer == IlluminaSequencer.HiSeq:
            # Hiseq 4000 has version 3 -- hopefully this is stable??
            app_str = re.search(r"[\d\.]+", self.application_version).group()
            main_app_ver = int(app_str.split(".")[0])
            if main_app_ver > 2:
                self.rc_i2 = True
            else:
                self.rc_i2 = False
        elif sequencer == IlluminaSequencer.NovaSeqX:
            self.rc_i2 = self._get_rci2_from_newer_xml()
        # according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-07.pdf
        # iSeq 100 = rc.  index has to be exactly 0, to avoid MiSeq conflict.
        elif sequencer == IlluminaSequencer.iSeq:
            self.rc_i2 = True
        elif sequencer == IlluminaSequencer.Unknown:
            martian.log_info(
                "Was not able to determine the sequencer used from the RunParameters.xml file."
            )
            self.rc_i2 = None
        else:
            raise RuntimeError("Invalid code path encountered.")

    def _get_rci2_from_newer_xml(self):
        """Determine whether we need to RC the I2 read by parsing the XML files in the run directory.

        The NovaSeqX and NextSeq 1000/2000 contain that IsReverseComplement tag.   Presumably any theoretical instruments released
         in the future will also have it, but none of the earlier instruments do.

         This method makes lots of assumptions about the file format, which we hope are true.

         We're looking to parse this from the RunInfo.xml file:
         <Reads>
            <Read Number="1" NumCycles="50" IsIndexedRead="N" IsReverseComplement="N"/>
            <Read Number="2" NumCycles="8" IsIndexedRead="Y" IsReverseComplement="N"/>
            <Read Number="3" NumCycles="24" IsIndexedRead="Y" IsReverseComplement="Y"/>
            <Read Number="4" NumCycles="49" IsIndexedRead="N" IsReverseComplement="N"/>
        </Reads>

        And this from the RunParameters.xml file:
          <PlannedReads>
            <Read ReadName="Read1" Cycles="50" />
            <Read ReadName="Index1" Cycles="8" />
            <Read ReadName="Index2" Cycles="24" />
            <Read ReadName="Read2" Cycles="49" />
          </PlannedReads>

        See CELLRANGER-7366
        Returns:
            Boolean
        """
        reads = self._root.find("PlannedReads")
        if reads is None:
            martian.exit("Expected a PlannedReads node in the XML file that wasn't found.")
            raise SystemExit()
        cycles = None
        read_index = None
        for read_index, child in enumerate(reads):
            read_name = child.attrib["ReadName"]
            if read_name == "Index2":
                cycles = child.attrib["Cycles"]
                break
        if cycles is None:
            # No read2 found? Log it and assume there is no I2 read.
            martian.log_info("Did not find a 'Index2' read in the RunParameters.xml file")
            return None
        run_info = os.path.join(self._input_path, "RunInfo.xml")
        if not os.path.exists(run_info):
            martian.exit(f"Expected to find file: {run_info} but it doesn't exist.")
            raise SystemExit()
        tree = etree.parse(run_info)
        reads = tree.getroot().find("Run").find("Reads")
        if len(reads) < read_index + 1:
            martian.exit(
                f"Expected to find a Read entry at index {read_index} of PlannedReads element in RunInfo.xml file, but it wasn't there."
            )
        index2_read = reads[read_index]
        # Assuming they're in the same order, checking that here with some guardrails
        if index2_read.attrib["NumCycles"] != cycles:
            martian.exit(
                "Read in XML position did not match expected cycles, contact support@10xgenomics.com"
            )
        if int(index2_read.attrib["Number"]) != read_index + 1:
            martian.exit(
                "Read in XML position did not match expected read number, contact support@10xgenomics.com"
            )
        if index2_read.attrib["IsIndexedRead"] != "Y":
            martian.exit(
                "Did not find an indexed read where expected, contact support@10xgenomics.com"
            )
        rc_txt = index2_read.attrib["IsReverseComplement"]
        return bool(rc_txt)

    def check_bcl2fastq(self, hostname: str) -> tuple[Bcl2FastqVersion | str]:
        """Determine the best available bcl2fastq version to use .

        Will call martian.exit() with an error message if
        there isn't a compatible version available.
        """
        # Try to get both, returns None  for v1/v2 if it doesn't work
        v1 = _get_bcl2fastq_v1(hostname)
        v2, v2_msg = _get_bcl2fastq_v2(hostname)
        if v1 is None and v2 is None:
            martian.exit(
                "No valid bcl2fastq found on path. Recommended version of bcl2fastq is v2.20.\n"
            )
            raise SystemExit()
        if v2 is not None:
            v2x, v2y = [int(v2_part) for v2_part in v2.split(b".")][0:2]

        not_novaseq_x = self.sequencer != IlluminaSequencer.NovaSeqX
        # The NovaSeq X does not have a RTAVersion to go with it
        if not_novaseq_x:
            x, y, z = [int(xi) for xi in self.rta_version.split(".")][0:3]
            if x == 1 and ((y < 18) or ((y == 18) and z < 54)):
                # RTA <1.18.54
                # must run 1.8.4
                if v1 is not None:
                    return Bcl2FastqVersion.V1, v1
                else:
                    msg = f"mkfastq requires bcl2fastq 1.8.4 for RTA version: {self.rta_version}"
                    martian.exit(msg)
                    raise SystemExit()
        # RTA >= 1.18.54 -> run 2.17 or higher
        # NovaSeq X -> run 2.20 or higher
        if v2 is None:
            msg = f"No valid bcl2fastq found on path. Recommended version of bcl2fastq is v2.20.\n\n{v2_msg}"
            martian.exit(msg)
            raise SystemExit()

        minv2y = 17 if not_novaseq_x else 20
        if v2x == 2 and v2y >= minv2y:
            return Bcl2FastqVersion.V2, v2
        else:
            martian.exit(
                "Incorrect bcl2fastq version found: %s. Recommended version of bcl2fastq is v2.20."
            )
            raise SystemExit()


def _get_bcl2fastq_v1(hostname: str) -> str | None:
    if shutil.which("configureBclToFastq.pl"):
        if not shutil.which("perl"):
            martian.log_info(
                f"On machine: {hostname}, perl not found on PATH. (Required for configureBclToFastq.pl).  "
                "This could indicate a newer bcl2fastq is being used."
            )
            return None
        return "1.8.4"
    else:
        martian.log_info(
            f"On machine: {hostname}, configureBclToFastq.pl not found on PATH.  "
            "This could indicate a newer bcl2fastq is being used."
        )
        return None


def _get_bcl2fastq_v2(hostname: str) -> tuple[bytes, None] | tuple[None, str]:
    if bcl2fastq := shutil.which("bcl2fastq"):
        # Restore the LD_LIBRARY_PATH set aside by sourceme.bash/shell10x.
        # Required for some installations of bcl2fastq.
        new_environ = dict(os.environ)
        new_environ["LD_LIBRARY_PATH"] = os.environ.get("_TENX_LD_LIBRARY_PATH", "")
        try:
            output = subprocess.check_output(
                [bcl2fastq, "--version"], env=new_environ, stderr=subprocess.STDOUT
            )
        except subprocess.CalledProcessError:
            msg = f"On machine: {hostname}, bcl2fastq does not work."
            return (None, msg)
        for l in output.split(b"\n"):
            match = re.match(b"bcl2fastq v([0-9.]+)", l)
            if match is not None:
                return match.groups()[0], None

        return (
            None,
            "bcl2fastq version not recognized -- please check the output of bcl2fastq --version",
        )
    else:
        msg = f"On machine: {hostname}, bcl2fastq not found on PATH."
        return (None, msg)


def _detect_rc_i2_via_recipe_novaseq(flowcell_path: str, application_version) -> bool | None:
    """Determine if I2 is RC by inspecting the Recipe xml file.

    Finds the recipe xml file given the flowcell_path. Return true if I2 is RC,
    false if RC is fwd-complement, and None if the recipe xml couldn't be
    found / read.
    """
    app_str = re.search(r"[\d\.]+", application_version).group()
    app_ver = [int(x) for x in app_str.split(".")]

    # We use a different detection scheme in ICS 1.8.0 and above.
    new_novaseq_detect_mode = False
    if app_ver >= [1, 8, 0]:
        new_novaseq_detect_mode = True

    run_info_xml = os.path.join(flowcell_path, "RunInfo.xml")
    if not os.path.exists(run_info_xml):
        return None
    (_, flowcell) = load_run_info(run_info_xml)

    recipe_file = None

    novaseq_recipe = os.path.join(flowcell_path, "Recipe", flowcell + ".xml")
    if os.path.exists(novaseq_recipe):
        recipe_file = novaseq_recipe

    if recipe_file is None:
        # Other sequencers (eg NextSeq) have a different convention
        # try and find any recipe file.
        xmls = glob.glob(os.path.join(flowcell_path, "Recipe", "*.xml"))
        if len(xmls) > 0:
            # We have no idea what >1 xml file here means
            # so we just pick one
            recipe_file = xmls[0]

    if recipe_file is not None:
        if new_novaseq_detect_mode:
            return _new_novaseq_detect_rc_i2_from_recipe_xml(recipe_file)
        else:
            return _old_novaseq_detect_rc_i2_from_recipe_xml(recipe_file)
    else:
        return None


def _old_novaseq_detect_rc_i2_from_recipe_xml(recipe_xml: str) -> bool:
    """Determine if I2 is RC by inspecting the Recipe xml file.

    Based on a scheme from Illumina, if the "IndexPreparation-i5" or "Index2Preparation" ChemistryStep
    exists and uses the "BP14" reagent, the the I2 read is RC.
    """
    tree = etree.parse(recipe_xml)
    chem_steps = [x for x in tree.iter() if x.tag == "ChemistryStep"]
    i5_steps = [x for x in chem_steps if x.attrib.get("Description").startswith("Index")]

    reagents = set()
    for step in i5_steps:
        for el in step.iter():
            r = el.attrib.get("ReagentName")
            if r:
                reagents.add(r)

    return "BP14" in reagents


def _new_novaseq_detect_rc_i2_from_recipe_xml(recipe_xml) -> bool | None:
    """New scheme for detecting workflow mode in NovaSeq SW version 1.8.0 and above.

    Parse if  `<ChemistryRef Description="PETurnaround" ChemistryName="PETurnaround" />`
    comes before or after `<ReadRef Description="IndexRead i5" ReadName="IndexRead2" />`.
    If PETuraround comes after, then you're in workflow A, if it comes before IndexRead2workflow B.
    """
    tree = etree.parse(recipe_xml)
    protocol = next(x for x in tree.iter() if x.tag == "Protocol")

    # find the 2 key steps in the protocol list
    index_read2 = [
        x
        for x in enumerate(protocol.iter())
        if x[1].tag == "ReadRef" and x[1].attrib.get("ReadName") == "IndexRead2"
    ]
    turnaround = [
        x
        for x in enumerate(protocol.iter())
        if x[1].tag == "ChemistryRef" and x[1].attrib.get("ChemistryName") == "PETurnaround"
    ]

    # Check which step comes first
    if len(index_read2) == 1 and len(turnaround) == 1:
        if turnaround > index_read2:
            # workflow A
            return False
        else:
            # workflow B
            return True

    else:
        # We don't recognize the steps in the protocol, so bail
        return None


class ReadInfo(TypedDict):
    read_length: int
    read_name: str
    index_read: bool
    original_read_name: str


def load_run_info(run_info_xml: str) -> tuple[list[ReadInfo], str | None]:
    """Get the read names and read lengths from the Illumina RunInfo.xml file."""
    tree = etree.parse(run_info_xml)
    reads_node = tree.getroot().find("Run").find("Reads")
    reads = reads_node.findall("Read")

    # Now we follow some hard-coded conventions on order in which reads appear.
    # R1 is always first
    # R2 is always last
    # Index reads (if they exist) are in the middle, in numerical order

    # NOTE -- if there is only one index read it is assumed to be I1.
    # BclToFastq should give us this
    # NOTE -- we assume paired end reads!
    read_info = [
        ReadInfo(
            read_length=x,
            index_read=(idx not in (0, len(reads) - 1)) or is_indexed_read,
            read_name=(
                "R1"
                if idx == 0 and not is_indexed_read
                else "R2" if idx == len(reads) - 1 and not is_indexed_read else "I" + str(idx)
            ),
            original_read_name="R" + str(idx + 1),
        )
        for idx, (x, is_indexed_read) in enumerate(
            (int(read.attrib["NumCycles"]), read.attrib.get("IsIndexedRead") == "Y")
            for read in reads
        )
    ]
    flowcell = tree.getroot().find("Run").find("Flowcell").text

    # NB: currently you have to comment out the next two lines to get
    # nosetests to run correctly outside of a stage.
    martian.log_info(f"Read Info: {read_info}")
    martian.log_info(f"Flowcell ID: {flowcell}")
    return (read_info, flowcell)


def make_bases_mask_val(
    read_info: Iterable[ReadInfo],
    barcode_read: str | None = None,
    sample_index_read: str | None = None,
    dual_indexed: bool = False,
    ignore_dual_index: bool = False,
) -> str:
    """Undocumented.

    Args:
        read_info: The ReadInfo block from RunInfo.xml
        barcode_read: The read to use as the barcode. Can be an index.
        sample_index_read: The ReadInfo read (I1, I2) to use as the sample index
        dual_indexed: If the input BCLs were dual-indexed intentionally, then
                      preserve the index status of the 2nd, non-sample index index
                      in the mask.  Too often, people inadvertently ran dual-indexed
                      values for the barcode read (Cell Ranger v1, GemCode), so the
                      default behavior is to treat the off-SI indexed read as a
                      normal, non-index read.
        ignore_dual_index: Stub out any dual index with Ns.
    """

    # We will emit all reads in RunInfo.xml
    def base_mask(read):
        if read["read_name"][0] == "R":
            return "Y" + str(read["read_length"])
        elif read["read_name"][0] == "I":
            if read["read_name"] == barcode_read:
                return "Y" + str(read["read_length"])
            elif read["read_name"] == sample_index_read:
                return "I" + str(read["read_length"])
            elif dual_indexed:
                if ignore_dual_index:
                    return "N" + str(read["read_length"])
                else:
                    return "I" + str(read["read_length"])
            else:
                return "Y" + str(read["read_length"])
        else:
            martian.throw("read name was not recognized: {}".format(read["read_name"]))
            raise SystemExit()

    # Special hack to convert the bases_mask
    # to only give 8bp in I1, when in dual indexing
    # mode but using `ignore_dual_index` (aka --filter-single-index)
    # This has the following effect on the a bases mask:
    # Y28,I10,I10,Y90 -> Y28,I8,N12,Y90.
    # See details in CELLRANGER-3909
    if dual_indexed and ignore_dual_index:
        print("original read_info: ", read_info)
        rr = collections.OrderedDict((r["read_name"], r) for r in read_info)
        old_i1_length = rr["I1"]["read_length"]
        rr["I1"]["read_length"] = 8
        rr["I2"]["read_length"] = rr["I2"]["read_length"] + (old_i1_length - 8)
        read_info = list(rr.values())

        print("edited read info for i1 only: ", read_info)

    masks = [base_mask(r) for r in read_info]
    return ",".join(masks)


def get_bcl2fastq_read_type_map(
    read_info: Iterable[ReadInfo],
    barcode_read: str | None = None,
    sample_index_read: str | None = None,
    dual_indexed: bool = False,
    ignore_dual_index: bool = False,
) -> dict[str, str]:
    """Get a mapping between read name and output file.

    Read name from ReadInfo (R1,I1,I2,R2) and bcl2fastq
    output file naming (R1/R2/R3/R4/I1)

    The guarantee here is that the 10X sample index will always be on I1,
    if generated.  If dual-indexing is specified, the secondary index will
    be on I2.  Upstream pipestances can expect read1 to be in R1, sample
    indexes to be on I1, and read2 to be on R2.

    Args:
        read_info: The ReadInfo block from RunInfo.xml
        barcode_read: The ReadInfo read to use as the barcode (can be I2)
        sample_index_read: The ReadInfo read (I1, I2) to use as the sample index
    """
    # read_names = [r["read_name"] for r in read_info]
    read_map = {}
    reads_counted = 0
    for r in read_info:
        read_name = r["read_name"]
        assert isinstance(read_name, str)
        if read_name == sample_index_read:
            read_map[read_name] = "I1"
        elif dual_indexed and r["index_read"]:
            # handle I2 == barcode_read case (ATAC)
            # ignore_dual_index will be True since ATAC sample
            # index is I7-only, so this goes first
            if read_name == barcode_read:
                reads_counted += 1
                read_map[read_name] = "R%d" % reads_counted
            # ignore I2 if ignore_dual_index specified
            elif ignore_dual_index:
                continue
            else:
                read_map[read_name] = "I2"
        else:
            reads_counted += 1
            read_map[read_name] = "R%d" % reads_counted
    return read_map
