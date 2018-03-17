#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import subprocess
import re
import os
import martian
import xml.etree.ElementTree as etree
import tenkit.constants as tk_constants

BCL2FASTQ_V1 = "OLD"
BCL2FASTQ_V2 = "NEW"

def get_bcl2fastq_v1(hostname):
    try:
        subprocess.check_call(["which", "configureBclToFastq.pl"])

        try:
            subprocess.check_call(["which", "perl"])
        except subprocess.CalledProcessError:
            msg = "On machine: %s, perl not found on PATH. (Required for configureBclToFastq.pl)" % hostname
            return (None, msg)

        return ("1.8.4", None)
    except subprocess.CalledProcessError:
        msg = "On machine: %s, configureBclToFastq.pl not found on PATH." % hostname
        return (None, msg)

def get_bcl2fastq_v2(hostname):
    try:
        subprocess.check_call(["which", "bcl2fastq"])
        # Restore the LD_LIBRARY_PATH set aside by sourceme.bash/shell10x.
        # Required for some installations of bcl2fastq.
        new_environ = dict(os.environ)
        new_environ['LD_LIBRARY_PATH'] = os.environ.get('_TENX_LD_LIBRARY_PATH', '')
        output = subprocess.check_output(["bcl2fastq", "--version"], env=new_environ, stderr=subprocess.STDOUT)
        match = None
        for l in output.split("\n"):
            match = re.match("bcl2fastq v([0-9.]+)", l)
            if match is not None:
                return (match.groups()[0], None)

        return (None, "bcl2fastq version not recognized -- please check the output of bcl2fastq --version")
    except subprocess.CalledProcessError:
        msg = "On machine: %s, bcl2fastq not found on PATH." % hostname
        return (None, msg)


def check_bcl2fastq(hostname, rta_str):
    x,y,z = [int(x) for x in rta_str.split(".")][0:3]

    v1,err1 = get_bcl2fastq_v1(hostname)
    v2,err2 = get_bcl2fastq_v2(hostname)

    # RTA <1.18.54
    # must run 1.8.4
    if x == 1 and ((y < 18) or ((y==18) and z < 54)):
        if v1 is not None:
            return (BCL2FASTQ_V1, v1)
        else:
            msg = "demux requires bcl2fastq 1.8.4 for RTA version: %s" % rta_str
            martian.exit(msg)


    # RTA >= 1.18.54
    # run 2.17 or higher
    else:

        if v2 is not None:
            v2x,v2y = [int(v2_part) for v2_part in v2.split(".")][0:2]
        else:
            msg = "No bcl2fastq found on path. demux requires bcl2fastq v2.17 or greater for RTA version: %s" % rta_str
            martian.exit(msg)

        if v2x == 2 and v2y >= 17:
            return (BCL2FASTQ_V2, v2)
        else:
            msg = "Incorrect bcl2fastq version found: %s. demux requires bcl2fastq v2.17 or greater for RTA version: %s" % (v2, rta_str)
            martian.exit(msg)


def get_rta_version(input_path):
    ''' Query the BCL folder for the RTA version of the
        run, and whether the I2 read needs to be reverse
        complemented. '''

    rp_nextseq = os.path.join(input_path, "RunParameters.xml")
    rp_other = os.path.join(input_path, "runParameters.xml")

    if os.path.exists(rp_nextseq):
        run_parameters_xml = rp_nextseq
    else:
        run_parameters_xml = rp_other

    tree = etree.parse(run_parameters_xml)

    # Do we need to RC the I2 read?
    # Our current understanding is that NextSeq and HiSeq X / 4000 require it
    rc_i2 = False

    # TENKIT-60 NovaSeq runParameters.xml doesn't have the "Setup" node
    setup_node = tree.getroot().find("Setup")
    if setup_node is not None:
        application = tree.getroot().find("Setup").find("ApplicationName").text
        application_version = tree.getroot().find("Setup").find("ApplicationVersion").text
    else:
        application = tree.getroot().find("Application").text
        application_version = tree.getroot().find("ApplicationVersion").text

    if application.find("NextSeq") >= 0:
        rc_i2 = True
    elif application.find("MiSeq") >= 0:
        rc_i2 = False
    # according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-03.pdf
    # NovaSeq follows MiSeq/HiSeq 2500 workflow for I5 index; this is effectively
    # a noop thanks to rc_i2 being False by default but let's make it explicit
    elif application.find("NovaSeq") >= 0:
        rc_i2 = False
    elif application.find("HiSeq") >= 0:
        # Hiseq 4000 has version 3 -- hopefully this is stable??
        app_str = re.search(r'[\d\.]+', application_version).group()
        main_app_ver = int(app_str.split(".")[0])
        if main_app_ver > 2:
            rc_i2 = True
        else:
            rc_i2 = False

    # Can we run bcl2fastq 2, or do we need to use the old version?
    if setup_node is not None:
        rta_tag = tree.getroot().find("Setup").find("RTAVersion")
        if rta_tag is None:
            rta_tag = tree.getroot().find("RTAVersion")
    # TENKIT-60 NovaSeq lowercases the tag
    else:
        rta_tag = tree.getroot().find("RtaVersion")

    rta_version = rta_tag.text
    if rta_version.startswith('v'):
        rta_version = rta_version[1:]

    # Also attempt to extract
    #    <ApplicationName>HiSeq Control Software</ApplicationName>
    #    <ApplicationVersion>3.3.20</ApplicationVersion>
    if setup_node is not None:
        application_name = setup_node.find("ApplicationName")
        if application_name is None:
            application_name = tree.getroot().find("ApplicationName")

        application_version = setup_node.find("ApplicationVersion")
        if application_version is None:
            application_version = tree.getroot().find("ApplicationVersion")
    else:
        application_name = tree.getroot().find("Application")
        application_version = tree.getroot().find("ApplicationVersion")

    params = { 'ApplicationName': application_name.text, 'ApplicationVersion': application_version.text }
    return (rta_version, rc_i2, params)


def get_sequencer_type(input_path):
    """
    Returns the sequencer type from runParameters.xml in the input path.

    TODO: Perhaps look at the flowcell ID instead (see Preyas'
    Illumina identification code to see if that's more robust)
    """
    _, _, params = get_rta_version(input_path)
    if 'ApplicationName' in params:
        return params['ApplicationName'].split()[0]
    else:
        return None


def load_run_info(run_info_xml):
    """
    Get the read names and read lengths from the Illumina RunInfo.xml file
    """
    tree = etree.parse(run_info_xml)
    reads_node = tree.getroot().find("Run").find("Reads")
    reads = reads_node.findall("Read")

    read_lengths = [ int(read.attrib["NumCycles"]) for read in reads ]
    read_info = [ { "read_length": x, "index_read": False } for x in read_lengths ]

    # Now we follow some hard-coded conventions on order in which reads appear.
    # R1 is always first
    # R2 is always last
    # Index reads (if they exist) are in the middle, in numerical order

    # NOTE -- if there is only one index read it is assumed to be I1.
    # BclToFastq should give us this
    # NOTE -- we assume paired end reads!
    read_info[0]["read_name"] = "R1"

    if len(read_info) > 1:
        read_info[-1]["read_name"] = "R2"

    index_read_count = 1
    for idx in range(1, len(read_info) - 1):
        index_read_name = "I" + str(index_read_count)
        read_info[idx]["read_name"] = index_read_name
        read_info[idx]["index_read"] = True
        index_read_count += 1

    # The original fastq files are labelled R1, R2, R3
    for idx in range(len(read_info)):
        read_info[idx]["original_read_name"] = "R" + str(idx + 1)

    flowcell = tree.getroot().find("Run").find("Flowcell").text

    # NB: currently you have to comment out the next two lines to get
    # nosetests to run correctly outside of a stage.
    martian.log_info("Read Info: %s" % read_info)
    martian.log_info("Flowcell ID: %s" % flowcell)
    return (read_info, flowcell)


def is_real_dual_index_flowcell(read_info):
    """
    Return whether we think this was really a dual-indexed flowcell.  The
    criteria is that there are two index reads with read length < 14.
    Anything above that, and we think it's a GemCode or SingleCell V1
    run with a barcode in one of the indices.

    :param read_info:
    :return:
    """
    real_index_count = 0
    for read in read_info:
        if read["index_read"] and read["read_length"] < tk_constants.DEMULTIPLEX_BARCODE_LENGTH:
            real_index_count += 1

    return real_index_count == 2


def make_bases_mask_val(read_info, sample_index_read=None, dual_indexed=False, ignore_dual_index=False):
    """
    :param read_info: The ReadInfo block from RunInfo.xml
    :param sample_index_read: The ReadInfo read (I1, I2) to use as the sample index
    :param dual_indexed: If the input BCLs were dual-indexed intentionally, then
                         preserve the index status of the 2nd, non-sample index index
                         in the mask.  Too often, people inadvertently ran dual-indexed
                         values for the barcode read (Cell Ranger v1, GemCode), so the
                         default behavior is to treat the off-SI indexed read as a
                         normal, non-index read.
    :param ignore_dual_index: Stub out the dual index with Ns, effectively allowing all values of
                              a dual index to pass into the corresponding reads.  Use this for
                              dual-indexed flowcells where a single sample is used on a single lane.
    """
    # We will emit all reads in RunInfo.xml
    def base_mask(read):
        if read["read_name"][0] == "R":
            return "Y" + str(read["read_length"])
        elif read["read_name"][0] == "I":
            if ignore_dual_index and read["read_name"] != sample_index_read:
                return "N" + str(read["read_length"])
            elif dual_indexed or read["read_name"] == sample_index_read:
                return "I" + str(read["read_length"])
            else:
                return "Y" + str(read["read_length"])
        else:
            martian.throw("read name was not recognized: %s" % read["read_name"])

    masks = [ base_mask(r) for r in read_info ]
    return ",".join(masks)


def get_bcl2fastq_read_type_map(read_info, sample_index_read=None, dual_indexed=False, ignore_dual_index=False):
    """
    Get a mapping between ReadInfo read name (R1,I1,I2,R2) and bcl2fastq
    output file naming (R1/R2/R3/R4/I1)

    The guarantee here is that the 10X sample index will always be on I1,
    if generated.  If dual-indexing is specified, the secondary index will
    be on I2.  Upstream pipestances can expect read1 to be in R1, sample
    indexes to be on I1, and read2 to be on R2.

    :param read_info: The ReadInfo block from RunInfo.xml
    :param sample_index_read: The ReadInfo read (I1, I2) to use as the sample index
    :param ignore_dual_index: Whether the dual index was ignored (and thus not sequenced)
    """
    #read_names = [r["read_name"] for r in read_info]
    read_map = {}
    reads_counted = 0
    for idx, r in enumerate(read_info):
        read_name = r["read_name"]
        if read_name == sample_index_read:
            read_map[read_name] = 'I1'
        elif ignore_dual_index and r["index_read"]:
            # we didn't read this index -- see make_bases_mask_val
            continue
        elif dual_indexed and r["index_read"]:
            read_map[read_name] = 'I2'
        else:
            reads_counted += 1
            read_map[read_name] = 'R%d' % reads_counted
    return read_map


if __name__ == '__main__':
    v = get_bcl2fastq_v2("host")
    print v
    check_bcl2fastq("host", "2.3.4")
