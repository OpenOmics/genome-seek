#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

from __future__ import print_function, division
import sys, gzip


# USAGE
# sys.argv[1] = sample_name.R1.fastq.gz
# sys.argv[2] = sample_name (name without PATH and .R?.fastq.gz extension)
# Example
# $ python3 flowcell_lane.py input.R1.fastq.gz input > flowcell_lanes.txt

# Input 1 (Normal FastQ from Casava > 1.8)
# @J00170:88:ANYVJBBXX:8:1101:1600:1244 1:N:0:ACTTGA
# GGGAAGTTGAAAGCTTCCAGTGCTCCCTGTCAATTCTAGTCCCTCCAGTCT
# +
# AAAFFJJFJJJJJJFJJJJJJJJJJFJAJJJJJFJJJJJFFJJAJJJJ7JJ

# Input 2 (SRA doesn't store FC ID, use intrument name instead)
# @SRR5351039.1 SN608:8:1101:31.20:96.50 length=51
# NTTTANNNNNNGNGCNCTGNNNNNNNNGNNNNNAAGGGNTNNNNNNNNNNN
# +SRR5351039.1 SN608:8:1101:31.20:96.50 length=51
# !0;??!!!!!!3!22!2:=!!!!!!!!1!!!!!00<=?!-!!!!!!!!!!!

# Input 3 (SRA download with no FC ID, instrument name, or lanes information)
# @SRR6755966.1 1 length=101
# GCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCA
# +SRR6755966.1 1 length=101
# CC@FFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJHIJJJJI

def usage(message = '', exitcode = 0):
    """Displays help and usage information. If provided invalid usage
    returns non-zero exit-code. Additional message can be displayed with
    the 'message' parameter.
    """
    print('Usage: python3 {} sampleName.R1.fastq.gz  sampleName > sampleName.flowcell_lanes.txt'.format(sys.argv[0]))
    if message:
        print(message)
    sys.exit(exitcode)


def reader(fname):
    """Returns correct file object handler or reader for gzipped
    or non-gzipped FastQ files based on the file extension. Assumes
    gzipped files endwith the '.gz' extension.
    """
    if fname.endswith('.gz'):
        # Opens up file with gzip handler
        return gzip.open
    else:
        # Opens up file normal, uncompressed handler
        return open


def get_flowcell_lane(sequence_identifer):
    """Returns flowcell and lane information for different fastq formats.
    FastQ files generated with older versions of Casava or downloaded from
    SRA have a different format than newer FastQ files generated with the
    current version of Casava. It is worth noting that FastQ files downloaded from SRA
    or FastQ files generated with Casava version < 1.8 do not have Flowcell
    IDs in its sequence indentifer.
    For more information visit: https://en.wikipedia.org/wiki/FASTQ_format
    """
    try:
        # Decode gzip byte-string object,
        # needed for gzipped input files
        id_list = sequence_identifer.decode('utf8').strip().split(':')
    except AttributeError:
        # Needed for non-gzipped inputs
        id_list = sequence_identifer.strip().split(':')
    if len(id_list) < 7:
        # No Flowcell IDs in this format
        # Return next instrument id instead (next best thing)
        if sequence_identifer.startswith('@SRR'):
            # SRA format or downloaded SRA FastQ file
            # SRA format 1: contains machine and lane information
            # @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
            # SRA format 2: contains nothing, grab SRR ID
            # @SRR6755966.1 1 length=101
            try:
                # SRA format 1
                id1 = id_list[0].split()[1]
                id2 = id_list[1]
            except IndexError:
                # SRA format 2
                id1 = id_list[0].split()[0].split(".")[0]
                id2 = id1.lstrip('@')
            return id1,id2
        else:
            # Casava < 1.8 (fastq format)
            # @HWUSI-EAS100R:6:73:941:1973#0/1
            return id_list[0],id_list[1]
    else:
        # Casava >= 1.8
        # Normal FastQ format
        # @J00170:88:HNYVJBBXX:8:1101:6390:1244 1:N:0:ACTTGA
        return id_list[2],id_list[3]


def md5sum(filename, blocksize = 65536):
    """Gets md5checksum of a file in memory-safe manner.
    The file is read in blocks defined by the blocksize parameter. This is a safer
    option to reading the entire file into memory if the file is very large.
    @param filename <str>:
        Input file on local filesystem to find md5 checksum
    @param blocksize <int>:
        Blocksize of reading N chunks of data to reduce memory profile
    @return hasher.hexdigest() <str>:
        MD5 checksum of the file's contents
    """
    import hashlib

    hasher = hashlib.md5()
    with open(filename, 'rb') as fh:
        buf = fh.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = fh.read(blocksize)

    return hasher.hexdigest()


if __name__ == '__main__':

    # Check Usage
    if '-h' in sys.argv or '--help' in sys.argv or '-help' in sys.argv:
        usage(exitcode = 0)
    elif len(sys.argv) != 3:
        usage(message = 'Error: failed to provide all required positional arguments!', exitcode = 1)

    # Get file name and sample name prefix
    filename = sys.argv[1]
    sample = sys.argv[2]
    # Get md5 checksum
    md5 = md5sum(filename)

    # Get Flowcell and Lane information
    handle = reader(filename)
    meta = {'flowcell': [], 'lane': [], 'flowcell_lane': []}
    i = 0  # keeps track of line number
    with handle(filename, 'r') as file:
        print('sample_name\ttotal_read_pairs\tflowcell_ids\tlanes\tflowcell_lanes\tmd5_checksum')
        for line in file:
            line = line.strip()
            if i%4 == 0: # read id or sequence identifer
                fc, lane = get_flowcell_lane(line)
                fc = fc.lstrip('@')
                fc_lane = "{}_{}".format(fc,lane)
                if fc not in meta['flowcell']:
                    meta['flowcell'].append(fc)
                if lane not in meta['lane']:
                    meta['lane'].append(lane)
                if fc_lane not in meta['flowcell_lane']:
                    meta['flowcell_lane'].append(fc_lane)
            i += 1

    print("{}\t{}\t{}\t{}\t{}\t{}".format(sample, int(i/4),",".join(sorted(meta['flowcell'])),",".join(sorted(meta['lane'])),",".join(sorted(meta['flowcell_lane'])), md5))
