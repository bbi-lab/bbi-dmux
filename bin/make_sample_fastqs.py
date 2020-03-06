#!/usr/bin/env python
# Make sample fastqs 2-level
# Andrew's barcode parser/fastq combiner modified for 2-level

import barcodeutils as bu
import argparse
import os
import json
from collections import OrderedDict
import re
import sys
import glob
import xml.etree.ElementTree as ET

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
P5_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p5.txt')
P7_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p7.txt')
RT_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/rt2.txt')
LIG_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/ligation.txt')
RT3_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/rt.txt')

def get_programmed_pcr_combos(p5_lookup, p7_lookup, p5_cols_used, p7_rows_used):
    """
    Assuming p5 are columns and p7 are rows, get the set of (p5, p7) wells that were programmed.
    These are the set that would have been physically combined.
    Args:
        p5_lookup (dict): p5_lookup dict mapping sequences to wells as passed to barcode specification
        p7_lookup (dict): p7_lookup dict mapping sequences to wells as passed to barcode specification
        p5_cols_used (list): A list of the cols used from P5 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. A B C for p7 and 1 2 3 for p5.
        p7_rows_used: A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. A B C for p7 and 1 2 3 for p5.
    Returns:
        set of (p5, p7): set of (p5, p7) well ID tuples that are valid.
    """

    valid_combos = set()

    for p5, p7 in zip(p5_cols_used, p7_rows_used):

        selected_p5 = [p5_well for p5_well in p5_lookup.values() if int(p5_well[1:]) == p5]
        selected_p7 = [p7_well for p7_well in p7_lookup.values() if p7_well[0] == p7]

        for selected_p5_i in selected_p5:
            for selected_p7_i in selected_p7:
                valid_combos.add((selected_p5_i, selected_p7_i))

    return valid_combos

def get_programmed_pcr_combos_wells(p5_wells_used, p7_wells_used):
    """
    Assuming p5 and p7 are wells, get the set of (p5, p7) wells that were programmed.
    These are the set that would have been physically combined.
    Args:
        p5_wells_used (list): A list of the wells used from P5 plate for PCR in same order as P7 to indicate the pairs of P7 and P5 used (e.g. A1 B1 C1 for p7 and C1 D2 E3 for p5.
        p7_wells_used: A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. A1 B1 C1 for p7 and C1 D2 E3 for p5.
    Returns:
        set of (p5, p7): set of (p5, p7) well ID tuples that are valid.
    """
    good_nums = {1:"01", 2:"02", 3:"03", 4:"04", 5:"05", 6:"06", 7:"07", 8:"08", 9:"09", 10:"10", 11:"11", 12:"12"}

    valid_combos = set()
    p5_wells_fixed = [p5_well[0] + good_nums[int(p5_well[1:])] for p5_well in p5_wells_used]
    p7_wells_fixed = [p7_well[0] + good_nums[int(p7_well[1:])] for p7_well in p7_wells_used]
    for selected_p5, selected_p7 in zip(p5_wells_fixed, p7_wells_fixed):
        valid_combos.add((selected_p5, selected_p7))

    return valid_combos


def quick_parse(file_path):
    # a copy of only the relevant lines from easygrid.read_delim
    fh = open(file_path)
    columns = next(fh).strip().split(",")

    # Parse file
    for line_number, line in enumerate(fh):
        entries = line.strip().split(",")

        if len(entries) != len(columns):
            raise ValueError('Length of entries: %s, not equal to length of columns %s, in file: %s, line number %s' % (len(entries), len(columns), file_path, line_number))

        entries_dict = dict(zip(columns, entries))
        yield entries_dict

def load_sample_layout(file_path):
    """
    Function that loads the sample layout file to an RT p5_lookup table.
    """

    lookup = {}
    for rt_well in quick_parse(file_path):
        lookup[rt_well['RT Barcode']] = rt_well['Sample ID'].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')
    return lookup

NEXTSEQ = 'NextSeq'
MISEQ = 'MiSeq'
NOVASEQ = 'NovaSeq'
HISEQ4000 = 'HiSeq4000'
HISEQ3000 = 'HiSeq3000'
HISEQ = 'HiSeq'
UNKNOWN_SEQUENCER = 'unknown'

SEQUENCERS_P5_RC_MAP = {
    NEXTSEQ: True,
    MISEQ: False,
    NOVASEQ: False,
    HISEQ4000: True,
    HISEQ3000: False
}

# Taken from barcodeutils, but excluding what's not needed and not always found
def reverse_complement_i5(name):
    """
    Take a BCL directory or instrument type (NextSeq, MiSeq, NovaSeq, HiSeq4000, HiSeq3000) and return whether or not i5 should be reverse complemented.
    This assumes that NextSeq instruments and other similar machines should be reverse complemeted whereas MiSeq should not.
    Args:
        name (str): BCL directory or one of the instrument types as mentioned above    
    
    Returns:
        bool: True if user would typically reverse complement i5 index and False otherwise.
    """
    
    if name in SEQUENCERS_P5_RC_MAP:
        sequencer_type = name
    elif os.path.exists(name):
        sequencer_type = get_run_info(name)['instrument_type']
        
        if sequencer_type not in SEQUENCERS_P5_RC_MAP:
            raise ValueError('Sequencer type detected from BCL is %s, which is not in our known list of which sequencers require P5 reverse complementing or not.' % sequencer_type)
    else:
        raise ValueError('Invalid input, could not detect BCL or instrument ID.')

    return SEQUENCERS_P5_RC_MAP[sequencer_type]

# Taken from barcodeutils, but excluding what's not needed and not always found
def get_run_info(flow_cell_path):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        flow_cell_path (str): Path to BCL directory for run.
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    bcl_run_info = os.path.join(flow_cell_path, 'RunParameters.xml*')
    bcl_run_info = glob.glob(bcl_run_info)
    if not bcl_run_info:
        raise ValueError('BCL RunParameters.xml not found for specified flowcell: %s' % bcl_run_info)
    else:
        bcl_run_info = bcl_run_info[0]

    # Set up a few nodes for parsing
    tree = ET.parse(bu.open_file(bcl_run_info))

    setup_node = tree.getroot().find("Setup")
    if setup_node is None:
        setup_node = tree.getroot()

    # Figure out instrument
    application = setup_node.find('ApplicationName')
    if application is None:
        application = setup_node.find('Application')

    application = application.text
    application_version = setup_node.find('ApplicationVersion')
    if NEXTSEQ in application:
        run_stats['instrument_type'] = NEXTSEQ
    elif MISEQ in application:
        run_stats['instrument_type'] = MISEQ
    elif NOVASEQ in application:
        run_stats['instrument_type'] = NOVASEQ
    elif HISEQ in application:
        app_string = re.search(r'[\d\.]+', application_version).group()
        app_major_version = int(app_string.split('.')[0])

        if app_major_version > 2:
            run_stats['instrument_type'] = HISEQ4000
        else:
            run_stats['instrument_type'] = HISEQ3000
    else:
        run_stats['instrument_type'] = UNKNOWN_SEQUENCER

    run_start_date_node = tree.getroot().find('RunStartDate')

    # Now actually populate various stats
    run_stats['date'] = run_start_date_node.text

    if run_stats['instrument_type'] == NOVASEQ:
        run_stats['p7_index_length'] = int(setup_node.find('PlannedIndex1ReadCycles').text)
        run_stats['p5_index_length'] = int(setup_node.find('PlannedIndex2ReadCycles').text)
    else:
        run_stats['p7_index_length'] = int(setup_node.find('Index1Read').text)
        run_stats['p5_index_length'] = int(setup_node.find('Index2Read').text)

    return run_stats

def write_to_undetermined(entry):
    output_name = entry['r1_name'] + "|" + entry['r1_seq']
    r2_seq = entry['r2_seq']
    r2_qual = entry['r2_qual']
    output_line = f'{output_name}\n{r2_seq}\n+\n{r2_qual}\n'
    undetermined.write(output_line)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to combine R1 and R2 file from sci run into a single file with barcodes in the fastq header.')
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    parser.add_argument('--read1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
    parser.add_argument('--read2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
    parser.add_argument('--file_name', required=True, help='The R1 file name.')
    parser.add_argument('--sample_layout', required=True, help='Text file containing the sample layout by RT well.')
    parser.add_argument('--p5_cols_used', nargs='+', type=int, required=True, help='A list of the columns used from P5 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p7_rows_used', nargs='+', required=True, help='A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p5_wells_used', nargs='+', required=True, help='A list of the wells used from P5 plate for PCR in same order as P7 to indicate the pairs of P7 and P5 used (e.g. --p7 A1 B1 C1 for p7 and --p5 A1 A2 A3 for p5. Alternative to p5_cols_used.')
    parser.add_argument('--p7_wells_used', nargs='+', required=True, help='A list of the wells used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A1 B1 C1 for p7 and --p5 A1 A2 A3 for p5. Alternative to p7_rows_used.')
    parser.add_argument('--output_dir', required=True, help='Output directory for files.')
    parser.add_argument('--p7_length', type=int, default=10, help='Expected P7 index length.')
    parser.add_argument('--p5_length', type=int, default=10, help='Expected P5 index length.')
    parser.add_argument('--level', required=True, help = "2 or 3 level sci?")
    parser.add_argument('--rt_barcode_file', required=True, help='Path to RT barcode file, or "default".')
    args = parser.parse_args()

    if args.rt_barcode_file == "default":
        if args.level == "3":
            rtfile = RT3_FILE
        else:
            rtfile = RT_FILE
    else:
        rtfile = args.rt_barcode_file
    lane_num = args.file_name
    lane_num = lane_num.replace("Undetermined_S0_L", "L")
    lane_num = lane_num.replace("_R1_001.fastq.gz", "")
    stats_file = os.path.join(args.output_dir, lane_num + ".stats.json")
    suffix = lane_num + ".fastq"

    reverse_complement_i5 = reverse_complement_i5(args.run_directory)

    if args.level == "3":
        ligation_lookup = bu.load_whitelist(LIG_FILE, variable_lengths=True)
        ligation_9_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 9}
        ligation_10_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 10}

    # Load barcodes
    p7_lookup = bu.load_whitelist(P7_FILE)
    p7_lookup = {sequence[0:args.p7_length]: well for sequence,well in p7_lookup.items()}

    p5_lookup = bu.load_whitelist(P5_FILE)
    if reverse_complement_i5:
        p5_lookup = {bu.reverse_complement(sequence): well for sequence,well in p5_lookup.items()}
    p5_lookup = {sequence[0:args.p5_length]: well for sequence,well in p5_lookup.items()}

    # Get the set of all valid PCR combos
    # Not very robust - to be improved
    if args.p5_cols_used != [0]:
        programmed_pcr_combos = get_programmed_pcr_combos(p5_lookup, p7_lookup, args.p5_cols_used, args.p7_rows_used)
    else:
        programmed_pcr_combos = get_programmed_pcr_combos_wells(args.p5_wells_used, args.p7_wells_used)
    # Define where all sequences are and what the whitelists are

    if args.level == "3":
        barcode_spec = {
            'p5': {
                'start': 1,
                'end':  args.p5_length,
                'read': 'i5',
                'whitelist': p5_lookup
            },
            'p7': {
                'start': 1,
                'end': args.p7_length,
                'read': 'i7',
                'whitelist': p7_lookup
            },
            'ligation_9': {
                'start': 1,
                'end': 9,
                'read': 'r1',
                'whitelist': ligation_9_lookup
            },
            'ligation_10': {
                'start': 1,
                'end': 10,
                'read': 'r1',
                'whitelist': ligation_10_lookup
            },
            'umi_9': {
                'start': 16,
                'end': 23,
                'read': 'r1'
            },
            'umi_10': {
                'start': 17,
                'end': 24,
                'read': 'r1'
            },
            'rt_9': {
                'start': 24,
                'end': 33,
                'read': 'r1',
                'whitelist': rtfile
            },
            'rt_10': {
                'start': 25,
                'end': 34,
                'read': 'r1',
                'whitelist': rtfile
            }  
        }
    else:
        barcode_spec = {
            'p5': {
                'start': 1,
                'end':  args.p5_length,
                'read': 'i5',
                'whitelist': p5_lookup
            },
            'p7': {
                'start': 1,
                'end': args.p7_length,
                'read': 'i7',
                'whitelist': p7_lookup
            },
            'umi': {
                'start': 1,
                'end': 8,
                'read': 'r1'
            },
            'rt': {
                'start': 9,
                'end': 18,
                'read': 'r1',
                'whitelist': rtfile
            }  
        }

  # Set up the output files
    sample_rt_lookup = load_sample_layout(args.sample_layout)
    undetermined = open(os.path.join(args.output_dir, '%s-%s' % ("Undetermined", suffix)), 'w')
    sample_to_output_filename_lookup = {sample: os.path.join(args.output_dir, '%s-%s' % (sample, suffix)) for well,sample in sample_rt_lookup.items()}
    sample_to_output_file_lookup = {sample: open(filename, 'w') for sample,filename in sample_to_output_filename_lookup.items()}
    print("Demuxing %s samples (%s total RT wells) into their own files..." % (len(sample_to_output_filename_lookup), len(sample_rt_lookup)))
     
    # Set up some basic tracking for each sample
    sample_read_counts = {}
    for sample in sample_to_output_file_lookup:
        sample_read_counts[sample] = 0

    if args.level == "3":
        total_reads = 0
        total_uncorrected = 0
        total_pcr_mismatch = 0
        total_ambiguous_ligation_length = 0
        total_unused_rt_well = 0
        total_corrected_9 = 0
        total_corrected_10 = 0
        read_pair_dict = {}
        rt_dict = {}
        lig_dict = {}        
        # Finally, process reads
        for read_number, entry in enumerate(bu.parse_fastq_barcodes(args.read1, args.read2, spec=barcode_spec, edit_distance=1)):

            total_reads += 1

            # Only allow the programmed PCR combos (helps clean things up a bit)
            p5 = entry['p5']
            p7 = entry['p7']
    
            ## Choose between the _8 and _9 options based on which correct properly
            corrected_9 = entry['ligation_9'] is not None and entry['rt_9'] is not None
            corrected_10 = entry['ligation_10'] is not None and entry['rt_10'] is not None
            corrected_p5_p7 = p5 is not None and p7 is not None

            if corrected_9 and corrected_10:
                total_ambiguous_ligation_length += 1
                write_to_undetermined(entry)
                continue

            if corrected_9 and corrected_p5_p7:
                total_corrected_9 += 1
                ligation_barcode = entry['ligation_9']
                rt_barcode = entry['rt_9']
                umi = entry['umi_9']
            elif corrected_10 and corrected_p5_p7:
                total_corrected_10 += 1
                ligation_barcode = entry['ligation_10']
                rt_barcode = entry['rt_10']
                umi = entry['umi_10']
            else:
                total_uncorrected += 1
                write_to_undetermined(entry)
                continue

            if not (p5, p7) in programmed_pcr_combos:
                total_pcr_mismatch += 1
                write_to_undetermined(entry)
                continue

            if rt_barcode not in sample_rt_lookup:
                total_unused_rt_well += 1
                write_to_undetermined(entry)
                continue

            sample = sample_rt_lookup[rt_barcode]
            sample_read_number = sample_read_counts[sample] + 1
            sample_read_counts[sample] += 1

            if (p5, p7) in read_pair_dict:
                read_pair_dict[(p5, p7)] += 1
            else:
                read_pair_dict[(p5, p7)] = 1

            if ligation_barcode in lig_dict:
                lig_dict[ligation_barcode] += 1
            else:
                lig_dict[ligation_barcode] = 1

            if rt_barcode in rt_dict:
                rt_dict[rt_barcode] += 1
            else:
                rt_dict[rt_barcode] = 1

            r2_qual = entry['r2_qual']
            r2_seq = entry['r2_seq']
            output_name = f'@{sample}-P5{p5}-P7{p7}_{sample_read_number}|{sample}|{p5}|{p7}|{rt_barcode}_{ligation_barcode}|{umi}'
            output_line = f'{output_name}\n{r2_seq}\n+\n{r2_qual}\n'
            sample_to_output_file_lookup[sample].write(output_line)

        # Close output files
        for f in sample_to_output_file_lookup.values():
            f.close()
        undetermined.close()

        # Output stats
        total_passed_reads = sum(list(sample_read_counts.values()))
        stats = OrderedDict()
        stats['total_input_reads'] = total_reads
        stats['total_passed_reads'] = total_passed_reads
        stats['fraction_passed_reads'] = total_passed_reads / total_reads
        stats['fraction_uncorrected_reads'] = total_uncorrected / total_reads
        stats['fraction_ambiguous_ligation_length'] = total_ambiguous_ligation_length / total_reads
        stats['fraction_invalid_rt_well'] = total_unused_rt_well / total_reads
        stats['fraction_pcr_mismatch'] = total_pcr_mismatch / total_reads
        stats['total_reads_corrected_when_9bp_ligation'] = total_corrected_9
        stats['total_reads_corrected_when_10bp_ligation'] = total_corrected_10
        stats['total_reads_passed_per_sample'] = sample_read_counts
        lig_dict_file = os.path.join(args.output_dir, lane_num + ".lig_counts.csv")
        # Output read_pair_dict
        with open(lig_dict_file, 'w') as f:
            for lig,val in lig_dict.items():
                if lig is not None:
                    f.write(lig + "," + str(val) + "\n")

    else:
        total_reads = 0
        total_uncorrected = 0
        total_pcr_mismatch = 0
        total_unused_rt_well = 0
        total_corrected = 0
        read_pair_dict = {}
        rt_dict = {}
        # Finally, process reads
        for read_number, entry in enumerate(bu.parse_fastq_barcodes(args.read1, args.read2, spec=barcode_spec, edit_distance=1)):

            total_reads += 1

            # Only allow the programmed PCR combos (helps clean things up a bit)
            p5 = entry['p5']
            p7 = entry['p7']

            if (p5, p7) in read_pair_dict:
                read_pair_dict[(p5, p7)] += 1
            else:
                read_pair_dict[(p5, p7)] = 1

            corrected = entry['rt'] is not None
            corrected_p5_p7 = p5 is not None and p7 is not None

            if corrected and corrected_p5_p7:
                total_corrected += 1
                rt_barcode = entry['rt']
                umi = entry['umi']
            else:
                total_uncorrected += 1
                write_to_undetermined(entry)
                continue

            if not (p5, p7) in programmed_pcr_combos:
                total_pcr_mismatch += 1
                write_to_undetermined(entry)
                continue

            if rt_barcode not in sample_rt_lookup:
                total_unused_rt_well += 1
                write_to_undetermined(entry)
                continue

            sample = sample_rt_lookup[rt_barcode]
            sample_read_number = sample_read_counts[sample] + 1
            sample_read_counts[sample] += 1

            if rt_barcode in rt_dict:
                rt_dict[rt_barcode] += 1
            else:
                rt_dict[rt_barcode] = 1

            r2_qual = entry['r2_qual']
            r2_seq = entry['r2_seq']
            output_name = f'@{sample}-P7{p5}-P5{p7}_{sample_read_number}|{sample}|{p5}|{p7}|{rt_barcode}|{umi}'
            output_line = f'{output_name}\n{r2_seq}\n+\n{r2_qual}\n'
            sample_to_output_file_lookup[sample].write(output_line)

        # Close output files
        for f in sample_to_output_file_lookup.values():
            f.close()
        undetermined.close()

        # Output stats
        total_passed_reads = sum(list(sample_read_counts.values()))
        stats = OrderedDict()
        stats['total_input_reads'] = total_reads
        stats['total_passed_reads'] = total_passed_reads
        stats['fraction_passed_reads'] = total_passed_reads / total_reads
        stats['fraction_uncorrected_reads'] = total_uncorrected / total_reads
        stats['fraction_invalid_rt_well'] = total_unused_rt_well / total_reads
        stats['fraction_pcr_mismatch'] = total_pcr_mismatch / total_reads
        stats['total_reads_corrected'] = total_corrected
        stats['total_reads_passed_per_sample'] = sample_read_counts
 
    dict_file = os.path.join(args.output_dir, lane_num + ".pcr_counts.csv")
    # Output read_pair_dict
    with open(dict_file, 'w') as f:
        for (p5, p7),val in read_pair_dict.items():
            if p5 is not None and p7 is not None:
                f.write(p5 + "," + p7 + "," + str(val) + "\n")
    
    rt_dict_file = os.path.join(args.output_dir, lane_num + ".rt_counts.csv")
    # Output rt_dict
    with open(rt_dict_file, 'w') as f:
        for rt,val in rt_dict.items():
            if rt is not None:
                f.write(rt + "," + str(val) + "\n")    
    with open(stats_file, 'w') as f:
        f.write(json.dumps(stats, indent=4))

