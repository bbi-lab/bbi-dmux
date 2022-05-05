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
import ast
import run_info
import gzip
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

    if p5_cols_used == ["none"]:
        p5_cols_used = p5_cols_used * len(p7_rows_used)
    if p7_rows_used == ["none"]:
        p7_rows_used = p7_rows_used * len(p5_cols_used)

    valid_combos = set()
    for p5, p7 in zip(p5_cols_used, p7_rows_used):
 
        if p7 == "none":
            selected_p7 = ["none"]
        else:
            selected_p7 = [p7_well for p7_well in p7_lookup.values() if p7_well[0] == p7 or p7_well == p7]

        if p5 == "none":
            selected_p5 = ["none"]
        else:   
            selected_p5 = [p5_well for p5_well in p5_lookup.values() if int(p5_well[1:]) == p5 or p5_well == p5]

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
    if p5_wells_used == ["none"]:
        p5_wells_fixed = p5_wells_used * len(p7_wells_used)
    else:
        p5_wells_fixed = [p5_well[0] + good_nums[int(p5_well[1:])] for p5_well in p5_wells_used]
    if p7_wells_used == ["none"]:
        p7_wells_fixed = p7_wells_used * len(p5_wells_used)
    else:
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

def load_sample_layout(file_path, multi_exp):
    """
    Function that loads the sample layout file to an RT p5_lookup table.
    """
    if multi_exp:
        lookup = {}
        for rt_well in quick_parse(file_path):
            lookup[(rt_well['RT Barcode'], rt_well['Experiment'])] = rt_well['Sample ID'].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')

    else:
        lookup = {}
        for rt_well in quick_parse(file_path):
            lookup[rt_well['RT Barcode']] = rt_well['Sample ID'].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')
    return lookup


def write_to_undetermined(entry):
    output_name = entry['r1_name'] + "|" + entry['r1_seq']
    r2_seq = entry['r2_seq']
    r2_qual = entry['r2_qual']
    output_line = f'@{output_name}\n{r2_seq}\n+\n{r2_qual}\n'
    undetermined.write(output_line)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to combine R1 and R2 file from sci run into a single file with barcodes in the fastq header.')
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    parser.add_argument('--read1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
    parser.add_argument('--read2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
    parser.add_argument('--file_name', required=True, help='The R1 file name.')
    parser.add_argument('--sample_layout', required=True, help='Text file containing the sample layout by RT well.')
    parser.add_argument('--p5_cols_used', nargs='+', required=True, help='A list of the columns used from P5 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p7_rows_used', nargs='+', required=True, help='A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p5_wells_used', nargs='+', required=True, help='A list of the wells used from P5 plate for PCR in same order as P7 to indicate the pairs of P7 and P5 used (e.g. --p7 A1 B1 C1 for p7 and --p5 A1 A2 A3 for p5. Alternative to p5_cols_used.')
    parser.add_argument('--p7_wells_used', nargs='+', required=True, help='A list of the wells used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A1 B1 C1 for p7 and --p5 A1 A2 A3 for p5. Alternative to p7_rows_used.')
    parser.add_argument('--output_dir', required=True, help='Output directory for files.')
    parser.add_argument('--p7_length', type=int, default=10, help='Expected P7 index length.')
    parser.add_argument('--p5_length', type=int, default=10, help='Expected P5 index length.')
    parser.add_argument('--multi_exp', type=ast.literal_eval)
    parser.add_argument('--level', required=True, help = "2 or 3 level sci?")
    parser.add_argument('--rt_barcode_file', required=True, help='Path to RT barcode file, or "default".')
    parser.add_argument('--p5_barcode_file', required=True, help='Path to p5 barcode file, or "default".')
    parser.add_argument('--p7_barcode_file', required=True, help='Path to p7 barcode file, or "default".')
    parser.add_argument('--lig_barcode_file', required=True, help='Path to ligation barcode file, or "default".')
    parser.add_argument('--buffer_blocks', required=True, help="Number of 8k buffer blocks to use for writing")
    args = parser.parse_args()


    buffer_size = int(args.buffer_blocks) * 8192
    run_info = run_info.get_run_info( args.run_directory, pipeline_type='RNA-seq' )

    if( run_info['paired_end'] == False ):
        raise ValueError('Single-end reads detected: paired-end reads required')

    reverse_complement_i5 = run_info['reverse_complement_i5']

    # Load and process PCR barcodes
    p7_lookup = bu.load_whitelist(P7_FILE)
    p7_lookup = {sequence[0:args.p7_length]: well for sequence,well in p7_lookup.items()}

    p5_lookup = bu.load_whitelist(P5_FILE)
    if reverse_complement_i5:
        p5_lookup = {bu.reverse_complement(sequence): well for sequence,well in p5_lookup.items()}
    p5_lookup = {sequence[0:args.p5_length]: well for sequence,well in p5_lookup.items()}

    multi_exp = True if ( args.multi_exp != 0 ) else False
    if multi_exp:
        all_p7 = ""
        all_p5 = ""
        exp_lookup = {}
        for key,value in args.multi_exp.items():
            all_p7 += value[0].strip() + " "
            all_p5 += value[1].strip() + " "
            if len(value[0].strip().split(" ")[0]) == 1:
                p5s = [int(x) for x in value[1].strip().split(" ")]
                combos = get_programmed_pcr_combos(p5_lookup, p7_lookup, p5s, value[0].strip().split(" "))
            else:
                combos = get_programmed_pcr_combos_wells(value[1].strip().split(" "), value[0].strip().split(" "))
            temp_dict = dict(zip(combos, [key]*len(combos)))
            exp_lookup.update(temp_dict)
        
        all_p7 = all_p7.strip().split(" ")
        all_p5 = all_p5.strip().split(" ")

        if len(all_p7[0]) == 1:
            args.p7_rows_used = all_p7
            args.p5_cols_used = all_p5
        else:
            args.p7_wells_used = all_p7
            args.p5_wells_used = all_p5

    if args.p5_cols_used == ["none"] or args.p5_wells_used == ["none"]:
        p5_none = True
    else:
        p5_none = False

    if args.p7_rows_used == ["none"] or args.p7_wells_used == ["none"]:
        p7_none = True
    else:
        p7_none = False

    if not p5_none:
        args.p5_cols_used = [int(x) for x in args.p5_cols_used]

    if args.p5_cols_used != [0]:
        programmed_pcr_combos = get_programmed_pcr_combos(p5_lookup, p7_lookup, args.p5_cols_used, args.p7_rows_used)
    else:
        programmed_pcr_combos = get_programmed_pcr_combos_wells(args.p5_wells_used, args.p7_wells_used)

    # Load and process RT barcodes
    if args.rt_barcode_file == "default":
        if args.level == "3":
            rtfile = RT3_FILE
        else:
            rtfile = RT_FILE
    else:
        rtfile = args.rt_barcode_file

    if args.p5_barcode_file == "default":
        p5file = P5_FILE
    else:
        p5file = args.p5_barcode_file

    if args.p7_barcode_file == "default":
        p7file = P7_FILE
    else:
        p7file = args.p7_barcode_file
    
    if args.lig_barcode_file == "default":
        ligfile = LIG_FILE
    else:
        ligfile = args.lig_barcode_file

    lane_num = args.file_name
    lane_num = lane_num.replace("Undetermined_S0_L", "L")
    lane_num = lane_num.replace("_R1_001.fastq.gz", "")
    stats_file = os.path.join(args.output_dir, lane_num + ".stats.json")
    suffix = lane_num + ".fastq"

    # Load and process ligation barcodes
    if args.level == "3":
        ligation_lookup = bu.load_whitelist(ligfile, variable_lengths=True)
        ligation_9_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 9}
        ligation_10_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 10}

    # Load barcodes
    p7_lookup = bu.load_whitelist(p7file)
    p7_lookup = {sequence[0:args.p7_length]: well for sequence,well in p7_lookup.items()}

    p5_lookup = bu.load_whitelist(p5file)
    if reverse_complement_i5:
        p5_lookup = {bu.reverse_complement(sequence): well for sequence,well in p5_lookup.items()}
    p5_lookup = {sequence[0:args.p5_length]: well for sequence,well in p5_lookup.items()}

    # Define where all sequences are and what the whitelists are
    if p5_none:
        pcr_spec = {
            'p7': {
                'start': 1,
                'end': args.p7_length,
                'read': 'i7',
                'whitelist': p7_lookup
             }
        }
    elif p7_none:
        pcr_spec = {
            'p5': {
                'start': 1,
                'end': args.p5_length,
                'read': 'i5',
                'whitelist': p5_lookup
             }
        }
    else:
        pcr_spec = {
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
            }
        }


    if args.level == "3":
        barcode_spec = {
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
    barcode_spec.update(pcr_spec)
    
    # Prep output files
    lane_num = args.file_name
    lane_num = lane_num.replace("Undetermined_S0_L", "L")
    lane_num = lane_num.replace("_R1_001.fastq.gz", "")
    stats_file = os.path.join(args.output_dir, lane_num + ".stats.json")
    suffix = lane_num + ".fastq.gz"

    if multi_exp: #, buffering=buffer_size
        sample_rt_exp_lookup = load_sample_layout(args.sample_layout, multi_exp)
        sample_to_output_filename_lookup = {sample: os.path.join(args.output_dir, '%s-%s' % (sample, suffix)) for well,sample in sample_rt_exp_lookup.items()}
        sample_to_output_file_lookup = {sample: open(filename, 'w') for sample,filename in sample_to_output_filename_lookup.items()}
        print("Demuxing %s samples (%s total RT wells) into their own files..." % (len(sample_to_output_filename_lookup), len(sample_rt_exp_lookup)))

     
    else:
        sample_rt_lookup = load_sample_layout(args.sample_layout, multi_exp)
        sample_to_output_filename_lookup = {sample: os.path.join(args.output_dir, '%s-%s' % (sample, suffix)) for well,sample in sample_rt_lookup.items()}
        sample_to_output_file_lookup = {sample: open(filename, 'w') for sample,filename in sample_to_output_filename_lookup.items()}
        print("Demuxing %s samples (%s total RT wells) into their own files..." % (len(sample_to_output_filename_lookup), len(sample_rt_lookup)))

    undetermined = open(os.path.join(args.output_dir, '%s-%s' % ("Undetermined", suffix)), 'w')
     
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

            if not p5_none:
                p5 = entry['p5']
            else:
                p5 = "none"
            if not p7_none:
                p7 = entry['p7']
            else:
                p7 = "none"
    
            ## Choose between the _8 and _9 options based on which correct properly
            corrected_9 = entry['ligation_9'] is not None and entry['rt_9'] is not None
            corrected_10 = entry['ligation_10'] is not None and entry['rt_10'] is not None
            corrected_p5_p7 = (p5_none or p5 is not None) and (p7_none or p7 is not None)

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

            # Only allow the programmed PCR combos (helps clean things up a bit)
            if not (p5, p7) in programmed_pcr_combos:
                total_pcr_mismatch += 1
                write_to_undetermined(entry)
                continue

            if multi_exp:
                experiment = exp_lookup[(p5, p7)]
                if (rt_barcode, experiment) not in sample_rt_exp_lookup:
                    total_unused_rt_well += 1
                    write_to_undetermined(entry)
                    continue
                sample = sample_rt_exp_lookup[(rt_barcode, experiment)]

            else:
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

#             if rt_barcode not in sample_rt_lookup:
#                 total_unused_rt_well += 1
#                 write_to_undetermined(entry)
#                 continue

            if multi_exp:
                experiment = exp_lookup[(p5, p7)]
                if (rt_barcode, experiment) not in sample_rt_exp_lookup:
                    total_unused_rt_well += 1
                    write_to_undetermined(entry)
                    continue
                sample = sample_rt_exp_lookup[(rt_barcode, experiment)]

            else:
                if rt_barcode not in sample_rt_lookup:
                    total_unused_rt_well += 1
                    write_to_undetermined(entry)
                    continue
                sample = sample_rt_lookup[rt_barcode]

#            sample = sample_rt_lookup[rt_barcode]

            sample_read_name = sample.split(".fq.part")[0]
            sample_read_number = sample_read_counts[sample] + 1
            sample_read_counts[sample] += 1

            if rt_barcode in rt_dict:
                rt_dict[rt_barcode] += 1
            else:
                rt_dict[rt_barcode] = 1

            r2_qual = entry['r2_qual']
            r2_seq = entry['r2_seq']
            output_name = f'@{sample_read_name}-P5{p5}-P7{p7}_{sample_read_number}|{sample_read_name}|{p5}|{p7}|{rt_barcode}|{umi}'
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

