#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import barcodeutils as bu
import argparse
import os
import glob
import xml.etree.ElementTree as ET
import operator

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
P5_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p5.txt')
P7_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p7.txt')
LIG_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/ligation.txt')

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
        selected_p5 = [p5_well for p5_well in p5_lookup.values() if int(p5_well[1:]) == int(p5)]
        selected_p7 = [p7_well for p7_well in p7_lookup.values() if p7_well[0] == p7]
        for selected_p5_i in selected_p5:
            for selected_p7_i in selected_p7:
                valid_combos.add((selected_p5_i, selected_p7_i))
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
    Function that loads the sample layout file to an RT lookup table.
    """
    lookup = {}
    for rt_well in quick_parse(file_path):
        lookup[rt_well['RT Barcode']] = rt_well['Sample ID'].replace('(', '.').replace(')', '.').replace(' ', '.').replace('-', '.').replace('_', '.').replace('/', '.')
    return lookup


def make_undetermined_dict_3lvl(read_file, out_file):
    bad_rt_dict = dict()
    bad_lig_dict = dict()
    not_in_samp_dict = dict()
    read_result = []
    pcr_combo_dict = dict()
    fo = open(out_file, "w")
    fo.write("read_type\trt_9_result\trt_9_value\trt_10_result\trt_10_value\tsample_assign_9\tsample_assign_10\tlig_9_result\tlig_9_value\tlig_10_result\tlig_10_value\tp5_result\tp7_result\tumi_9_value\tumi_10_value\tpcr_result\n")
    sum_dict = {'total_reads':0, 'amb_length':0, 'multi_wrong':0, 'bad_lig':0, 'bad_rt':0, 'pcr_mismatch':0, 'bad_pcr_comp':0, 'not_in_samp':0, 'unassigned':0}

    with open(read_file, 'rt') as f:
        ct = 0
        while True:
            ct += 1
            line1 = f.readline()
            if not line1:
                break
            read_res = {'read_type': "Unknown", 'sample_assign_9': "Unknown", 'sample_assign_10': "Unknown", 'rt_9_result': "Unknown", 'rt_9_value': "Unknown", 'rt_10_result': "Unknown", 'rt_10_value': "Unknown", 'lig_9_result': "Unknown", 'lig_9_value': "Unknown", 'lig_10_result': "Unknown", 'lig_10_value': "Unknown", 'p7_result': "Unknown", 'umi_9_value': "Unknown", 'umi_10_value': "Unknown", 'pcr_result': "Unknown"}
            barcs = line1.strip().split(" ")[1].split(":")[3]
            p7 = barcs.split("+")[0]
            p5 = barcs.split("+")[1].split("|")[0]
            r1 = barcs.split("|")[1]
            rt_9 = r1[23:33]
            rt_10 = r1[24:34]
            lig_9 = r1[0:9]
            lig_10 = r1[0:10]
            read_res['umi_9_value'] = r1[15:23]
            read_res['umi_10_value'] = r1[16:24]
            if rt_9 in rt3_whitelist[0].keys() or rt_9 in rt3_whitelist[1].keys():
                read_res["rt_9_result"] = "OK"
                if rt_9 in rt3_whitelist[1].keys():
                    rt_9 = rt3_whitelist[1][rt_9]
                rt_9 = rt3_lookup[rt_9]
                read_res["rt_9_value"] = rt_9
                if rt_9 in sample_rt_lookup.keys():
                    read_res["sample_assign_9"] = sample_rt_lookup[rt_9]
            else:
                read_res["rt_9_result"] = "Bad RT"
                read_res["rt_9_value"] = rt_9
            if rt_10 in rt3_whitelist[0].keys() or rt_10 in rt3_whitelist[1].keys():
                read_res["rt_10_result"] = "OK"
                if rt_10 in rt3_whitelist[1].keys():
                    rt_10 = rt3_whitelist[1][rt_10]
                rt_10 = rt3_lookup[rt_10]
                read_res["rt_10_value"] = rt_10
                if rt_10 in sample_rt_lookup.keys():
                    read_res["sample_assign_10"] = sample_rt_lookup[rt_10]
            else:
                read_res["rt_10_result"] = "Bad RT"
                read_res["rt_10_value"] = rt_10
            if lig_9 in lig_9_whitelist[0].keys() or lig_9 in lig_9_whitelist[1].keys():
                read_res["lig_9_result"] = "OK"
                if lig_9 in lig_9_whitelist[1].keys():
                    lig_9 = lig_9_whitelist[1][lig_9]
                lig_9 = lig_9_lookup[lig_9]
                read_res["lig_9_value"] = lig_9
            else:
                read_res["lig_9_result"] = "Bad Lig"
                read_res["lig_9_value"] = lig_9
            if lig_10 in lig_10_whitelist[0].keys() or lig_10 in lig_10_whitelist[1].keys():
                read_res["lig_10_result"] = "OK"
                if lig_10 in lig_10_whitelist[1].keys():
                    lig_10 = lig_10_whitelist[1][lig_10]
                lig_10 = lig_10_lookup[lig_10]
                read_res["lig_10_value"] = lig_10
            else:
                read_res["lig_10_result"] = "Bad Lig"
                read_res["lig_10_value"] = lig_10
            if read_res["rt_9_result"] == "OK" and not read_res["rt_10_result"] == "OK":
                rt_type = 9
            elif read_res["rt_10_result"] == "OK" and not read_res["rt_9_result"] == "OK":
                rt_type = 10
            else:
                rt_type = "Unknown"
            if read_res["lig_9_result"] == "OK" and not read_res["lig_10_result"] == "OK":
                lig_type = 9
            elif read_res["lig_10_result"] == "OK" and not read_res["lig_9_result"] == "OK":
                lig_type = 10
            else:
                lig_type = "Unknown"
            if lig_type == "Unknown" or rt_type == "Unknown":
                if rt_type == 9 or lig_type == 9:
                    read_type = 9
                elif rt_type == 10 or lig_type == 10:
                    read_type = 10
                else:
                    read_type = "Unknown"
            elif rt_type == 9 and lig_type == 9:
                read_type = 9
            elif rt_type == 10 and lig_type == 10:
                read_type = 10
            else:
                read_type = "Ambiguous"
            read_res['read_type'] = read_type
            if p5 in p5_whitelist[0].keys() or p5 in p5_whitelist[1].keys():
                read_res["p5_result"] = "OK"
                if p5 in p5_whitelist[1].keys():
                    p5 = p5_whitelist[1][p5]
                p5 = p5_lookup[p5]
            else:
                read_res["p5_result"] = "Bad p5"
                read_res["pcr_result"] = "Bad component"
            if p7 in p7_whitelist[0].keys() or p7 in p7_whitelist[1].keys():
                read_res["p7_result"] = "OK"
                if p7 in p7_whitelist[1].keys():
                    p7 = p7_whitelist[1][p7]
                p7 = p7_lookup[p7]
            else:
                read_res["p7_result"] = "Bad p7"
                read_res["pcr_result"] = "Bad component"
            if not read_res["pcr_result"] == "Bad component":
                if not (p5, p7) in programmed_pcr_combos:
                    read_res["pcr_result"] = "Barcode mismatch"
                    if not (p5, p7) in pcr_combo_dict.keys():
                        pcr_combo_dict[(p5, p7)] = 1
                    else:
                        pcr_combo_dict[(p5, p7)] += 1
                else:
                    read_res["pcr_result"] = "OK"
            r2 = f.readline()
            r3 = f.readline()
            r4 = f.readline()
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read_res['read_type'],\
            read_res['rt_9_result'],read_res['rt_9_value'],read_res['rt_10_result'],read_res['rt_10_value'],\
            read_res['sample_assign_9'],read_res['sample_assign_10'],read_res['lig_9_result'],read_res['lig_9_value'],\
            read_res['lig_10_result'],read_res['lig_10_value'],read_res['p5_result'],read_res['p7_result'],\
            read_res['umi_9_value'],read_res['umi_10_value'],read_res['pcr_result']))
            sum_dict['total_reads'] += 1
            if read_res['read_type'] == "Unknown" or read_res['read_type'] == "Ambiguous":
                sum_dict['amb_length'] += 1
            else:
                if read_res['read_type'] == 9:
                    if ((read_res['rt_9_result'] != 'OK') + (read_res['lig_9_result'] != 'OK') + (read_res['pcr_result'] != 'OK')) > 1:
                        sum_dict['multi_wrong'] += 1
                    elif read_res['lig_9_result'] != 'OK':
                        sum_dict['bad_lig'] += 1
                        if read_res['lig_9_value'] in bad_lig_dict.keys():
                            bad_lig_dict[read_res['lig_9_value']] += 1
                        else:
                            bad_lig_dict[read_res['lig_9_value']] = 1
                    elif read_res['rt_9_result'] != 'OK':
                        sum_dict['bad_rt'] += 1
                        if read_res['rt_9_value'] in bad_rt_dict.keys():
                            bad_rt_dict[read_res['rt_9_value']] += 1
                        else:
                            bad_rt_dict[read_res['rt_9_value']] = 1
                    elif read_res['pcr_result'] == "Barcode mismatch":
                        sum_dict['pcr_mismatch'] += 1
                    elif read_res['pcr_result'] == "Bad component":
                        sum_dict['bad_pcr_comp'] += 1
                    elif read_res['pcr_result'] == "OK" and read_res["sample_assign_9"] == "Unknown":
                        sum_dict['not_in_samp'] += 1
                        if read_res['rt_9_value'] in not_in_samp_dict.keys():
                            not_in_samp_dict[read_res['rt_9_value']] += 1
                        else:
                            not_in_samp_dict[read_res['rt_9_value']] = 1
                    else:
                        sum_dict['unassigned'] += 1
                if read_res['read_type'] == 10:
                    if ((read_res['rt_10_result'] != 'OK') + (read_res['lig_10_result'] != 'OK') + (read_res['pcr_result'] != 'OK')) > 1:
                        sum_dict['multi_wrong'] += 1
                    elif read_res['lig_10_result'] != 'OK':
                        sum_dict['bad_lig'] += 1
                        if read_res['lig_10_value'] in bad_lig_dict.keys():
                            bad_lig_dict[read_res['lig_10_value']] += 1
                        else:
                            bad_lig_dict[read_res['lig_10_value']] = 1
                    elif read_res['rt_10_result'] != 'OK':
                        sum_dict['bad_rt'] += 1
                        if read_res['rt_10_value'] in bad_rt_dict.keys():
                            bad_rt_dict[read_res['rt_10_value']] += 1
                        else:
                            bad_rt_dict[read_res['rt_10_value']] = 1
                    elif read_res['pcr_result'] == "Barcode mismatch":
                        sum_dict['pcr_mismatch'] += 1
                    elif read_res['pcr_result'] == "Bad component":
                        sum_dict['bad_pcr_comp'] += 1
                    elif read_res['pcr_result'] == "OK" and read_res["sample_assign_10"] == "Unknown":
                        sum_dict['not_in_samp'] += 1
                        if read_res['rt_10_value'] in not_in_samp_dict.keys():
                            not_in_samp_dict[read_res['rt_10_value']] += 1
                        else:
                            not_in_samp_dict[read_res['rt_10_value']] = 1
                    else:
                        sum_dict['unassigned'] += 1


    all_reads = sum_dict["bad_lig"] + sum_dict["amb_length"] + sum_dict["bad_rt"] + sum_dict["bad_pcr_comp"] + sum_dict["pcr_mismatch"] + sum_dict["multi_wrong"] + sum_dict["not_in_samp"] + sum_dict["unassigned"]
    out2 = out_file.replace(".txt", "-summary.txt")
    sf = open(out2, "w")
    sf.write("Read Recovery Summary File\n\n")

    rows = [(" ", "Count", "Percentage", "\n"),\
            ("Ambiguous length lig:", str(sum_dict["amb_length"]), str(round((sum_dict["amb_length"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Bad ligation barcode:", str(sum_dict["bad_lig"]), str(round((sum_dict["bad_lig"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Bad RT barcode:", str(sum_dict["bad_rt"]), str(round((sum_dict["bad_rt"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Bad PCR component:", str(sum_dict["bad_pcr_comp"]), str(round((sum_dict["bad_pcr_comp"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("PCR mismatch:", str(sum_dict["pcr_mismatch"]), str(round((sum_dict["pcr_mismatch"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Garbage reads*:", str(sum_dict["multi_wrong"]), str(round((sum_dict["multi_wrong"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("RT not in sample sheet:", str(sum_dict["not_in_samp"]), str(round((sum_dict["not_in_samp"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Unassigned issue:", str(sum_dict["unassigned"]), str(round((sum_dict["unassigned"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Total:", str(all_reads), str(round((all_reads*100.0)/sum_dict["total_reads"], 2)), "\n")]
    
    cols = zip(*rows)
    col_widths = [ max(len(value) for value in col) for col in cols ]
    format = ' '.join(['%%%ds' % width for width in col_widths ])
    for row in rows:
        sf.write(format % tuple(row))

    sf.write("\n\n* Garbage reads have multiple problems\n\n")

    sf.write("\n\nThe top RT barcodes not in the sample sheet are:\n")
    sorted_notrt = sorted(not_in_samp_dict.items(), key=operator.itemgetter(1), reverse=True)
    for item in sorted_notrt[:20]:
        sf.write('{}    {}\n'.format(item[0],item[1]))

    sf.write("\n\nThe top bad RT barcodes are:\n")
    sorted_rt = sorted(bad_rt_dict.items(), key=operator.itemgetter(1), reverse=True)
    for item in sorted_rt[:20]:
        sf.write('{}    {}\n'.format(item[0],item[1]))

    sf.write("\n\nThe top bad ligation barcodes are:\n")
    sorted_lig = sorted(bad_lig_dict.items(), key=operator.itemgetter(1), reverse=True)
    for item in sorted_lig[:20]:
        sf.write('{}    {}\n'.format(item[0],item[1]))

    sf.write("\n\nThe top bad pcr combos are:\n")
    sorted_pcr = sorted(pcr_combo_dict.items(), key=operator.itemgetter(1), reverse=True)
    for item in sorted_pcr[:20]:
        sf.write('{}    {}\n'.format(item[0],item[1]))

    sf.close()
    return pcr_combo_dict


def make_undetermined_dict_2lvl(read_file, out_file):
    bad_rt_dict = dict()
    not_in_samp_dict = dict()
    read_result = []
    pcr_combo_dict = dict()
    fo = open(out_file, "w")
    fo.write("rt_result\trt_value\tsample_assign\tp5_result\tp7_result\tumi_value\tpcr_result\n")
    sum_dict = {'total_reads':0, 'multi_wrong':0, 'bad_rt':0, 'pcr_mismatch':0, 'bad_pcr_comp':0, 'not_in_samp':0, 'unassigned':0}

    with open(read_file, 'rt') as f:
        ct = 0
        while True:
            ct += 1
            line1 = f.readline()
            if not line1:
                break
            read_res = {'sample_assign': "Unknown", 'rt_result': "Unknown", 'rt_value': "Unknown", 'p7_result': "Unknown", 'umi_value': "Unknown", 'pcr_result': "Unknown"}
            barcs = line1.strip().split(" ")[1].split(":")[3]
            p7 = barcs.split("+")[0]
            p5 = barcs.split("+")[1].split("|")[0]
            r1 = barcs.split("|")[1]
            rt = r1[8:19]
            read_res['umi_value'] = r1[0:8]
            if rt in rt2_whitelist[0].keys() or rt in rt2_whitelist[1].keys():
                read_res["rt_result"] = "OK"
                if rt in rt2_whitelist[1].keys():
                    rt = rt2_whitelist[1][rt]
                rt = rt2_lookup[rt]
                read_res["rt_value"] = rt
                if rt in sample_rt_lookup.keys():
                    read_res["sample_assign_9"] = sample_rt_lookup[rt]
            else:
                read_res["rt_result"] = "Bad RT"
                read_res["rt_value"] = rt
            if p5 in p5_whitelist[0].keys() or p5 in p5_whitelist[1].keys():
                read_res["p5_result"] = "OK"
                if p5 in p5_whitelist[1].keys():
                    p5 = p5_whitelist[1][p5]
                p5 = p5_lookup[p5]
            else:
                read_res["p5_result"] = "Bad p5"
                read_res["pcr_result"] = "Bad component"
            if p7 in p7_whitelist[0].keys() or p7 in p7_whitelist[1].keys():
                read_res["p7_result"] = "OK"
                if p7 in p7_whitelist[1].keys():
                    p7 = p7_whitelist[1][p7]
                p7 = p7_lookup[p7]
            else:
                read_res["p7_result"] = "Bad p7"
                read_res["pcr_result"] = "Bad component"
            if not read_res["pcr_result"] == "Bad component":
                if not (p5, p7) in programmed_pcr_combos:
                    read_res["pcr_result"] = "Barcode mismatch"
                    if not (p5, p7) in pcr_combo_dict.keys():
                        pcr_combo_dict[(p5, p7)] = 1
                    else:
                        pcr_combo_dict[(p5, p7)] += 1
                else:
                    read_res["pcr_result"] = "OK"
            r2 = f.readline()
            r3 = f.readline()
            r4 = f.readline()
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read_res['rt_result'],\
            read_res['rt_value'],read_res['sample_assign'],read_res['p5_result'],read_res['p7_result'],\
            read_res['umi_value'],read_res['pcr_result']))
            sum_dict['total_reads'] += 1

            if ((read_res['rt_result'] != 'OK') + (read_res['pcr_result'] != 'OK')) > 1:
                sum_dict['multi_wrong'] += 1
            elif read_res['rt_result'] != 'OK':
                sum_dict['bad_rt'] += 1
                if read_res['rt_value'] in bad_rt_dict.keys():
                    bad_rt_dict[read_res['rt_value']] += 1
                else:
                    bad_rt_dict[read_res['rt_value']] = 1
            elif read_res['pcr_result'] == "Barcode mismatch":
                sum_dict['pcr_mismatch'] += 1
            elif read_res['pcr_result'] == "Bad component":
                sum_dict['bad_pcr_comp'] += 1
            elif read_res['pcr_result'] == "OK" and read_res["sample_assign"] == "Unknown":
                sum_dict['not_in_samp'] += 1
                if read_res['rt_value'] in not_in_samp_dict.keys():
                    not_in_samp_dict[read_res['rt_value']] += 1
                else:
                    not_in_samp_dict[read_res['rt_value']] = 1
            else:
                sum_dict['unassigned'] += 1

    all_reads = sum_dict["bad_rt"] + sum_dict["bad_pcr_comp"] + sum_dict["pcr_mismatch"] + sum_dict["multi_wrong"] + sum_dict["not_in_samp"] + sum_dict["unassigned"]
    out2 = out_file.replace(".txt", "-summary.txt")
    sf = open(out2, "w")
    sf.write("Read Recovery Summary File\n\n")

    rows = [(" ", "Count", "Percentage", "\n"),\
            ("Bad RT barcode:", str(sum_dict["bad_rt"]), str(round((sum_dict["bad_rt"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Bad PCR component:", str(sum_dict["bad_pcr_comp"]), str(round((sum_dict["bad_pcr_comp"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("PCR mismatch:", str(sum_dict["pcr_mismatch"]), str(round((sum_dict["pcr_mismatch"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Garbage reads*:", str(sum_dict["multi_wrong"]), str(round((sum_dict["multi_wrong"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("RT not in sample sheet:", str(sum_dict["not_in_samp"]), str(round((sum_dict["not_in_samp"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Unassigned issue:", str(sum_dict["unassigned"]), str(round((sum_dict["unassigned"]*100.0)/sum_dict["total_reads"], 2)), "\n"),\
            ("Total:", str(all_reads), str(round((all_reads*100.0)/sum_dict["total_reads"], 2)), "\n")]
    
    cols = zip(*rows)
    col_widths = [ max(len(value) for value in col) for col in cols ]
    format = ' '.join(['%%%ds' % width for width in col_widths ])
    for row in rows:
        sf.write(format % tuple(row))

    sf.write("\n\n* Garbage reads have multiple problems\n\n")

    sf.write("\n\nThe top RT barcodes not in the sample sheet are:\n")
    sorted_notrt = sorted(not_in_samp_dict.items(), key=operator.itemgetter(1), reverse=True)
    for item in sorted_notrt[:20]:
        sf.write('{}    {}\n'.format(item[0],item[1]))

    sf.write("\n\nThe top bad RT barcodes are:\n")
    sorted_rt = sorted(bad_rt_dict.items(), key=operator.itemgetter(1), reverse=True)
    for item in sorted_rt[:20]:
        sf.write('{}    {}\n'.format(item[0],item[1]))

    sf.write("\n\nThe top bad pcr combos are:\n")
    sorted_pcr = sorted(pcr_combo_dict.items(), key=operator.itemgetter(1), reverse=True)
    for item in sorted_pcr[:20]:
        sf.write('{}    {}\n'.format(item[0],item[1]))

    sf.close()
    return pcr_combo_dict

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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to combine R1 and R2 file from sci run into a single file with barcodes in the fastq header.')
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    parser.add_argument('--input_file', required=True, help='File name for input file.')
    parser.add_argument('--output_file', required=True, help='File name for output file.')
    parser.add_argument('--rt_barcodes', default="default", help="file name for custom RT barcodes.")
#    parser.add_argument('--read2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
#    parser.add_argument('--file_name', required=True, help='The R1 file name.')
    parser.add_argument('--sample_layout', required=True, help='Text file containing the sample layout by RT well.')
    parser.add_argument('--p5_cols_used', nargs='+', type=int, required=True, help='A list of the columns used from P5 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p7_rows_used', nargs='+', required=True, help='A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
#    parser.add_argument('--output_dir', required=True, help='Output directory for files.')
    parser.add_argument('--p7_length', type=int, default=10, help='Expected P7 index length.')
    parser.add_argument('--p5_length', type=int, default=10, help='Expected P5 index length.')
    parser.add_argument('--level', required=True, help = "2 or 3 level sci?")
    args = parser.parse_args()

    rt_file = os.path.join(SCRIPT_DIR, 'barcode_files/rt2.txt')
    rt3_file = os.path.join(SCRIPT_DIR, 'barcode_files/rt.txt')
    if args.rt_barcodes != "default":
        rt_file = args.rt_barcodes
        rt3_file = args.rt_barcodes

    sample_rt_lookup = load_sample_layout(args.sample_layout)
    p7_lookup = bu.load_whitelist(P7_FILE)
    p7_lookup = {sequence[0:args.p7_length]: well for sequence,well in p7_lookup.items()}
    p7_whitelist = bu.construct_mismatch_to_whitelist_map(p7_lookup, edit_distance = 1)

    p5_lookup = bu.load_whitelist(P5_FILE)
    if reverse_complement_i5(args.run_directory):
        p5_lookup = {bu.reverse_complement(sequence): well for sequence,well in p5_lookup.items()}
    p5_lookup = {sequence[0:args.p5_length]: well for sequence,well in p5_lookup.items()}
    p5_whitelist = bu.construct_mismatch_to_whitelist_map(p5_lookup, edit_distance = 1)

    programmed_pcr_combos = get_programmed_pcr_combos(p5_lookup, p7_lookup, args.p5_cols_used, args.p7_rows_used)

    rt2_lookup = bu.load_whitelist(rt_file)
    rt2_whitelist = bu.construct_mismatch_to_whitelist_map(rt2_lookup, edit_distance = 1)

    rt3_lookup = bu.load_whitelist(rt3_file)
    rt3_whitelist = bu.construct_mismatch_to_whitelist_map(rt3_lookup, edit_distance = 1)


    ligation_lookup = bu.load_whitelist(LIG_FILE, variable_lengths=True)
    lig_9_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 9}
    lig_10_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 10}
            
    lig_10_whitelist = bu.construct_mismatch_to_whitelist_map(lig_10_lookup, edit_distance = 1)
    lig_9_whitelist = bu.construct_mismatch_to_whitelist_map(lig_9_lookup, edit_distance = 1)

    if args.level == "3":
        x = make_undetermined_dict_3lvl(args.input_file, args.output_file)
    else:
        x = make_undetermined_dict_2lvl(args.input_file, args.output_file)



