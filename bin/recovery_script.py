#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import barcodeutils as bu
import argparse
import os
import glob
import xml.etree.ElementTree as ET
import operator
import run_info

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

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
    fo.write("read_type\trt_9_result\trt_9_value\trt_10_result\trt_10_value\tsample_assign_9\tsample_assign_10\tlig_9_result\tlig_9_value\tlig_10_result\tlig_10_value\tp5_result\tp7_result\tumi_9_value\tumi_10_value\tpcr_result\tread_name\n")
    sum_dict = {'total_reads':0, 'amb_length':0, 'multi_wrong':0, 'bad_lig':0, 'bad_rt':0, 'pcr_mismatch':0, 'bad_pcr_comp':0, 'not_in_samp':0, 'unassigned':0}

    with open(read_file, 'rt') as f:
        ct = 0
        while True:
            ct += 1
            line1 = f.readline()
            if not line1:
                break
            read_name = line1.split(" ")[0]
            read_res = {'read_type': "Unknown", 'sample_assign_9': "Unknown", 'sample_assign_10': "Unknown", 'rt_9_result': "Unknown", 'rt_9_value': "Unknown", 'rt_10_result': "Unknown", 'rt_10_value': "Unknown", 'lig_9_result': "Unknown", 'lig_9_value': "Unknown", 'lig_10_result': "Unknown", 'lig_10_value': "Unknown", 'p7_result': "Unknown", 'umi_9_value': "Unknown", 'umi_10_value': "Unknown", 'pcr_result': "Unknown"}
            barcs = line1.strip().split(" ")[1].split(":")[3]
            if p5_none:
                p7 = barcs.split("|")[0]
                p5 = "none"
            elif p7_none:
                p7 = "none"
                p5 = barcs.split("|")[0]
            else:
                p7 = barcs.split("+")[0]
                p5 = barcs.split("+")[1].split("|")[0]
            r1 = barcs.split("|")[1]
            rt_9 = r1[23:33]
            rt_10 = r1[24:34]
            lig_9 = r1[0:9]
            lig_10 = r1[0:10]
            read_res['umi_9_value'] = r1[15:23]
            read_res['umi_10_value'] = r1[16:24]
            if rt_9 in rt3_whitelist[0] or rt_9 in rt3_whitelist[1]:
                read_res["rt_9_result"] = "OK"
                if rt_9 in rt3_whitelist[1]:
                    rt_9 = rt3_whitelist[1][rt_9]
                rt_9 = rt3_lookup[rt_9]
                read_res["rt_9_value"] = rt_9
                if rt_9 in sample_rt_lookup:
                    read_res["sample_assign_9"] = sample_rt_lookup[rt_9]
            else:
                read_res["rt_9_result"] = "Bad RT"
                read_res["rt_9_value"] = rt_9
            if rt_10 in rt3_whitelist[0] or rt_10 in rt3_whitelist[1]:
                read_res["rt_10_result"] = "OK"
                if rt_10 in rt3_whitelist[1]:
                    rt_10 = rt3_whitelist[1][rt_10]
                rt_10 = rt3_lookup[rt_10]
                read_res["rt_10_value"] = rt_10
                if rt_10 in sample_rt_lookup:
                    read_res["sample_assign_10"] = sample_rt_lookup[rt_10]
            else:
                read_res["rt_10_result"] = "Bad RT"
                read_res["rt_10_value"] = rt_10
            if lig_9 in lig_9_whitelist[0] or lig_9 in lig_9_whitelist[1]:
                read_res["lig_9_result"] = "OK"
                if lig_9 in lig_9_whitelist[1]:
                    lig_9 = lig_9_whitelist[1][lig_9]
                lig_9 = lig_9_lookup[lig_9]
                read_res["lig_9_value"] = lig_9
            else:
                read_res["lig_9_result"] = "Bad Lig"
                read_res["lig_9_value"] = lig_9
            if lig_10 in lig_10_whitelist[0] or lig_10 in lig_10_whitelist[1]:
                read_res["lig_10_result"] = "OK"
                if lig_10 in lig_10_whitelist[1]:
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
            if p5 in p5_whitelist[0] or p5 in p5_whitelist[1]:
                read_res["p5_result"] = "OK"
                if p5 in p5_whitelist[1]:
                    p5 = p5_whitelist[1][p5]
                p5 = p5_lookup[p5]
            else:
                read_res["p5_result"] = "Bad p5"
                read_res["pcr_result"] = "Bad component"
            if p7 in p7_whitelist[0] or p7 in p7_whitelist[1]:
                read_res["p7_result"] = "OK"
                if p7 in p7_whitelist[1]:
                    p7 = p7_whitelist[1][p7]
                p7 = p7_lookup[p7]
            else:
                read_res["p7_result"] = "Bad p7"
                read_res["pcr_result"] = "Bad component"
            if not read_res["pcr_result"] == "Bad component":
                if not (p5, p7) in programmed_pcr_combos:
                    read_res["pcr_result"] = "Barcode mismatch"
                    if not (p5, p7) in pcr_combo_dict:
                        pcr_combo_dict[(p5, p7)] = 1
                    else:
                        pcr_combo_dict[(p5, p7)] += 1
                else:
                    read_res["pcr_result"] = "OK"
            r2 = f.readline()
            r3 = f.readline()
            r4 = f.readline()
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read_res['read_type'],\
            read_res['rt_9_result'],read_res['rt_9_value'],read_res['rt_10_result'],read_res['rt_10_value'],\
            read_res['sample_assign_9'],read_res['sample_assign_10'],read_res['lig_9_result'],read_res['lig_9_value'],\
            read_res['lig_10_result'],read_res['lig_10_value'],read_res['p5_result'],read_res['p7_result'],\
            read_res['umi_9_value'],read_res['umi_10_value'],read_res['pcr_result'],read_name))
            sum_dict['total_reads'] += 1
            if read_res['read_type'] == "Unknown" or read_res['read_type'] == "Ambiguous":
                sum_dict['amb_length'] += 1
            else:
                if read_res['read_type'] == 9:
                    if ((read_res['rt_9_result'] != 'OK') + (read_res['lig_9_result'] != 'OK') + (read_res['pcr_result'] != 'OK')) > 1:
                        sum_dict['multi_wrong'] += 1
                    elif read_res['lig_9_result'] != 'OK':
                        sum_dict['bad_lig'] += 1
                        if read_res['lig_9_value'] in bad_lig_dict:
                            bad_lig_dict[read_res['lig_9_value']] += 1
                        else:
                            bad_lig_dict[read_res['lig_9_value']] = 1
                    elif read_res['rt_9_result'] != 'OK':
                        sum_dict['bad_rt'] += 1
                        if read_res['rt_9_value'] in bad_rt_dict:
                            bad_rt_dict[read_res['rt_9_value']] += 1
                        else:
                            bad_rt_dict[read_res['rt_9_value']] = 1
                    elif read_res['pcr_result'] == "Barcode mismatch":
                        sum_dict['pcr_mismatch'] += 1
                    elif read_res['pcr_result'] == "Bad component":
                        sum_dict['bad_pcr_comp'] += 1
                    elif read_res['pcr_result'] == "OK" and read_res["sample_assign_9"] == "Unknown":
                        sum_dict['not_in_samp'] += 1
                        if read_res['rt_9_value'] in not_in_samp_dict:
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
                        if read_res['lig_10_value'] in bad_lig_dict:
                            bad_lig_dict[read_res['lig_10_value']] += 1
                        else:
                            bad_lig_dict[read_res['lig_10_value']] = 1
                    elif read_res['rt_10_result'] != 'OK':
                        sum_dict['bad_rt'] += 1
                        if read_res['rt_10_value'] in bad_rt_dict:
                            bad_rt_dict[read_res['rt_10_value']] += 1
                        else:
                            bad_rt_dict[read_res['rt_10_value']] = 1
                    elif read_res['pcr_result'] == "Barcode mismatch":
                        sum_dict['pcr_mismatch'] += 1
                    elif read_res['pcr_result'] == "Bad component":
                        sum_dict['bad_pcr_comp'] += 1
                    elif read_res['pcr_result'] == "OK" and read_res["sample_assign_10"] == "Unknown":
                        sum_dict['not_in_samp'] += 1
                        if read_res['rt_10_value'] in not_in_samp_dict:
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
    fo.write("rt_result\trt_value\tsample_assign\tp5_result\tp7_result\tumi_value\tpcr_result\tread_name\n")
    sum_dict = {'total_reads':0, 'multi_wrong':0, 'bad_rt':0, 'pcr_mismatch':0, 'bad_pcr_comp':0, 'not_in_samp':0, 'unassigned':0}

    with open(read_file, 'rt') as f:
        ct = 0
        while True:
            ct += 1
            line1 = f.readline()
            if not line1:
                break
            read_name = line1.split(" ")[0]
            read_res = {'sample_assign': "Unknown", 'rt_result': "Unknown", 'rt_value': "Unknown", 'p7_result': "Unknown", 'umi_value': "Unknown", 'pcr_result': "Unknown"}
            barcs = line1.strip().split(" ")[1].split(":")[3]
            if p5_none:
                p7 = barcs.split("|")[0]
                p5 = "none"
            elif p7_none:
                p7 = "none"
                p5 = barcs.split("|")[0]
            else:
                p7 = barcs.split("+")[0]
                p5 = barcs.split("+")[1].split("|")[0]
            r1 = barcs.split("|")[1]
            rt = r1[8:18]
            read_res['umi_value'] = r1[0:8]
            if rt in rt2_whitelist[0] or rt in rt2_whitelist[1]:
                read_res["rt_result"] = "OK"
                if rt in rt2_whitelist[1]:
                    rt = rt2_whitelist[1][rt]
                rt = rt2_lookup[rt]
                read_res["rt_value"] = rt
                if rt in sample_rt_lookup:
                    read_res["sample_assign_9"] = sample_rt_lookup[rt]
            else:
                read_res["rt_result"] = "Bad RT"
                read_res["rt_value"] = rt
            if p5 in p5_whitelist[0] or p5 in p5_whitelist[1]:
                read_res["p5_result"] = "OK"
                if p5 in p5_whitelist[1]:
                    p5 = p5_whitelist[1][p5]
                p5 = p5_lookup[p5]
            else:
                read_res["p5_result"] = "Bad p5"
                read_res["pcr_result"] = "Bad component"
            if p7 in p7_whitelist[0] or p7 in p7_whitelist[1]:
                read_res["p7_result"] = "OK"
                if p7 in p7_whitelist[1]:
                    p7 = p7_whitelist[1][p7]
                p7 = p7_lookup[p7]
            else:
                read_res["p7_result"] = "Bad p7"
                read_res["pcr_result"] = "Bad component"
            if not read_res["pcr_result"] == "Bad component":
                if not (p5, p7) in programmed_pcr_combos:
                    read_res["pcr_result"] = "Barcode mismatch"
                    if not (p5, p7) in pcr_combo_dict:
                        pcr_combo_dict[(p5, p7)] = 1
                    else:
                        pcr_combo_dict[(p5, p7)] += 1
                else:
                    read_res["pcr_result"] = "OK"
            r2 = f.readline()
            r3 = f.readline()
            r4 = f.readline()
            fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(read_res['rt_result'],\
            read_res['rt_value'],read_res['sample_assign'],read_res['p5_result'],read_res['p7_result'],\
            read_res['umi_value'],read_res['pcr_result'],read_name))
            sum_dict['total_reads'] += 1

            if ((read_res['rt_result'] != 'OK') + (read_res['pcr_result'] != 'OK')) > 1:
                sum_dict['multi_wrong'] += 1
            elif read_res['rt_result'] != 'OK':
                sum_dict['bad_rt'] += 1
                if read_res['rt_value'] in bad_rt_dict:
                    bad_rt_dict[read_res['rt_value']] += 1
                else:
                    bad_rt_dict[read_res['rt_value']] = 1
            elif read_res['pcr_result'] == "Barcode mismatch":
                sum_dict['pcr_mismatch'] += 1
            elif read_res['pcr_result'] == "Bad component":
                sum_dict['bad_pcr_comp'] += 1
            elif read_res['pcr_result'] == "OK" and read_res["sample_assign"] == "Unknown":
                sum_dict['not_in_samp'] += 1
                if read_res['rt_value'] in not_in_samp_dict:
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


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Script to combine R1 and R2 file from sci run into a single file with barcodes in the fastq header.')
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    parser.add_argument('--input_file', required=True, help='File name for input file.')
    parser.add_argument('--output_file', required=True, help='File name for output file.')
    parser.add_argument('--rt_barcodes', default="default", help="file name for custom RT barcodes.")
    parser.add_argument('--sample_layout', required=True, help='Text file containing the sample layout by RT well.')
    parser.add_argument('--p5_cols_used', nargs='+', required=True, help='A list of the columns used from P5 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p7_rows_used', nargs='+', required=True, help='A list of the rows used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A B C for p7 and --p5 1 2 3 for p5.')
    parser.add_argument('--p5_wells_used', nargs='+', required=True, help='A list of the wells used from P5 plate for PCR in same order as P7 to indicate the pairs of P7 and P5 used (e.g. --p7 A1 B1 C1 for p7 and --p5 A1 A2 A3 for p5. Alternative to p5_cols_used.')
    parser.add_argument('--p7_wells_used', nargs='+', required=True, help='A list of the wells used from P7 plate for PCR in same order as P5 to indicate the pairs of P7 and P5 used (e.g. --p7 A1 B1 C1 for p7 and --p5 A1 A2 A3 for p5. Alternative to p7_rows_used.')
    parser.add_argument('--p7_length', type=int, default=10, help='Expected P7 index length.')
    parser.add_argument('--p5_length', type=int, default=10, help='Expected P5 index length.')
    parser.add_argument('--level', required=True, help = "2 or 3 level sci?")
    parser.add_argument('--p5_barcode_file', required=True, help='Path to p5 barcode file, or "default".')
    parser.add_argument('--p7_barcode_file', required=True, help='Path to p7 barcode file, or "default".')
    parser.add_argument('--lig_barcode_file', required=True, help='Path to ligation barcode file, or "default".')
    args = parser.parse_args()

    run_info = run_info.get_run_info( args.run_directory, pipeline_type='RNA-seq' )
    if( run_info['paired_end'] == False ):
        raise ValueError('Single-end reads detected: paired-end reads required')

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

    rt_file = os.path.join(SCRIPT_DIR, 'barcode_files/rt2.txt')
    rt3_file = os.path.join(SCRIPT_DIR, 'barcode_files/rt.txt')
    p5_file = os.path.join(SCRIPT_DIR, 'barcode_files/p5.txt')
    p7_file = os.path.join(SCRIPT_DIR, 'barcode_files/p7.txt') 
    lig_file = os.path.join(SCRIPT_DIR, 'barcode_files/ligation.txt')

    if args.rt_barcodes != "default":
        rt_file = args.rt_barcodes
        rt3_file = args.rt_barcodes

    if args.p5_barcode_file != "default":
        p5_file = args.p5_barcode_file

    if args.p7_barcode_file != "default":
        p7_file = args.p7_barcode_file

    if args.lig_barcode_file != "default":
        lig_file = args.lig_barcode_file

    sample_rt_lookup = load_sample_layout(args.sample_layout)

    if p7_none:
        p7_whitelist = ["none", "no_correct"]
        p7_lookup = {"none":"none"}
    else:
        p7_lookup = bu.load_whitelist(p7_file)
        p7_lookup = {sequence[0:args.p7_length]: well for sequence,well in p7_lookup.items()}
        p7_whitelist = bu.construct_mismatch_to_whitelist_map(p7_lookup, edit_distance = 1)

    if p5_none:
        p5_whitelist = ["none", "no_correct"]
        p5_lookup = {"none":"none"}
    else:
        p5_lookup = bu.load_whitelist(p5_file)
        if run_info['reverse_complement_i5']:
            p5_lookup = {bu.reverse_complement(sequence): well for sequence,well in p5_lookup.items()}
        p5_lookup = {sequence[0:args.p5_length]: well for sequence,well in p5_lookup.items()}
        p5_whitelist = bu.construct_mismatch_to_whitelist_map(p5_lookup, edit_distance = 1)

    if args.p5_cols_used != [0]:
        programmed_pcr_combos = get_programmed_pcr_combos(p5_lookup, p7_lookup, args.p5_cols_used, args.p7_rows_used)
    else:
        programmed_pcr_combos = get_programmed_pcr_combos_wells(args.p5_wells_used, args.p7_wells_used)
    rt2_lookup = bu.load_whitelist(rt_file)
    rt2_whitelist = bu.construct_mismatch_to_whitelist_map(rt2_lookup, edit_distance = 1)

    rt3_lookup = bu.load_whitelist(rt3_file)
    rt3_whitelist = bu.construct_mismatch_to_whitelist_map(rt3_lookup, edit_distance = 1)

    ligation_lookup = bu.load_whitelist(lig_file, variable_lengths=True)
    lig_9_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 9}
    lig_10_lookup = {barcode:well for barcode,well in ligation_lookup.items() if len(barcode) == 10}
            
    lig_10_whitelist = bu.construct_mismatch_to_whitelist_map(lig_10_lookup, edit_distance = 1)
    lig_9_whitelist = bu.construct_mismatch_to_whitelist_map(lig_9_lookup, edit_distance = 1)

    if args.level == "3":
        x = make_undetermined_dict_3lvl(args.input_file, args.output_file)
    else:
        x = make_undetermined_dict_2lvl(args.input_file, args.output_file)



