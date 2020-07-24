#!/usr/bin/env python
# Andrew's sample sheet creator

import os
import barcodeutils as bu
import argparse
import glob
import xml.etree.ElementTree as ET
import re

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
P5_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p5.txt')
P7_FILE = os.path.join(SCRIPT_DIR, 'barcode_files/p7.txt')

def get_sample_sheet_text(p7_index_length, p5_index_length):
    """
    Gets the sample sheet text that will demux cells into one set of files
    """

    sample_sheet_template = """[DATA]
Lane,Sample_ID,Sample_Name,index,index2
%s"""

    line = ',fake,fake,%s,%s' % ('N' * p7_index_length, 'N' * p5_index_length)
    return sample_sheet_template % line


def load_barcode_file(barcode_file, index_length=None):
    """
    Loads one of the sciRNA barcode files allowing it to be trimmed to a specified length.
    """
    barcode_dict = {}
    for line in open(barcode_file):
        index, barcode = line.strip().split('\t')

        if index_length:
            if len(barcode) < index_length:
                raise ValueError('The specified barcodes are not long enough to match the run index length: %s vs. %s' % (len(barcode), index_length))
            barcode = barcode[0:index_length]

        if index in barcode_dict and barcode_dict[index] != barcode:
            raise ValueError('Conflicting entries in barcode file for index: %s, %s' % (barcode_file, index))
        barcode_dict[index] = barcode
    return barcode_dict


NEXTSEQ = 'NextSeq'
NEXTSEQ2000 = 'NextSeq 1000/2000'
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
    HISEQ3000: False,
    NEXTSEQ2000: True
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
    if NEXTSEQ2000 in application:
        run_stats['instrument_type'] = NEXTSEQ2000
    elif NEXTSEQ in application:
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

    # Now actually populate various stats

    if run_stats['instrument_type'] == NOVASEQ:
        run_stats['p7_index_length'] = int(setup_node.find('PlannedIndex1ReadCycles').text)
        run_stats['p5_index_length'] = int(setup_node.find('PlannedIndex2ReadCycles').text)
    elif run_stats['instrument_type'] == NEXTSEQ2000:
        run_stats['p7_index_length'] = int(setup_node.find('PlannedCycles').find('Index1').text)
        run_stats['p5_index_length'] = int(setup_node.find('PlannedCycles').find('Index1').text)
    else:
        run_stats['p7_index_length'] = int(setup_node.find('Index1Read').text)
        run_stats['p5_index_length'] = int(setup_node.find('Index2Read').text)

    return run_stats


if __name__ == '__main__':
    # Script
    parser = argparse.ArgumentParser('Script make sample sheet')

    # Required args
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    args = parser.parse_args()

    #run_stats['date'] = run_start_date_node.text
    # Get simple things like index lengths and flow cell ID for the run
    run_info = get_run_info(args.run_directory)

    # Set up samplesheet for BCL2FASTQ
    p5_indices = load_barcode_file(P5_FILE)
    reverse_i5 = reverse_complement_i5(args.run_directory)
    if reverse_i5:
        p5_indices = {x:bu.reverse_complement(y) for x,y in p5_indices.items()}
    p5_indices = {x: y[0:run_info['p5_index_length']] for x,y in p5_indices.items()}

    p7_indices = load_barcode_file(P7_FILE, run_info['p7_index_length'])

    sample_sheet_text = get_sample_sheet_text(run_info['p7_index_length'], run_info['p5_index_length'])
    sample_sheet_path = os.path.join('SampleSheet.csv')    
    with open(sample_sheet_path, 'w') as sample_sheet:
        sample_sheet.write(sample_sheet_text)

