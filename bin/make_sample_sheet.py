#!/usr/bin/env python

import os
import barcodeutils as bu
import argparse
import glob
import xml.etree.ElementTree as ET
import re

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

def get_sample_sheet_text(p7_index_length, p5_index_length):
    """
    Gets the sample sheet text that will demux cells into one set of files
    """

    sample_sheet_template = """[DATA]
Lane,Sample_ID,Sample_Name,index,index2
%s"""

    line = ',fake,fake,%s,%s' % ('N' * p7_index_length, 'N' * p5_index_length)
    return sample_sheet_template % line


NEXTSEQ = 'NextSeq'
NEXTSEQ2000 = 'NextSeq 1000/2000'
MISEQ = 'MiSeq'
NOVASEQ = 'NovaSeq'
HISEQ4000 = 'HiSeq4000'
HISEQ3000 = 'HiSeq3000'
HISEQ = 'HiSeq'
UNKNOWN_SEQUENCER = 'unknown'

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

    sample_sheet_text = get_sample_sheet_text(run_info['p7_index_length'], run_info['p5_index_length'])
    sample_sheet_path = os.path.join('SampleSheet.csv')    
    with open(sample_sheet_path, 'w') as sample_sheet:
        sample_sheet.write(sample_sheet_text)

