#!/usr/bin/env python

import os
import argparse
import run_info

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

if __name__ == '__main__':
    # Script
    parser = argparse.ArgumentParser('Script make sample sheet')

    # Required args
    parser.add_argument('--run_directory', required=True, help='Path to BCL directory for sequencing run.')
    args = parser.parse_args()

    # Get simple things like index lengths and flow cell ID for the run
    run_info = run_info.get_run_info(args.run_directory)
    if( run_info['paired_end'] == False ):
        raise ValueError('Single-end reads detected: paired-end reads required')

    sample_sheet_text = get_sample_sheet_text(run_info['p7_index_length'], run_info['p5_index_length'])
    sample_sheet_path = os.path.join('SampleSheet.csv')    
    with open(sample_sheet_path, 'w') as sample_sheet:
        sample_sheet.write(sample_sheet_text)

