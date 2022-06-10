#!/usr/bin/env python3

import sys
import re
import csv


def check_pcr_rxn_header(row_one_list):
    """
    Check for a PCR index file header in the supplied list.
    Args:
        row_one_list (list): fields from a CSV file row
    Returns:
        test result (string): 'short_header', 'long_header', 'no_header'
        where a short header has 3 fields and a long header has 5 fields.
    """
    len_row_one = len(row_one_list)
    if(len_row_one == 3 and
       row_one_list[0] == 'pcr_rxn_name' and
       row_one_list[1] == 'p5_index' and
       row_one_list[2] == 'p7_index'):
           return('short_header')
    elif(len_row_one == 5 and
       row_one_list[0] == 'pcr_rxn_name' and
       row_one_list[1] == 'p5_well' and
       row_one_list[2] == 'p5_index' and
       row_one_list[3] == 'p7_well' and
       row_one_list[4] == 'p7_index'):
           return('long_header')
    else:
      return('no_header')


def check_pcr_rxn_row(file_path, i_row, row_length, row_list):
    """
    Check for valid PCR reaction values in CSV file row.
    Args:
        file_path (string): path to the CSV file. Used for diagnostic message
        i_row (integer): the row number under test
        row_length (integer): the expected number of elements in row_list
        row_list (list): the fields in the i_row-th row of the CSV file
    Returns:
        test result (boolean): True if the row passes, False if it fails.
    """
    if(len(row_list) != row_length):
        print('Unexpected field count in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
        return(False)

    if(row_length == 3):
        mobj = re.match('^[a-zA-Z0-9_:]+$', row_list[0])
        if(mobj == None):
             print('Bad reaction name in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
        mobj = re.match('^[acgtACGT]+$', row_list[1])
        if(mobj == None):
             print('Bad P5 index in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
        mobj = re.match('^[acgtACGT]+$', row_list[2])
        if(mobj == None):
             print('Bad P7 index in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
    elif(row_length == 5):
        mobj = re.match('^[a-zA-Z0-9_:]+$', row_list[0])
        if(mobj == None):
             print('Bad reaction name in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
#        mobj = re.match('^P[0-9][0-9][:][A-H]([0][1-9]|[1][0-2])$', row_list[1])
        mobj = re.match('^[A-H]([0][1-9]|[1][0-2])$', row_list[1])
        if(mobj == None):
             print('Bad P5 well name in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
        mobj = re.match('^[acgtACGT]+$', row_list[2])
        if(mobj == None):
             print('Bad P5 index in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
#        mobj = re.match('^P[0-9][0-9][:][A-H]([0][1-9]|[1][0-2])$', row_list[3])
        mobj = re.match('^[A-H]([0][1-9]|[1][0-2])$', row_list[3])
        if(mobj == None):
             print('Bad P7 well name in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
        mobj = re.match('^[acgtACGT]+$', row_list[4])
        if(mobj == None):
             print('Bad P7 index in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
             return(False)
    else:
        print('Bad number of columns in row %d of PCR reaction file %s' % (i_row, file_path), file=sys.stderr)
        return(False)

    return(True)


def load_pcr_indexlist(file_path):
    """
    Load PCR index list into xxx.
    PCR index list file headers:

      pcr_rxn_name,p5_index,p7_index
    or
      pcr_rxn_name,p5_well,p5_index,p7_well,p7_index

    In the latter case, the user can specify well ids for
    use in the cell read names rather than the index
    sequences.

    The pcr_rxn_name is printable text.

    The p5_well and p7_well name format is

      <plate_id>:<well_id>

    where the plate_id has the format

      Pnn

    where nn are decimal digits. The plate format is

      <row><column>

    where <row> is [A-H] and <column> is a decimal
    number in the range 01-12; that is, A01-H12.

    The index is a DNA sequence in lower or upper case
    but will is stored as upper case.

    The fields are separated by commas.
    Args:
        file_path (string): directory path to PCR reaction index file
    Returns:
        list of PCR reaction information. Each element is one
        PCR reaction.
    """
    i_row = 1
    pcr_rxn_list = []
    with open(file_path) as fh:
        csvreader = csv.reader(fh, delimiter=',', quotechar='"')

        row_list = csvreader.__next__()

        if(check_pcr_rxn_header(row_list) == 'long_header'):
            row_length = 5
            row_list[2].upper()
            row_list[4].upper()
        elif(check_pcr_rxn_header(row_list) == 'short_header'):
            row_length = 3
            row_list[1].upper()
            row_list[2].upper()
        else:
            raise ValueError('Missing header in file %s' % (file_path))

        for row_list in csvreader:
            i_row += 1
            for i in range(len(row_list)):
                row_list[i] = row_list[i].strip()

            if(check_pcr_rxn_row(file_path, i_row, row_length, row_list) == False):
                sys.exit(-1)

            if(row_length == 3):
                row_list[0] = row_list[0]
                row_list[1] = row_list[1].upper()
                row_list[2] = row_list[2].upper()
            else:
                row_list[0] = row_list[0]
                row_list[1] = row_list[1]
                row_list[2] = row_list[2].upper()
                row_list[3] = row_list[3]
                row_list[4] = row_list[4].upper()

            pcr_rxn_list.append(row_list)
    return(pcr_rxn_list)


def make_pcr_whitelist(pcr_rxn_list, pcr_index, index_length, well_ids=False):
    """
    Make a PCR sequence/well whitelist from a PCR reaction list.
    Args:
        pcr_rxn_list (list): PCR reactions from load_pcr_indexlist()
        pcr_index (string): the PCR index to use: either 'p5' or 'p7'
        well_ids (boolean): use wells to identify primer; otherwise, use sequence
    Return:
        dictionary of PCR primers. key=sequence, value=name
    """
    if(pcr_index != 'p5' and pcr_index != 'p7'):
        raise ValueError('pcr_index parameter must be either "p5" or "p7"')

    row_length = len(pcr_rxn_list[0])
    if(row_length == 3):
        if(pcr_index == 'p5'):
            i_name = 1
            i_seq = 1
        elif(pcr_index == 'p7'):
            i_name = 2
            i_seq = 2
    elif(row_length == 5):
        if(pcr_index == 'p5'):
            i_name = 1 if well_ids else 2
            i_seq  = 2
        elif(pcr_index == 'p7'):
            i_name = 3 if well_ids else 4
            i_seq  = 4

    whitelist = dict()
    for pcr_rxn in pcr_rxn_list:
        whitelist[pcr_rxn[i_seq][0:index_length]] = pcr_rxn[i_name]

    return(whitelist)


def get_programmed_pcr_combos_pcrlist(pcr_rxn_list, well_ids=False):
    """
    Make a set of valid PCR primer name pairs. The pairs are stored
    as tuples.
    Args:
        pcr_rxn_list (list): PCR reactions from load_pcr_indexlist()
        well_ids (boolean): use wells to identify primer; otherwise, use sequence
    Return:
        set of tuples of valid PCR primer pair names
    """
    row_length = len(pcr_rxn_list[0])
    if(row_length == 3):
        i_p5_name = 1 
        i_p7_name = 2 
    elif(row_length == 5):
        i_p5_name = 1 if well_ids else 2
        i_p7_name = 3 if well_ids else 4

    valid_combos = set()
    for pcr_rxn in pcr_rxn_list:
        valid_combos.add((pcr_rxn[i_p5_name], pcr_rxn[i_p7_name]))

    return(valid_combos)


# Stand-alone testing.
if __name__ == '__main__':

    well_ids = False
    pcr_rxn_list = load_pcr_indexlist('test_index_list.csv')
    print(pcr_rxn_list)
    print('p5 whitelist')
    print(make_pcr_whitelist(pcr_rxn_list, 'p5', 3, well_ids=well_ids))
    print('p7 whitelist')
    print(make_pcr_whitelist(pcr_rxn_list, 'p7', 3, well_ids=well_ids))
    print('pcr combinations')
    print(get_programmed_pcr_combos_pcrlist(pcr_rxn_list, well_ids=well_ids))


