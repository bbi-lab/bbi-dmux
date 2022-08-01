#!/usr/bin/env python3


# Version: 20220615a


import sys
import re
import csv
import string


if sys.version_info[0] >= 3:
    revcomp = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
else:
    from itertools import izip as zip
    revcomp = string.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')


def reverse_complement_sequence(seq):
    return seq.translate(revcomp)[::-1]


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
        print('Error: unexpected field count in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
        return(False)

    # Check for expected and consistent column counts.
    if(row_length == 3):
        # Check for valid characters in reaction names (first column).
        mobj = re.match('^[-a-zA-Z0-9_]+$', row_list[0])
        if(mobj == None):
             print('Error: bad reaction name in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
        # Check for valid bases in P5 index sequences.
        mobj = re.match('^([acgtACGT]+|none)$', row_list[1])
        if(mobj == None):
             print('Error: bad P5 index in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
        # Check for valid bases in P7 index sequences.
        mobj = re.match('^([acgtACGT]+|none)$', row_list[2])
        if(mobj == None):
             print('Error: bad P7 index in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
    elif(row_length == 5):
        # Check for valid characters in reaction names (first column).
        mobj = re.match('^[-a-zA-Z0-9_]+$', row_list[0])
        if(mobj == None):
             print('Error: bad reaction name in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
        mobj = re.match('^([pP][0-9]+[-])?[a-hA-H]([0][1-9]|[1][0-2])$', row_list[1])
        if(mobj == None):
             print('Error: bad P5 well name in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
        # Check for valid bases in P5 index sequences.
        mobj = re.match('^([acgtACGT]+|none)$', row_list[2])
        if(mobj == None):
             print('Error: bad P5 index in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
        mobj = re.match('^([pP][0-9]+[-])?[a-hA-H]([0][1-9]|[1][0-2])$', row_list[3])
        if(mobj == None):
             print('Error: bad P7 well name in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
        # Check for valid bases in P7 index sequences.
        mobj = re.match('^([acgtACGT]+|none)$', row_list[4])
        if(mobj == None):
             print('Error: bad P7 index in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
             return(False)
    else:
        print('Bad number of columns in row %d of PCR reaction file "%s".' % (i_row, file_path), file=sys.stderr)
        return(False)

    return(True)


def check_pcr_rxn_file(file_path, pcr_rxn_list, row_length):
    """
    Check PCR index file row data for various errors.
    Args:
        file_path (string): path to the CSV file. Used for diagnostic message
        pcr_rxn_list (list): PCR reactions from load_pcr_indexlist()
        row_length (integer): the expected number of elements in each row
    Return:
        test result (boolean): True if the file passes, False if it fails.
    """
    # The header is not stored in pcr_rxn_list.
    num_row = len(pcr_rxn_list)

    if(row_length ==  3):
        i_p5_index = 1
        i_p7_index = 2
    elif(row_length == 5):
        i_p5_name  = 1
        i_p5_index = 2
        i_p7_name  = 3
        i_p7_index = 4
    else:
        print('Error: unexpected column count in file "%s"' % (file_path), file=sys.stderr)
        return(False)

    errorFlag = 0

    # Check that each PCR rxn well name occurs no more than once.
    tmp_dict = dict()
    for i in range(num_row):
        key = pcr_rxn_list[i][0]
        tmp_dict[key] = tmp_dict.setdefault(key, 0) + 1
        if(tmp_dict[key] > 1):
            print('Error: PCR rxn well name (column 1) "%s" occurs more than once in file "%s"' % (pcr_rxn_list[i][0], file_path), file=sys.stderr)
            errorFlag = 1

    # Check that if one P5 index entry is 'none', all are 'none'.
    num_none = 0
    for i in range(num_row):
        if(pcr_rxn_list[i][i_p5_index] == 'none'):
            num_none += 1
    if(num_none != 0 and num_none != num_row):
        print('Error: P5 index must not be a mix of "none" and not "none" in file "%s".' % (file_path), file=sys.stderr)
        errorFlag = 1

    # Check that if one P7 index entry is 'none', all are 'none'.
    num_none = 0
    for i in range(num_row):
        if(pcr_rxn_list[i][i_p7_index] == 'none'):
            num_none += 1
    if(num_none != 0 and num_none != num_row):
        print('Error: P7 index must not be a mix of "none" and not "none" in file "%s".' % (file_path), file=sys.stderr)
        errorFlag = 1

    # Check that a P5 and P7 sequence pair is not repeated.
    tmp_dict = dict()
    for i in range(num_row):
        if(pcr_rxn_list[i][i_p5_index] == 'none' or pcr_rxn_list[i][i_p7_index] == 'none'):
          continue
        key = (pcr_rxn_list[i][i_p5_index], pcr_rxn_list[i][i_p7_index])
        tmp_dict[key] = tmp_dict.setdefault(key, 0) + 1
        if(tmp_dict[key] > 1):
            print('Error: P5 and P7 index sequence pair %s occurs more than once in file "%s".' % (key, file_path), file=sys.stderr)
            errorFlag = 1

    if(row_length == 5):
        # Check that each P5 well id matches one sequence value.
        tmp_dict = dict()
        for i in range(num_row):
            key = pcr_rxn_list[i][i_p5_name]
            value = pcr_rxn_list[i][i_p5_index]
            if(tmp_dict.get(key) != None and tmp_dict[key] != value):
                print('Error: P5 well id "%s" matches more than one sequence in file "%s".' % (key, file_path), file=sys.stderr)
            else:
                tmp_dict[key] = value

        # Check that each P5 sequence matches one well id.
        tmp_dict = dict()
        for i in range(num_row):
            key = pcr_rxn_list[i][i_p5_index]
            if(key == 'none'):
              continue
            value = pcr_rxn_list[i][i_p5_name]
            if(tmp_dict.get(key) != None and tmp_dict[key] != value):
                print('Error: P5 index sequence "%s" matches more than one well id in file "%s".' % (key, file_path), file=sys.stderr)
            else:
                tmp_dict[key] = value

        # Check that each P7 well id matches one sequence value.
        tmp_dict = dict()
        for i in range(num_row):
            key = pcr_rxn_list[i][i_p7_name]
            value = pcr_rxn_list[i][i_p7_index]
            if(tmp_dict.get(key) != None and tmp_dict[key] != value):
                print('Error: P7 well id "%s" matches more than one sequence in file "%s".' % (key, file_path), file=sys.stderr)
            else:
                tmp_dict[key] = value

        # Check that each P7 sequence matches one well id.
        tmp_dict = dict()
        for i in range(num_row):
            key = pcr_rxn_list[i][i_p7_index]
            if(key == 'none'):
              continue
            value = pcr_rxn_list[i][i_p7_name]
            if(tmp_dict.get(key) != None and tmp_dict[key] != value):
                print('Error: P7 index sequence "%s" matches more than one well id in file "%s".' % (key, file_path), file=sys.stderr)
            else:
                tmp_dict[key] = value

        # Check that P5 and P7 well id pairs are not repeated.
        tmp_dict = dict()
        for i in range(num_row):
            key = (pcr_rxn_list[i][i_p5_name], pcr_rxn_list[i][i_p7_name])
            tmp_dict[key] = tmp_dict.setdefault(key, 0) + 1

        for i in range(num_row):
            key = (pcr_rxn_list[i][i_p5_name], pcr_rxn_list[i][i_p7_name])
            if(tmp_dict[key] > 1):
                print('Error: P5 and P7 well id pair %s occurs more than once in file %s.' % (key, file_path), file=sys.stderr)
                errorFlag = 1

    return(False if (errorFlag > 0) else True)


def load_pcr_indexlist(file_path):
    """
    Load PCR index list into a (Python) list.
    PCR index list file headers:

      pcr_rxn_name,p5_index,p7_index
    or
      pcr_rxn_name,p5_well,p5_index,p7_well,p7_index

    In the latter case, the user can specify well ids for
    use in the cell read names rather than the index
    sequences.

    The pcr_rxn_name is printable text restricted to the
    characters in the regex '^[-a-zA-Z0-9_]+$'.

    The p5_well and p7_well name format is

      <plate_id>:<well_id>

    where the plate_id has the format

      Pnn

    where nn are decimal digits. The plate format is

      <row><column>

    where <row> is [A-H] and <column> is a decimal
    number in the range 01-12; that is, A01-H12.

    The index is a DNA sequence in lower or upper case
    but will be stored as upper case. The sequence can
    also be given as the string 'none'. If any entry,
    in the column is 'none', all must be.

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

        for i in range(len(row_list)):
            row_list[i] = row_list[i].strip()

        if(check_pcr_rxn_header(row_list) == 'long_header'):
            row_length = 5
        elif(check_pcr_rxn_header(row_list) == 'short_header'):
            row_length = 3
        else:
            raise ValueError('Missing header in file %s' % (file_path))

        for row_list in csvreader:
            i_row += 1
            for i in range(len(row_list)):
                row_list[i] = row_list[i].strip()

            if(check_pcr_rxn_row(file_path, i_row, row_length, row_list) == False):
                sys.exit(-1)

            pcr_rxn_list.append(row_list)

    # Check file for various types of errors.
    if(check_pcr_rxn_file(file_path, pcr_rxn_list, row_length) == False):
        sys.exit(-1)

    # Set all alphabetic characters to upper case.
    for i_row in range(len(pcr_rxn_list)):
        for i_col in range(row_length):
            pcr_rxn_list[i_row][i_col] = pcr_rxn_list[i_row][i_col].upper()

    # Set 'none' values to lower case.
    if(row_length == 3):
        i_p5_index = 1
        i_p7_index = 2
    else:
        i_p5_index = 2
        i_p7_index = 4

    for i in range(len(pcr_rxn_list)):
        if(pcr_rxn_list[i][i_p5_index] == 'NONE'):
            pcr_rxn_list[i][i_p5_index] = 'none'
        if(pcr_rxn_list[i][i_p7_index] == 'NONE'):
            pcr_rxn_list[i][i_p7_index] = 'none'

    return(pcr_rxn_list)


def make_pcr_whitelist(pcr_rxn_list, pcr_index, index_length, reverse_complement=False, well_ids=False):
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

    if(not isinstance(reverse_complement, bool)):
        raise ValueError('reverse_complement argument must be boolean.')

    if(not isinstance(well_ids, bool)):
        raise ValueError('well_ids argument must be boolean.')

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
    pcr_names = dict() # Check for loss of identity below.
    for pcr_rxn in pcr_rxn_list:
        if(not reverse_complement):
            key = pcr_rxn[i_seq][0:index_length]
            whitelist[key] = pcr_rxn[i_name]
            pcr_names.setdefault(key, []).append(pcr_rxn[i_seq])
        else:
            tmp_index = reverse_complement_sequence(pcr_rxn[i_seq])
            key = tmp_index[0:index_length]
            whitelist[key] = pcr_rxn[i_name]
            pcr_names.setdefault(key, []).append(pcr_rxn[i_seq])

    # Check that PCR sequence trimming does not result in
    # sequences that are not distinct.
    errorFlag = False
    for key in pcr_names.keys():
        if(len(set(pcr_names[key])) != 1):
            print('Error: trimmed %s PCR primers "%s" are not distinct after' % (pcr_index.upper(), set(pcr_names[key])), file=sys.stderr)
            print('       %strimming to %d bases.' % ('reverse complementing and ' if reverse_complement else '', index_length), file=sys.stderr)
            errorFlag = True
    if(errorFlag):
       sys.exit(-1)

    return(whitelist)


def get_programmed_pcr_combos_pcrlist(pcr_rxn_list, well_ids=False):
    """
    Make a set of valid PCR primer name pairs. The pairs are stored
    as a set.
    Args:
        pcr_rxn_list (list): PCR reactions from load_pcr_indexlist()
        well_ids (boolean): use wells to identify primer; otherwise, use sequence
    Return:
        set of tuples of valid PCR primer pair names
    """
    row_length = len(pcr_rxn_list[0])
    if(row_length == 3):
        i_p5_name  = 1 
        i_p5_index = 1 
        i_p7_name =  2 
        i_p7_index = 2 
    elif(row_length == 5):
        i_p5_name  = 1 if well_ids else 2
        i_p5_index = 2 
        i_p7_name =  3 if well_ids else 4
        i_p7_index = 4

    valid_combos = set()
    for pcr_rxn in pcr_rxn_list:
        p5_name = pcr_rxn[i_p5_name] if pcr_rxn[i_p5_index] != 'none' else 'none'
        p7_name = pcr_rxn[i_p7_name] if pcr_rxn[i_p7_index] != 'none' else 'none'
        valid_combos.add((p5_name, p7_name))

    return(valid_combos)


def pcr_index_list_is_none(file_path, pcr_rxn_list, p5_p7='p5'):
    """
    Return whether the P5 or P7 indices are 'none' in the
    pcr_rxn_list.
    Args:
        pcr_rxn_list (list): PCR reactions from load_pcr_indexlist()
        p5_p7 (string): Check either the 'p5' or 'p7' values for
                        'none'. We test earlier for the required
                        condition that if one index is 'none', all
                        must be 'none'.
    Return
        A boolean that reflects whether the values in the specified
        column are all 'none' (True is returned) or none of the values
        are 'none' (False is returned).
    """
    row_length = len(pcr_rxn_list[0])
    if(row_length == 3):
        i_p5_index = 1
        i_p7_index = 2
    elif(row_length == 5):
        i_p5_index = 2
        i_p7_index = 4
    else:
        print('Error: there must be either 3 or 5 columns in "%s"' % (file_path), file=sys.stderr)
        sys.exit(-1)

    if(p5_p7 == 'p5' and pcr_rxn_list[0][i_p5_index] == 'none'):
      return(True)
    elif(p5_p7 == 'p7' and pcr_rxn_list[0][i_p7_index] == 'none'):
      return(True)

    return(False)


# Stand-alone testing.
if __name__ == '__main__':

    pcr_rxn_list = load_pcr_indexlist('pcr_index.csv')
    print('load succeeds')
    print(pcr_rxn_list)

    well_ids = True
    index_length=5
    reverse_complement_p5=True

    print('p5 whitelist')
    print(make_pcr_whitelist(pcr_rxn_list, 'p5', index_length, reverse_complement=reverse_complement_p5, well_ids=well_ids))
    print('p7 whitelist')
    print(make_pcr_whitelist(pcr_rxn_list, 'p7', index_length, reverse_complement=False, well_ids=well_ids))
    print('pcr combinations')
    print(get_programmed_pcr_combos_pcrlist(pcr_rxn_list, well_ids=well_ids))

