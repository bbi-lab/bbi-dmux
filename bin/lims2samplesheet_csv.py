#!/usr/bin/env python3

#
# Convert a LIMS CSV manifest file to a bbi-dmux samplesheet file.
# The manifest file format is set on 10-2025. It may change in
# the future.
#
# This a minimal script at this time.
#


import sys
import argparse
import csv
import re
import operator


# "name","experiment","assay","plate","coordinates","RT_block","BBI_ID","investigatorSpecimenId","investigator","organism","tissue","genome"
# "RNA3-074-P01-A1-24.0262","RNA3-074","sci-RNA-seq","P01","A1","RNA3-074_24.0262","24.0262","RA-gastruloid-2day","Hamazaki","Human","Gastruloid","Human"
# "RNA3-074-P01-A10-24.0270","RNA3-074","sci-RNA-seq","P01","A10","RNA3-074_24.0270","24.0270","RA-gastruloid-36h","Hamazaki","Human","Gastruloid","Human"
# "RNA3-074-P01-A11-24.0270","RNA3-074","sci-RNA-seq","P01","A11","RNA3-074_24.0270","24.0270","RA-gastruloid-36h","Hamazaki","Human","Gastruloid","Human"



# Check the LIMS manifest CSV file header.
def check_header(row):
  if(row[0] !=  'name' or
     row[1] !=  'experiment' or
     row[2] !=  'assay' or
     row[3] !=  'plate' or
     row[4] !=  'coordinates' or
     row[5] !=  'RT_block' or
     row[6] !=  'BBI_ID' or
     row[7] !=  'investigatorSpecimenId' or
     row[8] !=  'investigator' or
     row[9] !=  'organism' or
     row[10] != 'tissue' or
     row[11] != 'genome'):
    print('Error: unexpected header token in row: ', row, file=sys.stderr)
    sys.exit(-1)



# Read the LIMS manifest CSV file.
def read_lims_manifest(rows):
  rows = []
  with open(infile, newline='') as ifp:
    reader = csv.reader(ifp, dialect='unix')
    # Read rows from input file.
    for irow, row in enumerate(reader):
      if(irow == 0):
        check_header(row)
        continue
      rows.append(row)
  return(rows)


# Format fields in inrows and store in outrows.
# Return outrows.
def format_rows(inrows):
  pobj1 = re.compile('^P([0-9]+)$')
  pobj2 = re.compile('^([A-H])([0-9]+)$')
  outrows = []
  for inrow in inrows:
    mobj1 = re.match(pobj1, inrow[3])
    if(mobj1 == None):
      print('Error: no match to plate \'%s\'' % (inrow[3]), file=sys.stderr)
      sys.exit(-1)
    mobj2 = re.match(pobj2, inrow[4])
    if(mobj2 == None):
      print('Error: no match to well \'%s\'' % (inrow[5]), file=sys.stderr)
      sys.exit(-1)

#    plate_str              = inrow[2]
    plate_str              = 'P' + '%d' % (int(mobj1.group(1)))
    well_str               = '%s%02d' % (mobj2.group(1), int(mobj2.group(2)))
    sample_str             = str(inrow[5])
    investigator_sample_id = inrow[7]
    tstr                   = inrow[8].split(':')
    if(len(tstr) > 2):
      print('More than one semicolon in investigator string \'%s\'' % (inrow[6]))
    investigator_str       = inrow[8].split(':')[0]
    tissue_str             = inrow[10]
    genome_str             = inrow[11]
    outrow = ['%s%s' % (plate_str, well_str), sample_str, investigator_str, genome_str, tissue_str, '%s-%s' % (plate_str, well_str), investigator_sample_id]
#   0    sample
#   1    investigator
#   2    genome
#   3    tissue
#   4    plate-well
#   5    samples
    outrows.append(outrow)
  return(outrows)


# Write samplesheet CSV file to outfile.
def write_sample_sheet(outrows, outfile):
  with open(outfile, 'w+', newline='') as ofp:
    writer = csv.writer(ofp, dialect='unix', quoting=csv.QUOTE_NONE)
#    writer.writerow(['RT Barcode','Sample ID','Reference Genome'])
    for outrow in outrows:
      writer.writerow([outrow[1], outrow[2], outrow[3], outrow[4], outrow[5], outrow[6]])


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to make a samplesheet CSV file for the bbi-dmux pipeline from a LIMS CSV experiment manifest file.')
  parser.add_argument('-i', '--in_file', required=True, help='Input LIMS manifest CSV filename.')
  parser.add_argument('-o', '--out_file', required=False, default='samplesheet_spreadsheet.csv', help='Output CSV filename (default is SampleSheet.csv).')
  args = parser.parse_args()

  infile = args.in_file
  outfile = args.out_file

  print('Input file: %s' % (infile))
  print('Output file: %s' % (outfile))

  # Read input file.
  inrows = read_lims_manifest(infile)

  # Set up out rows.
  outrows = format_rows(inrows)

  # Sort rows by key in outrows[i][0]
  outrows.sort(key=operator.itemgetter(0))

  # Write samplesheet file.
  write_sample_sheet(outrows, outfile)

