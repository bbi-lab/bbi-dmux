#!/usr/bin/env python3

#
# Convert a LIMS CSV manifest file to a bbi-dmux samplesheet file.
# The manifest file format is set on 20231025. It may change in
# the future.
#
# This a minimal script at this time.
#


import sys
import argparse
import csv
import re



def check_header(row):
  if(row[0] != 'name' or
     row[1] != 'experiment' or
     row[2] != 'plate' or
     row[3] != 'coordinates' or
     row[4] != 'specimen' or
     row[5] != 'organism' or
     row[6] != 'genome'):
    print('Error: unexpected header token in row: ', row, file=sys.stderr)
    sys.exit(-1)


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to make a samplesheet CSV file for the bbi-dmux pipeline from a LIMS CSV experiment manifest file.')
  parser.add_argument('-i', '--in_file', required=True, help='Input LIMS manifest CSV filename.')
  parser.add_argument('-o', '--out_file', required=False, default='SampleSheet.csv', help='Output CSV filename (default is SampleSheet.csv).')
  args = parser.parse_args()

  infile = args.in_file
  outfile = args.out_file

  print('Input file: %s' % (infile))
  print('Output file: %s' % (outfile))

  pobj = re.compile('^([A-H])([0-9]+)$')

  ofp = open(outfile, 'w+', newline='')
  writer = csv.writer(ofp, dialect='unix', quoting=csv.QUOTE_NONE) 

  # Write header
  writer.writerow(['RT Barcode','Sample ID','Reference Genome'])

  # Read and write well rows.
  with open(infile, newline='') as ifp:
    reader = csv.reader(ifp, dialect='unix')
    # Read rows from input file.
    for irow, row in enumerate(reader):
      if(irow == 0):
        check_header(row)
        continue
      # Convert well A1 to well A01, etc.
      mobj = re.match(pobj, row[3])
      if(mobj == None):
        print('Error: no match to well \'%s\'' % (row[3]), file=sys.stderr)
        sys.exit(-1)
#      print('%s-%s%02d,%s,%s' % (row[2], mobj.group(1), int(mobj.group(2)), row[4], row[6]))
      # Write row to output file.
      writer.writerow(['%s-%s' % (row[2], mobj.group(1)+'%02d' % int(mobj.group(2))), row[4], row[6]])

  ofp.close()
