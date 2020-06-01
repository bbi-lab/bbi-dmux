#!/usr/bin/env python    

# Script for generating the required sample sheets for the BBI sci-RNA pipelines
# from a single input sheet

import argparse
import sys
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script for generating the required sample sheets for the BBI sci-RNA pipelines from a single input csv.')
    parser.add_argument('input_csv', help='Path to input sheet. Input sheet should contain the following columns in order, with no header: sample name, project name (for wrap), species, tissue (for Garnett), RT barcode.')
    parser.add_argument('--out_path', required=False, default = ".", help='Folder to output to. Default is current directory.')
    parser.add_argument('--no_map', action='store_true', default=False, required=False, help='If included, no SampleMap.csv will be generated.')
    parser.add_argument('--garnett_map', default="garnett_map.csv", required=False, help='Map of Garnett models to apply given a species and tissue. csv with three columns and no header - species, tissue, path to classifier. Default is provided map.')
    args = parser.parse_args()

    # Load data
    sample_dict = dict()
    garnett = False
    with open(args.input_csv, 'r') as incsv:
        for line in incsv:
            line = line.strip().split(',')
            samp = {'sample_id' : line[0], 'project_id' : line[1], 'species' : line[2].capitalize(), 
                    'tissue' : line[3].capitalize(), 'rt_barc' : line[4]}
            if samp['tissue'] != "":
                garnett = True
            sample_dict[line[0]] = samp
    
    # Generate SampleSheet.csv
    sampsheet = open(os.path.join(args.out_path, "SampleSheet.csv"), 'w')
    sampsheet.write('RT Barcode,Sample ID,Reference Genome\n')
    for key, value in sample_dict.items():
        ent = value['rt_barc'] + "," + key + "," + value['species'] + "\n"
        sampsheet.write(ent)
    sampsheet.close()

    # Generate SampleMap.csv
    if not args.no_map:
        sampmap = open(os.path.join(args.out_path, "SampleMap.csv"), 'w')
        sampmap.write('BAT Sample ID,Group,Organism,Fastqs\n')
        for key, value in sample_dict.items():
            if value['project_id'] != "":
                ent = key + "," + value['project_id'] + "," + value['species'] + ",\n"
                sampmap.write(ent)
        sampmap.close()
    
    # Generate GarnettSheet.csv
    garnettsheet = open(os.path.join(args.out_path, "GarnettSheet.csv"), 'w')
    if garnett:
        garnett_dict = dict()
        with open(args.garnett_map, 'r') as garmap:
            for line in garmap:
                line = line.strip().split(',')
                tissue_key = (line[0].capitalize(), line[1].capitalize())
                if tissue_key in garnett_dict:
                    garnett_dict[tissue_key].append(line[2])
                else:
                    garnett_dict[tissue_key] = [line[2]]

        for key, value in sample_dict.items():
            if (value['species'], value['tissue']) in garnett_dict:
                for item in garnett_dict[(value['species'], value['tissue'])]:
                    ent = value['sample_id'] + "," + item + "\n"
                    garnettsheet.write(ent)
    garnettsheet.close()









