#!/usr/bin/env python    

# Script for generating the required sample sheets for the BBI sci-RNA pipelines
# from a single input sheet

import argparse
import sys
import os
import codecs


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script for generating the required sample sheets for the BBI sci-RNA pipelines from a single input csv.')
    parser.add_argument('input_csv', help='Path to input sheet. Input sheet should contain the following columns in order, with no header: sample name, project name (for wrap), species, tissue (for Garnett), RT barcode, investigator ID.')
    parser.add_argument('--out_path', required=False, default = ".", help='Folder to output to. Default is current directory.')
    parser.add_argument('--no_map', action='store_true', default=False, required=False, help='If included, no SampleMap.csv will be generated.')
    parser.add_argument('--garnett_map', default="garnett_classifiers_map.csv", required=False, help='Map of Garnett models to apply given a species and tissue. csv with three columns and no header - species, tissue, path to classifier. Default is provided map.')
    args = parser.parse_args()

    # Load data
    sample_dict = dict()
    garnett = False
    with codecs.open(args.input_csv, "r", encoding="utf-8-sig") as incsv:
        for line in incsv:
            line = line.strip().split(',')
            tiss = line[3].split(":")
            tiss = [x.capitalize() for x in tiss]
            if line[0] in sample_dict:
                sample_dict[line[0]]['rt_barc'].add(line[4])
            else:
                samp = {'sample_id' : line[0], 'project_id' : line[1].replace('/', '.'), 'species' : line[2].capitalize(), 
                        'tissue' : tiss, 'rt_barc' : {line[4]}, 'inv_id' : line[5]}
                if samp['tissue'] != [""]:
                    garnett = True
                sample_dict[line[0]] = samp
    
    # Generate SampleSheet.csv
    sampsheet = open(os.path.join(args.out_path, "SampleSheet.csv"), 'w')
    sampsheet.write('RT Barcode,Sample ID,Reference Genome\n')
    for key, value in sample_dict.items():
        for rt in value['rt_barc']:
            ent = rt + "," + key + "," + value['species'] + "\n"
            sampsheet.write(ent)
    sampsheet.close()
    # Generate SampleMap.csv and SampleIDMap.csv
    if not args.no_map:
        proj_dict = {}
        sampmap = open(os.path.join(args.out_path, "SampleMap.csv"), 'w')
        sampmap.write('BAT Sample ID,Group,Organism,Fastqs\n')
        for key, value in sample_dict.items():
            if value['project_id'] != "":
                ent = key + "," + value['project_id'] + "," + value['species'] + ",\n"
                sampmap.write(ent)
                if value['project_id'] in proj_dict:
                    proj_dict[value['project_id']].append((value['sample_id'], value['inv_id']))
                else:
                    proj_dict[value['project_id']] = [(value['sample_id'], value['inv_id'])]
        sampmap.close()
        for key, value in proj_dict.items():
            projmap = open(os.path.join(args.out_path, key + "_SampleIDMap.csv"), 'w')
            projmap.write("BAT Sample ID,Sample ID\n")
            for entry in value:
                projmap.write(entry[0] + "," + entry[1] + "\n")
            projmap.close()
    
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
            for tiss in value['tissue']:
                if (value['species'], tiss) in garnett_dict:
                    for item in garnett_dict[(value['species'], tiss)]:
                        ent = value['sample_id'] + "," + item + "\n"
                        garnettsheet.write(ent)
    garnettsheet.close()









