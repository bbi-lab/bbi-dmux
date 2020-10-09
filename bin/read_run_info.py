#!/usr/bin/python3

#
# Notes:
#   o  the master file is bbi-dmux-data/illumina_run_parameters/src/run_info.py
#      Make changes to the master file and distribute
#      to all other pipelines that use it in order to
#      maintain consistency.
#   o  this script is used in
#        o  bbi-dmux-data/illumina_run_parameters/src
#        o  bbi-dmux/main.nf
#        o  bbi-sciatac-demux/main.nf
#

import sys
import json
import run_info

if __name__ == '__main__':
    pathname = sys.argv[1]
    run_info_json = run_info.get_run_info( pathname )
    # at least some Java JsonSlurper.parseText() method versions refuse JSON files made with 'indent=4')
    json.dump( run_info_json, sys.stdout )
