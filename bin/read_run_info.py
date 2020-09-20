#!/usr/bin/python3

#
# This script is used in
#    o  bbi-dmux/data/illumina_run_parameters/
#    o  bbi-sciatac-demux/main.nf
#

import sys
import json
import run_info

if __name__ == '__main__':
    pathname = sys.argv[1]
    run_info_json = run_info.get_run_info( pathname )
    # at least some Java JsonSlurper.parseText() method versions refuse JSON files made with 'indent=4')
    json.dump( run_info_json, sys.stdout )
