#!/bin/bash

#
# Path to the Nextflow processing run configuration file.
#
PWD=`pwd`
CONFIG_FILE="$PWD/experiment.config"

#
# Nextflow executable and pipeline script locations.
# Note: set the paths in the three variables below.
#
NEXTFLOW="$HOME/bin/nextflow"
NF_HOME="$HOME/git/bbi-dmux"
NF_MAIN="${NF_HOME}/main.nf"

#
# Current date and time.
#
NOW=`date '+%Y%m%d_%H%M%S'`

#
# Get the path to the demux output directory from
# the configuration file and set the Nextflow work
# directory to be in the demux output directory.
#
OUTPUT_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.output_dir"){print$2}}' | sed 's/"//g'`
DEMUX_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.demux_out"){print$2}}' | sed 's/"//g'`
WORK_DIR="$DEMUX_DIR/work_demux"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
# REPORT_FIL=$DEMUX_DIR/demux.report.${NOW}.html
TRACE_FIL=$DEMUX_DIR/demux.trace.${NOW}.tsv
# TIMELINE_FIL=$DEMUX_DIR/demux.timeline.${NOW}.html

#
# Nextflow run parameters.
#
# PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-trace $TRACE_FIL -resume"

mkdir -p $DEMUX_DIR
pushd $DEMUX_DIR

date > ./run_start.${NOW}.txt

#
# Run Nextflow sci-RNA demux pipeline.
#
$NEXTFLOW run $NF_MAIN $PARS

date > ./run_finish.${NOW}.txt

popd

#
# Set run directory file and directory permissions.
#
${NF_HOME}/scripts/set_run_permissions.sh ${OUTPUT_DIR}

