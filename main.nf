import groovy.json.JsonSlurper

/*
** Where to find scripts.
** Note: script_dir needs to be visible within Groovy functions
**       so there is no 'def', which makes it global.
*/
pipeline_path="$workflow.projectDir"
script_dir="${pipeline_path}/bin"


/*
** Check that Nextflow version meets minimum version requirements.
*/
def minMajorVersion = 20
def minMinorVersion = 07
checkNextflowVersion( minMajorVersion, minMinorVersion )


/*
** Check OS version.
** Notes:
**   o  works only for Linux systems
**   o  used to distinguish between CentOS 6 and CentOS 7
*/
( osName, osDistribution, osRelease ) = getOSInfo()


// Parse input parameters
params.help = false
params.rerun = false
params.star_file = "$baseDir/bin/star_file.txt"
params.level = 3
params.bcl_max_mem = 40
params.fastq_chunk_size = 100000000
params.run_recovery = false
params.rt_barcode_file="default"
params.large = false
params.generate_samplesheets = 'no_input'
params.max_cores = 16
params.max_wells_per_sample = 20

params.p5_cols = 0
params.p7_rows = 0
params.p5_wells = 0
params.p7_wells = 0
// Note the p5 and p7 index lengths must have default value 'false' in order
// to use the values in RunParameters.xml, which is the typical case.
params.p5_index_length = false
params.p7_index_length = false


//print usage
if (params.help) {
    log.info ''
    log.info 'BBI sci-RNA-seq Demultiplexer'
    log.info '--------------------------------'
    log.info ''
    log.info 'For reproducibility, please specify all parameters to a config file'
    log.info 'by specifying -c CONFIG_FILE.config.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run bbi-dmux -c CONFIG_FILE'
    log.info ''
    log.info 'Help: '
    log.info '    --help                                     Show this message and exit.'
    log.info ''
    log.info 'Required parameters (specify in your config file):'
    log.info '    params.run_dir = RUN_DIRECTORY             Path to the sequencer output.'
    log.info '    params.output_dir OUTPUT DIRECTORY         Output directory.'
    log.info '    params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.'
    log.info '    params.level = 3                           Level of run - either 2 or 3.'
    log.info ''
    log.info 'Required parameters (one of the pairs below is required - p7_wells and p5_wells or p7_rows and p5_cols):'
    log.info '    params.p7_wells = "A1 B1 C1"               Alternative to p7_rows and p5_cols - specify specific PCR wells instead of full rows/columns. Must match order of params.p5_wells.'
    log.info '    params.p5_wells = "A1 A2 A3"               Alternative to p7_rows and p5_cols - specify specific PCR wells instead of full rows/columns. Must match order of params.p7_wells.'
    log.info '    params.p7_rows = "A B C"                   The PCR rows used - must match order of params.p5_cols.'
    log.info '    params.p5_cols = "1 2 3"                   The PCR columns used - must match order of params.p7_rows.'
    log.info ''
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.rt_barcode_file = "default"         The path to a custom RT barcode file. If "default", default BBI barcodes will be used.'
    log.info '    params.p5_index_length = false             The P5 index length if it is not the value in the Illumina RunParameters.xml file.'
    log.info '    params.p7_index_length = false             The P7 index length if it is not the value in the Illumina RunParameters.xml file.'
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info '    params.star_file = PATH/TO/FILE            File with the genome to star maps, similar to the one included with the package.'
    log.info '    params.fastq_chunk_size = 100000000        The number of reads that should be processed together for demultiplexing.'
    log.info '    params.bcl_max_mem = 40                    The maximum number of GB of RAM to assign for bcl2fastq'
    log.info '    params.large = false                       Is this a very large run? If true, the fastqs will be split - note that for smaller runs this will make the pipeline run more slowly.'
    log.info '    params.max_wells_per_sample = 20           The maximum number of wells per sample - if a sample is in more wells, the fastqs will be split then reassembled.'
    log.info '    --run_recovery true                        Add this to run the recovery script AFTER running the normal pipeline.'
    log.info '    --generate_samplesheets input_csv          Add this to generate the necessary samplesheet from the BBI universal input sheet.'    
    log.info ''
    log.info 'Issues? Contact hpliner@uw.edu'
    exit 1
}

/*
** Read Illumina run information.
*/
illuminaRunInfoMap = readIlluminaRunInfo( params )

/*
** Set p5_index_length and p7_index_length.
*/
if( params.p5_index_length == false ) {
  p5_index_length = illuminaRunInfoMap['p5_index_length']
}
else {
  p5_index_length = params.p5_index_length
}

if( params.p7_index_length == false ) {
  p7_index_length = illuminaRunInfoMap['p7_index_length']
}
else {
  p7_index_length = params.p7_index_length
}


process generate_sheets {
    publishDir path: "${params.output_dir}", pattern: "SampleSheet.csv", mode: 'copy'
    publishDir path: "${params.output_dir}", pattern: "SampleMap.csv", mode: 'copy'
    publishDir path: "${params.output_dir}", pattern: "GarnettSheet.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/sample_id_maps", pattern: "*_SampleIDMap.csv", mode: 'copy'

    input:
        file insamp from Channel.fromPath(params.generate_samplesheets)

    output:
        file "*Sheet.csv"
        file "SampleMap.csv" optional true
        file "*_SampleIDMap.csv" optional true

    when:
        params.generate_samplesheets != 'no_input'


    """
    generate_sample_sheets.py $params.generate_samplesheets
    """
}


// check required options
if (!params.run_dir || !params.output_dir || !params.sample_sheet ) {
    exit 1, "Must include config file using -c CONFIG_FILE.config that includes output_dir, sample_sheet and run_dir."
}

//if (!(params.p7_rows && params.p5_cols) && !(params.p7_wells && params.p5_wells)) {
//    exit 1, "Must include config file using -c CONFIG_FILE.config that includes p7_rows and p5_cols or p5_wells and p7_wells"
//}

star_file = file(params.star_file)

// check sample sheet
process check_sample_sheet {
    input:
	val params.sample_sheet
    file star_file

    output:
        file "*.csv" into good_sample_sheet

    when:
        params.generate_samplesheets == "no_input"

    """
    check_sample_sheet.py --sample_sheet $params.sample_sheet --star_file $star_file --level $params.level --rt_barcode_file $params.rt_barcode_file --max_wells_per_samp $params.max_wells_per_sample
    """
}

sample_sheet_file1 = good_sample_sheet
sample_sheet_file2 = good_sample_sheet
sample_sheet_file3 = good_sample_sheet
sample_sheet_file4 = good_sample_sheet
sample_sheet_file5 = good_sample_sheet

process make_sample_sheet {
    cache 'lenient'

    input:
        val params.run_dir
        file good_sample_sheet

    output:
        file "SampleSheet.csv" into bcl_samp_sheet

    when:
        !params.run_recovery

    """
    make_sample_sheet.py --run_directory $params.run_dir --p5_length $p5_index_length --p7_length $p7_index_length
    """    
}

// Run bcl2fastq
if (params.max_cores > 16) {
    max_cores_bcl = 16
    bcl_mem = params.bcl_max_mem/16
} else {
    max_cores_bcl = params.max_cores
    bcl_mem = params.bcl_max_mem/max_cores_bcl
}

process bcl2fastq {
    cache 'lenient'
    cpus max_cores_bcl
    memory "${bcl_mem} GB"

    input:
        file bcl_samp_sheet

    output:
        file "lane_fastqs" into bcl2fastq_output
        set file("lane_fastqs/Undetermined_S0_*_R1_001.fastq.gz"), file("lane_fastqs/Undetermined_S0_*_R2_001.fastq.gz") into fastqs mode flatten
        file "lane_fastqs/fake*.gz" optional true into fakes mode flatten

    """
    min_threads=\$((($max_cores_bcl/2)<4 ? ($max_cores_bcl/2):4))

    bcl2fastq -R $params.run_dir --output-dir ./lane_fastqs \
        --sample-sheet $bcl_samp_sheet \
        --loading-threads \$min_threads \
        --processing-threads $max_cores_bcl  \
        --writing-threads \$min_threads \
        --barcode-mismatches 1 \
        --ignore-missing-positions \
        --ignore-missing-controls \
        --ignore-missing-filter \
        --ignore-missing-bcls \
        --minimum-trimmed-read-length 15 \
        --mask-short-adapter-reads 15

    """
}

fastqs.into { fastqs_path1; fastqs_path2 }

/*
** ================================================================================
** PATH 1 - For small datasets, no fastq splitting
** ================================================================================
*/

process seg_sample_fastqs1 {
    cache 'lenient'

    publishDir path: "${params.output_dir}/", pattern: "demux_out/*fastq.gz", mode: 'link'     
    publishDir  path: "${params.output_dir}/demux_out/", pattern: "*.csv", mode: 'copy'
    publishDir  path: "${params.output_dir}/demux_out/", pattern: "*.json", mode: 'copy'

    input:
        set file(R1), file(R2) from fastqs_path1
        file sample_sheet_file1

    output:
        file "demux_out/*" into seg_output1
        file "demux_out/*.fastq.gz" into samp_fastqs_check1
        file "demux_out/*.stats.json" into json_stats1 mode flatten
        file "demux_out/*.csv" into csv_stats1
    
    when:
        !params.large
 
    """
    mkdir demux_out
    make_sample_fastqs.py --run_directory $params.run_dir \
        --read1 <(zcat $R1) --read2 <(zcat $R2) \
        --file_name $R1 --sample_layout $sample_sheet_file1 \
        --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
        --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
        --rt_barcode_file $params.rt_barcode_file \
        --p5_length $p5_index_length \
        --p7_length $p7_index_length \
        --output_dir ./demux_out --level $params.level
    pigz -p 8 demux_out/*.fastq
    """    
}


out_dir_str = params.output_dir.replaceAll("/\\z", "");
project_name = out_dir_str.substring(out_dir_str.lastIndexOf("/")+1);

process demux_dash1 {
    publishDir path: "${params.output_dir}/", pattern: "demux_dash", mode: 'copy'

    input:
        file demux_stats_csvs from csv_stats1.collect()
        file jsons from json_stats1.collect()
        file sample_sheet_file2
    output:
        file "demux_dash" into demux_dash1

    """
    mkdir demux_dash
    cp -R $baseDir/bin/skeleton_dash/* demux_dash/
    generate_html.R \
        "." --p7_rows "$params.p7_rows" --p5_cols "$params.p5_cols" --p7_wells "$params.p7_wells" --p5_wells "$params.p5_wells" --level "$params.level" --project_name "${project_name}" --sample_sheet "$sample_sheet_file2"

    """

}


/*
** ================================================================================
** PATH 2 - For large datasets, with fastq splitting
** ================================================================================
*/


process seg_sample_fastqs2 {
    cache 'lenient'

    input:
        set file(R1), file(R2) from fastqs_path2.splitFastq(by: params.fastq_chunk_size, file: true, pe: true)
        file sample_sheet_file3

    output:
        file "demux_out/*" into seg_output
        file "demux_out/*.fastq.gz" into samp_fastqs_check mode flatten
        file "demux_out/*.stats.json" into json_stats mode flatten
        file "demux_out/*.csv" into csv_stats mode flatten

    when:
        params.large

    """
    mkdir demux_out
    make_sample_fastqs.py --run_directory $params.run_dir \
        --read1 $R1 --read2 $R2 \
        --file_name $R1 --sample_layout $sample_sheet_file3 \
        --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
        --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
        --rt_barcode_file $params.rt_barcode_file \
        --p5_length $p5_index_length \
        --p7_length $p7_index_length \
        --output_dir ./demux_out --level $params.level
    pigz -p 8 demux_out/*.fastq    
    """    
}

get_prefix = { fname ->
    (fname - ~/_R1_001\.[0-9]+\.fastq.fastq.gz/)
}

samp_fastqs_check
    .map { file -> tuple(get_prefix(file.name), file) }
    .groupTuple()
    .set { grouped_files }

process recombine_fastqs {
    cache 'lenient'
    publishDir  path: "${params.output_dir}/demux_out", pattern: "*.fastq.gz", mode: 'move'

    input:
        set prefix, file(all_fqs) from grouped_files 

    output:
        file "*.gz" into gz_fqs 

    """
    cat $all_fqs > ${prefix}.fastq.gz 
    """
}

csv_prefix = { fname ->
    (fname - ~/_R1_001\.[0-9]+\.fastq\.[a-z]+_[a-z]+\.csv/)
}

csv_stats
    .map { file -> tuple(csv_prefix(file.name), file) }
    .groupTuple()
    .set { grouped_csvs }

process recombine_csvs {
    cache 'lenient'
    publishDir  path: "${params.output_dir}/demux_out/", pattern: "*.csv", mode: 'copy'

    input:
        set prefix, file(all_csvs) from grouped_csvs

    output:
        file "*.csv" into all_csv

    """
    csvs="$all_csvs"    
    arr=(\$csvs)
    if [ "$params.level" = "3" ]; then
        cat \$(IFS=\$'\n'; echo "\${arr[*]}" | grep lig_counts) | awk -F ',' 'BEGIN {OFS = ","} {a[\$1] += \$2} END {for (i in a) print i, a[i]}' > ${prefix}.lig_counts.csv    
    fi
    cat \$(IFS=\$'\n'; echo "\${arr[*]}" | grep rt_counts) | awk -F ',' 'BEGIN {OFS = ","} {a[\$1] += \$2} END {for (i in a) print i, a[i]}' > ${prefix}.rt_counts.csv
    cat \$(IFS=\$'\n'; echo "\${arr[*]}" | grep pcr_counts) | awk -F ',' 'BEGIN {OFS = ","; SUBSEP = OFS = FS} {a[\$1,\$2] += \$3} END {for (i in a) print i, a[i]}' > ${prefix}.pcr_counts.csv

    """
}

json_prefix = { fname ->
    (fname - ~/_R1_001\.[0-9]+\.fastq\.stats\.json/)
}

json_stats
    .map { file -> tuple(json_prefix(file.name), file) }
    .groupTuple()
    .set { grouped_jsons }

process recombine_jsons {
    cache 'lenient'
    publishDir  path: "${params.output_dir}/demux_out/", pattern: "*.json", mode: 'copy'

    input:
        set prefix, file(all_jsons) from grouped_jsons

    output:
        file "*.json" into all_json

    """
#!/usr/bin/env python
import json
import glob
import os
from collections import OrderedDict

if $params.level == 3:
    total_input_reads = 0
    total_passed_reads = 0
    total_uncorrected = 0
    total_ambiguous_ligation_length = 0
    total_unused_rt_well = 0
    total_pcr_mismatch = 0
    total_corrected_9 = 0
    total_corrected_10 = 0
    sample_read_counts = {}

    for file_name in glob.glob("*.json"):
        with open (file_name, "r") as read_file:
            data = json.load(read_file)
        total_input_reads += data['total_input_reads']
        total_passed_reads += data['total_passed_reads']
        total_uncorrected += (data['fraction_uncorrected_reads'] * data['total_input_reads'])
        total_ambiguous_ligation_length += (data['fraction_ambiguous_ligation_length'] * data['total_input_reads'])
        total_unused_rt_well += (data['fraction_invalid_rt_well'] * data['total_input_reads'])
        total_pcr_mismatch += (data['fraction_pcr_mismatch'] * data['total_input_reads'])
        total_corrected_9 += data['total_reads_corrected_when_9bp_ligation'] 
        total_corrected_10 += data['total_reads_corrected_when_10bp_ligation']
        sample_read_counts = { k: sample_read_counts.get(k, 0) + data['total_reads_passed_per_sample'].get(k, 0) for k in set(sample_read_counts) | set(data['total_reads_passed_per_sample']) }

    stats = OrderedDict()
    stats['total_input_reads'] = total_input_reads
    stats['total_passed_reads'] = total_passed_reads
    stats['fraction_passed_reads'] = total_passed_reads / total_input_reads
    stats['fraction_uncorrected_reads'] = total_uncorrected / total_input_reads
    stats['fraction_ambiguous_ligation_length'] = total_ambiguous_ligation_length / total_input_reads
    stats['fraction_invalid_rt_well'] = total_unused_rt_well / total_input_reads
    stats['fraction_pcr_mismatch'] = total_pcr_mismatch / total_input_reads
    stats['total_reads_corrected_when_9bp_ligation'] = total_corrected_9
    stats['total_reads_corrected_when_10bp_ligation'] = total_corrected_10
    stats['total_reads_passed_per_sample'] = sample_read_counts

if $params.level == 2:
    total_input_reads = 0
    total_passed_reads = 0
    total_uncorrected = 0
    total_ambiguous_ligation_length = 0
    total_unused_rt_well = 0
    total_pcr_mismatch = 0
    total_corrected = 0
    sample_read_counts = {}

    for file_name in glob.glob("*.json"):
        with open (file_name, "r") as read_file:
            data = json.load(read_file)
        total_input_reads += data['total_input_reads']
        total_passed_reads += data['total_passed_reads']
        total_uncorrected += (data['fraction_uncorrected_reads'] * data['total_input_reads'])
        total_unused_rt_well += (data['fraction_invalid_rt_well'] * data['total_input_reads'])
        total_pcr_mismatch += (data['fraction_pcr_mismatch'] * data['total_input_reads'])
        total_corrected += data['total_reads_corrected'] 
        sample_read_counts = { k: sample_read_counts.get(k, 0) + data['total_reads_passed_per_sample'].get(k, 0) for k in set(sample_read_counts) | set(data['total_reads_passed_per_sample']) }

    stats = OrderedDict()
    stats['total_input_reads'] = total_input_reads
    stats['total_passed_reads'] = total_passed_reads
    stats['fraction_passed_reads'] = total_passed_reads / total_input_reads
    stats['fraction_uncorrected_reads'] = total_uncorrected / total_input_reads
    stats['fraction_invalid_rt_well'] = total_unused_rt_well / total_input_reads
    stats['fraction_pcr_mismatch'] = total_pcr_mismatch / total_input_reads
    stats['total_reads_corrected'] = total_corrected
    stats['total_reads_passed_per_sample'] = sample_read_counts

with open("${prefix}.stats.json", 'w') as f:
    f.write(json.dumps(stats, indent=4))
    """
}

out_dir_str = params.output_dir.replaceAll("/\\z", "");
project_name = out_dir_str.substring(out_dir_str.lastIndexOf("/")+1);

process demux_dash {
    publishDir path: "${params.output_dir}/", pattern: "demux_dash", mode: 'copy'

    input:
        file demux_stats_csvs from all_csv.collect()
        file jsons from all_json.collect()
        file sample_sheet_file4
    output:
        file demux_dash

    """
    mkdir demux_dash
    cp -R $baseDir/bin/skeleton_dash/* demux_dash/
    generate_html.R \
        "." --p7_rows "$params.p7_rows" --p5_cols "$params.p5_cols" --p7_wells "$params.p7_wells" --p5_wells "$params.p5_wells" --level "$params.level" --project_name "${project_name}" --sample_sheet "$sample_sheet_file4"

    """

}


/*
** ================================================================================
** End path 2 - Start recovery
** ================================================================================
*/



save_recovery2 = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz-summary.txt/ + "-recovery_summary.txt"}
save_recovery = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz.txt/ + "-recovery_table.txt"}
process run_recovery {
    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery, pattern: "*.gz.txt", mode: 'link'
    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery2, pattern: "*-summary.txt", mode: 'link'

    input:
        file input from Channel.fromPath("${params.demux_out}/Undetermined*")
        file sample_sheet_file5

    output:
        file "*.txt"
        file "*summary.txt" into summaries
    when:
        params.run_recovery


    """
    recovery_script.py --input_file <(zcat $input) --output_file ${input}.txt \
        --run_directory $params.run_dir \
        --sample_layout $sample_sheet_file5 \
        --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
        --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
        --level $params.level \
        --rt_barcodes $params.rt_barcode_file

    """

}

process sum_recovery {
    publishDir path: "${params.output_dir}/demux_dash/js/", pattern: "recovery_summary.js", mode: 'move'

    input:
        file summary from summaries.collect()


    output: 
        file "*summary.js"

    """
   echo "const log_data = {" > recovery_summary.js
   for file in $summary
   do
     filename=\$(basename \$file); 
     part=\${filename/Undetermined-L00/};
     lane=\${part/.fastq.gz-summary.txt/};   
     printf "\$lane : \\`" >> recovery_summary.js;
     cat \$file >> recovery_summary.js;
     printf "\\`," >> recovery_summary.js;
   done
   echo "}" >> recovery_summary.js   


    """


}




workflow.onComplete {
    println "bbi-dmux completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}


/*************
Groovy functions
*************/

/*
** checkNextflowVersion
**
** Purpose: check Nextflow version information to minimum version values.
**
** Returns:
**   exits when Nextflow version is unacceptable
*/
def checkNextflowVersion( Integer minMajorVersion, Integer minMinorVersion )
{
  def sVersion = nextflow.version.toString()
  def aVersion = sVersion.split( /[.]/ )
  def majorVersion = aVersion[0].toInteger()
  def minorVersion = aVersion[1].toInteger()
  if( majorVersion < minMajorVersion || minorVersion < minMinorVersion )
  {
    def serr = "This pipeline requires Nextflow version at least %s.%s: you have version %s."
    println()
    println( '****  ' + String.format( serr, minMajorVersion, minMinorVersion, sVersion ) + '  ****' )
    println()
    System.exit( -1 )
    /*
    ** An exception produces an exceptionally verbose block of confusing text. I leave
    ** the command here in case the println() output is obscured by fancy Nextflow tables.
    **
    ** throw new Exception( String.format( serr, minMajorVersion, minMinorVersion, sVersion ) )
    */
  }
  return( 0 )
}


/*
** getOSInfo()
**
** Purpose: get information about the operating system.
**
** Returns:
**    list of strings with OS name, OS distribution, OS distribution release
**
** Notes:
**   o  limited to Linux operating systems at this time
*/
def getOSInfo()
{
  def osName = System.properties['os.name']
  def osDistribution
  def osRelease
  if( osName == 'Linux' )
  {
    def proc
    proc = "lsb_release -a".execute() | ['awk', 'BEGIN{FS=":"}{if($1=="Distributor ID"){print($2)}}'].execute()
    proc.waitFor()
    osDistribution = proc.text.trim()
    proc = "lsb_release -a".execute() | ['awk', 'BEGIN{FS=":"}{if($1=="Release"){print($2)}}'].execute()
    proc.waitFor()
    osRelease = proc.text.trim()
  }
  return( [ osName, osDistribution, osRelease ] )
}


/*
** readIlluminaRunInfo()
**
** Purpose: read Illumina RunParameters.xml file and return selected parameters.
**
** Returns:
**   Groovy 'map' that contains run info.
*/
def readIlluminaRunInfo( params ) {
    def command = "${script_dir}/read_run_info.py ${params.run_dir}"
    def strOut = new StringBuffer()
    def strErr = new StringBuffer()
    def proc = command.execute()
    def jsonSlurper = new JsonSlurper()

    proc.consumeProcessOutput(strOut, strErr)
    proc.waitForProcessOutput()
    if( proc.exitValue() != 0 ) {
        System.err << strErr.toString()
        System.exit( -1 )
    }
    illuminaRunInfoMap = jsonSlurper.parseText(strOut.toString())
    return( illuminaRunInfoMap )
}


