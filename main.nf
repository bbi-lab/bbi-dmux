
// Parse input parameters
params.help = false
params.rerun = false
params.star_file = "$baseDir/bin/star_file.txt"
params.level = 3
params.bcl_max_mem = 40
params.fastq_chunk_size = 100000000
params.run_recovery = false
params.rt_barcode_file="default"

params.p5_cols = 0
params.p7_rows = 0
params.p5_wells = 0
params.p7_wells = 0

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
    log.info '    --help                              Show this message and exit.'
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
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info '    params.rerun = [sample1, sample2]          Add to only rerun certain samples from trimming on.'
    log.info '    params.star_file = PATH/TO/FILE            File with the genome to star maps, similar to the one included with the package.'
    log.info '    params.fastq_chunk_size = 100000000        The number of reads that should be processed together for demultiplexing.'
    log.info '    params.bcl_max_mem = 40                    The maximum number of GB of RAM to assign for bcl2fastq'
    log.info '    --run_recovery true                        Add this to run the recovery script AFTER running the normal pipeline.'
    log.info ''
    log.info 'Issues? Contact hpliner@uw.edu'
    exit 1
}

// check required options
if (!params.run_dir || !params.output_dir || !params.sample_sheet ) {
    exit 1, "Must include config file using -c CONFIG_FILE.config that includes output_dir, sample_sheet and run_dir."
}

// check required options
if (!(params.p7_rows && params.p5_cols) && !(params.p7_wells && params.p5_wells)) {
    exit 1, "Must include config file using -c CONFIG_FILE.config that includes p7_rows and p5_cols or p5_wells and p7_wells"
}

star_file = file(params.star_file)

//check sample sheet
process check_sample_sheet {
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'

    input:
	val params.sample_sheet
        file star_file

    output:
        file "*.csv" into good_sample_sheet

    """
    check_sample_sheet.py --sample_sheet $params.sample_sheet --star_file $star_file --level $params.level --rt_barcode_file $params.rt_barcode_file
    """
}

sample_sheet_file = good_sample_sheet
sample_sheet_file2 = good_sample_sheet
sample_sheet_file3 = good_sample_sheet

process make_sample_sheet {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'

    input:
        val params.run_dir
        file good_sample_sheet

    output:
        file "SampleSheet.csv" into samp_sheet

    when:
        !params.run_recovery

    """
    make_sample_sheet.py --run_directory $params.run_dir

    """    
}

if (params.max_cores > 16) {
    max_cores_bcl = 16
    bcl_mem = params.bcl_max_mem/16
} else {
    max_cores_bcl = params.max_cores
    bcl_mem = params.bcl_max_mem/max_cores_bcl
}

process bcl2fastq {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:gmp/5.0.2'
    module 'mpfr/3.1.0:mpc/0.8.2:gcc/4.9.1:bcl2fastq/2.20'
//    publishDir path: "$params.output_dir", pattern: "lane_fastqs/Undetermined_S0_*.fastq.gz", mode: 'copy'
    penv 'serial'
    cpus max_cores_bcl
    memory "$bcl_mem" + " GB"    

    input:
        file samp_sheet

    output:
        file "lane_fastqs" into bcl2fastq_output
        set file("lane_fastqs/Undetermined_S0_*_R1_001.fastq.gz"), file("lane_fastqs/Undetermined_S0_*_R2_001.fastq.gz") into fastqs mode flatten
        file "lane_fastqs/fake*.gz" optional true into fakes mode flatten

    """
    min_threads=\$((($max_cores_bcl/2)<4 ? ($max_cores_bcl/2):4))

    bcl2fastq -R $params.run_dir --output-dir ./lane_fastqs \
        --sample-sheet $samp_sheet \
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

process seg_sample_fastqs {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
    memory '1 GB'    

    input:
        set file(R1), file(R2) from fastqs.splitFastq(by: params.fastq_chunk_size, file: true, pe: true)
        file sample_sheet_file

    output:
        file "demux_out/*" into seg_output
        file "demux_out/*.fastq.gz" into samp_fastqs_check mode flatten
        file "demux_out/*.stats.json" into json_stats mode flatten
        file "demux_out/*.csv" into csv_stats mode flatten

    """
    mkdir demux_out
    make_sample_fastqs.py --run_directory $params.run_dir \
        --read1 $R1 --read2 $R2 \
        --file_name $R1 --sample_layout $sample_sheet_file \
        --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
        --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
        --rt_barcode_file $params.rt_barcode_file \
        --output_dir ./demux_out --level $params.level
    gzip demux_out/*.fastq
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
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
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
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
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
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'
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
    module 'java/latest:modules:modules-init:modules-gs:gcc/8.1.0:R/3.6.1'
    memory '8 GB'    

    publishDir path: "${params.output_dir}/", pattern: "demux_dash", mode: 'copy'


    input:
        file demux_stats_csvs from all_csv.collect()
        file jsons from all_json.collect()
        file sample_sheet_file2
    output:
        file demux_dash

    """
    mkdir demux_dash
    cp -R $baseDir/bin/skeleton_dash/* demux_dash/
    generate_html.R \
        "." --p7_rows "$params.p7_rows" --p5_cols "$params.p5_cols" --p7_wells "$params.p7_wells" --p5_wells "$params.p5_wells" --level "$params.level" --project_name "${project_name}" --sample_sheet "$sample_sheet_file2"

    """

}

save_recovery2 = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz-summary.txt/ + "-recovery_summary.txt"}
save_recovery = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz.txt/ + "-recovery_table.txt"}
process run_recovery {
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'
    memory '4 GB'

    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery, pattern: "*.gz.txt", mode: 'move'
    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery2, pattern: "*-summary.txt", mode: 'move'
    input:
        file input from Channel.fromPath("${params.demux_out}/Undetermined*")
        file sample_sheet_file3

    output:
        file "*.txt"

    when:
        params.run_recovery


    """
    
    recovery_script.py --input_file <(zcat $input) --output_file ${input}.txt \
        --run_directory $params.run_dir \
        --sample_layout $sample_sheet_file3 \
        --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
        --level $params.level \
        --rt_barcodes $params.rt_barcode_file



    """



}





