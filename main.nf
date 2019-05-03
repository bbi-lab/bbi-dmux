
// Parse input parameters
params.help = false
params.rerun = false

//print usage
if (params.help) {
    log.info ''
    log.info 'BBI 2-level sci-RNA-seq Demultiplexer'
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
    log.info '    params.p7_rows = "A B C"                   The PCR rows used - must match order of params.p5_cols.'
    log.info '    params.p5_cols = "1 2 3"                   The PCR columns used - must match order of params.p7_rows.'
    log.info ''
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info '    params.rerun = [sample1, sample2]          Add to only rerun certain samples from trimming on.'
    log.info ''
    log.info 'Issues? Contact hpliner@uw.edu'
    exit 1
}

// check required options
if (!params.run_dir || !params.output_dir || !params.sample_sheet || !params.p7_rows || !params.p5_cols) {
    exit 1, "Must include config file using -c CONFIG_FILE.config that includes output_dir, sample_sheet, run_dir, p7_rows and p5_cols"
}

//rt2_file = Channel.fromPath('barcode_files/rt2.txt')
//check sample sheet
process check_sample_sheet {
    module 'modules:java/latest:modules-init:modules-gs:python/3.6.4'

    input:
	val params.sample_sheet

    output:
        file "*.csv" into good_sample_sheet

    """
    check_sample_sheet.py --sample_sheet $params.sample_sheet
    """
}

sample_sheet_file = good_sample_sheet

process make_sample_sheet {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:python/3.6.4'

    input:
        val params.run_dir
        file good_sample_sheet

    output:
        file "SampleSheet.csv" into samp_sheet

    """
    make_sample_sheet.py --run_directory $params.run_dir

    """    
}

if (params.max_cores > 16) {
    max_cores_bcl = 16
    bcl_mem = 2.5
} else {
    max_cores_bcl = params.max_cores
    bcl_mem = 40/max_cores_bcl
}

process bcl2fastq {
    cache 'lenient'
    module 'java/latest:modules:modules-init:modules-gs:gmp/5.0.2'
    module 'mpfr/3.1.0:mpc/0.8.2:gcc/4.9.1:bcl2fastq/2.20'
    publishDir path: "$params.output_dir", pattern: "lane_fastqs/Undetermined_S0_*.fastq.gz", mode: 'copy'
    clusterOptions "-pe serial $max_cores_bcl -l mfree=$bcl_mem" + "G"

    input:
        file samp_sheet

    output:
        file "lane_fastqs" into bcl2fastq_output
        set file("lane_fastqs/Undetermined_S0_*_R1_001.fastq.gz"), file("lane_fastqs/Undetermined_S0_*_R2_001.fastq.gz") into fastqs mode flatten
        file "lane_fastqs/fake*.gz" into fakes mode flatten

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
    clusterOptions "-l mfree=1G"
    publishDir = [path: "${params.output_dir}/", pattern: "demux_out/*.stats.json", mode: 'copy']
    publishDir = [path: "${params.output_dir}/", pattern: "demux_out/*.fastq", mode: 'copy']
    publishDir = [path: "${params.output_dir}/", pattern: "demux_out/*.csv", mode: 'copy'] 

    input:
        set file(R1), file(R2) from fastqs
        file sample_sheet_file

    output:
        file "demux_out/*" into seg_output
        file "demux_out/*.fastq" into samp_fastqs_check mode flatten
        file "demux_out/*.stats.json" into stats
        file "demux_out/*.csv" into csv_stats

    """
    mkdir demux_stats
    make_sample_fastqs.py --run_directory $params.run_dir \
        --read1 <(zcat $R1) --read2 <(zcat $R2) \
        --file_name $R1 --sample_layout $sample_sheet_file \
        --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
        --output_dir ./demux_out

    """    
}


process demux_dash {
    module 'java/latest:modules:modules-init:modules-gs:gcc/8.1.0:R/3.5.2'
    clusterOptions "-l mfree=8G"

    publishDir path: "${params.output_dir}/", pattern: "demux_dash", mode: 'copy'


    input:
        file demux_stats_files from seg_output.collect()
        file icon from Channel.fromPath('$baseDir/bin/bbi_icon.png')
    output:
        file demux_dash

    """
    mkdir demux_dash
    mkdir demux_dash/img
    cp $icon demux_dash/img/
    generate_html.R \
        "." --p7_rows "$params.p7_rows" --p5_cols "$params.p5_cols"

    """



}
