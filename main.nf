/*
** Check that Nextflow version meets minimum version requirements.
*/
def minMajorVersion = 20
def minMinorVersion = 07
checkNextflowVersion( minMajorVersion, minMinorVersion )
nextflow.enable.dsl=2

DEFAULT = "default"
default_rt2_barcode_file = "$baseDir/bin/barcode_files/rt2.txt"
default_rt3_barcode_file = "$baseDir/bin/barcode_files/rt.txt"
default_p5_barcode_file = "$baseDir/bin/barcode_files/p5.txt"
default_p7_barcode_file = "$baseDir/bin/barcode_files/p7.txt"
default_lig_barcode_file = "$baseDir/bin/barcode_files/ligation.txt"
default_star_file = "$baseDir/bin/star_file.txt"

// Parse input parameters
params.help = false
params.rerun = false
params.star_file = DEFAULT
params.level = 3
params.bcl_max_mem = 40
params.fastq_chunk_size = 100000000
params.run_recovery = false
params.rt_barcode_file = DEFAULT
params.p5_barcode_file = DEFAULT
params.p7_barcode_file = DEFAULT
params.lig_barcode_file = DEFAULT
params.large = false
params.generate_samplesheets = false
params.max_cores = 16
params.max_wells_per_sample = 20
params.demux_buffer_blocks = 16

params.multi_exp = 0
params.p5_cols = 0
params.p7_rows = 0
params.p5_wells = 0
params.p7_wells = 0

include {
    generate_sheets
    check_sample_sheet
    make_sample_sheet
} from './module.samplesheets'

include {
    bcl2fastq
} from './module.bcl2fastq'

include {
    seg_sample_fastqs1
} from './module.process-small'

include {
    seg_sample_fastqs2
} from './module.process-large'

include {
    demux_dash
} from './module.dashboard'

def validate(params) {
    if (!params.run_dir || !params.output_dir || !params.sample_sheet ) {
        exit 1, "Must include config file using -c CONFIG_FILE.config that includes output_dir, sample_sheet and run_dir."
    }

    if (!params.help) {
        return;
    }

    log.info"""
        BBI sci-RNA-seq Demultiplexer'
        --------------------------------'
        '
        For reproducibility, please specify all parameters to a config file'
        by specifying -c CONFIG_FILE.config.'
        '
        Usage: '
            nextflow run bbi-dmux -c CONFIG_FILE'
        '
        Help: '
            --help                                     Show this message and exit.'
        '
        Required parameters (specify in your config file):'
            params.run_dir = RUN_DIRECTORY             Path to the sequencer output.'
            params.output_dir OUTPUT DIRECTORY         Output directory.'
            params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.'
            params.level = 3                           Level of run - either 2 or 3.'
        '
        Required parameters (one of the pairs below is required - p7_wells and p5_wells or p7_rows and p5_cols or mult-exp):'
            params.p7_wells = "A1 B1 C1"               Alternative to p7_rows and p5_cols - specify specific PCR wells instead of full rows/columns. Must match order of params.p5_wells.'
            params.p5_wells = "A1 A2 A3"               Alternative to p7_rows and p5_cols - specify specific PCR wells instead of full rows/columns. Must match order of params.p7_wells.'
            params.p7_rows = "A B C"                   The PCR rows used - must match order of params.p5_cols.'
            params.p5_cols = "1 2 3"                   The PCR columns used - must match order of params.p7_rows.'
            params.multi_exp = "see config"            The PCR columns used for each experiment in map format - see example.config.'
        '
        '
        Optional parameters (specify in your config file):'
            params.rt_barcode_file = "default"         The path to a custom RT barcode file. If "default", default BBI barcodes will be used.'
            params.p7_barcode_file = "default"         The path to a custom p7 barcode file. If "default", default BBI barcodes will be used.'
            params.p5_barcode_file = "default"         The path to a custom p5 barcode file. If "default", default BBI barcodes will be used.'
            params.lig_barcode_file = "default"        The path to a custom ligation barcode file. If "default", default BBI barcodes will be used.'
            params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
            process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
            process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
            params.star_file = PATH/TO/FILE            File with the genome to star maps, similar to the one included with the package.'
            params.fastq_chunk_size = 100000000        The number of reads that should be processed together for demultiplexing.'
            params.bcl_max_mem = 40                    The maximum number of GB of RAM to assign for bcl2fastq'
            params.large = false                       Is this a very large run? If true, the fastqs will be split - note that for smaller runs this will make the pipeline run more slowly.'
            params.max_wells_per_sample = 20           The maximum number of wells per sample - if a sample is in more wells, the fastqs will be split then reassembled.'
            params.demux_buffer_blocks = 16            The number of 8K blocks to use for demux output buffer.'
            --run_recovery true                        Add this to run the recovery script AFTER running the normal pipeline.'
            --generate_samplesheets input_csv          Add this to generate the necessary samplesheet from the BBI universal input sheet.'    
        '
        Issues? Contact hpliner@uw.edu'
    """
    exit 1
}

sample_sheet_file = file(params.sample_sheet)
run_dir = Channel.fromPath(params.run_dir)
run_parameters_file = file("${params.run_dir}/[Rr]unparameters.xml*")

if (params.rt_barcode_file == DEFAULT) {
    if (params.level == "2") {
        rt_barcode_file = file(default_rt2_barcode_file)
    }
    if (params.level == "3") {
        rt_barcode_file = file(default_rt3_barcode_file)
    }
} else {
    rt_barcode_file = file(params.rt_barcode_file)
}

star_file = params.star_file == DEFAULT ? file(params.default_star_file) : file(params.star_file)
lig_barcode_file = params.lig_barcode_file == DEFAULT ? file(default_lig_barcode_file) : file(params.lig_barcode_file)
p5_barcode_file = params.p5_barcode_file == DEFAULT ? file(default_p5_barcode_file) : file(params.p5_barcode_file)
p7_barcode_file = params.p7_barcode_file == DEFAULT ? file(default_p7_barcode_file) : file(params.p7_barcode_file)

if (params.max_cores > 16) {
    max_cores_bcl = 16
    bcl_mem = params.bcl_max_mem/16
} else {
    max_cores_bcl = params.max_cores
    bcl_mem = params.bcl_max_mem/max_cores_bcl

workflow {
    main:

    validate(params)

    if ( params.generate_samplesheets ) {
        bbi_universal_sheet_file = Channel.fromPath(params.generate_sample_sheets)
        generate_sheets( bbi_universal_sheet_file )

        sample_sheet_file = generate_sheets.out.sample_sheet
    } else {
        check_sample_sheet( sample_sheet_file, star_file, rt_barcode_file, params.level, params.max_wells_per_sample )
        
        sample_sheet_file = check_sample_sheet.out.good_sample_sheet
    }

    if ( !params.run_recovery ) {
        make_sample_sheet( run_parameters_file, good_sample_sheet )

        bcl_sample_sheet = make_sample_sheet.out.bcl_sample_sheet
    }

    bcl2fastq( run_dir, bcl_sample_sheet, max_cores_bcl, bcl_mem )

    fastqs = bcl2fastq.out.fastqs

    /*
    ** ================================================================================
    ** PATH 1 - For small datasets, no fastq splitting
    ** ================================================================================
    */
    if (!params.large) {
        seg_sample_fastqs1( fastqs, sample_sheet_file )

        csv_stats = seg_sample_fastqs1.out.csv_stats
        json_stats = seg_sample_fastqs1.out.json_stats
    }

    /*
    ** ================================================================================
    ** PATH 2 - For large datasets, with fastq splitting
    ** ================================================================================
    */
    else {
        fastq_chunks = fastqs.splitFastq(by: params.fastq_chunk_size, file: true, pe: true)
        seg_sample_fastqs2( fastq_chunks, sample_sheet_file )

        all_csv_stats = seg_sample_fastqs2.out.csv_stats
        all_json_stats = seg_sample_fastqs2.out.json_stats
        sample_fastqs_check = seg_sample_fastqs2.out.sample_fastqs_check

        // Combine Fastqs
        get_prefix = { fname ->
            (fname - ~/_R1_001\.[0-9]+\.fastq.fastq.gz/)
        }

        grouped_fastqs = samp_fastqs_check
            .map { file -> tuple(get_prefix(file.name), file) }
            .groupTuple()
        
        recombine_fastqs( grouped_fastqs )

        // Combine CSVs
        csv_prefix = { fname ->
            (fname - ~/_R1_001\.[0-9]+\.fastq\.[a-z]+_[a-z]+\.csv/)
        }

        grouped_csvs = csv_stats
            .map { file -> tuple(csv_prefix(file.name), file) }
            .groupTuple()

        recombine_csvs( grouped_csvs )
        
        // Combine JSON
        json_prefix = { fname ->
            (fname - ~/_R1_001\.[0-9]+\.fastq\.stats\.json/)
        }

        grouped_jsons = json_stats
            .map { file -> tuple(json_prefix(file.name), file) }
            .groupTuple()
        
        recombine_jsons( grouped_jsons )

        csv_stats = recombine_csvs.out
        json_stats = recombine_jsons.out
    }
    
    demux_dash( csv_stats.collect(), json_stats.collect(), sample_sheet_file )

    if (params.run_recovery) {
        // run_recovery
        // sum_recovery
    }

    println "bbi-dmux completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}



fastqs.into { fastqs_path1; fastqs_path2 }


/*
** ================================================================================
** End path 2 - Start recovery
** ================================================================================
*/



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
  if( majorVersion < minMajorVersion || ( majorVersion == minMajorVersion && minorVersion < minMinorVersion ) )
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
