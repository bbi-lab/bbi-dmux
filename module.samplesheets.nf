
process generate_sheets {
    publishDir path: "${params.output_dir}", pattern: "SampleSheet.csv", mode: 'copy'
    publishDir path: "${params.output_dir}", pattern: "SampleMap.csv", mode: 'copy'
    publishDir path: "${params.output_dir}", pattern: "GarnettSheet.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/sample_id_maps", pattern: "*_SampleIDMap.csv", mode: 'copy'

    input:
        file( sample_sheet )
    output:
        file ("*Sheet.csv"), emit: sample_sheet
        file ("SampleMap.csv") optional true
        file ("*_SampleIDMap.csv") optional true

"""/bin/bash
set -Eeuo pipefail

generate_sample_sheets.py $sample_sheet
"""

}

process check_sample_sheet {
    input:
        file( sample_sheet )
        file( star_file )
        file( rt_barcode_file )
        val( level )
        val( max_wells_per_sample )

    output:
        file ("*.csv" ), emit: good_sample_sheet

    when:
        params.generate_samplesheets == "no_input"

"""/bin/bash
set -Eeuo pipefail

check_sample_sheet.py \
    --sample_sheet $sample_sheet \ 
    --star_file $star_file \
    --level $level \
    --rt_barcode_file $rt_barcode_file \
    --max_wells_per_samp $max_wells_per_sample    
"""

}

/**
 * Generates sample sheet used for bcl2fastq
 */
process make_sample_sheet {
    cache 'lenient'

    input:
        file( run_parameters_file )
        file( good_sample_sheet )

    output:
        file( "SampleSheet.csv" ), emit: bcl_sample_sheet

    when:
        !params.run_recovery

"""/bin/bash
set -Eeuo pipefail

make_sample_sheet.py --run_directory .

"""    
}
