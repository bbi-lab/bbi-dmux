process seg_sample_fastqs1 {
    cache 'lenient'

    publishDir path: "${params.output_dir}/", pattern: "demux_out/*fastq.gz", mode: 'link'     
    publishDir path: "${params.output_dir}/demux_out/", pattern: "*.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/demux_out/", pattern: "*.json", mode: 'copy'

    input:
        set file(R1), file(R2) from fastqs
        file( run_parameters_file )
        file( sample_sheet_file )
        file( rt_barcode_file )
        file( p5_barcode_file )
        file( p7_barcode_file )
        file( lig_barcode_file )

    output:
        file "demux_out/*", emit: seg_output
        file "demux_out/*.fastq.gz", emit: samp_fastqs_check
        file "demux_out/*.stats.json" mode flatten, emit: json_stats 
        file "demux_out/*.csv", emit: csv_stats
 
"""/bin/bash
set -Eeuo pipefail

mkdir demux_out
make_sample_fastqs.py --run_directory . \
    --read1 <(zcat $R1) --read2 <(zcat $R2) \
    --file_name $R1 --sample_layout $sample_sheet_file \
    --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
    --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
    --rt_barcode_file $rt_barcode_file \
    --p5_barcode_file $p5_barcode_file \
    --p7_barcode_file $p7_barcode_file \
    --lig_barcode_file $lig_barcode_file \
    --multi_exp "$params.multi_exp" \
    --buffer_blocks $params.demux_buffer_blocks \
    --output_dir ./demux_out --level $params.level
pigz -p 8 demux_out/*.fastq
"""
}
