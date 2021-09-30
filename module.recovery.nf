// TODO


save_recovery2 = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz-summary.txt/ + "-recovery_summary.txt"}
save_recovery = {params.output_dir + "/recovery_output/" +  it - ~/.fastq.gz.txt.gz/ + "-recovery_table.txt.gz"}
process run_recovery {
    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery, pattern: "*.gz.txt.gz", mode: 'link'
    publishDir path: "${params.output_dir}/recovery_output", saveAs: save_recovery2, pattern: "*-summary.txt", mode: 'link'

    input:
        file input from Channel.fromPath("${params.demux_out}/Undetermined*")
        file sample_sheet_file5

    output:
        file "*gz.txt.gz"
        file "*summary.txt" into summaries
    when:
        params.run_recovery


    """
    recovery_script.py --input_file <(zcat $input) --output_file ${input}.txt \
        --run_directory $params.run_dir \
        --sample_layout $sample_sheet_file5 \
        --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
        --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
        --p7_barcode_file $params.p7_barcode_file \
        --p5_barcode_file $params.p5_barcode_file \
        --lig_barcode_file $params.lig_barcode_file \
        --level $params.level \
        --rt_barcodes $params.rt_barcode_file
     pigz -p 1 *.fastq.gz.txt
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

