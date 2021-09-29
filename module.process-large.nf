process seg_sample_fastqs2 {
    cache 'lenient'

    input:
        set file(R1), file(R2) from fastqs
        file sample_sheet_file

    output:
        file "demux_out/*", emit: seg_output
        file "demux_out/*.fastq.gz" mode flatten, emit: samp_fastqs_check
        file "demux_out/*.stats.json" mode flatten, emit: json_stats
        file "demux_out/*.csv" mode flatten, emit: csv_stats

"""/bin/bash
set -Eeuo pipefail

mkdir demux_out
make_sample_fastqs.py --run_directory $params.run_dir \
    --read1 $R1 --read2 $R2 \
    --file_name $R1 --sample_layout $sample_sheet_file \
    --p5_cols_used $params.p5_cols --p7_rows_used $params.p7_rows \
    --p5_wells_used $params.p5_wells --p7_wells_used $params.p7_wells \
    --rt_barcode_file $params.rt_barcode_file \
    --p5_barcode_file $params.p5_barcode_file \
    --p7_barcode_file $params.p7_barcode_file \
    --lig_barcode_file $params.lig_barcode_file \
    --multi_exp "$params.multi_exp" \
    --buffer_blocks $params.demux_buffer_blocks \
    --output_dir ./demux_out --level $params.level
pigz -p 8 demux_out/*.fastq    
"""
}

process recombine_fastqs {
    cache 'lenient'
    publishDir  path: "${params.output_dir}/demux_out", pattern: "*.fastq.gz", mode: 'move'

    input:
        set prefix, file(all_fqs) from grouped_files 

    output:
        file "*.gz"

"""/bin/bash
set -Eeuo pipefail

cat $all_fqs > ${prefix}.fastq.gz 
"""
}

process recombine_csvs {
    cache 'lenient'
    publishDir  path: "${params.output_dir}/demux_out/", pattern: "*.csv", mode: 'copy'

    input:
        set prefix, file(all_csvs) from grouped_csvs

    output:
        file "*.csv"

"""/bin/bash
set -Eeuo pipefail

csvs="$all_csvs"    
arr=(\$csvs)
if [ "$params.level" = "3" ]; then
    cat \$(IFS=\$'\n'; echo "\${arr[*]}" | grep lig_counts) | awk -F ',' 'BEGIN {OFS = ","} {a[\$1] += \$2} END {for (i in a) print i, a[i]}' > ${prefix}.lig_counts.csv    
fi
cat \$(IFS=\$'\n'; echo "\${arr[*]}" | grep rt_counts) | awk -F ',' 'BEGIN {OFS = ","} {a[\$1] += \$2} END {for (i in a) print i, a[i]}' > ${prefix}.rt_counts.csv
cat \$(IFS=\$'\n'; echo "\${arr[*]}" | grep pcr_counts) | awk -F ',' 'BEGIN {OFS = ","; SUBSEP = OFS = FS} {a[\$1,\$2] += \$3} END {for (i in a) print i, a[i]}' > ${prefix}.pcr_counts.csv

"""
}

process recombine_jsons {
    cache 'lenient'
    publishDir  path: "${params.output_dir}/demux_out/", pattern: "*.json", mode: 'copy'

    input:
        set prefix, file(all_jsons) from grouped_jsons

    output:
        file "*.json"

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
