process bcl2fastq {
    cache 'lenient'
    cpus max_cores_bcl
    memory "${bcl_mem} GB"

    input:
        path( run_dir )
        file( bcl_samp_sheet )
        val( max_cores_bcl )
        val( bcl_mem )

    output:
        path "lane_fastqs" emit: bcl2fastq_output
        file( "lane_fastqs/Undetermined_S0_*_R1_001.fastq.gz" ), file("lane_fastqs/Undetermined_S0_*_R2_001.fastq.gz") mode flatten, emit: fastqs
        file( "lane_fastqs/fake*.gz" ) optional true mode flatten, emit: fakes

    """/bin/bash
set -Eeuo pipefail

min_threads=\$((($max_cores_bcl/2)<4 ? ($max_cores_bcl/2):4))

bcl2fastq -R $run_dir --output-dir ./lane_fastqs \
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
