
process demux_dash {
    publishDir path: "${params.output_dir}/", pattern: "demux_dash", mode: 'copy'

    input:
        file( demux_stats )
        file( jsons )
        file( sample_sheet )
    output:
        file( demux_dash )

    script:
    
    out_dir_str = params.output_dir.replaceAll("/\\z", "");
    project_name = out_dir_str.substring(out_dir_str.lastIndexOf("/")+1);

"""/bin/bash
set -Eeuo pipefail

mkdir demux_dash
cp -R $baseDir/bin/skeleton_dash/* demux_dash/
generate_html.R \
    "." --p7_rows "$params.p7_rows" --p5_cols "$params.p5_cols" --p7_wells "$params.p7_wells" \
    --p5_wells "$params.p5_wells" --level "$params.level" --project_name "${project_name}" \
    --sample_sheet "$sample_sheet"
"""

}
