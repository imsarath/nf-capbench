
process DOWNSAMPLING_FACTOR {
    tag   "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8b/8b36354c45846a274872a174e1a351e571ae6a377696c2f24b23bdda478ca68d/data' :
        'community.wave.seqera.io/library/awk_pandas_python:1d39e9f44134204d' }"


    input:
    tuple val(meta), path(duplex_metrics)
    val  targets_str

    output:
    path "downsampling_factors.tsv",   emit: factors_tsv
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def samples_tsv = "all_summary_dups.tsv"
    def targets_arg = targets_str ? "--targets ${targets_str}" : ''

    """

    awk -F "\\t" -v OFS="\\t" 'NR==1 {print "sample", "totalReadPairs"} \$1 == 1 {print FILENAME, \$2}' *duplex_yield_metrics.txt | \\
        sed 's/\\.duplex_yield_metrics\\.txt//g' > ${samples_tsv}

    calculate_downsampling_factors.py  $samples_tsv \
        $args \
        $targets_arg  \
        --out downsampling_factors.tsv


    # Stamp versions
    PYVER=\$(python --version 2>&1 | awk '{print \$2}')
    PDVER=\$(python - <<'PY'
    import pandas
    print(pandas.__version__)
    PY
    )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$PYVER"
        pandas: "\$PDVER"
        wrapper: "DOWNSAMPLING_FACTOR/1.0"
    END_VERSIONS

    """

    stub:

    """

    # Take the first two targets from the passed list (defaults to 1 and 5 if missing)
    T1=\$(echo "${targets_str}" | awk -F',' '{gsub(/ /,""); print \$1}')
    T2=\$(echo "${targets_str}" | awk -F',' '{gsub(/ /,""); print \$2}')
    [ -z "\$T1" ] && T1=1
    [ -z "\$T2" ] && T2=5

    # Convert target (in millions) to read pairs
    T1RP=\$(( T1 * 1000000 ))
    T2RP=\$(( T2 * 1000000 ))

    samples=(\$(find . -name "*.duplex_yield_metrics.txt" -exec basename {} \\; | sed 's/.duplex_yield_metrics.txt//g'))

    echo -e "sample_id\\ttotal_readpairs\\ttarget_readpairs\\ttarget_M\\tfraction\\testimated_pairs_kept" > downsampling_factors.tsv
    for sample in "\${samples[@]}"; do
        echo -e "\${sample}\\t42000000\\t\${T1RP}\\t\${T1}.0\\t0.25000000\\t10500000" >> downsampling_factors.tsv
        echo -e "\${sample}\\t42000000\\t\${T2RP}\\t\${T2}.0\\t0.50000000\\t21000000" >> downsampling_factors.tsv
    done

    # Stamp versions
    PYVER=\$(python --version 2>&1 | awk '{print \$2}')
    PDVER=\$(python - <<'PY'
    import pandas
    print(pandas.__version__)
    PY
    )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$PYVER"
        pandas: "\$PDVER"
        wrapper: "DOWNSAMPLING_FACTOR/1.0"
    END_VERSIONS

    """
}
