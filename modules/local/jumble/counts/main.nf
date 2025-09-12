
process JUMBLE_COUNTS {
    tag "$meta.id"
    label 'process_medium'

    /*
     * Either use the provided Conda environment (recommended for portability)
     * or set a container image that contains R and jumble-run.R.
     */
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rocker-r-ver:4.3.2' :
        'rocker/r-ver:4.3.2' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(targets), optional: true

    output:
    tuple val(meta), path("*.counts.RDS"), emit: outdir
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    /*
     * Variables used in both output and script blocks.
     * - prefix: used for naming logs and as default for outdir
     * - outdir: what gets passed to -o; defaults to meta.id
     * - args: extra CLI flags (e.g. task.ext.args = '--some-flag')
     *
     * These variables are now set inside the script block below.
     */

    /*
     * You can compute memory hints if jumble-count.R supports them.
     * For example, to pass RAM (MB) or threads if supported by the tool:
     *
     * def memory_in_mb = task.memory ? task.memory.toUnit('MB') : null
     * def threads      = task.cpus ?: 1
     */

    script:
    def args = task.ext.args ?: ''
    def targets_bed = targets ? "-t ${targets} " : ''
    """

    # Run the tool
    jumble-count.R \\
        $args \\
        -c $task.cpus \\
        ${targets_bed} \\
        -b ${bam} \\
        -o "./"

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jumble-count.R: \$( jumble-count.R --version 2>&1 | head -n 1 || echo "unknown" )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.counts.RDS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jumble-run.R: "stub"
        R: \$( R --version | sed -n '1p' )
    END_VERSIONS
    """
}
