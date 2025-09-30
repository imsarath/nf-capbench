
process SAMTOOLS_SUBSAMPLE {
    tag   "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37998e74f9455bc772603426441a58036a4d4a6e236cc5c0b80aceab72738bf4/data' :
        'community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509' }"


    input:
    tuple val(meta), path(bam), path(bai), val(fraction), val(targetM)
    val  seed

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"),  emit: bam
    path "versions.yml",               emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${targetM.toString()}M"
    def perc = (fraction * 100).round(2).toString()

    """
    set -euo pipefail

    samtools view \
        $args \
        -hs $seed.$perc \
        -Sbo ${prefix}.bam $bam

    samtools index ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """


    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${targetM.toString()}M"

    """
    set -euo pipefail

    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS

    """

}
