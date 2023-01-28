#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPARE_GENOME as PREPARE_GENOME_BWA } from '../../subworkflows/local/prepare_genome.nf'
include { PREPARE_GENOME as PREPARE_GENOME_BOWTIE } from '../../subworkflows/local/prepare_genome.nf' addParams( gene_bed: "*" )
include { PREPARE_GENOME as PREPARE_GENOME_CHROMAP } from '../../subworkflows/local/prepare_genome.nf' addParams( gene_bed: "*" )
include { PREPARE_GENOME as PREPARE_GENOME_STAR } from '../../subworkflows/local/prepare_genome.nf' addParams( gene_bed: "*" )

workflow {
    PREPARE_GENOME_BWA ( "bwa" )
    PREPARE_GENOME_BOWTIE ( "bowtie2" )
    PREPARE_GENOME_CHROMAP ( "chromap" )
    PREPARE_GENOME_STAR ( "star" )
}