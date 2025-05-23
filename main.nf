#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fasta = "$baseDir/genomeInput/GRCh38.p13.genome.fa"
params.gtf = "$baseDir/genomeInput/gencode.v39.annotation.gtf"
//This pipeline expects paired read fastq files in their own folder under this dir
params.fastq = "$baseDir/samples"

fasta_channel = Channel.fromPath( params.fasta, checkIfExists: true )
gtf_channel = Channel.fromPath( params.gtf, checkIfExists: true )
fastq_channel = Channel.fromPath( params.fastq, checkIfExists: true)

workflow {
    genomeGenerate1(fasta_channel, gtf_channel).view()
    alignPass1(fastq_channel, genomeGenerate1.out).view()
    genomeGenerate2(fasta_channel, gtf_channel, alignPass1.out).view()
    alignPass2(fastq_channel, genomeGenerate2.out).view()
}


process genomeGenerate1 {
    /*
     * This conda environment is used with the "waveDynamic" profile to dynamically
     *  fetch a singularity OCI-SIF container with the dependencies in this xml.
     * By defauly this is ignored and the static container specified in the
     *  "standard" profile is used.
     */
    conda 'conda/env.yml'

    publishDir "STEP1GENOME", mode: 'copy'
    
    input:
    path fasta
    path gtf

    output:
    val true

    script:
    """
    STAR --runMode genomeGenerate --runThreadN 14 --genomeDir step1_Genome --genomeFastaFiles $fasta --sjdbGTFfile $gtf;
    """

}

process alignPass1 {
    /*
     * This conda environment is used with the "waveDynamic" profile to dynamically
     *  fetch a singularity OCI-SIF container with the dependencies in this xml.
     * By defauly this is ignored and the static container specified in the
     *  "standard" profile is used.
     */
    conda 'conda/env.yml'

    publishDir "ALIGN1", mode: 'copy'
    
    input:
    path fastqs
    val ready

    output:
    val true

    script:
    """
    find $baseDir/samples -type d -mindepth 1 | while read sampleDir; do 
        sampleName=\$(basename "\$sampleDir");
        STAR --runMode alignReads --runThreadN 14 --genomeDir $baseDir/STEP1GENOME/step1_Genome --readFilesIn "\$sampleDir"/*.fastq;
        mkdir "\$sampleName"; mv ./* "\$sampleName"
        done;
    """

}


process genomeGenerate2 {
    /*
     * This conda environment is used with the "waveDynamic" profile to dynamically
     *  fetch a singularity OCI-SIF container with the dependencies in this xml.
     * By defauly this is ignored and the static container specified in the
     *  "standard" profile is used.
     */
    conda 'conda/env.yml'

    publishDir "STEP2GENOME", mode: 'copy'
    
    input:
    path fasta
    path gtf
    val ready

    output:
    val true

    /*
     * This step performs a second alignment and passes all SJ.out files for all samples 
     *  from the first alignReads step. This is recommended in the latest STAR manual and 
     *  is not currently performed by LENS, the most complete existing pipeline.
     */
    script:
    """
    STAR --runMode genomeGenerate --runThreadN 14 --genomeDir step2_Genome --genomeFastaFiles $fasta --sjdbGTFfile $gtf --sjdbFileChrStartEnd $baseDir/ALIGN1/*/SJ.out;
    """

}

process alignPass2 {
    /*
     * This conda environment is used with the "waveDynamic" profile to dynamically
     *  fetch a singularity OCI-SIF container with the dependencies in this xml.
     * By defauly this is ignored and the static container specified in the
     *  "standard" profile is used.
     */
    conda 'conda/env.yml'

    publishDir "ALIGN2", mode: 'copy'
    
    input:
    path fastqs
    val ready

    output:
    val true

    script:
    """
    find $baesDir/samples -type d -mindepth 1 | while read sampleDir; do 
        sampleName=\$(basename "\$sampleDir");
        STAR --runMode alignReads --runThreadN 14 --genomeDir $baseDir/STEP2GENOME/step2_Genome --readFilesIn "\$sampleDir"/*.fastq;
        mkdir "\$sampleName"; mv ./* "\$sampleName"
        done;
    """

}
