#!/usr/bin/env python3
import peppy

#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    input_keys = ['raw_r1', 'raw_r2', 'demuxed_r1', 'demuxed_r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

c5_samples = ["AS1F24_1", "AS1F24_2", "AS1F24_3", "AS1F73_1", "AS1F73_2", "AS1F73_3", "AS2F24_1", "AS2F24_2", "AS2F24_3", "AS2F71_1", "AS2F71_2", "AS2F71_3", "AS4F26_1", "AS4F26_2", "AS4F26_3", "AS4F71_1", "AS4F71_2", "AS4F71_3", "AS6F36_1", "AS6F36_2", "AS6F46_1", "AS6F72_1", "AS6F72_2", "AS6F72_3", "AS9F19_1", "AS9F20_1", "AS9F21_1", "AS9F72_1", "AS9F72_2", "AS9F66_1", "AS28F26_1", "AS28F26_2", "AS28F26_3", "AS28F72_1", "AS28F72_2", "AS28F72_3", "BA16F16_1", "BA16F17_1", "BA16F65_1", "BA6F31_1", "BA6F31_2", "BA6F31_3", "BA6F62_1", "BA6F62_2", "BA6F62_3", "BA12F20_1", "BA12F20_2", "BA12F20_3", "BA12F71_1", "BA12F71_2", "BA12F71_3", "BA13F22_1", "BA13F22_2", "BA13F22_3", "BA13F69_1", "BA13F69_2", "BA13F69_3", "BA14F23_1", "BA14F23_2", "BA14F23_3", "BA14F68_1", "BA14F68_2", "BA14F68_3", "BA17F30_1", "BA17F30_2", "BA17F67_1", "BA17F67_2", "BA17F67_3", "BR5F27_1", "BR5F27_2", "BR5F27_3", "BR5F74_1", "BR5F74_2", "BR5F74_3", "BR14F35_1", "BR14F43_1", "BR14F72_1", "BR14F72_2", "BR14F70_1", "BR23F27_1", "BR23F27_2", "BR23F27_3", "BR23F74_1", "BR23F74_2", "BR24F25_1", "BR24F25_2", "BR24F25_3", "BR24F73_2", "BR24F73_3", "BR25F35_3", "BR25F72_1", "BR25F72_2", "BR25F72_3", "BR27F35_1", "BR27F35_2", "BR27F73_1", "BR27F73_2", "BR27F73_3", "BR48F28_1", "BR48F28_2", "BR48F28_3", "BR48F73_1", "BR48F73_2", "BR48F73_3", "BR12F22_1", "BR12F22_2", "BR12F22_3", "BR12F75_1", "BR12F75_2", "BR12F75_3", "MN1F23_1", "MN1F23_3", "MN1F68_1", "MN1F68_2", "MN1F68_3", "LS15F36_1", "LS15F36_2", "LS15F36_3", "LS15F69_1", "LS15F69_2", "LS15F69_3", "LS17F32_1", "LS17F32_2", "LS17F32_3", "LS17F69_1", "LS17F69_2", "LS17F69_3", "LS34A6F16_1", "LS34A6F16_2", "LS34A6F16_3", "LS34F61_1", "LS34F61_2", "LS36F29_1", "LS36F67_1", "LS36F67_2", "LS36F67_3", "LS38F22_1", "LS38F22_2", "LS38F22_3", "LS38F68_1", "LS48F69_1", "LS48F69_2", "LS48F69_3", "RN2F34_1", "RN2F34_2", "RN2F34_3", "RN2F68_1", "RN2F68_2", "RN2F68_3", "RN8F35_1", "RN8F35_2", "RN8F35_3", "RN8F72_1", "RN8F72_2", "RN8F72_3", "RN16F34_1", "RN16F73_1", "RN16F73_2", "RN16F73_3", "RN19F35_1", "RN19F35_2", "RN19F35_3", "RN19F77_1", "RN19F77_2", "RN19F77_3", "RN21F17_1", "RN21F17_2", "RN21F17_3", "RN21F72_1", "RN21F72_2", "RN21F72_3", "RN27F24_1", "RN27F24_2", "RN27F24_3", "RN27F69_1", "RN27F69_2", "RN27F69_3", "RN30F30_1", "RN30F30_2", "RN30F31_1", "RN30F69_1", "RN30F69_2", "RN30F69_3", "RN36F18_1", "RN36F19_1", "RN36F19_2", "RN36F69_1", "RN36F69_2", "SC2F36_1", "SC2F36_2", "SC2F36_3", "SC2F67_1", "SC2F67_3", "SC3F29_1", "SC3F30_1", "SC3F30_2", "SC3F65_1", "SC3F65_2", "SC3F65_3", "SC7F35_1", "SC7F35_2", "SC7F35_3", "SC7F68_1", "SC7F68_2", "SC7F68_3", "SC8F17_1", "SC8F17_2", "SC8F17_3", "SC8F68_1", "SC8F68_2", "SC8F68_3", "SC11F36_1", "SC11F36_2", "SC11F42_1", "SC11F67_1", "SC11F67_2", "SC11F67_3", "UR6F18_1", "UR6F18_2", "UR6F18_3", "UR6F70_1", "UR6F70_2", "UR6F70_3", "UR12F22_1", "UR12F22_2", "UR12F22_3", "UR12F68_1", "UR12F68_2", "UR12F68_3", "UR15F21_1", "UR15F21_2", "UR15F21_3", "UR15F70_1", "UR15F70_2", "UR15F70_3", "UR20F31_1", "UR20F31_2", "UR20F31_3", "UR20F41_1", "UR20F41_2", "UR20F41_3", "UR21F68_1", "UR21F68_2", "UR21F68_3", "UR32F21_1", "UR32F21_2", "UR32F21_3", "UR32F21_4", "UR32F21_5", "UR32F21_6", "UR32F70_1", "UR32F70_2", "UR32F70_3", "UR47F21_1", "UR47F21_2", "UR47F21_3", "UR47F68_1", "UR47F68_2", "UR47F68_3", "BA17F29_1", "BR14F35_2", "BA16F17_2", "BR23F74_3", "RN36F69_3", "BR24F73_1", "BR25F35_1", "BR25F35_2", "BR25F35_3", "LS34F61_1", "LS36F31_1", "LS36F31_2", "LS38F68_1", "LS38F68_2", "LS43F23_1", "LS43F33_1", "LS43F60_1", "LS43F60_2", "LS43F60_3", "LS45F23_1", "LS45F23_2", "LS45F23_3", "LS45F67_1", "LS45F67_1", "LS45F67_3", "LS48F35_1", "LS48F35_2", "MN1F23_1", "SC2F67_1", "RN16F34_2", "RN16F34_3", "UR21F21_1", "UR21F21_2", "UR21F21_3", "19RU_1", "19RU_1993_2", "20RU_1993_1", "20RU_1993_2", "21RU_1993_1", "21RU_1993_2", "10RE_1993_1", "10RE_1993_2", "10RE_1993_3", "10RE_1993_4", "11RE_1993_1", "12RE_1993_1", "12RE_1993_2", "13RE_1993_1", "13RE_1993_2", "13RE_1993_3", "Reporoa_1996_1", "Reporoa_1996_2", "Reporoa_1996_3", "Reporoa_1996_4", "08WL_1", "08WL_2", "08WL_3", "08WL_4", "08WL_5", "Wellsford_1996_1", "Wellsford_1996_2", "Wellsford_1996_3", "Wellsford_1996_4", "Linc_farm_1993_40", "Linc_farm_1993_41", "Linc_farm_1993_42", "Linc_farm_1993_42", "Linc_farm_1993_43", "Linc_farm_1993_44", "Linc_farm_1993_45", "Linc_farm_1993_46", "Linc_farm_1993_46", "Linc_farm_1993_47", "Lincoln_dispersal_study_01_10_1996_1", "Lincoln_dispersal_study_01_10_1996_2", "Lincoln_dispersal_study_01_10_1996_3", "Lincoln_dispersal_study_01_10_1996_4", "Lincoln_dispersal_study_01_10_1996_5", "Lincoln_dispersal_study_01_10_1996_6", "Lincoln_dispersal_study_01_10_1996_7", "Lincoln_dispersal_study_01_10_1996_8", "122_Lincoln_1996_1", "122_Lincoln_1996_2", "Hawkes_bay_Mokotutu_south_2004_1", "Hawkes_bay_Mokotutu_south_2004_2", "Hawkes_bay_Mokotutu_south_2004_3", "Hawkes_bay_Mokotutu_south_2004_4", "Hawkes_bay_Mokotutu_south_2004_5", "Hawkes_bay_Mokotutu_south_2004_6", "Hawkes_bay_Mokotutu_south_2004_7", "Hawkes_bay_Mokotutu_south_2004_8", "Hawkes_bay_Mokotutu_south_2004_9", "Hawkes_bay_Mokotutu_south_2004_10"]
c5_fix_samples = ["Linc_farm_1993_42_1", "Linc_farm_1993_42_2", "Linc_farm_1993_46_1", "Linc_farm_1993_46_2", "SC2F67_2", "MN1F23_2", "LS45F67_2", "LS38F68_3", "LS34F61_3", "BR27F35_3", "SC2F67_1", "MN1F23_1", "LS45F67_1", "LS38F68_1", "LS34F61_1", "BR25F35_3"]

#containers
bbduk_container = 'https://github.com/deardenlab/container-bbmap/releases/download/0.0.3/container-bbmap.bbmap_38.90.sif'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'
kraken_container = 'docker://staphb/kraken2:2.1.2-no-db'
multiqc_container = 'docker://ewels/multiqc:1.9'
bwa_container = 'docker://staphb/bwa:0.7.17'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.14:0.0.1'
seqkit_container= 'docker://nanozoo/seqkit:2.3.1--401ce9e'
samtools_container = 'docker://staphb/samtools:1.10'

#########
# RULES #
#########

rule target:
    input:
        ## trim, fastQC --> multiQC
        'output/03_multiqc_Mh_MhFV/multiqc_report.html',
        ## mapped - indexed bams & QC
        expand('output/05_bwa_{genome}/{sample}_sorted.bam', genome=["Mh_MhmtDNA_MhFV"], sample=all_samples),
        'output/07_mapping_qc_analysis/mapping_QC_res.csv',

################
## mapping QC ##
################

## QC analysis identifies any sample as having an issue if they have:
    # contamination:
        # more than 10% of reads classified as bacterial
        # less than 70% reads unclassified
    # are or may be M. aeth rather than M. hyp (better mapping to M. aeth)
    # have poor coverage/depth/low reads
        # less than 70% coverage of first M. hyp Hi-C scaffold
        # less than 9x mean depth of first M. hyp Hi-C scaffold
        # less than 2 million reads mapped

rule mapping_qc_analysis:
    input:
        kraken = expand('output/04_kraken/reports/kraken_{sample}_report.txt', sample=all_samples),
        mapping = expand('output/05_bwa_{genome}/samtools_flagstat/{sample}.out', sample=all_samples, genome=["Mh_MhmtDNA_MhFV", "MaIR"]),
        coverage = expand('output/06_samtools_Mh_MhFV/coverage/{sample}_coverage.out', sample=all_samples),
        sample_table_file = 'data/sample_table.csv',
        bad_samples = 'data/low_read_contaminated_samples.txt'
    output:
        mapping_qc_res = 'output/07_mapping_qc_analysis/mapping_QC_res.csv',
        concern_sample_ids = 'output/07_mapping_qc_analysis/lowread_contam_potMaeth_ids.txt'
    params:
        kraken_dir = 'output/04_kraken/reports',
        bwa_mh_dir = 'output/05_bwa_Mh_MhmtDNA_MhFV/samtools_flagstat',
        bwa_mair_dir = 'output/05_bwa_MaIR/samtools_flagstat',
        mh_coverage_dir = 'output/06_samtools_Mh_MhFV/coverage'
    log:
        'output/logs/mapping_qc_analysis.log'
    singularity:
        bioconductor_container
    script:
        'src/mapping_QC_analysis.R'

# count number of reads mapped to each contig
rule seqkit_count_reads:
    input:
        bam = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam',
        index = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam.bai'
    output:
        out = temp('output/06_seqkit_Mh_MhFV/{sample}.out')
    container:
        seqkit_container
    shell:
        'nice seqkit bam -C '
        '{input.bam} '
        '2> {output.out} '

rule samtools_coverage:
    input:
        bam = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam'
    output:
        coverage_out = temp('output/06_samtools_Mh_MhFV/coverage/{sample}_coverage.out')
    log:
        'output/logs/samtools_Mh_MhFV_coverage/samtools_Mh_MhFV_coverage_{sample}.log'
    shell:
        'nice samtools coverage '
        '{input.bam} '
        '-o {output.coverage_out} '
        '2> {log}'

rule filter_depth_file:
    input:
        depth_out = 'output/06_samtools_Mh_MhFV/depth/{sample}_depth.out',
        viral_hic_ids_list = 'data/contig_ids/viral_and_hic_contig_ids.txt'
    output:
        filtered_depth = temp('output/06_samtools_Mh_MhFV/depth_filtered/filtered_{sample}_depth.out')
    threads:
        10
    shell:
        'egrep -wf {input.viral_hic_ids_list} {input.depth_out} > {output.filtered_depth}'

rule samtools_depth:
    input:
        sorted_bam = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam'
    output:
        depth_out = temp('output/06_samtools_Mh_MhFV/depth/{sample}_depth.out')
    log:
        'output/logs/samtools_depth/samtools_depth_Mh_MhFV_{sample}.log'
    shell:
        'nice samtools depth '
        '{input.sorted_bam} '
        '-aa ' # write positions/contigs even when depth = 0
        '> {output.depth_out} '
        '2> {log}'

#############
## mapping ##
#############

# mapping rates - bwa doesn't output
rule samtools_flagstat:
    input:
        'output/05_bwa_{genome}/{sample}_sorted.bam'
    output:
        temp('output/05_bwa_{genome}/samtools_flagstat/{sample}.out')
    log:
        'output/logs/samtools_flagstat/samtools_flagstat_{sample}_{genome}.log'
    shell:
        'nice samtools flagstat '
        '{input} > {output} '
        '2> {log}'

rule samtools_index_MaIR:
    input:
        bam = 'output/05_bwa_MaIR/{sample}_sorted.bam'
    output:
        index = temp('output/05_bwa_MaIR/{sample}_sorted.bam.bai')
    log:
        'output/logs/samtools_index/samtools_index_MaIR_{sample}.log'
    shell:
        'nice samtools index '
        '{input.bam} '
        '2> {log}'

rule samtools_sort_MaIR:
    input:
        sam = 'output/05_bwa_MaIR/{sample}_bwa_mem.sam'
    output:
        sorted_bam = temp('output/05_bwa_MaIR/{sample}_sorted.bam')
    log:
        'output/logs/samtools_sort/samtools_sort_MaIR_{sample}.log'
    shell:
        'nice samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '2> {log}'

rule samtools_index_Mh:
    input:
        bam = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam'
    output:
        index = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam.bai'
    log:
        'output/logs/samtools_index/samtools_index_Mh_MhmtDNA_MhFV_{sample}.log'
    shell:
        'nice samtools index '
        '{input.bam} '
        '2> {log}'

rule samtools_sort_Mh:
    input:
        sam = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_bwa_mem.sam'
    output:
        sorted_bam = 'output/05_bwa_Mh_MhmtDNA_MhFV/{sample}_sorted.bam'
    log:
        'output/logs/samtools_sort/samtools_sort_Mh_MhmtDNA_MhFV_{sample}.log'
    shell:
        'nice samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '2> {log}'

rule bwa_mem:
    input:
        index = 'output/05_bwa_{genome}/index/index.bwt',
        r1 = 'output/01_bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/01_bbduk_trim/{sample}_r2.fq.gz',
    output:
        sam = temp('output/05_bwa_{genome}/{sample}_bwa_mem.sam')
    params:
        index_dir = 'output/05_bwa_{genome}/index/index'
    threads:
        20
    log:
        'output/logs/bwa_mem/{sample}_{genome}.log'
    singularity:
        bwa_container
    shell:
        'nice bwa mem '
        '-t {threads} '
        '{params.index_dir} '
        '{input.r1} {input.r2} '
        '> {output.sam} '
        '2> {log}'

rule bwa_index:
    input:
        genome = 'data/genomes/{genome}.fa'
    output:
        index = 'output/05_bwa_{genome}/index/index.bwt'
    params:
        outdir = 'output/05_bwa_{genome}/index/index'
    threads:
        20
    log:
        'output/logs/bwa_index_{genome}.log'
    singularity:
        bwa_container
    shell:
        'nice bwa index '
        '{input.genome} '
        '-p {params.outdir} '
        '2> {log} '

################
## initial QC ##
################

##kraken db from https://benlangmead.github.io/aws-indexes/k2 as build std was failing
rule kraken:
    input:
        r1 = 'output/01_bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/01_bbduk_trim/{sample}_r2.fq.gz',
        db = 'bin/kraken_std_PF'
    output:
        report = 'output/04_kraken/reports/kraken_{sample}_report.txt'
    log:
        'output/logs/04_kraken/kraken_{sample}.log'
    threads:
        20
    singularity:
        kraken_container
    shell:
        'nice kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--report {output.report} '
        '--output - '
        '--use-names '
        '{input.r1} {input.r2} '
        '2> {log}'

rule multiqc_Mh_MhFV:
    input:
        fastqc = expand('output/02_fastqc/{sample}_r{n}_fastqc.zip', sample=all_samples, n=[1,2])
    output:
        'output/03_multiqc_Mh_MhFV/multiqc_report.html'
    params:
        outdir = 'output/03_multiqc_Mh_MhFV',
        indirs = ['output/02_fastqc']
    log:
        'output/logs/03_multiqc_Mh_MhFV.log'
    container:
        multiqc_container
    shell:
        'nice multiqc '
        '-f '
        '-o {params.outdir} '
        '{params.indirs} '
        '2> {log}'

rule fastqc:
    input:
        'output/01_bbduk_trim/{sample}_r{n}.fq.gz'
    output:
        'output/02_fastqc/{sample}_r{n}_fastqc.zip'
    params:
        outdir = directory('output/02_fastqc')
    log:
        'output/logs/fastqc_{sample}_r{n}.log'
    singularity:
        fastqc_container
    threads:
        10
    shell:
        'mkdir -p {params.outdir} ; '
        'nice fastqc '
        '--outdir {params.outdir} '
        '{input} '
        '-t {threads} '
        '2> {log}'

################
## prep reads ##
################

##trim adapters and low quality seq.
rule bbduk_trim:
    input:
        r1 = 'output/00_joined/{sample}_r1.fq.gz',
        r2 = 'output/00_joined/{sample}_r2.fq.gz'
    output:
        r1 = 'output/01_bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/01_bbduk_trim/{sample}_r2.fq.gz',
        stats = 'output/01_bbduk_trim/stats/{sample}_stats.txt'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/01_bbduk_trim/{sample}.log'
    singularity:
        bbduk_container
    threads:
        20
    shell:
        'nice bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo forcetrimmod=5 qtrim=r trimq=15 '
        'stats={output.stats} '
        'threads={threads} '
        '&> {log}'
