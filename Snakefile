shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())
BWA_INDEX = config['BWA_INDEX']
chromsizes = config['chromsizes']

TARGETS = []
genome = config['genome']

## test


## constructe the target if the inputs are fastqs
bam = expand("01_bam/{sample}.bam", sample = SAMPLES)
# raw_pairsam = expand("pairs-hg38/{sample}/{sample}.raw.pairsam.gz", sample = SAMPLES)
selected_pairsam = expand("pairs-{genome}/{sample}/{sample}.selected.pairsam.gz", sample = SAMPLES, genome = genome)

TARGETS.extend(bam) ##append all list to 
TARGETS.extend(selected_pairsam) ## check later




localrules: all

rule all:
    input: TARGETS




rule bwa_align:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: "01_bam/{sample}.bam"
    threads: 12
    message: "bwa {input}: {threads} threads"
    log:
         "00_log/{sample}.bwa"
    shell:
        """
        module load bwa
        module load samtools
        bwa mem  -SP -t {threads} {BWA_INDEX} {input} | samtools view -bS - > {output}  2> {log}
        """


rule prase_bam:
    input:  "01_bam/{sample}.bam"
    output: "pairs-hg38/{sample}/{sample}.raw.pairsam.gz"
    message: "prase bam {input} "
    threads: 5
    shell:
        """
        pairtools parse -c {chromsizes}  \
        --assembly {genome} --min-mapq 1 \
        --max-molecule-size 2000 --max-inter-align-gap 20 \
        --walks-policy mask --no-flip --drop-seq --drop-sam  -o {output} {input}
        """

rule sort_pairsam:
    input:  "pairs-hg38/{sample}/{sample}.raw.pairsam.gz"
    output: "pairs-hg38/{sample}/{sample}.sort.pairsam.gz"
    message: "sort pairsam {input} "
    threads: 8
    shell:
        """
        pairtools sort  --nproc 8  --memory 20G  -o {output} {input}  
        """

rule seleted_pairsam:
    input:  "pairs-hg38/{sample}/{sample}.sort.pairsam.gz"
    output: "pairs-hg38/{sample}/{sample}.selected.pairsam.gz"
    message: "selected pairsam {input} "
    threads: 5
    shell:
        """
        pairtools  select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' -o  {output} {input}  
        """


