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
frag_path = config['frag_path']

## test


## constructe the target if the inputs are fastqs
bam = expand("01_bam/{sample}.bam", sample = SAMPLES)
# raw_pairsam = expand("pairs-hg38/{sample}/{sample}.raw.pairsam.gz", sample = SAMPLES)
selected_pairsam = expand("pairs-{genome}/{sample}/{sample}.selected.pairsam.gz", sample = SAMPLES, genome = genome)
cool = expand("coolers-hg38/{sample}.cool", sample = SAMPLES)
TARGETS.extend(bam) ##append all list to 
TARGETS.extend(selected_pairsam) ## check later
TARGETS.extend(cool)



localrules: all

rule all:
    input: TARGETS




rule bwa_align:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: "01_bam/{sample}.bam"
    threads: 10
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
    output: temp("pairs-hg38/{sample}/{sample}.raw.pairsam.gz")
    message: "prase bam {input} "
    threads: 5
    shell:
        """
        pairtools parse -c {chromsizes}  \
        --assembly {genome} --min-mapq 1 \
        --max-molecule-size 2000 --max-inter-align-gap 20 --drop-readid \
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

rule restrict_pairsam:
    input:  "pairs-hg38/{sample}/{sample}.selected.pairsam.gz"
    output: "pairs-hg38/{sample}/{sample}.pairsam.gz", "pairs-hg38/{sample}/{sample}.select.samefrag.pairsam.gz"
    message: "selected pairsam {input} "
    threads: 5
    shell:
        """
        pairtools restrict -f {frag_path} {input} |  \
        pairtools select '(COLS[-6]==COLS[-3]) and (chrom1==chrom2)' \
        --output-rest {output[0]} -o {output[1]} 
        """

rule deduplicated_pairsam:
    input:  "pairs-hg38/{sample}/{sample}.pairsam.gz"
    output: "filtered-hg38/{sample}.filtered.pairsam.gz"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
         pairtools dedup --max-mismatch 1 --method max    -o {output[0]} {input}
        """

rule flip_pairsam:
    input:  "filtered-hg38/{sample}.filtered.pairsam.gz"
    output: "filtered-hg38/{sample}.filtered.flip.pairs.gz"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
         pairtools flip -c {chromsizes}   -o {output[0]} {input}
        """
        
rule sort_pairsam:
    input:  "filtered-hg38/{sample}.filtered.flip.pairs.gz"
    output: "filtered-hg38/{sample}.filtered.flip.sorted.pairs.gz"
    message: "sort pairsam {input} "
    threads: 8
    shell:
        """
        pairtools sort  --nproc 8  --memory 20G  -o {output} {input}  
        """
        
rule index_pairsam:
    input:  "filtered-hg38/{sample}.filtered.flip.sorted.pairs.gz"
    output: "filtered-hg38/{sample}.filtered.flip.sorted.pairs.gz.px2"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
         pairix -p pairs     {input}
        """
        
rule load_cooler:
    input:  "filtered-hg38/{sample}.filtered.flip.sorted.pairs.gz", "filtered-hg38/{sample}.filtered.flip.sorted.pairs.gz.px2"
    output: "coolers-hg38/{sample}.cool"
    message: "cooler {input} "
    threads: 10
    shell:
        """
        cooler cload pairix --assembly hg38 --nproc {threads} \
                   --max-split 2 {chromsizes}:10000 {input[0]} {output}
        """

