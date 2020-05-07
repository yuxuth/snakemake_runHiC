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
cool_bin = config['cool_bin']

## test


## constructe the target if the inputs are fastqs
# bam = expand("01_bam/{sample}.bam", sample = SAMPLES)
# raw_pairsam = expand("pairs-hg38/{sample}/{sample}.raw.pairsam.gz", sample = SAMPLES)

# peak_pairs = expand("peaks-{genome}/{sample}.final.pairs.gz", sample = SAMPLES, genome = genome)

# all_pairs = expand("pairs-{genome}/{sample}.valid.pairs.gz", sample = SAMPLES, genome = genome)
deup_pairs = expand("filtered-{genome}/{sample}.dedup.pairs.gz", sample = SAMPLES, genome = genome)
# cool = expand("coolers-{genome}/{sample}.{cool_bin}.cool", sample = SAMPLES, cool_bin = cool_bin, genome = genome)

# TARGETS.extend(bam) ##append all list to 
# TARGETS.extend(peak_pairs) ## check later
# cool40k = expand("coolers-hg38/{sample}.40k.cool", sample = SAMPLES) ## will lead to keyerror
# cool40k = expand("coolers-hg38_40k/{sample}.cool", sample = SAMPLES)


# TARGETS.extend(all_pairs)
# TARGETS.extend(cool)
TARGETS.extend(deup_pairs)

stat1 = expand("pairs-{genome}/{sample}.raw.pairsam.stat", sample = SAMPLES, genome = genome)
stat2 = expand("filtered-{genome}/{sample}.dedup.pairs.stat", sample = SAMPLES, genome = genome)
TARGETS.extend(stat1)
TARGETS.extend(stat2)

localrules: all

rule all:
    input: TARGETS




rule bwa_align:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: "01_bam/{sample}.bam"
    threads: 24
    message: "bwa {input}: {threads} threads"
    log:
         "00_log/{sample}.bwa"
    shell:
        """
        module load bwa
        module load samtools
        bwa mem  -SP -t {threads} {BWA_INDEX} {input} | samtools view -bS - > {output}  2> {log}
        """


rule prase_bam_no_flip: ## no flip to makesure the R1 R2 position for the peak calling
    input:  "01_bam/{sample}.bam"
    output: ("pairs-{genome}/{sample}.raw.pairsam.gz")
    message: "prase bam {input} "
    threads: 5
    shell:
        """
        pairtools parse -c {chromsizes}  \
        --assembly {genome} --min-mapq 1 \
        --max-molecule-size 2000 --max-inter-align-gap 20 --drop-readid \
        --walks-policy mask --no-flip --drop-seq --drop-sam  -o {output} {input}  
        """

rule stat1 :
    input:  "pairs-{genome}/{sample}.raw.pairsam.gz"
    output: "pairs-{genome}/{sample}.raw.pairsam.stat"
    message: "stat1 bam {input} "
    threads: 5
    shell:
        """
        pairtools stats -o {output} {input}  
        """


rule flip_pairsam_sort:
    input:  "pairs-{genome}/{sample}.raw.pairsam.gz"
    output: "pairs-{genome}/{sample}.sorted.pairs.gz"
    message: "flip to filted {input} "
    threads: 8
    shell:
        """
         pairtools flip -c {chromsizes} {input} | \
         pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' | \
         pairtools sort  --nproc 8  --memory 15G  -o {output}
        """

rule dedup_pairs:
    input:  "pairs-{genome}/{sample}.sorted.pairs.gz"
    output: "filtered-{genome}/{sample}.dedup.pairs.gz"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
         pairtools dedup --max-mismatch 1 --method max -o {output[0]} {input}
        """

rule stat2 :
    input:  "filtered-{genome}/{sample}.dedup.pairs.gz"
    output: "filtered-{genome}/{sample}.dedup.pairs.stat"
    message: "stat1 bam {input} "
    threads: 5
    shell:
        """
        pairtools stats -o {output} {input}  
        """


rule enzyme_fragment_detection:
    input:  "pairs-{genome}/{sample}.selected.pairs.gz"
    output: valid = "pairs-{genome}/{sample}.valid.pairs.gz", same_f = "pairs-{genome}/{sample}.samefrag.pairs.gz"
    message: "selected pairsam {input} "
    threads: 5
    shell:
        """
        pairtools restrict -f {frag_path} {input} |  \
        pairtools select '(COLS[-6]==COLS[-3]) and (chrom1==chrom2)' \
        --output-rest {output[valid]} -o {output[same_f]} 
        """







rule R2peak_pairs_enzyme_fragment_dedup:
    input:  "peaks-{genome}/{sample}.sorted.pairs.gz"
    output: "peaks-{genome}/{sample}.final.pairs.gz"
    message: "flip to filted {input} "
    threads: 8
    shell:
        """
         pairtools dedup --max-mismatch 1 --method max {input} | \
         pairtools restrict -f {frag_path} -o {output}
        """

rule R2peak_pairs_sort_mapping_filter:
    input:  "pairs-{genome}/{sample}.raw.pairsam.gz"
    output: "peaks-{genome}/{sample}.sorted.pairs.gz"
    message: "flip to filted {input} "
    threads: 8
    shell:
        """
         pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' {input} | \
         pairtools sort  --nproc 8  --memory 15G  -o {output}
        """

# rule sort_pairsam:
#     input:  "pairs-{genome}/{sample}/{sample}.raw.pairsam.gz"
#     output: "pairs-{genome}/{sample}/{sample}.sort.pairsam.gz"
#     message: "sort pairsam {input} "
#     threads: 8
#     shell:
#         """
#         pairtools sort  --nproc 8  --memory 20G  -o {output} {input}  
#         """



        
rule index_pairs:
    input:  "filtered-{genome}/{sample}.valid.pairs.gz"
    output: "filtered-{genome}/{sample}.valid.pairs.gz.px2"
    message: "dedup to filted {input} "
    threads: 5
    shell:
        """
         pairix -p pairs  {input}
        """
        
rule load_cooler:
    input:  "filtered-{genome}/{sample}.valid.pairs.gz", "filtered-{genome}/{sample}.valid.pairs.gz.px2"
    output: "coolers-{genome}/{sample}.{cool_bin}.cool"
    message: "cooler {input} "
    params: res = {cool_bin}
    threads: 10
    shell:
        """
        cooler cload pairix --assembly hg38 --nproc {threads} \
                   --max-split 2 {chromsizes}:{params.res} {input[0]} {output}
        """
      
                 
        
# rule load_cooler_40k:
#     input:  "filtered-hg38/{sample}.filtered.flip.sorted.pairs.gz", "filtered-hg38/{sample}.filtered.flip.sorted.pairs.gz.px2"
#     output: "coolers-hg38_40k/{sample}.cool"
#     message: "cooler {input} "
#     threads: 10
#     shell:
#         """
#         cooler cload pairix --assembly hg38 --nproc {threads} \
#                    --max-split 2 {chromsizes}:40000 {input[0]} {output}
#         """

