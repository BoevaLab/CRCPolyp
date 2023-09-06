IDS, = glob_wildcards("fastq_files/{id}_R2_001.fastq.gz")

rule all:
    input:
        expand("raw_counts/{id}/quant.sf", id=IDS)
    script:
        "echo {input}"

rule test_salmon_decoy:
    input:
        transcriptome="gencode.v44.transcripts.fa.gz",
        genome="GRCh38.primary_assembly.genome.fa.gz",
    output:
        gentrome="gentrome.fasta.gz",
        decoys="decoys.txt",
    threads: 2
    log:
        "decoys.log"
    wrapper:
        "file://biowrapper/salmon/decoys"

rule salmon_index:
    input:
        sequences="gentrome.fasta.gz",
    output:
        multiext(
            "salmon/transcriptome_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
    log:
        "logs/salmon/transcriptome_index.log",
    threads: 2
    params:
        # optional parameters
        extra="",
    wrapper:
        "file://biowrapper/salmon/index"

rule trimmomatic_pe:
    input:
        r1="fastq_files/{sample}_R1_001.fastq.gz",
        r2="fastq_files/{sample}_R2_001.fastq.gz",
    output:
        r1="trimmed/{sample}.1.fastq.gz",
        r2="trimmed/{sample}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired="trimmed/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="trimmed/{sample}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:TrueSeq3-PE.fa:2:30:10","ILLUMINACLIP:NexteraPE-PE.fa:2:30:10","LEADING:3","TRAILING:3","SLIDINGWINDOW:4:15","MINLEN:30"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        32
    wrapper:
        "file://biowrapper/trimmomatic"

rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1="trimmed/{sample}.1.fastq.gz",
        r2="trimmed/{sample}.2.fastq.gz",
        index="salmon/transcriptome_index",
    output:
        quant="raw_counts/{sample}/quant.sf",
        lib="raw_counts/{sample}/lib_format_counts.json",
    log:
        "logs/salmon/{sample}.log",
    resources:
        mem_mb=16000
    params:
        # optional parameters
        libtype="A",
        extra="",
    threads: 2
    wrapper:
        "file://biowrapper/salmon/quant"
