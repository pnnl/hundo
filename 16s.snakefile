import os
from subprocess import check_output


def get_samples(eid):
    samples = set()
    for f in os.listdir(os.path.join("results", eid, "demux")):
        samples.add(f.partition(".")[0].partition("_")[0])
    return samples


def fix_tax_entry(tax, kingdom=""):
    """
    >>> t = "p:Basidiomycota,c:Tremellomycetes,o:Tremellales,f:Tremellales_fam_Incertae_sedis,g:Cryptococcus"
    >>> fix_tax_entry(t)
    'k__,p__Basidiomycota,c__Tremellomycetes,o__Tremellales,f__Tremellales_fam_Incertae_sedis,g__Cryptococcus,s__'
    """

    if tax == "" or tax == "*":
        taxonomy = dict()
    else:
        taxonomy = dict(x.split(":") for x in tax.split(","))
    if "d" in taxonomy:
        taxonomy["k"] = taxonomy["d"]
    else:
        taxonomy["k"] = kingdom

    # add the underscores
    underscores = dict()
    for k, v in taxonomy.items():
        v = v.strip('"')
        underscores[k] = "%s__%s" % (k, v)

    # fill in missing
    for idx in "kpcofgs":
        if idx not in underscores:
            underscores[idx] = "%s__" % idx

    # build needed string
    new_taxonomy = ""
    for idx in "kpcofgs":
        new_taxonomy += underscores[idx] + ","

    return new_taxonomy.strip(",")


def fix_fasta_tax_entry(tax, kingdom=""):
    """
    >>> t = ">OTU_7;tax=p:Basidiomycota,c:Microbotryomycetes,o:Sporidiobolales,f:Sporidiobolales_fam_Incertae_sedis,g:Rhodotorula;"
    >>> fix_fasta_tax_entry(t)
    '>OTU_7;tax=k__,p__Basidiomycota,c__Microbotryomycetes,o__Sporidiobolales,f__Sporidiobolales_fam_Incertae_sedis,g__Rhodotorula,s__;'
    """
    toks = tax.split(";")
    otu = toks[0]
    tax_piece = toks[1]
    if not tax_piece.startswith("tax"):
        raise ValueError
    sequence_tax = tax_piece.split("=")[1]
    new_tax = fix_tax_entry(sequence_tax, kingdom)
    return "%s;tax=%s;" % (toks[0], new_tax)


configfile: "config.yaml"
ruleorder: utax > fix_utax_taxonomy
USEARCH_VERSION = check_output("usearch --version", shell=True).strip()
CLUSTALO_VERSION = check_output("/people/brow015/apps/cbb/clustalo/1.2.0/clustalo --version", shell=True).strip()
EID = config['eid']
SAMPLES = get_samples(EID)
# name output folder appropriately
CLUSTER_THRESHOLD = 100 - config['clustering']['percent_of_allowable_difference']


rule all:
    input:
        expand("results/{eid}/logs/{pid}/quality_filtering_stats.txt", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/demux/{sample}_R1.fastq", eid=EID, sample=SAMPLES),
        expand("results/{eid}/demux/{sample}_R2.fastq", eid=EID, sample=SAMPLES),
        expand("results/{eid}/{pid}/OTU.txt", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/OTU.biom", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/OTU_tax.txt", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/OTU.tree", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/unfiltered/OTU.txt", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/unfiltered/OTU.biom", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/unfiltered/OTU_tax.txt", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/unfiltered/OTU.tree", eid=EID, pid=CLUSTER_THRESHOLD)


rule make_reference_database:
    input:
        fasta = "ref/rdp_16s_trainset15/fasta/refdb.fa",
        trained_parameters = "ref/rdp_16s_trainset15/taxconfs/250.tc"
    output:
        "ref/rdp_16s_trainset15/udb/250.udb"
    version: USEARCH_VERSION
    message: "Creating a UTAX database based on RDP trainset version 15 trained on 250 nt sequence length"
    shell:
        '''
        usearch -makeudb_utax {input.fasta} -taxconfsin {input.trained_parameters} -output {output}
        '''


rule make_uchime_reference:
    input: "ref/chimera/rdp_gold.fa"
    output: "ref/chimera/udb/rdp_gold.udb"
    version: USEARCH_VERSION
    message: "Creating chimera reference index based on RDP classifier training database (v9) which contains 10,049 reference sequences."
    shell: "usearch -makeudb_usearch {input} -output {output}"


rule quality_filter_reads:
    input:
        r1 = "results/{eid}/demux/{sample}_R1.fastq",
        r2 = "results/{eid}/demux/{sample}_R2.fastq"
    output:
        r1 = temp("results/{eid}/filtered/{sample}_R1.fastq"),
        r2 = temp("results/{eid}/filtered/{sample}_R2.fastq"),
        stats = temp("results/{eid}/{sample}_quality_filtering_stats.txt")
    message: "Filtering reads using BBDuk2 to remove adapters and phiX with matching kmer length of {params.k} at a hamming distance of {params.hdist} and quality trim both ends to Q{params.quality}. Reads shorter than {params.minlength} were discarded."
    params:
        adapters = "ref/phix174_ill.ref.fa.gz,ref/adapters.fa",
        quality = config['filtering']['minimum_base_quality'],
        hdist = config['filtering']['allowable_kmer_mismatches'],
        k = config['filtering']['reference_kmer_match_length'],
        qtrim = "rl",
        ktrim = "l",
        minlength = config['filtering']['minimum_passing_read_length']
    threads: 4
    shell:
        '''
        bbduk2.sh -Xmx8g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} \
            fref={params.adapters} stats={output.stats} hdist={params.hdist} k={params.k} \
            trimq={params.quality} qtrim={params.qtrim} threads={threads} ktrim={params.ktrim} \
            minlength={params.minlength}
        '''


rule combine_filtering_stats:
    input: expand("results/{eid}/{sample}_quality_filtering_stats.txt", eid=EID, sample=SAMPLES)
    output: "results/{eid}/logs/quality_filtering_stats.txt"
    shell: "cat {input} > {output}"


rule merge_reads:
    input:
        r1 = "results/{eid}/filtered/{sample}_R1.fastq",
        r2 = "results/{eid}/filtered/{sample}_R2.fastq"
    output: temp("results/{eid}/merged/{sample}.fastq")
    version: USEARCH_VERSION
    message: "Merging paired-end reads with USEARCH at a minimum merge length of {params.minimum_merge_length}"
    params:
        minimum_merge_length = config['merging']['minimum_merge_length']
    shell:
        '''
        samplename={wildcards.sample}
        usearch -fastq_mergepairs {input.r1} -relabel ${{samplename//-/}}. \
            -fastq_minmergelen {params.minimum_merge_length} \
            -fastqout {output}
        '''


rule combine_merged_reads:
    input: expand("results/{eid}/merged/{sample}.fastq", eid=EID, sample=SAMPLES)
    output: "results/{eid}/merged.fastq"
    message: "Concatenating the merged reads into a single file"
    shell: "cat {input} > {output}"


rule fastq_filter:
    input: "results/{eid}/merged.fastq"
    output: temp("results/{eid}/merged_filtered.fasta")
    version: USEARCH_VERSION
    message: "Filtering FASTQ with USEARCH with an expected maximum error rate of {params.maxee}"
    params:
        maxee = config['filtering']['maximum_expected_error']
    shell: "usearch -fastq_filter {input} -fastq_maxee {params.maxee} -fastaout {output} -relabel Filt"


rule dereplicate_sequences:
    input: "results/{eid}/merged_filtered.fasta"
    output: temp("results/{eid}/uniques.fasta")
    version: USEARCH_VERSION
    message: "Dereplicating with USEARCH"
    shell: "usearch -derep_fulllength {input} -sizeout -relabel Uniq -fastaout {output}"


rule cluster_sequences:
    input: "results/{eid}/uniques.fasta"
    output: "results/{eid}/{pid}/unfiltered/OTU.fasta"
    version: USEARCH_VERSION
    message: "Clustering sequences with USEARCH where OTUs have a minimum size of {params.minsize} and where the maximum difference between an OTU member sequence and the representative sequence of that OTU is {params.otu_radius_pct}%"
    params:
        minsize = config['filtering']['minimum_sequence_abundance'],
        otu_radius_pct = config['clustering']['percent_of_allowable_difference']
    shell:
        '''
        usearch -cluster_otus {input} -minsize {params.minsize} -otus {output} -relabel OTU_ \
            -otu_radius_pct {params.otu_radius_pct}
        '''


rule remove_chimeric_otus:
    input:
        fasta = "results/{eid}/{pid}/unfiltered/OTU.fasta",
        reference = "ref/chimera/udb/rdp_gold.udb"
    output: "results/{eid}/{pid}/OTU.fasta"
    version: USEARCH_VERSION
    message: "Chimera filtering OTU seed sequences against RDP classifier training database (v9), canonically known as rdp_gold.fa"
    threads: 12
    log: "results/{eid}/{pid}/logs/uchime_ref.log"
    shell:
        '''usearch -uchime_ref {input.fasta} -db {input.reference} -uchimeout {output} \
            -strand plus -threads {threads} 2> {log}
        '''


rule utax:
    input:
        fasta = "results/{eid}/{pid}/OTU.fasta",
        db = "ref/rdp_16s_trainset15/udb/250.udb"
    output:
        fasta = "results/{eid}/{pid}/utax/OTU_tax.fasta",
        txt = "results/{eid}/{pid}/utax/OTU_tax.txt"
    version: USEARCH_VERSION
    message: "Assigning taxonomies with UTAX algorithm using USEARCH with a confidence cutoff of {params.utax_cutoff}"
    params:
        utax_cutoff = config['taxonomy']['prediction_confidence_cutoff']
    threads: 12
    log: "results/{eid}/{pid}/logs/utax.log"
    shell:
        '''
        usearch -utax {input.fasta} -db {input.db} -strand both -threads {threads} \
            -fastaout {output.fasta} -utax_cutoff {params.utax_cutoff} \
            -utaxout {output.txt} 2> {log}
        '''


rule utax_unfiltered:
    input:
        fasta = "results/{eid}/{pid}/unfiltered/OTU.fasta",
        db = "ref/rdp_16s_trainset15/udb/250.udb"
    output:
        fasta = "results/{eid}/{pid}/unfiltered/utax/OTU_tax.fasta",
        txt = "results/{eid}/{pid}/unfiltered/utax/OTU_tax.txt"
    version: USEARCH_VERSION
    message: "Assigning taxonomies with UTAX algorithm using USEARCH with a confidence cutoff of {params.utax_cutoff}"
    params:
        utax_cutoff = config['taxonomy']['prediction_confidence_cutoff']
    threads: 12
    log: "results/{eid}/{pid}/logs/utax.log"
    shell:
        '''
        usearch -utax {input.fasta} -db {input.db} -strand both -threads {threads} \
            -fastaout {output.fasta} -utax_cutoff {params.utax_cutoff} \
            -utaxout {output.txt} 2> {log}
        '''


rule fix_utax_taxonomy:
    input:
        fasta = "results/{eid}/{pid}/utax/OTU_tax.fasta",
        txt = "results/{eid}/{pid}/utax/OTU_tax.txt"
    output:
        fasta = "results/{eid}/{pid}/OTU_tax.fasta",
        txt = "results/{eid}/{pid}/OTU_tax.txt"
    message: "Altering taxa to reflect QIIME style annotation"
    run:
        with open(input.fasta) as ifh, open(output.fasta, 'w') as ofh:
            for line in ifh:
                line = line.strip()
                if not line.startswith(">OTU_"):
                    print(line, file=ofh)
                else:
                    print(fix_fasta_tax_entry(line), file=ofh)
        with open(input.txt) as ifh, open(output.txt, 'w') as ofh:
            for line in ifh:
                toks = line.strip().split("\t")
                print(toks[0], fix_tax_entry(toks[1]), fix_tax_entry(toks[2]), toks[3], sep="\t", file=ofh)


rule fix_utax_taxonomy_unfiltered:
    input:
        fasta = "results/{eid}/{pid}/unfiltered/utax/OTU_tax.fasta",
        txt = "results/{eid}/{pid}/unfiltered/utax/OTU_tax.txt"
    output:
        fasta = "results/{eid}/{pid}/unfiltered/OTU_tax.fasta",
        txt = "results/{eid}/{pid}/unfiltered/OTU_tax.txt"
    message: "Altering taxa to reflect QIIME style annotation"
    run:
        with open(input.fasta) as ifh, open(output.fasta, 'w') as ofh:
            for line in ifh:
                line = line.strip()
                if not line.startswith(">OTU_"):
                    print(line, file=ofh)
                else:
                    print(fix_fasta_tax_entry(line), file=ofh)
        with open(input.txt) as ifh, open(output.txt, 'w') as ofh:
            for line in ifh:
                toks = line.strip().split("\t")
                print(toks[0], fix_tax_entry(toks[1]), fix_tax_entry(toks[2]), toks[3], sep="\t", file=ofh)


rule compile_counts:
    input:
        fastq = "results/{eid}/merged.fastq",
        db = "results/{eid}/{pid}/OTU_tax.fasta",
    output:
        txt = "results/{eid}/{pid}/OTU.txt",
    params:
        threshold = config['mapping_to_otus']['read_identity_requirement']
    threads: 12
    shell:
        '''
        usearch -usearch_global {input.fastq} -db {input.db} -strand plus \
            -id {params.threshold} -otutabout {output.txt} \
            -threads {threads}
        '''


rule compile_counts_unfiltered:
    input:
        fastq = "results/{eid}/merged.fastq",
        db = "results/{eid}/{pid}/unfiltered/OTU_tax.fasta",
    output:
        txt = "results/{eid}/{pid}/unfiltered/OTU.txt",
    params:
        threshold = config['mapping_to_otus']['read_identity_requirement']
    threads: 12
    shell:
        '''
        usearch -usearch_global {input.fastq} -db {input.db} -strand plus \
            -id {params.threshold} -otutabout {output.txt} \
            -threads {threads}
        '''


rule biom:
    input:
        "results/{eid}/{pid}/OTU.txt"
    output:
        "results/{eid}/{pid}/OTU.biom"
    shell:
        '''
        sed 's|\"||g' {input} | sed 's|\,|\;|g' > results/{wildcards.eid}/{wildcards.pid}/OTU_converted.txt
        biom convert -i results/{wildcards.eid}/{wildcards.pid}/OTU_converted.txt \
            -o results/{wildcards.eid}/{wildcards.pid}/OTU_converted.biom --to-json \
            --process-obs-metadata sc_separated
        sed 's|\"type\": \"Table\"|\"type\": \"OTU table\"|g' results/{wildcards.eid}/{wildcards.pid}/OTU_converted.biom > {output}
        rm results/{wildcards.eid}/{wildcards.pid}/OTU_converted.txt results/{wildcards.eid}/{wildcards.pid}/OTU_converted.biom
        '''


rule biom_unfiltered:
    input:
        "results/{eid}/{pid}/unfiltered/OTU.txt"
    output:
        "results/{eid}/{pid}/unfiltered/OTU.biom"
    shell:
        '''
        sed 's|\"||g' {input} | sed 's|\,|\;|g' > results/{wildcards.eid}/{wildcards.pid}/unfiltered/OTU_converted.txt
        biom convert -i results/{wildcards.eid}/{wildcards.pid}/unfiltered/OTU_converted.txt \
            -o results/{wildcards.eid}/{wildcards.pid}/unfiltered/OTU_converted.biom --to-json \
            --process-obs-metadata sc_separated
        sed 's|\"type\": \"Table\"|\"type\": \"OTU table\"|g' results/{wildcards.eid}/{wildcards.pid}/unfiltered/OTU_converted.biom > {output}
        rm results/{wildcards.eid}/{wildcards.pid}/unfiltered/OTU_converted.txt results/{wildcards.eid}/{wildcards.pid}/unfiltered/OTU_converted.biom
        '''


rule multiple_align:
    input: "results/{eid}/{pid}/OTU.fasta"
    output: "results/{eid}/{pid}/OTU_aligned.fasta"
    message: "Multiple alignment of samples using Clustal Omega"
    version: CLUSTALO_VERSION
    threads: 12
    shell:
        '''
        /people/brow015/apps/cbb/clustalo/1.2.0/clustalo -i {input} -o {output} \
            --outfmt=fasta --threads {threads} --force
        '''


rule multiple_align_unfiltered:
    input: "results/{eid}/{pid}/unfiltered/OTU.fasta"
    output: "results/{eid}/{pid}/unfiltered/OTU_aligned.fasta"
    message: "Multiple alignment of samples using Clustal Omega"
    version: CLUSTALO_VERSION
    threads: 12
    shell:
        '''
        /people/brow015/apps/cbb/clustalo/1.2.0/clustalo -i {input} -o {output} \
            --outfmt=fasta --threads {threads} --force
        '''


rule newick_tree:
    input: "results/{eid}/{pid}/OTU_aligned.fasta"
    output: "results/{eid}/{pid}/OTU.tree"
    message: "Building tree from aligned OTU sequences with FastTree2"
    log: "results/{eid}/{pid}/logs/fasttree.log"
    shell: "FastTree -nt -gamma -spr 4 -log {log} -quiet {input} > {output}"


rule newick_tree_unfiltered:
    input: "results/{eid}/{pid}/unfiltered/OTU_aligned.fasta"
    output: "results/{eid}/{pid}/unfiltered/OTU.tree"
    message: "Building tree from aligned OTU sequences with FastTree2"
    log: "results/{eid}/{pid}/logs/fasttree.log"
    shell: "FastTree -nt -gamma -spr 4 -log {log} -quiet {input} > {output}"


"""

"""
