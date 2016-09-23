"""
Notes pertaining to eventual addition of Phinch compatible biom file...
otu.txt , commas to semicolons
biom convert -i otu.nocommas.txt -o some.biom --to-json
biom add-metadata -i some.biom -o phinch.biom --sample-metadata-fp tsv --output-as-json
#SampleID phinchID
"""
import os
from snakemake.utils import report
from subprocess import check_output


def read_count(fastq):
    total = 0
    count_file = fastq + '.count'
    if os.path.exists(fastq) and os.path.getsize(fastq) > 100:
        if not os.path.exists(count_file):
            check_output("awk '{n++}END{print n/4}' %s > %s" % (fastq, fastq + '.count'), shell=True)
        with open(count_file) as fh:
            for line in fh:
                total = int(line.strip())
                break
    return total


def get_samples(eid):
    samples = set()
    omitted = set()
    input_dir = os.path.join("results", eid, "demux")
    for f in os.listdir(input_dir):
        if (f.endswith("fastq") or f.endswith("fq")) and ("_r1" in f or "_R1" in f):
            if read_count(os.path.join(input_dir, f)) > 1000:
                samples.add(f.partition(".")[0].partition("_")[0])
            else:
                omitted.add(f.partition(".")[0].partition("_")[0])
    return samples, omitted


def fix_tax_entry(tax, kingdom="?"):
    """
    >>> t = "p:Basidiomycota,c:Tremellomycetes,o:Tremellales,f:Tremellales_fam_Incertae_sedis,g:Cryptococcus"
    >>> fix_tax_entry(t, "Fungi")
    'k__Fungi,p__Basidiomycota,c__Tremellomycetes,o__Tremellales,f__Tremellales_fam_Incertae_sedis,g__Cryptococcus,s__?'
    """
    if tax == "" or tax == "*":
        taxonomy = dict()
    else:
        taxonomy = dict(x.split(":") for x in tax.split(","))
    if "d" in taxonomy and not "k" in taxonomy:
        taxonomy["k"] = taxonomy["d"]
    else:
        taxonomy["k"] = kingdom

    new_taxonomy = []
    for idx in "kpcofgs":
        new_taxonomy.append("%s__%s" % (idx, taxonomy.get(idx, "?")))
    return ",".join(new_taxonomy)


def fix_fasta_tax_entry(tax, kingdom="?"):
    """
    >>> t = ">OTU_7;tax=p:Basidiomycota,c:Microbotryomycetes,o:Sporidiobolales,f:Sporidiobolales_fam_Incertae_sedis,g:Rhodotorula;"
    >>> fix_fasta_tax_entry(t)
    '>OTU_7;tax=k__?,p__Basidiomycota,c__Microbotryomycetes,o__Sporidiobolales,f__Sporidiobolales_fam_Incertae_sedis,g__Rhodotorula,s__?;'
    """
    toks = tax.split(";")
    otu = toks[0]
    tax_piece = toks[1]
    if not tax_piece.startswith("tax"):
        raise ValueError
    sequence_tax = tax_piece.split("=")[1]
    new_tax = fix_tax_entry(sequence_tax, kingdom)
    return "%s;tax=%s;" % (toks[0], new_tax)


USEARCH_VERSION = check_output("usearch --version", shell=True).strip()
CLUSTALO_VERSION = check_output("clustalo --version", shell=True).strip()
EID = config['eid']
SAMPLES, OMITTED = get_samples(EID)
# name output folder appropriately
CLUSTER_THRESHOLD = 100 - config['clustering']['percent_of_allowable_difference']


rule all:
    input:
        expand("results/{eid}/logs/quality_filtering_stats.txt", eid=EID),
        expand("results/{eid}/demux/{sample}_R1.fastq", eid=EID, sample=SAMPLES),
        expand("results/{eid}/demux/{sample}_R2.fastq", eid=EID, sample=SAMPLES),
        expand("results/{eid}/logs/{sample}_R1.fastq.count", eid=EID, sample=SAMPLES),
        expand("results/{eid}/logs/{sample}_filtered_R1.fastq.count", eid=EID, sample=SAMPLES),
        expand("results/{eid}/logs/{sample}_merged.fastq.count", eid=EID, sample=SAMPLES),
        expand("results/{eid}/{pid}/OTU.biom", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/OTU.tree", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/utax/OTU.biom", eid=EID, pid=CLUSTER_THRESHOLD),
        expand("results/{eid}/{pid}/README.html", eid=EID, pid=CLUSTER_THRESHOLD)


rule make_tax_database:
    input:
        fasta = config['taxonomy_database']['fasta'],
        trained_parameters = config['taxonomy_database']['trained_parameters']
    output:
        os.path.splitext(config['taxonomy_database']['trained_parameters'])[0] + '.udb'
    version: USEARCH_VERSION
    message: "Creating a UTAX database trained on {input.fasta} using {input.trained_parameters}"
    shell:
        '''
        usearch -makeudb_utax {input.fasta} -taxconfsin {input.trained_parameters} -output {output}
        '''


rule make_uchime_database:
    input: config['chimera_database']['fasta']
    output: os.path.splitext(config['chimera_database']['fasta'])[0] + '.udb'
    version: USEARCH_VERSION
    message: "Creating chimera reference index based on {input}"
    shell: "usearch -makeudb_usearch {input} -output {output}"


rule make_blast_db:
    input: config['blast_database']['fasta']
    output: expand(config['blast_database']['fasta'] + ".{idx}", idx=['nhr', 'nin', 'nsq'])
    message: "Formatting BLAST database"
    shell: "makeblastdb -in {input} -dbtype nucl"


rule count_raw_reads:
    input:
        "results/{eid}/demux/{sample}_R1.fastq"
    output:
        "results/{eid}/logs/{sample}_R1.fastq.count"
    shell:
        "awk '{{n++}}END{{print n/4}}' {input} > {output}"


rule quality_filter_reads:
    input:
        r1 = "results/{eid}/demux/{sample}_R1.fastq",
        r2 = "results/{eid}/demux/{sample}_R2.fastq"
    output:
        r1 = temp("results/{eid}/{sample}_filtered_R1.fastq"),
        r2 = temp("results/{eid}/{sample}_filtered_R2.fastq"),
        stats = temp("results/{eid}/{sample}_quality_filtering_stats.txt")
    message: "Filtering reads using BBDuk2 to remove adapters and phiX with matching kmer length of {params.k} at a hamming distance of {params.hdist} and quality trim both ends to Q{params.quality}. Reads shorter than {params.minlength} were discarded."
    params:
        adapters = config['filtering']['adapters'],
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
            minlength={params.minlength} overwrite=true
        '''


rule count_filtered_reads:
    input:
        "results/{eid}/{sample}_filtered_R1.fastq"
    output:
        "results/{eid}/logs/{sample}_filtered_R1.fastq.count"
    shell:
        "awk '{{n++}}END{{print n/4}}' {input} > {output}"


rule combine_filtering_stats:
    input: expand("results/{eid}/{sample}_quality_filtering_stats.txt", eid=EID, sample=SAMPLES)
    output: "results/{eid}/logs/quality_filtering_stats.txt".format(eid=EID)
    shell: "cat {input} > {output}"


rule merge_reads:
    input:
        r1 = "results/{eid}/{sample}_filtered_R1.fastq",
        r2 = "results/{eid}/{sample}_filtered_R2.fastq"
    output: temp("results/{eid}/{sample}_merged.fastq")
    version: USEARCH_VERSION
    message: "Merging paired-end reads with USEARCH at a minimum merge length of {params.minimum_merge_length}"
    params:
        minimum_merge_length = config['merging']['minimum_merge_length']
    log: "results/{eid}/{pid}/logs/fastq_mergepairs.log"
    shell:
        '''
        usearch -fastq_mergepairs {input.r1} -relabel @ -sample {wildcards.sample} \
            -fastq_minmergelen {params.minimum_merge_length} \
            -fastqout {output} -log {log}
        '''


rule count_joined_reads:
    input:
        "results/{eid}/{sample}_merged.fastq"
    output:
        "results/{eid}/logs/{sample}_merged.fastq.count"
    shell:
        "awk '{{n++}}END{{print n/4}}' {input} > {output}"


rule combine_merged_reads:
    input: expand("results/{eid}/{sample}_merged.fastq", eid=EID, sample=SAMPLES)
    output: "results/{eid}/merged.fastq"
    message: "Concatenating the merged reads into a single file"
    shell: "cat {input} > {output}"


rule fastq_filter:
    input: "results/{eid}/merged.fastq"
    output: "results/{eid}/merged_%s.fasta" % str(config['filtering']['maximum_expected_error'])
    version: USEARCH_VERSION
    message: "Filtering FASTQ with USEARCH with an expected maximum error rate of {params.maxee}"
    params:
        maxee = config['filtering']['maximum_expected_error']
    log: "results/{eid}/{pid}/logs/fastq_filter.log"
    shell: "usearch -fastq_filter {input} -fastq_maxee {params.maxee} -fastaout {output} -relabel Filt -log {log}"


rule dereplicate_sequences:
    input: rules.fastq_filter.output
    output: temp("results/{eid}/uniques.fasta")
    version: USEARCH_VERSION
    message: "Dereplicating with USEARCH"
    threads: 22
    log: "results/{eid}/{pid}/logs/uniques.log"
    shell: "usearch -fastx_uniques {input} -fastaout {output} -sizeout -threads {threads} -log {log}"


rule cluster_sequences:
    input: "results/{eid}/uniques.fasta"
    output: temp("results/{eid}/{pid}/OTU_unfiltered.fasta")
    version: USEARCH_VERSION
    message: "Clustering sequences with USEARCH where OTUs have a minimum size of {params.minsize} and where the maximum difference between an OTU member sequence and the representative sequence of that OTU is {params.otu_radius_pct}%"
    params:
        minsize = config['clustering']['minimum_sequence_abundance'],
        otu_radius_pct = config['clustering']['percent_of_allowable_difference']
    log: "results/{eid}/{pid}/logs/cluster_sequences.log"
    shell:
        '''
        usearch -cluster_otus {input} -minsize {params.minsize} -otus {output} -relabel OTU_ \
            -otu_radius_pct {params.otu_radius_pct} -log {log}
        '''


# cluster
# -cluster_otus $out/uniques.fa -minsize 2 -otus $out/otus.fa -relabel Otu -log $out/cluster_otus.log
# followed by error correction
# -unoise $out/uniques.fa -fastaout $out/denoised.fa -relabel Den -log $out/unoise.log -minampsize 4


rule remove_chimeric_otus:
    input:
        fasta = rules.cluster_sequences.output,
        reference = rules.make_uchime_database.output
    output: "results/{eid}/{pid}/OTU.fasta"
    version: USEARCH_VERSION
    message: "Chimera filtering OTU seed sequences against %s" % config['chimera_database']['metadata']
    threads: 22
    log: "results/{eid}/{pid}/logs/uchime_ref.log"
    shell:
        '''usearch -uchime_ref {input.fasta} -db {input.reference} -nonchimeras	{output} \
            -strand plus -threads {threads} 2> {log}
        '''


rule utax:
    input:
        fasta = rules.remove_chimeric_otus.output,
        db = rules.make_tax_database.output
    output:
        fasta = temp("results/{eid}/{pid}/utax/OTU_tax_utax.fasta"),
        txt = temp("results/{eid}/{pid}/utax/OTU_tax_tax.txt")
    version: USEARCH_VERSION
    message: "Assigning taxonomies with UTAX algorithm using USEARCH with a confidence cutoff of {params.utax_cutoff}"
    params:
        utax_cutoff = config['taxonomy']['prediction_confidence_cutoff']
    threads: 22
    log: "results/{eid}/{pid}/logs/utax.log"
    shell:
        '''
        usearch -utax {input.fasta} -db {input.db} -strand both -threads {threads} \
            -fastaout {output.fasta} -utax_cutoff {params.utax_cutoff} \
            -utaxout {output.txt} 2> {log}
        '''


rule fix_utax_taxonomy:
    input:
        fasta = rules.utax.output.fasta,
        txt = rules.utax.output.txt
    output:
        fasta = "results/{eid}/{pid}/utax/OTU_tax.fasta",
        txt = "results/{eid}/{pid}/utax/OTU_tax.txt"
    params:
        kingdom = config['kingdom']
    message: "Altering taxa to reflect QIIME style annotation"
    run:
        with open(input.fasta) as ifh, open(output.fasta, 'w') as ofh:
            for line in ifh:
                line = line.strip()
                if not line.startswith(">OTU_"):
                    print(line, file=ofh)
                else:
                    print(fix_fasta_tax_entry(line, params.kingdom), file=ofh)
        with open(input.txt) as ifh, open(output.txt, 'w') as ofh:
            for line in ifh:
                toks = line.strip().split("\t")
                print(toks[0], fix_tax_entry(toks[1], params.kingdom),
                      fix_tax_entry(toks[2], params.kingdom), toks[3], sep="\t", file=ofh)


rule blast:
    input:
        fasta = rules.remove_chimeric_otus.output,
        db = rules.make_blast_db.output
    output:
        "results/{eid}/{pid}/blast/blast_hits.txt"
    params:
        P = config['taxonomy']['lca_cutoffs'],
        L = config['taxonomy']['prediction_confidence_cutoff'],
        db = config['blast_database']['fasta']
    threads: 22
    shell:
        '''
        blastn -query {input.fasta} -db {params.db} -num_alignments 200 \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -out {output} -num_threads {threads}
        '''


rule lca:
    input: rules.blast.output
    output: "results/{eid}/{pid}/blast/lca_assignments.txt"
    params:
        tax = config['blast_database']['taxonomy'],
        P = config['taxonomy']['lca_cutoffs'],
        L = config['taxonomy']['prediction_confidence_cutoff']
    shell: "resources/lca_src/lca -L {params.L} -P {params.P} -i {input} -r {params.tax} -o {output}"


rule assignments_from_lca:
    input:
        tsv = rules.lca.output,
        fasta = rules.remove_chimeric_otus.output
    output:
        "results/{eid}/{pid}/OTU_tax.fasta"
    run:
        lca_results = {}
        with open(input.tsv[0]) as fh:
            for line in fh:
                line = line.strip().split("\t")
                # file has a header, but doesn't matter
                tax = []
                for i, level in enumerate('kpcofgs'):
                    try:
                        tax.append("%s__%s" % (level, line[i + 1]))
                    except IndexError:
                        # ensure unassigned
                        for t in line[1:]:
                            assert t == "?"
                        tax.append("%s__?" % level)
                lca_results[line[0]] = tax
        with open(input.fasta[0]) as fasta_file, open(str(output), "w") as outfile:
            for line in fasta_file:
                line = line.strip()
                if not line.startswith(">"):
                    print(line, file=outfile)
                else:
                    try:
                        print(">%s;tax=%s" % (line[1:], ",".join(lca_results[line[1:]])), file=outfile)
                    except KeyError:
                        print(">%s;tax=k__?,p__?,c__?,o__?,f__?,g__?,s__?" % line[1:], file=outfile)


rule compile_counts:
    input:
        fastq = rules.combine_merged_reads.output,
        utax_db = rules.fix_utax_taxonomy.output.fasta,
        lca_db = rules.assignments_from_lca.output
    output:
        utax_txt = "results/{eid}/{pid}/utax/OTU.txt",
        lca_txt = "results/{eid}/{pid}/OTU.txt"
    params:
        threshold = config['mapping_to_otus']['read_identity_requirement']
    threads: 22
    shell:
        '''
        usearch -usearch_global {input.fastq} -db {input.utax_db} -strand plus \
            -id {params.threshold} -otutabout {output.utax_txt} \
            -threads {threads}
        usearch -usearch_global {input.fastq} -db {input.lca_db} -strand plus \
            -id {params.threshold} -otutabout {output.lca_txt} \
            -threads {threads}
        '''


rule biom:
    input:
        utax_txt = rules.compile_counts.output.utax_txt,
        lca_txt = rules.compile_counts.output.lca_txt
    output:
        utax_biom = "results/{eid}/{pid}/utax/OTU.biom",
        lca_biom = "results/{eid}/{pid}/OTU.biom"
    shadow: "shallow"
    shell:
        '''
        sed 's|\"||g' {input.utax_txt} | sed 's|\,|\;|g' > OTU_converted.txt
        biom convert -i OTU_converted.txt -o {output.utax_biom} --to-json --process-obs-metadata sc_separated --table-type "OTU table"
        sed 's|\"||g' {input.lca_txt} | sed 's|\,|\;|g' > OTU_converted.txt
        biom convert -i OTU_converted.txt -o {output.lca_biom} --to-json --process-obs-metadata sc_separated --table-type "OTU table"
        '''


rule multiple_align:
    input: rules.remove_chimeric_otus.output
    output: "results/{eid}/{pid}/OTU_aligned.fasta"
    message: "Multiple alignment of samples using Clustal Omega"
    version: CLUSTALO_VERSION
    threads: 1
    shell: "clustalo -i {input} -o {output} --outfmt=fasta --threads {threads} --force"


rule newick_tree:
    input: "results/{eid}/{pid}/OTU_aligned.fasta"
    output: "results/{eid}/{pid}/OTU.tree"
    message: "Building tree from aligned OTU sequences with FastTree2"
    log: "results/{eid}/{pid}/logs/fasttree.log"
    shell: "FastTree -nt -gamma -spr 4 -log {log} -quiet {input} > {output}"


rule report:
    input:
        file1 = "results/{eid}/{pid}/OTU.biom",
        file2 = "results/{eid}/{pid}/OTU.fasta",
        file3 = "results/{eid}/{pid}/OTU.tree",
        raw_counts = expand("results/{eid}/logs/{sample}_R1.fastq.count", eid=EID, sample=SAMPLES),
        filtered_counts = expand("results/{eid}/logs/{sample}_filtered_R1.fastq.count", eid=EID, sample=SAMPLES),
        merged_counts = expand("results/{eid}/logs/{sample}_merged.fastq.count", eid=EID, sample=SAMPLES),
        css = "resources/report.css"
    shadow: "shallow"
    params:
        kmer_len = config['filtering']['reference_kmer_match_length'],
        ham_dist = config['filtering']['allowable_kmer_mismatches'],
        min_read_len = config['filtering']['minimum_passing_read_length'],
        min_merge_len = config['merging']['minimum_merge_length'],
        max_ee = config['filtering']['maximum_expected_error'],
        tax_cutoff = config['taxonomy']['prediction_confidence_cutoff'],
        tax_levels = config['taxonomy']['lca_cutoffs'],
        min_seq_abundance = config['clustering']['minimum_sequence_abundance'],
        tax_metadata = config['blast_database']['metadata'],
        tax_citation = config['blast_database']['citation'],
        alt_tax_metadata = config['taxonomy_database']['metadata'],
        chimera_metadata = config['chimera_database']['metadata'],
        chimera_citation = config['chimera_database']['citation'],
        samples = SAMPLES
    output:
        html = "results/{eid}/{pid}/README.html"
    run:
        from biom import parse_table
        from biom.util import compute_counts_per_sample_stats
        from operator import itemgetter
        from numpy import std

        summary_csv = "stats.csv"
        sample_summary_csv = "samplesummary.csv"
        samples_csv = "samples.csv"
        biom_per_sample_counts = {}
        with open(input.file1) as fh, open(summary_csv, 'w') as sumout, open(samples_csv, 'w') as samout, open(sample_summary_csv, 'w') as samplesum:
            bt = parse_table(fh)
            stats = compute_counts_per_sample_stats(bt)
            biom_per_sample_counts = stats[4]
            sample_counts = list(stats[4].values())

            # summary
            print("Samples", len(bt.ids()), sep=",", file=sumout)
            print("OTUs", len(bt.ids(axis='observation')), sep=",", file=sumout)
            print("OTU Total Count", sum(sample_counts), sep=",", file=sumout)
            print("OTU Table Density", bt.get_table_density(), sep=",", file=sumout)

            # sample summary within OTU table
            print("Minimum Count", stats[0], sep=",", file=samplesum)
            print("Maximum Count", stats[1], sep=",", file=samplesum)
            print("Median", stats[2], sep=",", file=samplesum)
            print("Mean", stats[3], sep=",", file=samplesum)
            print("Standard Deviation", std(sample_counts), sep=",", file=samplesum)

            for k, v in sorted(stats[4].items(), key=itemgetter(1)):
                print(k, '%1.1f' % v, sep=",", file=samout)

        sample_counts = {}

        for sample in params.samples:
            sample_counts[sample] = {}
            # get raw count
            for f in input.raw_counts:
                if "%s_R1.fastq.count" % sample in f:
                    with open(f) as fh:
                        for line in fh:
                            sample_counts[sample]['raw_counts'] = int(line.strip())
                            break
            # filtered count
            for f in input.filtered_counts:
                if "%s_filtered_R1.fastq.count" % sample in f:
                    with open(f) as fh:
                        for line in fh:
                            sample_counts[sample]['filtered_counts'] = int(line.strip())
                            break
            # merged count
            for f in input.merged_counts:
                if "%s_merged.fastq.count" % sample in f:
                    with open(f) as fh:
                        for line in fh:
                            sample_counts[sample]['merged_counts'] = int(line.strip())
                            break

        raw_counts = []
        filtered_counts = []
        merged_counts = []
        biom_counts = []
        samps = []
        # sort this by the raw counts total and get the strings for the report
        for s in sorted(sample_counts.items(), key=lambda k_v: k_v[1]['raw_counts']):
            samps.append(s[0])
            raw_counts.append(s[1]['raw_counts'])
            filtered_counts.append(s[1]['filtered_counts'])
            merged_counts.append(s[1]['merged_counts'])
            # read count contribution to OTUs
            biom_counts.append(biom_per_sample_counts[s[0]])

        # quoted strings within brackets
        samples_str = "['%s']" % "', '".join(map(str, samps))
        # non-quoted ints or floats within brackets
        raw_counts_str = "[%s]" % ", ".join(map(str, raw_counts))
        filtered_counts_str = "[%s]" % ", ".join(map(str, filtered_counts))
        merged_counts_str = "[%s]" % ", ".join(map(str, merged_counts))
        biom_counts_str = "[%s]" % ", ".join(map(str, biom_counts))

        report("""
        =============================================================
        README - {wildcards.eid}
        =============================================================

        .. raw:: html

            <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
            <script src="https://code.highcharts.com/highcharts.js"></script>
            <script src="https://code.highcharts.com/modules/exporting.js"></script>
            <script type="text/javascript">
            $(function () {{
                $('#raw-count-plot').highcharts({{
                    chart: {{
                        type: 'column'
                    }},
                    title: {{
                        text: 'Sequence Counts'
                    }},
                    xAxis: {{
                        categories: {samples_str},
                        crosshair: true
                    }},
                    yAxis: {{
                        min: 0,
                        title: {{
                            text: 'Count'
                        }}
                    }},
                    tooltip: {{
                        headerFormat: '<span style="font-size:10px">{{point.key}}</span><table>',
                        pointFormat: '<tr><td style="color:{{series.color}};padding:0">{{series.name}}: </td>' +
                            '<td style="padding:0"><b>{{point.y:.1f}}</b></td></tr>',
                        footerFormat: '</table>',
                        shared: true,
                        useHTML: true
                    }},
                    credits: {{
                        enabled: false
                    }},
                    plotOptions: {{
                        column: {{
                            pointPadding: 0.2,
                            borderWidth: 0
                        }}
                    }},
                    series: [{{
                                name: 'Raw',
                                data: {raw_counts_str}
                            }},
                            {{
                                name: 'Filtered',
                                data: {filtered_counts_str},
                                visible: false
                            }},
                            {{
                                name: 'Merged',
                                data: {merged_counts_str},
                                visible: false
                            }},
                            {{
                                name: 'Assigned to OTUs',
                                data: {biom_counts_str},
                                visible: false
                            }}]
                    }});
            }});
            </script>

        .. contents:: Contents
            :backlinks: none

        Summary
        -------

        .. csv-table::
            :file: {summary_csv}

        .. raw:: html

            <div id="raw-count-plot" style="min-width: 310px; height: 500px; margin: 0 auto"></div>


        Output
        ------

        Chimera Removal
        ***************

        This occurs de novo during clustering and reference-based post clustering on the OTU seed
        sequences.

        Chimera database - {params.chimera_metadata}

        Biom Table
        **********

        Counts observed per sample as represented in the biom file (file1_). This count is
        representative of quality filtered reads that were assigned per sample to OTU seed
        sequences.

        .. csv-table::
            :file: {sample_summary_csv}

        Taxonomy was assigned to the OTU sequences at an overall cutoff of {params.tax_cutoff}%
        and per taxonomic level from Species to Kingdom at: {params.tax_levels}.

        Taxonomy database - {params.tax_metadata}

        OTU Sequences
        *************

        The OTU sequences are available in FASTA format (file2_) and aligned as newick tree
        (file3_).

        To build the tree, sequences were aligned using Clustalo [1] and FastTree2 [2] was used
        to generate the phylogenetic tree.


        Methods
        -------

        Raw sequence reads were demultiplexed with using EA-Utils [3] with zero
        mismatches allowed in the barcode sequence. Reads were quality filtered with BBDuk2 [4]
        to remove adapter sequences and PhiX with matching kmer length of {params.kmer_len}
        bp at a hamming distance of {params.ham_dist}. Reads shorter than {params.min_read_len} bp
        were discarded. Reads were merged using USEARCH [5] with a minimum length
        threshold of {params.min_merge_len} bp and maximum error rate of {params.max_ee}%. Sequences
        were dereplicated (minimum sequence abundance of {params.min_seq_abundance}) and clustered
        using the distance-based, greedy clustering method of USEARCH [6] at
        {wildcards.pid}% pairwise sequence identity among operational taxonomic unit (OTU) member
        sequences. De novo prediction of chimeric sequences was performed using USEARCH during
        clustering. Taxonomy was assigned to OTU sequences using BLAST [7] alignments followed by
        least common ancestor assignments across {params.tax_metadata} [8]. OTU seed sequences were
        filtered against {params.chimera_metadata} [9] to identify chimeric OTUs using USEARCH.


        References
        ----------

        1. Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, et al. 2011. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol 7: 539
        2. Price MN, Dehal PS, Arkin AP. 2010. FastTree 2--approximately maximum-likelihood trees for large alignments. ed. A.F.Y. Poon. PLoS One 5: e9490
        3. Erik Aronesty (2013). TOBioiJ : "Comparison of Sequencing Utility Programs", DOI:10.2174/1875036201307010001
        4. Bushnell, B. (2014). BBMap: A Fast, Accurate, Splice-Aware Aligner. URL https://sourceforge.net/projects/bbmap/
        5. Edgar, RC (2010). Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461. doi: 10.1093/bioinformatics/btq461
        6. Edgar, RC (2013). UPARSE: highly accurate OTU sequences from microbial amplicon reads. Nat Methods.
        7. Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., & Madden T.L. (2008) "BLAST+: architecture and applications." BMC Bioinformatics 10:421.
        8. {params.tax_citation}
        9. {params.chimera_citation}


        All Files
        ---------

        More files are available in relation to this analysis than are presented here. They can
        be accessed from the results directory and are organized by your experiment ID
        ({wildcards.eid})::

            {wildcards.eid}/
            ├── {wildcards.pid}                                 # clustering pairwise identity threshold
            │   ├── blast
            │   │   ├── blast_hits.txt             # raw blast hits per OTU seed seq
            │   │   └── lca_assignments.txt        # raw lca results TSV from blast hits
            │   ├── logs                           # output from these applications
            │   │   ├── fasttree.log
            │   │   ├── uchime_ref.log
            │   │   └── utax.log
            │   ├── OTU_aligned.fasta              # multiple alignment file of otu seed seqs
            │   ├── OTU.biom                       # tax annotated biom (no metadata, no normalization, SILVA v123)
            │   ├── OTU.fasta                      # otu seqs
            │   ├── OTU_tax.fasta                  # otu seqs with tax in FASTA header
            │   ├── OTU.tree                       # newick tree of multiple alignment
            │   ├── OTU.txt                        # tab delimited otu table
            │   ├── README.html                    # results report
            │   └── utax                           # utax results using {params.alt_tax_metadata}
            │       ├── OTU.biom
            │       ├── OTU_tax.fasta
            │       ├── OTU_tax.txt
            │       └── OTU.txt
            ├── demux
            │   └── *.fastq                        # all of the samples that were analyzed
            ├── logs
            │   ├── quality_filtering_stats.txt
            │   └── *.counts
            ├── merged_1.fasta                     # error corrected FASTA prior to clustering into OTU seqs
            └── merged.fastq                       # all sample reads merged into single file with updated headers

        """, output.html, metadata="Author: Joe Brown (joe.brown@pnnl.gov)",
        stylesheet=input.css, file1=input.file1, file2=input.file2, file3=input.file3)
