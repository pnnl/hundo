import logging
import multiprocessing
import os
import subprocess
import sys
from collections import OrderedDict

import click

from hundo import __version__

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M",
    format="[%(asctime)s %(levelname)s] %(message)s",
)


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """Documentation is available at:

        \b
        https://hundo.rtfd.io

    Issues can be submitted to:

        \b
        https://github.com/pnnl/hundo/issues
    """


@cli.command("lca", short_help="runs LCA across aligned hits")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("aligned_hits", type=click.File("r"))
@click.argument("mapfile", type=click.File("r"))
@click.argument("trefile", type=click.File("r"))
@click.argument("outfasta", type=click.File("w"))
@click.argument("outtab", type=click.File("w"))
@click.option(
    "--min-score",
    type=int,
    default=125,
    show_default=True,
    help="minimum allowable bitscore",
)
@click.option(
    "--min-pid", type=float, default=0.85, show_default=True, help="min percent id"
)
@click.option(
    "--top-fraction",
    type=float,
    default=0.99,
    show_default=True,
    help="calculate LCA based on HSPS within this fraction of highest scoring HSP. Vsearch accuracy decreases with a lower number than this",
)
@click.option(
    "--aligner",
    type=str,
    default="blast",
    show_default=True,
    help="pick which aligner that you would like to use, blast or vsearch. default:blast",
)
def run_lca(
    fasta,
    aligned_hits,
    mapfile,
    trefile,
    outfasta,
    outtab,
    min_score,
    min_pid,
    top_fraction,
    aligner,
):
    """Classifies BLAST or VSEARCH HSPs using associated newick tree with corresponding
    names and map.
    """
    import statistics
    import hundo.unite_classifier as unite_lca
    import hundo.crest_classifier as crest_lca
    from hundo.blast import parse_blasthits
    from hundo.fasta import read_fasta, format_fasta_record

    def print_unite(name, seq, tx, fa, tsv):
        full_name = "{name};tax={taxonomy}".format(
            name=name.strip(";"), taxonomy=",".join(tx)
        )
        print(format_fasta_record(full_name, seq), file=fa)
        print(name, ";".join(tx), sep="\t", file=tsv)

    def unite_namemap(mapfile):
        m = dict()
        for line in mapfile:
            toks = line.strip().split("\t")
            m[toks[0]] = toks[1]
        return m

    protocol = "16S" if not "unite" in os.path.basename(mapfile.name) else "ITS"
    logging.info("Parsing aligned hits")
    if aligner == "vsearch":
        hsps = parse_vsearchhits(aligned_hits, min_pid, top_fraction)
    else:
        hsps = parse_blasthits(aligned_hits, min_score, top_fraction)
    unknown_taxonomy = ["%s__?" % i for i in "kpcofgs"]
    if protocol == "ITS":
        tree = unite_lca.Tree(trefile)
        namemap = unite_namemap(mapfile)
        logging.info("Performing LCA per OTU")
        with open(fasta) as fh:
            for name, seq in read_fasta(fh):
                if name in hsps:
                    hits = hsps[name]
                    # need to translate fasta names of blast hits to
                    # species names of the taxonomy map
                    translated_hits = [namemap[i] for i in hits.names]
                    # only slightly more stringent
                    taxonomy = tree.lca(
                        translated_hits, statistics.mean(hits.percent_ids)
                    )
                    print_unite(name, seq, taxonomy, outfasta, outtab)
                else:
                    # unknown
                    print_unite(name, seq, unknown_taxonomy, outfasta, outtab)
    # silva/gg
    else:
        otus = OrderedDict()
        with open(fasta) as fh:
            for name, seq in read_fasta(fh):
                otus[name] = crest_lca.OTU(name, seq)
        tree = crest_lca.Tree(mapfile, trefile)
        # update OTUs with a classification
        for otu_name, hits in hsps.items():
            otu = otus[otu_name]
            lca_node = tree.get_common_ancestor(hits.names)
            if not lca_node:
                logging.debug("No LCA -- assigning to No Hits: %s" % hits.names)
                otu.classification = tree.no_hits
                continue

            while (
                lca_node.name in tree.assignment_min
                and hits.percent_ids[-1] < tree.assignment_min[lca_node.name]
                and lca_node is not tree.root
            ):
                lca_node = tree.get_parent(lca_node)
            otu.classification = lca_node
        for otu_id, otu in otus.items():
            taxonomy = tree.get_taxonomy(otu.classification)
            full_name = ("{name};" "tax={taxonomy}").format(
                name=otu.name.strip(";"),
                taxonomy=",".join(
                    ["%s__%s" % (abb, tax) for abb, tax in taxonomy.items()]
                ),
            )
            print(format_fasta_record(full_name, otu.sequence), file=outfasta)
            print(
                otu_id,
                ";".join(["%s__%s" % (abb, tax) for abb, tax in taxonomy.items()]),
                sep="\t",
                file=outtab,
            )


def get_snakefile():
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Snakefile")
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


@cli.command(
    "download",
    context_settings=dict(ignore_unknown_options=True),
    short_help="download reference data (Optional)",
)
@click.option(
    "-d",
    "--database-dir",
    default="references",
    show_default=True,
    help="directory containing reference data or new directory into which to download reference data",
)
@click.option(
    "-j",
    "--jobs",
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help="use at most this many cores in parallel",
)
@click.option(
    "--rd",
    "--reference-database",
    default="silva",
    type=click.Choice(["silva", "greengenes", "unite"]),
    show_default=True,
    help="two 16S databases are supported along with Unite for ITS; only the reference specified will be downloaded",
)
@click.option(
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="do not execute anything, just show the commands that would be executed",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_download(database_dir, jobs, reference_database, dryrun, snakemake_args):
    """
    Download the reference databases, but do not execute any of the 16S or ITS
    annotation protocol. Running this prior to `hundo annotate` is not
    required, but may be useful in instances where compute nodes do not have
    internet access.
    """
    database_dir = os.path.realpath(database_dir)
    targets = []
    if reference_database == "greengenes":
        targets = " ".join(
            [
                os.path.join(database_dir, i)
                for i in [
                    "greengenes.fasta.nhr",
                    "greengenes.fasta.nin",
                    "greengenes.fasta.nsq",
                    "greengenes.map",
                    "greengenes.tre",
                ]
            ]
        )
    elif reference_database == "silva":
        targets = " ".join(
            [
                os.path.join(database_dir, i)
                for i in [
                    "silvamod128.fasta.nhr",
                    "silvamod128.fasta.nin",
                    "silvamod128.fasta.nsq",
                    "silvamod128.map",
                    "silvamod128.tre",
                ]
            ]
        )
    else:
        targets = " ".join(
            [
                os.path.join(database_dir, i)
                for i in [
                    "unite.fasta.nhr",
                    "unite.fasta.nin",
                    "unite.fasta.nsq",
                    "unite.map",
                    "unite.tre",
                ]
            ]
        )
    cmd = (
        "snakemake --snakefile {snakefile} --printshellcmds "
        "--jobs {jobs} --rerun-incomplete "
        "--nolock {dryrun} "
        "--config database_dir={database_dir} workflow=download "
        "reference_database={reference_database} "
        "{add_args} {args} {targets}"
    ).format(
        snakefile=get_snakefile(),
        jobs=jobs,
        dryrun="--dryrun" if dryrun else "",
        database_dir=database_dir,
        reference_database=reference_database,
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
        targets=targets,
    )
    logging.info("Executing: " + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)


@cli.command(
    "annotate",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run annotation protocol",
)
@click.argument("fastq-dir", type=click.Path(exists=True))
@click.option(
    "-i",
    "--input-dir",
    multiple=True,
    help=(
        "add directories in which to search for sample input file pairs "
        "in addition to FASTQ_DIR; may be specified multiple times"
    ),
)
@click.option(
    "--prefilter-file-size",
    default=100000,
    type=int,
    show_default=True,
    help="any file smaller than this in bytes is omitted from being processed",
)
@click.option(
    "-j",
    "--jobs",
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help="use at most this many cores in parallel; total running tasks at any given time will be jobs/threads",
)
@click.option(
    "-o",
    "--out-dir",
    default=os.path.realpath("."),
    show_default=True,
    help="results output directory",
)
@click.option(
    "--no-conda",
    is_flag=True,
    default=False,
    show_default=True,
    help="do not use conda environments",
)
@click.option(
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="do not execute anything, just show the commands that would be executed",
)
@click.option(
    "-a",
    "--author",
    default=subprocess.check_output(["whoami"]).decode("utf-8").strip(),
    show_default=True,
    help="will show in footer of summary HTML document",
)
@click.option(
    "--aligner",
    type=click.Choice(choices=["blast", "vsearch"]),
    default="blast",
    show_default=True,
    help="local aligner; `blast` is more sensitive while `vsearch` is much faster",
)
@click.option(
    "-t",
    "--threads",
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help="when a step is multi-threaded, use this many threads; this all or a subset of --jobs",
)
@click.option(
    "-d",
    "--database-dir",
    default="references",
    show_default=True,
    help="directory containing reference data or new directory into which to download reference data",
)
@click.option(
    "-fa",
    "--filter-adapters",
    default="",
    show_default=True,
    help="file path to adapters FASTA to use for trimming read ends",
)
@click.option(
    "-fc",
    "--filter-contaminants",
    default="",
    show_default=True,
    help="file path to FASTA to use for filtering reads",
)
@click.option(
    "-ak",
    "--allowable-kmer-mismatches",
    type=int,
    default=1,
    show_default=True,
    help="kmer mismatches allowed during adapter trim process",
)
@click.option(
    "-kl",
    "--reference-kmer-match-length",
    type=int,
    default=27,
    show_default=True,
    help="length of kmer to search against contaminant sequences",
)
@click.option(
    "-km",
    "--reduced-kmer-min",
    type=int,
    default=8,
    show_default=True,
    help="look for shorter kmers at read tips down to this length; 0 disables",
)
@click.option(
    "-rl",
    "--minimum-passing-read-length",
    type=int,
    default=100,
    show_default=True,
    help="passing single-end read length prior to merging",
)
@click.option(
    "-bq",
    "--minimum-base-quality",
    type=int,
    default=10,
    show_default=True,
    help="regions with average quality below this will be trimmed",
)
@click.option(
    "-ml",
    "--minimum-merge-length",
    type=int,
    default=150,
    show_default=True,
    help="minimum allowable read length after merging",
)
@click.option(
    "-am",
    "--allow-merge-stagger",
    is_flag=True,
    default=False,
    show_default=True,
    help="allow merging of staggered reads",
)
@click.option(
    "-md",
    "--max-diffs",
    type=int,
    default=5,
    show_default=True,
    help="maximum number of different bases in overlap",
)
@click.option(
    "-mo",
    "--min-overlap",
    type=int,
    default=16,
    show_default=True,
    help="minimum length of overlap between reads",
)
@click.option(
    "-ee",
    "--maximum-expected-error",
    type=float,
    default=1,
    show_default=True,
    help="after merging; the allowable limit of erroneous bases; accepts fractions as well",
)
@click.option(
    "-cf",
    "--reference-chimera-filter",
    default=True,
    show_default=True,
    help="define a file path or set to true to use BLAST reference database",
)
@click.option(
    "-sa",
    "--minimum-sequence-abundance",
    default=2,
    type=int,
    show_default=True,
    help="when clustering, don't create any clusters with fewer than this many representative sequences",
)
@click.option(
    "-pd",
    "--percent-of-allowable-difference",
    default=3,
    type=float,
    show_default=True,
    help="maximum difference between an OTU member sequence and the representative sequence of that OTU",
)
@click.option(
    "-rd",
    "--reference-database",
    default="silva",
    type=click.Choice(["silva", "greengenes", "unite"]),
    show_default=True,
    help="two 16S databases are supported along with Unite for ITS; references will be downloaded as needed",
)
@click.option(
    "-mb",
    "--blast-minimum-bitscore",
    default=125,
    type=int,
    show_default=True,
    help="filter out alignments below this bitscore threshold",
)
@click.option(
    "-tf",
    "--blast-top-fraction",
    default=0.95,
    type=float,
    show_default=True,
    help="when calculating LCA, only use this fraction of HSPs from the best scoring alignment",
)
@click.option(
    "-ir",
    "--read-identity-requirement",
    default=0.97,
    type=float,
    show_default=True,
    help="reflects the difference between OTU clusters to reduce ambiguous assignment",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_annotate(
    fastq_dir,
    input_dir,
    prefilter_file_size,
    jobs,
    out_dir,
    no_conda,
    dryrun,
    author,
    aligner,
    threads,
    database_dir,
    filter_adapters,
    filter_contaminants,
    allowable_kmer_mismatches,
    reference_kmer_match_length,
    reduced_kmer_min,
    minimum_passing_read_length,
    minimum_base_quality,
    minimum_merge_length,
    allow_merge_stagger,
    max_diffs,
    min_overlap,
    maximum_expected_error,
    reference_chimera_filter,
    minimum_sequence_abundance,
    percent_of_allowable_difference,
    reference_database,
    blast_minimum_bitscore,
    blast_top_fraction,
    read_identity_requirement,
    snakemake_args,
):
    """
    Run annotation across paired-end sequence data contained in FASTQ_DIR.
    Both R1 and R2 are expected to be present in the same directory and have
    the same name except for the index ID (R1 and R2).

    FASTQ_DIR may be a comma separated list of directories or additional
    input directories may be added using --input-dir multiple times.

    By using SILVA, you agree to their license terms which are available at:

        \b
        https://www.arb-silva.de/silva-license-information

    For complete documentation and parameter definitions, please see:

        \b
        https://hundo.rtfd.io
    """
    fq_dir = list()
    # combine single or comma separated list with input_dir
    all_input_paths = fastq_dir.replace(" ", "").split(",") + list(input_dir)
    for input_path in all_input_paths:
        fq_dir.append(os.path.realpath(input_path))
    # format the input paths in order to send to Snakemake command
    fq_dir = ",".join(fq_dir)
    database_dir = os.path.realpath(database_dir)
    filter_adapters = os.path.realpath(filter_adapters) if filter_adapters else ""
    filter_contaminants = (
        os.path.realpath(filter_contaminants) if filter_contaminants else ""
    )
    no_temp_declared = False
    for sa in snakemake_args:
        if sa == "--nt" or sa == "--notemp":
            no_temp_declared = True
    cmd = (
        "snakemake --snakefile {snakefile} --directory {out_dir} "
        "--printshellcmds --jobs {jobs} --rerun-incomplete "
        "--nolock {conda} {dryrun} "
        "--config fastq_dir={fq_dir} author='{author}' threads={threads} "
        "database_dir={database_dir} filter_adapters={filter_adapters} "
        "filter_contaminants={filter_contaminants} "
        "allowable_kmer_mismatches={allowable_kmer_mismatches} "
        "reference_kmer_match_length={reference_kmer_match_length} "
        "reduced_kmer_min={reduced_kmer_min} "
        "minimum_passing_read_length={minimum_passing_read_length} "
        "minimum_base_quality={minimum_base_quality} "
        "minimum_merge_length={minimum_merge_length} "
        "fastq_allowmergestagger={allow_merge_stagger} "
        "fastq_maxdiffs={max_diffs} "
        "fastq_minovlen={min_overlap} "
        "maximum_expected_error={maximum_expected_error} "
        "reference_chimera_filter={reference_chimera_filter} "
        "minimum_sequence_abundance={minimum_sequence_abundance} "
        "percent_of_allowable_difference={percent_of_allowable_difference} "
        "reference_database={reference_database} "
        "blast_minimum_bitscore={blast_minimum_bitscore} "
        "blast_top_fraction={blast_top_fraction} "
        "read_identity_requirement={read_identity_requirement} "
        "prefilter_file_size={prefilter_file_size} "
        "aligner={aligner} "
        "no_temp_declared={no_temp_declared} {add_args} "
        "{args}"
    ).format(
        snakefile=get_snakefile(),
        out_dir=os.path.realpath(out_dir),
        jobs=jobs,
        conda="" if no_conda else "--use-conda",
        dryrun="--dryrun" if dryrun else "",
        fq_dir=os.path.realpath(fastq_dir),
        author=author,
        threads=threads,
        database_dir=database_dir,
        filter_adapters=filter_adapters,
        filter_contaminants=filter_contaminants,
        allowable_kmer_mismatches=allowable_kmer_mismatches,
        reference_kmer_match_length=reference_kmer_match_length,
        reduced_kmer_min=reduced_kmer_min,
        minimum_passing_read_length=minimum_passing_read_length,
        minimum_base_quality=minimum_base_quality,
        minimum_merge_length=minimum_merge_length,
        allow_merge_stagger=allow_merge_stagger,
        max_diffs=max_diffs,
        min_overlap=min_overlap,
        maximum_expected_error=maximum_expected_error,
        reference_chimera_filter=reference_chimera_filter,
        minimum_sequence_abundance=minimum_sequence_abundance,
        percent_of_allowable_difference=percent_of_allowable_difference,
        reference_database=reference_database,
        blast_minimum_bitscore=blast_minimum_bitscore,
        blast_top_fraction=blast_top_fraction,
        read_identity_requirement=read_identity_requirement,
        prefilter_file_size=prefilter_file_size,
        aligner=aligner,
        no_temp_declared=no_temp_declared,
        add_args="" if snakemake_args and snakemake_args[0].startswith("-") else "--",
        args=" ".join(snakemake_args),
    )
    logging.info("Executing: " + cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
