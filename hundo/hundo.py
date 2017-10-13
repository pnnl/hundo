import click
import logging
import multiprocessing
import os
import subprocess
import sys
from collections import OrderedDict
from hundo import __version__
from hundo.crest_classifier import run_crest_classifier

logging.basicConfig(level=logging.INFO,
                    datefmt="%Y-%m-%d %H:%M",
                    format="[%(asctime)s %(levelname)s] %(message)s")


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """hundo"""


@cli.command("lca", short_help="runs LCA across BLAST hits")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("blasthits", type=click.Path(exists=True))
@click.argument("mapfile", type=click.Path(exists=True))
@click.argument("trefile", type=click.Path(exists=True))
@click.argument("outfasta", type=click.Path())
@click.argument("outtab", type=click.Path())
@click.option("--min-score",
              type=int,
              default=100,
              show_default=True,
              help="minimum allowable bitscore")
@click.option(
    "--top-fraction",
    type=float,
    default=0.98,
    show_default=True,
    help=
    "calculate LCA based on HSPS within this fraction of highest scoring HSP")
@click.option("--filter-euks",
              is_flag=True,
              default=False,
              show_default=True,
              help="filter eukaryotic assigned OTUs")
def run_lca(fasta,
            blasthits,
            mapfile,
            trefile,
            outfasta,
            outtab,
            min_score=100,
            top_fraction=0.98,
            filter_euks=False):
    """Classifies BLAST HSPs using associated newick tree with corresponding
    names and map.
    """
    run_crest_classifier(fasta, blasthits, mapfile, trefile, outfasta, outtab,
                         min_score, top_fraction, filter_euks)


def get_snakefile():
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Snakefile")
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


@cli.command("annotate",
             context_settings=dict(ignore_unknown_options=True),
             short_help="run annotation protocol")
@click.argument("fastq-dir", type=click.Path(exists=True))
@click.option("--prefilter-file-size",
    default=100000,
    type=int,
    show_default=True,
    help="any file smaller than this in bytes is omitted from being processed")
@click.option(
    "-j",
    "--jobs",
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help=
    "use at most this many cores in parallel; total running tasks at any given time will be jobs/threads"
)
@click.option("-o",
              "--out-dir",
              default=os.path.realpath("."),
              show_default=True,
              help="results output directory")
@click.option("--no-conda",
              is_flag=True,
              default=False,
              show_default=True,
              help="do not use conda environments")
@click.option(
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help=
    "do not execute anything, just show the commands that would be executed")
@click.option(
    "-a",
    "--author",
    default=subprocess.check_output(["whoami"]).decode("utf-8").strip(),
    show_default=True,
    help="will show in footer of summary HTML document")
@click.option(
    "-t",
    "--threads",
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help=
    "when a step is multi-threaded, use this many threads; this all or a subset of --jobs"
)
@click.option(
    "-d",
    "--database-dir",
    default="references",
    show_default=True,
    help=
    "directory containing reference data or new directory into which to download reference data"
)
@click.option("-fa",
              "--filter-adapters",
              default="",
              show_default=True,
              help="file path to adapters FASTA to use for trimming read ends")
@click.option("-fl",
              "--filter-contaminants",
              default="",
              show_default=True,
              help="file path to FASTA to use for filtering reads")
@click.option("-ms",
              "--allowable-kmer-mismatches",
              type=int,
              default=1,
              show_default=True,
              help="kmer mismatches allowed during adapter trim process")
@click.option("-kl",
              "--reference-kmer-match-length",
              type=int,
              default=27,
              show_default=True,
              help="length of kmer to search against contaminant sequences")
@click.option(
    "-km",
    "--reduced-kmer-min",
    type=int,
    default=8,
    show_default=True,
    help="look for shorter kmers at read tips down to this length; 0 disables")
@click.option("-rl",
              "--minimum-passing-read-length",
              type=int,
              default=100,
              show_default=True,
              help="passing single-end read length prior to merging")
@click.option("-bq",
              "--minimum-base-quality",
              type=int,
              default=10,
              show_default=True,
              help="regions with average quality below this will be trimmed")
@click.option("-ml",
              "--minimum-merge-length",
              type=int,
              default=150,
              show_default=True,
              help="minimum allowable read length after merging")
@click.option(
    "-ee",
    "--maximum-expected-error",
    type=float,
    default=1,
    show_default=True,
    help=
    "after merging; the allowable limit of erroneous bases; accepts fractions as well"
)
@click.option(
    "-cf",
    "--reference-chimera-filter",
    default=True,
    show_default=True,
    help="define a file path or set to true to use BLAST reference database")
@click.option(
    "-sa",
    "--minimum-sequence-abundance",
    default=2,
    type=int,
    show_default=True,
    help=
    "when clustering, don't create any clusters with fewer than this many representative sequences"
)
@click.option(
    "-pd",
    "--percent-of-allowable-difference",
    default=3,
    type=float,
    show_default=True,
    help=
    "maximum difference between an OTU member sequence and the representative sequence of that OTU"
)
@click.option(
    "-rd",
    "--reference-database",
    default="silva",
    type=click.Choice(["silva", "greengenes", "unite"]),
    show_default=True,
    help=
    "two 16S databases are supported along with Unite for ITS; references will be downloaded as needed"
)
@click.option("-mb",
              "--blast-minimum-bitscore",
              default=100,
              type=int,
              show_default=True,
              help="filter out alignments below this bitscore threshold")
@click.option(
    "-tf",
    "--blast-top-fraction",
    default=0.95,
    type=float,
    show_default=True,
    help=
    "when calculating LCA, only use this fraction of HSPs from the best scoring alignment"
)
@click.option(
    "-ir",
    "--read-identity-requirement",
    default=0.97,
    type=float,
    show_default=True,
    help=
    "reflects the difference between OTU clusters to reduce ambiguous assignment"
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_annotate(
        fastq_dir, prefilter_file_size, jobs, out_dir, no_conda, dryrun,
        author, threads, database_dir, filter_adapters, filter_contaminants,
        allowable_kmer_mismatches, reference_kmer_match_length,
        reduced_kmer_min, minimum_passing_read_length, minimum_base_quality,
        minimum_merge_length, maximum_expected_error, reference_chimera_filter,
        minimum_sequence_abundance, percent_of_allowable_difference,
        reference_database, blast_minimum_bitscore, blast_top_fraction,
        read_identity_requirement, snakemake_args):
    """
    For complete documentation and parameter definitions, please see:
    https://hundo.rtfd.io
    """
    database_dir = os.path.realpath(database_dir)
    filter_adapters = os.path.realpath(
        filter_adapters) if filter_adapters else ""
    filter_contaminants = os.path.realpath(
        filter_contaminants) if filter_contaminants else ""

    cmd = ("snakemake --snakefile {snakefile} --directory {out_dir} "
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
           "maximum_expected_error={maximum_expected_error} "
           "reference_chimera_filter={reference_chimera_filter} "
           "minimum_sequence_abundance={minimum_sequence_abundance} "
           "percent_of_allowable_difference={percent_of_allowable_difference} "
           "reference_database={reference_database} "
           "blast_minimum_bitscore={blast_minimum_bitscore} "
           "blast_top_fraction={blast_top_fraction} "
           "read_identity_requirement={read_identity_requirement} "
           "prefilter_file_size={prefilter_file_size} {add_args} "
           "{args}").format(
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
               maximum_expected_error=maximum_expected_error,
               reference_chimera_filter=reference_chimera_filter,
               minimum_sequence_abundance=minimum_sequence_abundance,
               percent_of_allowable_difference=percent_of_allowable_difference,
               reference_database=reference_database,
               blast_minimum_bitscore=blast_minimum_bitscore,
               blast_top_fraction=blast_top_fraction,
               read_identity_requirement=read_identity_requirement,
               prefilter_file_size=prefilter_file_size,
               add_args="" if snakemake_args and snakemake_args[
                   0].startswith("-") else "--",
               args=" ".join(snakemake_args))
    logging.info("Executing: " + cmd)
    subprocess.check_call(cmd, shell=True)
