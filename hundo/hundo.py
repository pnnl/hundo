import click
import logging
import multiprocessing
import os
import subprocess
import yaml
from collections import OrderedDict
from hundo import __version__
from hundo.classifier import Classifier


logging.basicConfig(level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s %(levelname)s] %(message)s")


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """hundo"""


@cli.command("make-config", short_help="prepopulate a configuration file with defaults")
@click.argument("config")
@click.option("--database-dir", default="databases", show_default=True, help="location to save reference databases")
@click.option("--threads", default=multiprocessing.cpu_count(), type=int, show_default=True, help="number of threads to use per multi-threaded job")
def make_config(config, database_dir, threads):
    represent_dict_order = lambda self, data: self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)

    conf = OrderedDict()

    conf["author"] = subprocess.check_output(["whoami"]).decode("utf-8").strip()
    conf["threads"] = threads
    conf["threads_small"] = max(int(conf["threads"] / 4), 1)
    conf["database_dir"] = database_dir
    conf["filter_adapters"] = ""
    conf["filter_contaminants"] = ""
    conf["allowable_kmer_mismatches"] = 1
    conf["reference_kmer_match_length"] = 31
    conf["reduced_kmer_min"] = 8
    conf["minimum_passing_read_length"] = 100
    conf["minimum_base_quality"] = 10
    conf["minimum_merge_length"] = 150
    conf["perform_error_correction"] = True
    conf["maximum_expected_error"] = 1
    conf["denovo_chimera_filter"] = True
    conf["reference_chimera_filter"] = True
    conf["minimum_sequence_abundance"] = 2
    conf["percent_of_allowable_difference"] = 3
    conf["reference_database"] = "silva"
    conf["blast_minimum_bitscore"] = 100
    conf["blast_top_fraction"] = 0.95
    conf["read_identity_requirement"] = 0.97
    with open(config, "w") as f:
        print(yaml.dump(conf, default_flow_style=False), file=f)
    logging.info("Configuration file written to %s" % config)
    logging.info("For parameter definitions, please see our documentation at ")


@cli.command("lca", short_help="runs LCA across BLAST hits")
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("blasthits", type=click.Path(exists=True))
@click.argument("mapfile", type=click.Path(exists=True))
@click.argument("trefile", type=click.Path(exists=True))
@click.argument("outfasta", type=click.Path())
@click.argument("outtab", type=click.Path())
@click.option("--min-score", type=int, default=100, show_default=True, help="minimum allowable bitscore")
@click.option("--top-fraction", type=float, default=0.95, show_default=True, help="calculate LCA based on HSPS within this fraction of highest scoring HSP")
def run_lca(fasta, blasthits, mapfile, trefile, outfasta, outtab, min_score=100, top_fraction=0.95):
    """Classifies BLAST HSPs using associated newick tree with corresponding
    names and map.
    """
    lca_classifier = Classifier(fasta)
    lca_classifier.read_tree(mapfile, trefile)
    lca_classifier.parse_blast(blasthits, top_fraction=top_fraction, min_bitscore=min_score)
    lca_classifier.prune_unassigned()
    lca_classifier.print_sequences(outfasta)
    lca_classifier.print_table(outtab)


def get_snakefile():
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Snakefile")
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


@cli.command("annotate", context_settings=dict(ignore_unknown_options=True), short_help="run annotation protocol")
@click.argument("config", type=click.Path(exists=True))
@click.option("-j", "--jobs", default=multiprocessing.cpu_count(), type=int, show_default=True, help="use at most this many cores in parallel; total running tasks at any given time will be jobs/threads")
@click.option("-o", "--out-dir", default=os.path.realpath("."), show_default=True, help="results output directory")
@click.option("--dryrun", is_flag=True, default=False, show_default=True, help="do not execute anything, just show the commands that would be executed")
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_annotate(config, jobs, out_dir, dryrun, snakemake_args):
    cmd = ("snakemake --snakefile {snakefile} --directory {out_dir} --printshellcmds "
           "--jobs {jobs} --rerun-incomplete --configfile '{config}' --nolock "
           "--config {args} --{dryrun}").format(snakefile=get_snakefile(),
                                                                    out_dir=out_dir,
                                                                    jobs=jobs,
                                                                    config=config,
                                                                    dryrun="dryrun" if dryrun else "",
                                                                    args=" ".join(snakemake_args))
    logging.info("Executing: " + cmd)
    check_call(cmd, shell=True)