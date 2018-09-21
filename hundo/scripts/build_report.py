#!/usr/bin/env python
# coding=utf-8

import argparse
import base64
import csv
import datetime
import io
import mimetypes
import os
import textwrap
from collections import defaultdict
from glob import glob

import numpy as np
import pandas as pd
import plotly.figure_factory as ff
import plotly.graph_objs as go
from biom import parse_table
from docutils.core import publish_file, publish_parts
from docutils.parsers.rst import directives
from plotly import offline

import relatively

STYLE = """
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap.min.css"/>
    <style type="text/css">
    body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#summary::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:22%;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:20%;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}h1{line-height:1.6}.simple{padding-left:20px}.docutils.container{width:100%}
    </style>
    """
SCRIPT = """
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/dataTables.bootstrap.min.js"></script>
    <script>
    $(document).ready( function () {
        $('#summaryTable').DataTable();
        $('#otuTable').DataTable( {
            "searching": false,
            "paging": false,
            "ordering": false,
            "info": false,
        } );
    } );
    </script>
    """
BLUE = "rgba(31,119,180,0.3)"
RED = "rgba(214,39,40,0.3)"


def build_uri(file, file_type=None, default_encoding="utf8"):
    fname = os.path.basename(file)
    type, encoding = mimetypes.guess_type(file)
    if file_type is None:
        file_type = type
    if encoding is not None:
        default_encoding = encoding
    with open(file, mode="rb") as fh:
        fdata = base64.b64encode(fh.read())
    return f'data:{file_type};charset={default_encoding};filename={fname};base64,{fdata.decode("utf-8")}'


def parse_biom(biom):
    with open(biom) as fh:
        bt = parse_table(fh)
        df = bt.to_dataframe()
        otu_df = [
            ["Samples", "OTUs", "OTU Total Count", "OTU Table Density"],
            [
                bt.length("sample"),
                bt.length("observation"),
                bt.sum(),
                "{:f}".format(bt.get_table_density()),
            ],
        ]
        otu_df = pd.DataFrame(otu_df)
        # make the first row the column labels
        otu_df.rename(columns=otu_df.iloc[0], inplace=True)
        # drop the first row
        otu_df = otu_df.reindex(otu_df.index.drop(0))
        # build the table
        summary_table = otu_df.to_html(
            index=False,
            bold_rows=False,
            classes=["table", "table-bordered"],
            table_id="otuTable",
        )
        # fix for rst
        summary_table = summary_table.replace("\n", "\n" + 10 * " ")
    return df, summary_table


def make_div(figure_or_data, include_plotlyjs=False, show_link=False):
    div = offline.plot(
        figure_or_data,
        include_plotlyjs=include_plotlyjs,
        show_link=show_link,
        output_type="div",
    )
    if ".then(function ()" in div:
        div = f"""{div.partition(".then(function ()")[0]}</script>"""
    return div


def build_summary_table(raw, biom_df, div_df, omitted=None):
    # sample summary table
    header = ["Sample", "Raw", "Filtered", "Merged", "In OTUs"]
    count_data = []
    omit = []
    if omitted:
        for sample_count in omitted.split(","):
            sample, count = sample_count.split(":")
            omit.append(sample)
            count_data.append([sample, count])
    for count_file in raw:
        sample_name = count_file.partition("logs/raw_counts/")[-1].partition(
            "_R1.count"
        )[0]
        sd = [sample_name]
        with open(count_file) as fh:
            for line in fh:
                sd.append(int(float(line.strip())))
        with open(count_file.replace("raw_counts/", "filtered_counts/")) as fh:
            for line in fh:
                sd.append(int(float(line.strip())))
        with open(count_file.replace("raw_counts/", "merged_counts/")) as fh:
            for line in fh:
                sd.append(int(float(line.strip())))
        sd.append(int(biom_df[sample_name].sum()))
        count_data.append(sd)
    df = pd.DataFrame(count_data, columns=header)
    df = df.merge(div_df, how="outer", left_on="Sample", right_index=True)
    df.sort_values(by="shannon", inplace=True)
    # summary table
    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["table", "table-bordered"],
        table_id="summaryTable",
    )
    df_str = df_str.replace("\n", "\n" + 10 * " ")
    df = df.reset_index()
    # pull out the omitted samples
    keep = [idx for idx, sample in df.Sample.items() if sample not in omit]
    df = df.iloc[keep]
    df.sort_values(by="Raw", inplace=True)
    plot_data = []
    visible = [True] * 4
    if len(df.Sample) > 10:
        visible = [True] + ["legendonly"] * 3
    for i, metric in enumerate(header[1:]):
        plot_data.append(
            go.Bar(x=df.Sample, y=df[metric], name=metric, visible=visible[i])
        )
    layout = dict(
        title="Counts Per Sample By Stage",
        barmode="group",
        height=500,
        yaxis={"title": "Count"},
        showlegend=True,
        hovermode="x",
        bargroupgap=0.1,
        legend={"orientation": "h", "x": 0, "y": 1.06},
    )
    fig = go.Figure(data=plot_data, layout=layout)
    return df_str, make_div(fig)


def build_quality_plot(r1_quals):
    # quality plot
    raw_qual_stats = defaultdict(lambda: defaultdict(list))
    for r1_ee_file in r1_quals:
        sample_name = r1_ee_file.partition("logs/")[-1].partition("_R1")[0]
        r2_ee_file = "_R2".join(r1_ee_file.rsplit("_R1", 1))
        with open(r1_ee_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                raw_qual_stats["R1"][sample_name].append(float(row["Mean_Q"]))
        with open(r2_ee_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                q = float(row["Mean_Q"])
                raw_qual_stats["R2"][sample_name].append(q)
    data = []
    for read_index, sample_data in raw_qual_stats.items():
        color = BLUE if read_index == "R1" else RED
        for sample_name, sample_stats in sample_data.items():
            data.append(
                go.Scatter(
                    x=list(range(1, len(sample_stats))),
                    y=sample_stats,
                    name=sample_name,
                    text=sample_name,
                    hoverinfo="text+x+y",
                    legendgroup=read_index,
                    mode="lines",
                    line=dict(color=color, dash="solid"),
                )
            )
    layout = go.Layout(
        title="Mean Quality Scores for R1 and R2",
        xaxis={"title": "Position"},
        yaxis={"title": "Quality (Phred score)"},
        hovermode="closest",
        showlegend=False,
        autosize=True,
        annotations=[
            dict(
                x=0,
                y=1.1,
                xref="paper",
                yref="paper",
                text="Forward",
                showarrow=False,
                font=dict(size=16, color="#ffffff"),
                align="left",
                borderpad=4,
                bgcolor="#1f77b4",
            ),
            dict(
                x=0.15,
                y=1.1,
                xref="paper",
                yref="paper",
                text="Reverse",
                showarrow=False,
                font=dict(size=16, color="#ffffff"),
                align="left",
                borderpad=4,
                bgcolor="#d62728",
            ),
        ],
    )
    fig = go.Figure(data=data, layout=layout)
    return make_div(fig)


def build_taxonomy_plot(txt, value_cols, height=900):
    df = pd.read_table(txt)
    levels = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    df[levels] = df["taxonomy"].str.split(",", expand=True)
    # remove "k__", etc., from taxonomies
    for level in levels:
        df[level] = df[level].str.replace(f"{level[0]}__", "")
    hierarchy = ["phylum", "class", "order"]
    df[hierarchy] = df[hierarchy].fillna("NA")
    df[value_cols] = df[value_cols].fillna(0)
    df = df[hierarchy + value_cols]
    dfs = relatively.get_dfs_across_hierarchy(
        df, hierarchy, value_cols, reorder=None, dependent="; "
    )
    fig = relatively.get_abundance_figure_from_dfs(
        dfs, hierarchy, "Assigned Taxonomy Per Sample", height=height
    )
    return make_div(fig)


def build_ee_plot(eeplot):
    # merged reads quality plot
    ee_stats = defaultdict(list)
    for stats_file in eeplot:
        sample_name = stats_file.partition("logs/")[-1].partition("_merged")[0]
        with open(stats_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                ee_stats[sample_name].append(row["Mean_EE"])
    data = []
    for sample_name, sample_stats in ee_stats.items():
        data.append(
            go.Scatter(
                x=list(range(1, len(sample_stats))),
                y=sample_stats,
                name=sample_name,
                text=sample_name,
                hoverinfo="text+x+y",
                mode="lines",
                line=dict(color=BLUE, dash="solid"),
            )
        )
    layout = go.Layout(
        title="Mean Expected Error Rate Across Merged Reads",
        yaxis={"title": "Expected Error"},
        xaxis={"title": "Base Position"},
        hovermode="closest",
        showlegend=False,
    )
    fig = go.Figure(data=data, layout=layout)
    return make_div(fig)


def format_from_file(txt):
    ret = ""
    with open(txt) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                ret += line
            else:
                ret += " " * 4 + line
    return ret


def get_metadata(author=None, version=None, sep=" // "):
    md = [datetime.date.today().isoformat()]
    if version:
        md.insert(0, version)
    if author:
        md.insert(0, author)
    return sep.join(md)


def build_attachment(zip):
    data = build_uri(zip)
    name = os.path.basename(zip)
    return f"""<a href="{data}" download="{name}" draggable="true">{name}</a>"""


def main(
    biom,
    txt,
    zip,
    env,
    params,
    eeplot,
    r1quals,
    raw,
    omitted,
    html,
    author,
    version,
    samples,
):
    # resolve the file patterns
    eeplot = glob(eeplot, recursive=False)
    r1quals = glob(r1quals, recursive=False)
    raw = glob(raw, recursive=False)

    biom_df, otu_summary_table = parse_biom(biom)
    div_idx = relatively.calculate_diversity(
        biom_df,
        biom_df.index,
        biom_df.columns,
        index=["shannon", "simpson", "invsimpson"],
    )
    sample_summary_table, count_summary_plot = build_summary_table(
        raw, biom_df, div_idx, omitted
    )
    sample_order = div_idx.sort_values(by="shannon").index.tolist()
    quality_plot = build_quality_plot(r1quals)
    taxonomy_plot = build_taxonomy_plot(txt, sample_order, 800)
    ee_plot = build_ee_plot(eeplot)
    params_str = format_from_file(params)
    conda_str = format_from_file(env)
    metadata = get_metadata(author, version)
    attachment = build_attachment(zip)
    samples_attachment = build_attachment(samples)
    report_str = f"""

.. raw:: html

    {STYLE}
    {SCRIPT}

=============================================================
hundo_ - Experiment Summary
=============================================================

.. contents::
    :backlinks: none
    :depth: 2

Summary
-------

Sequence Counts
***************

.. raw:: html

    <div style="overflow-x:auto;">
    {sample_summary_table}
    </div>
    {count_summary_plot}

Sequence Quality
****************

.. raw:: html

    {quality_plot}
    {ee_plot}

OTUs
****

.. raw:: html

    <div style="overflow-x:auto;">
    {otu_summary_table}
    </div>

Taxonomy
********

Samples are ordered from least diverse to most.

.. raw:: html

    {taxonomy_plot}

Methods
-------

``hundo`` wraps a number of functions and its exact steps are
dependent upon a given user's configuration.

Paired-end sequence reads are quality trimmed and can be trimmed of
adapters and filtered of contaminant sequences using BBDuk2 of the
BBTools_ package. Passing reads are merged using
VSEARCH_ then aggregated into a single FASTA file with headers
describing the origin and count of the sequence.

Prior to creating clusters, reads are filtered again based on expected
error rates:

To create clusters, the aggregated, merged reads are dereplicated to
remove singletons by default using VSEARCH. Sequences are preclustered
into centroids using VSEARCH to accelerate chimera filtering. Chimera
filtering is completed in two steps: *de novo* and then reference
based. The reference by default is the entire annotation database.
Following chimera filtering, sequences are placed into clusters using
distance-based, greedy cluster with VSEARCH based on the allowable
percent difference of the configuration.

After OTU sequences have been determined, BLAST_ or VSEARCH is used to align
sequences to the reference database. Reference databases for 16S were curated
by the CREST_ team and hundo incorporates the CREST LCA method. ITS databases
are maintained by UNITE_.

Counts are assigned to OTUs using the global alignment method of
VSEARCH, which outputs the final OTU table as a tab-delimited text
file. The Biom_ command line tool is used to convert the tab table
to biom.

Multiple alignment of sequences is completed using MAFFT_.
A tree based on the aligned sequences is built using FastTree2_.

This workflow is built using Snakemake_ and makes use of Bioconda_
to install its dependencies.

.. _BBTools: https://sourceforge.net/projects/bbmap/
.. _VSEARCH: https://github.com/torognes/vsearch
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi
.. _CREST: https://github.com/lanzen/CREST
.. _UNITE: https://unite.ut.ee/
.. _Biom: http://biom-format.org/
.. _MAFFT: https://mafft.cbrc.jp/alignment/software/
.. _FastTree2: http://www.microbesonline.org/fasttree/
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/
.. _Bioconda: https://bioconda.github.io/

Configuration
-------------

::

    {params_str}


Execution Environment
---------------------

::

    {conda_str}


Output Files
------------

Not all files written by the workflow are contained within the
Downloads section of this page to minimize the size of the this
document. Other output files are described and are written to the
results directory.


Attached
********

The zip archive contains the following files:

OTU.biom
````````

Biom table with raw counts per sample and their associated
taxonomic assignment formatted to be compatible with downstream tools
like phyloseq_.

OTU.fasta
`````````

Representative DNA sequences of each OTU.

OTU.tree
````````

Newick tree representation of aligned OTU sequences.

OTU.txt
```````

Tab-delimited text table with columns OTU ID, a column for each sample,
and taxonomy assignment in the final column as a comma delimited list.


Other Result Files
******************

Other files that may be needed, but that are not included in the
attached archive include:

OTU_aligned.fasta
`````````````````

OTU sequences after alignment using MAFFT.

all-sequences.fasta
```````````````````

Quality-controlled, dereplicated DNA sequences of all samples. The
header of each record identifies the sample of origin and the count
resulting from dereplication.

blast-hits.txt
``````````````

The BLAST assignments per OTU sequence.

Downloads
---------

.. _hundo: https://github.com/pnnl/hundo
.. _phyloseq: https://joey711.github.io/phyloseq/


.. container::
   :name: attachments


        .. container::
            :name: archive

            archive:
                .. raw:: html

                    {attachment}

        .. container::
            :name: samples

            samples:
                .. raw:: html

                    {samples_attachment}


.. container::
   :name: metadata

   {metadata}

"""
    with open(html, "w") as fh:
        publish_file(
            source=io.StringIO(report_str),
            destination=fh,
            writer_name="html",
            settings_overrides={"stylesheet_path": ""},
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--biom")
    p.add_argument("--txt")
    p.add_argument("--zip")
    p.add_argument("--env")
    p.add_argument("--params")
    p.add_argument("--eeplot")
    p.add_argument("--r1quals")
    p.add_argument("--raw")
    p.add_argument("--omitted")
    p.add_argument("--html")
    p.add_argument("--author")
    p.add_argument("--version")
    p.add_argument("--samples")
    args = p.parse_args()
    main(
        args.biom,
        args.txt,
        args.zip,
        args.env,
        args.params,
        args.eeplot,
        args.r1quals,
        args.raw,
        args.omitted,
        args.html,
        args.author,
        args.version,
        args.samples,
    )
