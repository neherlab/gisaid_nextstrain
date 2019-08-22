from datetime import date
import pandas as pd
from treetime.utils import numeric_date

def reference_strain(v):
    references = {'h3n2':"A/Beijing/32/1992",
                  'h1n1pdm':"A/California/07/2009",
                  'vic':"B/HongKong/02/1993",
                  'yam':"B/Singapore/11/1994"
                  }
    return references[v.lineage]

genes_to_translate = {'ha':['SigPep', 'HA1', 'HA2'], 'na':['NA'],
                      'pb1':['PB1-F2'], 'pb2':['PB2'], 'pa':['PA'],
                      'np':['NP'], 'm':['M'], 'ns':['NEP']}

def gene_names(w):
    return genes_to_translate[w.segment]


def clock_rate(w):
    rate = {
        ('h3n2', 'ha'): 0.0043, ('h3n2', 'na'):0.0029,
        ('h1n1pdm', 'ha'): 0.0040, ('h1n1pdm', 'na'):0.0032,
        ('vic', 'ha'): 0.0024, ('vic', 'na'):0.0015,
        ('yam', 'ha'): 0.0019, ('yam', 'na'):0.0013
    }
    return rate.get((w.lineage, w.segment), 0.001)


def clock_std_dev(w):
    return 0.2*clock_rate(w)


# file names
outliers = "config/outliers_{lineage}.txt",
exclude_sites = "config/exclude-sites_{lineage}.txt",
references = "config/references_{lineage}.txt",
reference = "config/reference_{lineage}_{segment}.gb",
colors = "config/colors.tsv",
auspice_config = "config/auspice_config_{lineage}.json",
input_data = "data/{lineage}_{segment}.fasta"


rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = input_data
    output:
        sequences = "results/sequences_{lineage}_{segment}.fasta",
        metadata = "results/metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus accession date region country division location passage authors"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule filter:
    message:
        """
        Filtering {wildcards.lineage} {wildcards.segment} {wildcards.passage} sequences:
          - less than {params.min_length} bases
          - outliers
          - samples with missing region and country metadata
          - samples that are egg-passaged if cell build
        """
    input:
        metadata = rules.parse.output.metadata,
        sequences = rules.parse.output.sequences,
        exclude = outliers
    output:
        sequences = 'results/filtered_{lineage}_{segment}.fasta'
    params:
        min_length = 900
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --exclude {input.exclude} \
            --output {output}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = reference
    output:
        alignment = "results/aligned_{lineage}_{segment}.fasta"
    shell:
        """
        python3 scripts/codon_align.py \
            --sequences {input.sequences} \
            --reference {input.reference} \
            --output {output.alignment}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment,
        exclude_sites = exclude_sites
    output:
        tree = "results/tree-raw_{lineage}_{segment}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads 1 \
            --exclude-sites {input.exclude_sites}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree_{lineage}_{segment}.nwk",
        node_data = "results/branch-lengths_{lineage}_{segment}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = clock_rate,
        clock_std_dev = clock_std_dev
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --no-covariance \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt-muts_{lineage}_{segment}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa-muts_{lineage}_{segment}.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/traits_{lineage}_{segment}.json",
    params:
        columns = "region"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """


def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.traits.output.node_data,
    ]
    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs


rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = auspice_config,
        node_data = _get_node_data_for_export
    output:
        auspice_tree = "auspice/{lineage}_{segment}_tree.json",
        auspice_meta = "auspice/{lineage}_{segment}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta} \
            --minify-json
        """


rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice ",
        "logs"
    shell:
        "rm -rfv {params}"
