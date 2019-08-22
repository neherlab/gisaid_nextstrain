# get the influenza/gisaid repository
git clone https://github.com/neherlab/gisaid_nextstrain

cd gisaid_nextstrain
mkdir data && cp example_data/*fasta data
conda activate nextstrain
snakemake auspice/h1n1pdm_ha_tree.json
auspice view


