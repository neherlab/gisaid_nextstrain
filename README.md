# Nextstrain for GISAID downloads

This is a simple [nextstrain](https://nextstrain.org) pipeline for seasonal influenza virus data sets download from [GISAID](https://gisaid.org).
It is derived from the main nextstrain seasonal flu analysis workflow available at [nextstrain/seasonal-flu](https://github.com/nextstrain/seasonal-flu).

To install nextstrain, you need a linux-like environment (for example Windows Subsystem for Linux) and execute the script [`install_conda.sh`](install_conda.sh):
```
# This script installs miniconda and augur. It can be run without
# administrator privileges. It has the following prerequisites
#  - curl
#  - gcc
# You need to install these via your package manager (e.g. apt-get)

# get miniconda
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o conda_installer.sh

# install miniconda -- accepting defaults should be fine
bash conda_installer.sh

## install nextstrain
curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml
conda env create -f nextstrain.yml
```
These commands install miniconda and create an environment in which you can run nextstrain.
Running these commands will prompt you to accept terms of conditions and specify the location of the installation -- the default answers are fine here.
The conda installer will also ask you whether you want to set-up the conda installation -- choose `yes`.

Next, we need to install the visualization component.
Open a new terminal and type
```
conda activate nextstrain
npm install --global auspice
```
These commands are also provided as a script [`install_auspice.sh`](install_auspice.sh).

Once this is done, you can download this repository and run a test analysis using the following commands (provided as [`test_installation.sh`](test_installation.sh)).
```
# get the influenza/gisaid repository
git clone https://github.com/neherlab/gisaid_nextstrain

cd gisaid_nextstrain
mkdir data && cp example_data/*fasta data
conda activate nextstrain
snakemake auspice/h3n2_ha_tree.json
auspice view
```
This will run for a few minutes.
Once it completed, you should be able view an analysis in your browser at [`http://localhost:4000/h3n2_ha`](http://localhost:4000/h3n2_ha)

[Nextstrain]: https://nextstrain.org
[fauna]: https://github.com/nextstrain/fauna
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice
[snakemake cli]: https://snakemake.readthedocs.io/en/stable/executable.html#all-options
[nextstrain-cli]: https://github.com/nextstrain/cli
[nextstrain-cli README]: https://github.com/nextstrain/cli/blob/master/README.md
