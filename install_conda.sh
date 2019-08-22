# This script installs augur and auspice and can be run without
# administrator privileges. It has the following prerequisits
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

source ~/.bashrc

conda activate nextstrain
npm install --global auspice



