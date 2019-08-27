# This script installs miniconda and augur. It can be run without
# administrator privileges. It has the following prerequisites
#  - curl
#  - gcc
# You need to install these via your package manager (e.g. apt-get)
# In case of macOS, you will need to install the developer tools

# get miniconda for linux (use for WSL and Linux, for MacOS, comment out this line and use the line below)
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o conda_installer.sh

## use the following line for macOS
#curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o conda_installer.sh


# install miniconda -- accepting defaults should be fine
bash conda_installer.sh

## install nextstrain
curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml
conda env create -f nextstrain.yml
