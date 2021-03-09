#!/bin/bash

conda create -n boost-hic python=3.8
conda activate boost-hic
pip install -e .

BOOST_HIC_HOME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo "You can activate conda environment like this:"
echo "conda activate boost-hic"
echo "You can also add the following line to your .bashrc, or .zshrc, or .bash_profile:"
echo "export PATH=\$PATH:$BOOST_HIC_HOME"
echo "You can check your active shell like this:"
echo "echo \$SHELL"
echo "You can also change and update your start shell script boost-hic.sh command."
