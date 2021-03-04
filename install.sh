#!/bin/bash

conda create -n boost-hic python=3.8
conda activate boost-hic
pip install -e .
