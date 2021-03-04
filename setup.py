from setuptools import setup, find_packages

## to create conda environment call
#conda create -n boost-hic python=3.8
#conda activate boost-hic
#pip install -e .

setup(name='Boost-HiC',
      version='0.0.1',
      description='Boost-HiC: Computational enhancement of long-range contacts in chromosomal contact maps, '
                  'https://doi.org/10.1093/bioinformatics/bty1059',
      url='https://github.com/CellFateNucOrg/Boost-HiC',
      packages=find_packages(),
      python_requires='>=3.8',
      install_requires=['numpy', 'scipy', 'h5py', 'scikit-learn', 'scikit-image'],
      zip_safe=False)
