import os
from setuptools import setup, find_packages

# Read the contents of your README file
PACKAGE_DIR = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(PACKAGE_DIR, 'README.md'), encoding='utf-8') as f:
    LONG_DESCRIPTION = f.read()

setup(name='pmd',
      version='0.1.0',
      author='Kuan-Hsuan Shen',
      author_email='kshen64@gatech.edu',
      description='Automated generation of LAMMPS data and input files for polymer molecular dynamics simulations',
      long_description=LONG_DESCRIPTION,
      long_description_content_type='text/markdown',
      keywords=['LAMMPS', 'polymer', 'SMILES'],
      url='https://github.com/Ramprasad-Group/High-throughput-MD-simulations',
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
      packages=find_packages(),
      zip_safe=False
      )