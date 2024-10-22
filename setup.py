from setuptools import setup, find_packages

setup(
    name='DeepExonas',
    version='0.1.0',
    
    packages=find_packages(),
    install_requires=[
    'pytorch-lightning==1.7.7',
    'torch==1.12.1+cu113',
    'torch-geometric==2.2.0',
    'torch-scatter==2.1.0+pt112cu113',
    'torch-sparse==0.6.16+pt112cu113',
    'torchaudio==0.12.1+cu113',
    'torchmetrics==0.10.3',
    'torchvision==0.13.1+cu113',
    'pyro-api==0.1.2',
    'pyro-ppl==1.8.3',
    'numpyro==0.10.1',
    'numpy==1.23.5',
    'optuna==3.0.5',
    'PyYAML==6.0',
    'matplotlib==3.6.2',
    'scanpy==1.9.1'
    ],
    author='Kailu Song',
    author_email='kailu.song@mail.mcgill.ca',
    description='DeepExonas: Advancing Cell Representation Beyond Gene-Level by Integrating Exon-Level Quantification and Junction Reads with Deep Neural Networks',
    url='https://github.com/mcgilldinglab/DeepExonas.git',
)