from setuptools import setup, find_packages

setup(
    name='DOLPHIN',
    version='0.1.0',
    
    packages=find_packages(),
    install_requires=[
    'matplotlib>=3.0.0,<4.0.0', 
    'optuna>=2.0.0,<3.0.0',
    'torch-geometric==2.2.0',
    'torch-sparse==0.6.16',
    'torch-scatter==2.1.0',
    'h5py>=3.9.0',
    'igraph>=0.10.7',
    'igraph<=0.11',
    'networkx>=3.1',
    'numba>=0.57.1',
    'pandas>=2.1.0',
    'parso>=0.8.3',
    'pillow>=9.4.0',
    'pip>=23.2.1',
    'psutil>=5.9.5',
    'ptyprocess>=0.7.0',
    'pure_eval>=0.2.2',
    'pycparser>=2.21',
    'pygments>=2.16.1',
    'pynndescent>=0.5.10',
    'pyopenssl>=23.2.0',
    'pyparsing>=3.1.1',
    'pysocks>=1.7.1',
    'sphinx>=5.0',
    'python-dateutil>=2.8.2',
    'python-igraph>=0.10.8',
    'pytz>=2023.3.post1',
    'pyzmq>=25.1.1',
    'lightning==2.0.8',
    'numpyro==0.13.0',
    'requests>=2.31.0',
    'anndata==0.8.0',
    'scanpy>=1.9.5',
    'jax==0.4.20',
    'jaxlib==0.4.20',
    'scikit-learn>=1.3.0',
    'scipy>=1.9.3',
    'seaborn>=0.12.2',
    'setuptools>=68.0.0',
    'six>=1.16.0',
    'numpy==1.24.4',
    'stack_data>=0.6.2',
    'statsmodels>=0.14.0',
    'stdlib-list>=0.8.0',
    'sympy>=1.11.1',
    'tbb>=2021.8.0',
    'texttable>=1.6.7',
    'threadpoolctl>=3.2.0',
    'torchaudio>=2.0.0',
    'torchvision>=0.15.0',
    'tornado>=6.3.3',
    'tqdm>=4.66.1',
    'traitlets>=5.9.0',
    'typing_extensions>=4.7.1',
    'tzdata>=2023c',
    'umap-learn>=0.5.3',
    'wcwidth>=0.2.6',
    'wheel>=0.38.4',
    'pyro-ppl>=1.8.6',
    'scvi-tools>=1.0.3',
    'leidenalg>=0.10.1'
    ],
    author='Kailu Song',
    author_email='kailu.song@mail.mcgill.ca',
    description='DOLPHIN: Advancing Single-cell RNA-seq Analysis Beyond Gene-Level by Integrating Exon-Level Quantification and Junction Reads with Deep Neural Networks',
    url='https://github.com/mcgilldinglab/DeepExonas.git',
)