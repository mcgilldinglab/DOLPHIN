from setuptools import find_packages, setup


setup(
    name="DOLPHIN",
    version="0.1.0",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "DOLPHIN.EDEG": ["*.R", "README.md"],
        "DOLPHIN.alternative_splicing": ["README.md"],
        "DOLPHIN.alternative_splicing.runtime_support": ["README.md"],
        "DOLPHIN.alternative_splicing.runtime_support.pandas_compat": ["*.py"],
    },
    install_requires=[],
    entry_points={
        "console_scripts": [
            "dolphin-prepare-reference=DOLPHIN.preprocess.cli:main",
            "dolphin-preprocess=DOLPHIN.preprocess.cli:main",
            "dolphin-edeg=DOLPHIN.EDEG.cli:main",
            "dolphin-as=DOLPHIN.alternative_splicing.pipeline:main",
        ],
    },
)
