from setuptools import setup, find_packages

setup(
    name="beam_visualization",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "ete3",
        "pandas",
        "dendropy",
        "numpy",
        "matplotlib",
        "seaborn",
    ],
    extras_require={
        "test": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
        ],
    },
    author="Stephen Staklinski",
    author_email="staklins@cshl.edu",
    description="A package for visualizing and analyzing BEAM phylogenetic tree outputs",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/StephenStaklinski/beam_visualization",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
) 