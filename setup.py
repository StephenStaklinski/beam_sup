from setuptools import setup, find_packages
import sys

# Check Python version
if sys.version_info >= (3, 13):
    raise RuntimeError(
        "beam_visualization requires Python < 3.13 due to ete3 dependency. "
        "Please use Python 3.6-3.12 instead."
    )

setup(
    packages=find_packages(),
    python_requires=">=3.6,<3.13",
) 