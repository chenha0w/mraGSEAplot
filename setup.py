from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="mraGSEAplot",
    version="0.0.1",
    description='''
    Plot the GSEA-like plot for the results of master regulator analysis (MRA). 
    Flexibly sub-plot arrangement and color selection is allowed.
    Support both R-version and python-version of regulon.
    ''',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Chenhao Wu",
    author_email="tierwuchen@gmail.com"
    packages=find_packages(include=["mraGSEAplot"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=required,
    license="MIT",
    url="https://github.com/chenha0w/mraGSEAplot",
)