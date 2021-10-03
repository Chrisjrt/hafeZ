from setuptools import setup
from os import path
import setuptools

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), "r") as fh:
    long_description = fh.read()

setup(
    name='hafeZ',
    description='A tool for identifying active prophage elements through read mapping',
    long_description=long_description,
    long_description_content_type="text/markdown",
    version="1.0.2",
    license='GPL-3',
    author='Christopher J. R. Turkington',
    author_email='chrisjrt1@gmail.com',
    url='https://github.com/Chrisjrt/hafeZ',
     packages=[
        'hZ'
    ],
    install_requires=[
        'pyrodigal>=0.4.7',
        'biopython>=1.78',
        'matplotlib>=3.3.4',
        'numpy>=1.20.1',
        'pandas>=1.2.4',
        'pysam>=0.16.0.1',
        'scipy>=1.6.2',
        'seaborn>=0.11.1'
    ],
)
