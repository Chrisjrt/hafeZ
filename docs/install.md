# Installation

## Conda

The easiest way to install `hafeZ` is via conda. For inexperienced command line users, this method is highly recommended.

```
conda install -c bioconda hafeZ
```

This will install all the dependencies along with `hafeZ`. The dependencies are listed in environment.yml.

If conda is taking a long time to solve the environment, try using mamba:

```
conda install mamba
mamba install -c bioconda pharokka
```

## Pip

As of v2.0.0, you can also install the python components of `hafeZ` with pip.

```
pip install hafeZ
```

You will still need to install the non-python dependencies manually.

## Source

Alternatively, the development version of `hafeZ` (which may include new, untested features) can be installed manually via github. 

```
git clone "https://github.com/Chrisjrt/hafeZ.git"
cd hafeZ
pip install -e .
hafeZ -h
```

The dependencies found in environment.yml will then need to be installed manually.

For example using conda to install the required dependencies:

```
conda env create -f environment.yml
conda activate hafeZ
# assuming you are in the hafeZ directory 
# installs hafeZ from source
pip install -e .
hafeZ -h
```

# Database Installation

To install the hafez database to the default directory:

`hafeZ database`

If you would like to specify a different database directory (recommended), that can be achieved as follows:

`hafeZ database -d <path/to/databse_dir>`

If this does not work, you an alternatively download the databases from Zenodo at https://zenodo.org/record/8402631/files/hafeZ_v2.0.0_databases.tar.gz  and untar the directory in a location of your choice.

If you prefer to use the command line:

```
wget "https://zenodo.org/record/8402631/files/hafeZ_v2.0.0_databases.tar.gz"
tar -xzf hafeZ_v2.0.0_databases.tar.gz
```

which will create a directory called "hafeZ_v2.0.0_databases" containing the databases.

# Beginner Conda and Mamba Installation

If you are new to using the command-line, please install mamba or conda using either option detailed below.

I would prefer and recommend the miniforge option as it comes with mamba!

##  Install miniforge

1. Install [miniforge](https://github.com/conda-forge/miniforge). This will automatically install conda.

Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate [here](https://github.com/conda-forge/miniforge)

`curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh`

Install miniforge and follow the prompts.

`sh Miniforge3-Linux-x86_64.sh`


2. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Finally, I would recommend installing hafeZ into a fresh environment. For example to create an environment called hafeZ_env with hafeZ installed:

```
mamba create -n hafeZ_env hafeZ
conda activate hafeZ_env
hafeZ -h
```

## Install Miniconda

1. Install [Anaconda](https://www.anaconda.com/products/distribution). I would recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniconda](https://docs.conda.io/en/latest/miniconda.html) website).

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

Install miniconda and follow the prompts.

`sh Miniconda3-latest-Linux-x86_64.sh`


2. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. After this, conda should be installed (you may need to restart your terminal). It is recommended that mamba is also installed, as it will solve the enviroment quicker than conda:

```
conda install mamba
```

4. Finally, I would recommend installing hafeZ into a fresh environment. For example to create an environment called hafeZ_env with hafeZ installed:

```
mamba create -n hafeZ_env hafeZ
conda activate hafeZ_env
hafeZ -h
```

