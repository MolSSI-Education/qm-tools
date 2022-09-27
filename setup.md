---
title: Setup
---

## Using your personal computer

### Installing Python through Anaconda
[Python](https://python.org/) is a popular language for scientific computing, and great for general-purpose programming as well. Installing all of its scientific packages individually can be a bit difficult, however, so we recommend the all-in-one installer Anaconda.

1. Navigate to the [download page](https://www.anaconda.com/products/distribution) for Anaconda.
2. Double click the installer icon and follow the set-up instructions, keeping most of the default options. If you are Windows, make sure to choose to choose the option **Make Anaconda the default Python** during installation.

## Installing a Text Editor

You will need a text editor for this workshop. If you do not have a preferred text editor for writing code, we recommend [Visual Studio Code](https://code.visualstudio.com/). Download Visual Studio Code at the link and install on your computer.

## Installing Python Packages
We recommend using a conda environment for software installations. 

1.  <a href="https://education.molssi.org/qm-tools/files/conda-env.yaml" download>Download the conda environment file.</a> 

2. After you have Anaconda installed and have downloaded the file, you can install the software for this workshop the environment yaml. We recommend using the `mamba` package manager (installable with `conda`) to create this environment. 
The environment can be created with conda, but will take several minutes to install. Using mamba is much faster.

~~~
conda install mamba
mamba env create -f conda-env.yaml
~~~
{: .bash}

## Start a Jupyter notebook
Open the Anaconda Navigator, and choose your `qm-tooos` environment. After switching to your environment, find the Jupyter notebook button on your Home tab. Click Launch.

You can also activate your environment and start a Jupyter notebook from the the command line:

~~~
conda activate qm-tools
jupyter-notebook
~~~
{: .bash}

It may take a few seconds to load the page, especially if it is the first time you have ever used the jupyter notebook, so don't panic if nothing loads for a few seconds.  Then a new window should open in your default internet browser. Use the file navigation window to navigate to the `mm-tools` folder.  In the upper right hand corner, click New, then choose Python 3 from the dropdown list.  You're ready to go!


{% include links.md %}
