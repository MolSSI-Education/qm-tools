---
title: Setup
---
## Installing Python through Anaconda
[Python](https://python.org/) is a popular language for scientific computing, and great for general-purpose programming as well. Installing all of its scientific packages individually can be a bit difficult, however, so we recommend the all-in-one installer Anaconda.

1. Navigate to the [download page](https://www.anaconda.com/products/individual) for Anaconda.
2. Download the appropriate installer for your operating system. **Make sure you get the installer listed under Python 3** (not 2.7).
3. Double click the installer icon and follow the set-up instructions, keeping most of the default options. If you are Windows, make sure to choose to choose the option **Make Anaconda the default Python** during installation.

## Installing a Text Editor

You will need a text editor for this workshop. If you do not have a preferred text editor for writing code, we recommend [Atom](https://atom.io). Download Atom at the link and install on your computer.

## Installing OpenMM (computational chemistry software)
We recommend using a conda environment for software installations. If you know how to use a conda environment, you should create one for this tutorial and install OpenMM there. If you do not know how to create an environment, do not create one and just follow the instructions below.

After you have Anaconda installed, you can install the software using `conda`. Open a Terminal if you are on Mac or Linux. If you are on Windows, you should open a program called Anaconda Prompt. Type the following command into the terminal/anaconda prompt and press Enter. When prompted, answer 'yes' and wait for the install to finish:

## Windows command
~~~
conda install -c raimis -c psi4 -c conda-forge psi4=1.3.2
~~~
{: .language-bash}

## MacOS or Linux
~~~
conda install psi4 psi4-rt python=3.7 -c psi4
~~~
{: .language-bash}


## Start a Jupyter notebook
Open the Anaconda Navigator, and find the Jupyter notebook button. Click Launch.

It may take a few seconds to load the page, especially if it is the first time you have ever used the jupyter notebook, so don't panic if nothing loads for a few seconds.  Then a new window should open in your default internet browser. Use the file navigation window to navigate to the `qm-tools` folder.  In the upper right hand corner, click New, then choose Python 3 from the dropdown list.  You're ready to go!




{% include links.md %}
