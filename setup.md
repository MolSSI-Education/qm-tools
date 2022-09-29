---
title: Setup
---

You can either install software on your computer, or you can use a website called [ChemCompute](https://chemcompute.org/) which will allow you to access computing resources and a pre-configured environment for this workshop without installing any software on your own computer. 
ChemCompute is only accessible to users with an academic email address. **Warning** If you are an academic user outside of the United States, there is a chance your email won't be properly recognized. If this happens, you may have to install software locally.

## Using a cloud computing environment

1. [Register](https://chemcompute.org/register/) for a ChemCompute account. The easiest thing to do is to choose "Sign in with your University Login" on the right. You can login using your university email. You will need to confirm account registration from an email chemcompute will send you.
2. After confirming your account, you will be able to access a Jupyter notebook using the menu on the top. Select "Jupyter", then click "Use Psi4 or JupyterHub".
3. A Jupyter notebook environment will start for you. We recommend creating a new folder for and creating notebooks in that folder for this workshop.

## Using your personal computer

### Installing Python through Anaconda
[Python](https://python.org/) is a popular language for scientific computing, and great for general-purpose programming as well. We recommend using environments and managing packages using the `conda` package manager. If you have [Anaconda]() installed, you have the `conda` package manager and can use it through the command line.

If you don't have Anaconda installed already, we recommend installing Python and conda using `miniconda`.

1. Navigate to the [download page for miniconda](https://docs.conda.io/en/latest/miniconda.html). Choose an installer appropriate for your operating system.
2. 

## Installing a Text Editor

You will need a text editor for this workshop. If you do not have a preferred text editor for writing code, we recommend [Visual Studio Code](https://code.visualstudio.com/). Download Visual Studio Code at the link and install on your computer.

## Installing Python Packages
We recommend using a conda environment for software installations. 
We provide a  <a href="https://education.molssi.org/qm-tools/files/conda-env.yaml" download>conda environment yaml</a> file that you can use and download if you're familiar with creating environments from a file.
If you use the environment yaml, we recommend the [mamba package manager](https://mamba.readthedocs.io/en/latest/) to create the environment.

Otherwise, type the following commands into your terminal
~~~ 
conda env create -n qm-tools python="3.10"
conda activate qm-tools
conda install -c psi4 psi4
conda install -c conda-forge notebook matplotlib numpy nglview
~~~
{: .bash}

You now have an environment called `qm-tools` with Psi4, Jupyter Notebook, Matplotlib, NumPy, and NGLView installed. 
Any time you want to use Psi4, you should do

~~~
$ conda activate qm-tools
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
