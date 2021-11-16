.. _install:

Installation
============

``Nemaktis`` is a mixture of C++ and python codes, and has been successfully tested both on
Windows 10 and Linux. The recommended way of installing and using Nematkis is through the
package manager ``conda`` (see next subsection).

If you are an experienced Mac user and want to become maintainer of a MacOS version of
Nematkis (using conda), please contact me
(https://nemaktis.readthedocs.io/en/latest/intro/overview.html#contact). Otherwise, it is
probably possible to compile and use the software by hand on other OSs than the ones
mentioned above (windows 7, macOS...), provided you know what you are doing :)

With Conda (windows or linux)
-----------------------------

a. Install Anaconda
...................

The simplest way of installing the ``Nemaktis`` package is through the package manager
``conda`` (no compilation and no dependency to be installed by hand). I provide precompiled
packages for both **Linux** and **Windows 10**, which means the following method will work
on these two operating systems.

The first step is to install the software ``miniconda``, which contains a barebone version of
the package manager ``conda``. If you already have the full python distribution ``Anaconda``
installed on your machine, this step is not necessary. The installation files for Windows/Linux
are available at this address: https://docs.conda.io/en/latest/miniconda.html

b1. (Windows) Install Nemaktis automatically
............................................

If you are a Windows 10 user and do not want to copy-paste commands in a terminal, the next
step is as simple as running the following installation script 

https://github.com/warthan07/Nemaktis/releases/download/v1.4.3/Install_Nemaktis-1.4.3.cmd

This script will create a special environment for ``Nemaktis`` named *nm* and will install
everything needed in it. It will also install the python editor ``Spyder`` and create a
shortcut named *Spyder (Nemakis environment)* for it on your Desktop (this is necessary even
if you already installed Spyder, since it has to be run from inside the conda environment
*nm*).

b2. (Windows/Linux) Install Nemaktis on the command line
........................................................

Alternatively, if you are a Linux user or want to type the installation commands yourself
(they are not very complicated after all!), open a terminal (Windows: application "Conda
terminal" installed with miniconda, Linux: any terminal) and type the following command: ::

  conda create -n nm -c conda-forge -c warthan07 -y nemaktis=1.4.3

(Optional) If you want to use your favourite python editor when using ``Nemaktis``, you have
to install and run it from the same conda environment. You can search https://anaconda.org/
to find the associated package and installation command. For example, to install ``Spyder``
you just need to type: ::

  conda activate nm
  conda install spyder

Note that when you want to run python scripts using nemaktis, the installed python editor
should always be run from inside the *nm* environment. For example, to run ``Spyder``, you
should type in the terminal: ::

  conda activate nm
  spyder

The advantage of step b1 is that it creates a shortcut for Spyder which automatically does
this activation step for you. 


c. (Windows/Linux) How to update
................................

I do not recommend updating nematkis using the command ``conda update``, since I do not
compile Nematkis sufficiently often for a correct update of all dependencies. In other
words, running ``conda update`` has the risk of breaking your conda environment! I
apologize for this current limitation, which mostly stems from my inexperience at conda
packaging with complex dependencies. 

In this case, how can you safely update Nemaktis? The simplest way is to fully remove the
conda environment *nm*, either by removing the folder "envs/nm" inside the root folder of
Anaconda or by opening a terminal ("Conda terminal" for windows, any terminal for linux) and
typing: ::

  conda remove --name nm --all 

Then, simply repeat the installation step b1 or b2 above, eventually adjusting the version
number of Nemaktis (it can be obtained from https://anaconda.org/warthan07/nemaktis/files)
inside the commands or script if I forgot to do it :)

Another possibility is to repeat the installation step b1/b2 with a different environment
name than *nm* (for the windows script method, you need to manually edit the script file),
for example *nm[version number]*. Although this method is probably fine to test new features
of the software without removing the old version , it is probably not very good in the long
term since each conda environment takes a non-negligible portion of disk space. 

Developper method (only linux)
------------------------------

If you want to be able to modify the code, but still want to enjoy the simplicity of conda
packages (no relative paths to manage, everything works as with a system package), you can build
yourselves the nemaktis package for Linux:

1. Get the code of the ``Nemaktis`` repository: ::

     git clone git@github.com:warthan07/Nemaktis.git

   And implements the changes that you want. For the C++ codes, compilation steps are provided
   in the subfolders bpm-solver and rt-solver if you want to test them locally (in which case
   you will have to install yourselves the dependencies listed in each CMakeLists).

2. In a terminal, go to the subfolder **conda_recipe** of the ``Nemaktis`` repository and activate
   any conda environment in which you have write access. If you don't have any conda environment
   yet, you can type: ::
     
     conda create -n build
     conda activate build

   If necessary, install the conda-build tools: ::

     conda install conda-build conda-verify

3. Run the following command, which will create a sub-environment, install all dependencies
   listed in meta.yaml, and compile/package everything (it should take between 5 and 10
   minutes): ::

     conda-build . -c conda-forge

4. Once the package is built, you can install it in your current environment by typing: ::

     conda install -c anaconda -c conda-forge -c ${CONDA_PREFIX}/conda-bld/ nemaktis


