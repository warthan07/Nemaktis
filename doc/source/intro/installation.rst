.. _install:

Installation
============

``Nemaktis`` is a mixture of C++ and python codes, and has been successfully tested both on
Windows 10 and Linux. The recommended way of installing and using Nematkis is through the
package manager ``mamba`` (see next subsection).

If you are an experienced Mac user and want to become maintainer of a MacOS version of
Nematkis (using mamba), please contact me
(https://nemaktis.readthedocs.io/en/latest/intro/overview.html#contact). Otherwise, it is
probably possible to compile and use the software by hand on other OSs than the ones
mentioned above (windows 7, macOS...), provided you know what you are doing :)

With Mamba (windows or linux)
-----------------------------

a. Install Miniforge
...................

The simplest way of installing the ``Nemaktis`` package is through the package manager
``mamba`` (no compilation and no dependency to be installed by hand). I provide precompiled
packages for both **Linux** and **Windows 10**, which means the following method will work
on these two operating systems. 

``mamba`` is a package manager based on the ``conda`` python package manager, but with far better
performance since it is written in C++. Previous versions of Nemaktis were using ``conda``
(or ``mamba`` installed on top of ``conda``) but this proved to be too difficult to maintain because
of these poor performances. For this reason, any Nemaktis version >=1.4.7 requires the 
``Miniforge`` distribution to be installed locally on your computer (i.e. on your user folder), 
and will not necessarily work with the ``Anaconda`` or ``Miniconda`` distributions. The installation 
files of ``Miniforge`` for Windows/Linux are available at this address (be careful to choose the
``Miniforge`` distribution, not the ``Miniforge-pypy3`` distribution):
https://github.com/conda-forge/miniforge#miniforge3

b. Install Nemaktis on the command line
.......................................

Open a terminal (Windows: application "Miniforge prompt" installed with Miniforge, Linux: any
terminal) and type the following command: ::
  
  mamba create -n nm -c warthan07 -y nemaktis

(Optional) If you want a python editor (e.g. ``spyder``), you have to also type the following
commands: ::

  mamba activate nm
  conda install spyder

The above command should have automatically created a shortcut for ``Spyder`` (there is currently
a bug in mamba for creating it, which is why the above command uses conda), usually named
``Spyder (nm)`` to make it clear that it will automatically activate the newly created ``nm``
environment. Always use this shortcut when using nemaktis (it can be searched quickly on windows
with the shortcut ``Win`` + ``S``, or ``Win`` on Ubuntu) or alternatively type the following
command in a terminal: ::

  mamba activate nm
  spyder

You can also, of course, run directly your python scripts in the terminal after activating the ``nm``
environment.

c. How to update
................

I do not recommend updating nematkis using the command ``mamba update``, since I do not
compile Nematkis sufficiently often for a correct update of all dependencies. In other
words, running ``mamba update`` has the risk of breaking your mamba environment! I
apologize for this current limitation, which mostly stems from my inexperience at mamba
packaging with complex dependencies. 

In this case, how can you safely update Nemaktis? The simplest way is to fully remove the
mamba environment ``nm``, either by removing the folder "envs/nm" inside the root folder of
Miniforge or by opening a terminal ("Miniforge prompt" for windows, any terminal for linux)
and typing: ::

  mamba remove --name nm --all 

Then, simply repeat the installation step b above.

Another possibility is to repeat the installation step b with a different environment name
than ``nm``, for example ``nm[version number]``. Although this method is probably fine to
test new features of the software without removing the old version , it is probably not very
good in the long term since each mamba environment takes a non-negligible portion of disk space. 

Developper method
-----------------

If you want to be able to modify the code, but still want to enjoy the simplicity of mamba
packages (no relative paths to manage, everything works as with a system package), you can build
yourselves the nemaktis package for Linux:

1. Get the code of the ``Nemaktis`` repository: ::

     git clone git@github.com:warthan07/Nemaktis.git

   And implements the changes that you want. For the C++ codes, compilation steps are provided
   in the subfolders bpm-solver and rt-solver if you want to test them locally (in which case
   you will have to install yourselves the dependencies listed in each CMakeLists).

2. In a terminal, go to the subfolder **conda_recipe** of the ``Nemaktis`` repository and type
   the following (the first line is only needed if the build environment does not exist yet): ::
     
     mamba create -n build -y rattler-build
     mamba activate build

3. Run the following command, which will create an ``output`` folder in which the conda package
   is built: ::

     rattler-build build

4. Once the package is built, you can install it in a new ``nm`` environment by typing: ::

     mamba create -n nm -c file://{ full path to output folder } nemaktis

