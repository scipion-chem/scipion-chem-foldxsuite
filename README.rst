===========================
Scipion FoldX Suite plugin
===========================

This is a plugin for **scipion** wrapping some algorithms for analyzing the structure 
of proteins and protein complexes developed by the group of Luis Serrano at the `Centre 
for Genomic Regulation (CRG) <https://www.crg.eu>`_. 

These tools are available in the `FoldX Suite <https://foldxsuite.crg.eu>`_ and provides 
a fast and quantitative estimation of the importance of the interactions contributing to 
the stability of proteins and protein complexes. 

You need to download the FoldX Suite files before installing the plugin, see section "Binary 
Files" for details.


===========================
Install this plugin
===========================

You will need to use `Scipion3 <https://scipion-em.github.io/docs/docs/scipion
-modes/how-to-install.html>`_ to run these protocols.

1. **Binary Files**

**FoldX** binaries will **NOT** be downloaded automatically with the plugin.

The independent download of FoldX software suite by the user is required before running 
the programs. 

The installation path or any other of your preference has to be set in *FOLDX_HOME* in 
*scipion.conf* file.

The steps to download the Rosetta files are as follows:

    - Go to the FoldX software page <https://foldxsuite.crg.eu>`_.
    - FoldX is available through academic and commercial licenses. To download the software, 
      you need a FoldX Suite license. You can view and download the license types here 
      <https://foldxsuite.crg.eu/licensing-and-services>`_.
    - With the license (a link is sent to your email) you can download the FoldX files.
    - Save the package wherever you prefer.
    - Write the installation path in the *FOLDX_HOME* variable in the *scipion.conf* file.


2. **Install the plugin in Scipion**

- **Install the stable version (Not available yet)**

    Through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

    or by running command:

.. code-block::

    scipion3 installp -p scipion-chem-foldxsuite


- **Developer's version**

    1. **Download repository**:

    .. code-block::

        git clone https://github.com/scipion-chem/scipion-chem-foldxsuite.git

    2. **Install**:

    .. code-block::

        scipion3 installp -p path_to_scipion-chem-foldxsuite --devel


