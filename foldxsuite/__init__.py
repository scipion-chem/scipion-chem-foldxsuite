# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Natl. Center of Biotechnology CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# General imports
import os

# Scipion em imports
import pyworkflow.utils as pwutils
import pwem

# Plugin imports
from .constants import *

__version__ = "0.1"  # plugin version
_logo = "foldX_icon.png"
_references = ['schymkowitz2005']


class Plugin(pwem.Plugin):
    """
    Definition of class variables. For each package, a variable will be created.
    _<packageNameInLowercase>Home will contain the full path of the package, ending with a folder whose name will be <packageNameFirstLetterLowercase>-<defaultPackageVersion> variable.
    
    Inside that package, for each binary, there will also be another variable.
    _<binaryNameInLowercase>Binary will be a folder inside _<packageNameInLowercase>Home and its name will be <binaryName>.
    """
    _url = "https://github.com/scipion-chem/scipion-chem-foldxsuite"
    _supportedVersions = V1_0  # binary version


    @classmethod
    def _defineVariables(cls):
        """
        Return and write a home and conda enviroment variable in the config file.
        Each package will have a variable called <packageNameInUppercase>_HOME, and another called <packageNameInUppercase>_ENV
        <packageNameInUppercase>_HOME will contain the path to the package installation."
        <packageNameInUppercase>_ENV will contain the name of the conda enviroment for that package."
        """
        cls._defineVar(FOLDX_HOME, "FoldX")
        cls._defineEmVar(FOLDX_HOME, f"FoldX-{V1_0}")


    @classmethod
    def defineBinaries(cls, env):
        """
        This function defines the binaries for each protocol.
        """
        pass


    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch my program. """
        environ = pwutils.Environ(os.environ)
        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. """
        neededProgs = []
        return neededProgs

    @classmethod
    def runFOLDX(cls, protocol, args):
        """ 
        Run FoldX command from a given protocol. 
        """        
        args += ' --rotabaseLocation="%s"'%(os.path.join(cls.getVar(FOLDX_HOME),"rotabase.txt"))
        protocol.runJob(os.path.join(cls.getVar(FOLDX_HOME),"foldx_20241231"), args)
