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

from pwchem.wizards.wizard_add_mutations import AA_THREE_TO_ONE

FOLDX_SM_FILE = 'FoldX_SM.tsv'
FOLDX_ZSCORE_FILE = 'FoldX_zscore.tsv'

# ------------------------------------ INSTALLATION VARIABLES ------------------------------------
PLUGIN_NAME = 'FOLDXSUITE'
FOLDX_HOME = 'FOLDX_HOME'

# Supported versions

# THE VERSION OF FOLDX CAN CHANGE IN A FUTURE!
V5_1 = '5.1'

# Plugin version
FOLDXSUITE_VERSION = '0.1'

# Protocol versions 
FOLDX_DEFAULT_VERSION = V5_1

# Protocol repo versions
FOLDX_REPO_DEFAULT_VERSION = V5_1
