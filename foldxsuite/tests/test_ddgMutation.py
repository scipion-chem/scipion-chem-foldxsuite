# **************************************************************************
# *
# * Authors:     Judith Maestro
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

# Scipion em imports
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols import ProtImportPdb

# Scipion chem imports
from ..protocols import ProtocolDDGFoldX
from pwchem.tests import TestImportBase

class TestDDGFoldX(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    def _runDDGFoldX(self):
        protDDG = self.newProtocol(
            ProtocolDDGFoldX,
            inputAtomStruct=self.protImportPDB.outputPdb,
            toMutateList='TA108A')
        
        self.launchProtocol(protDDG)
        setOut = getattr(protDDG, 'outputStats', None)
        self.assertIsNotNone(setOut)
        self.assertTrue(setOut.getSize() > 0)

    def test_DDGFoldX(self):
        self._runDDGFoldX()
