# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Natalia del Rey
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


from foldxsuite.protocols import ProtocolFoldX
from foldxsuite.constants import *

from pwem.wizards import EmWizard
import pwem.convert as emconv


class AddMutationsFoldX(EmWizard):
    _targets = [(ProtocolFoldX, ['addMutation'])]    
    
    def getPositions(self, form):
        protocol = form.protocol
        ROIOrigin = protocol.ROIOrigin.get()

        if ROIOrigin == 0:
            allRanPos = protocol.RangPositions.get().split(", ")
        else:
            structROI = protocol.inputStructROI.get()
            allRanPos = []
            for item in structROI:
                chain_res = item.getDecodedCResidues()
                for roi in chain_res:
                    res = roi.split("_")
                    if f"{res[1]}-{res[1]}" not in allRanPos:
                        allRanPos.append(f"{res[1]}-{res[1]}")
        return allRanPos

    def getaaTo(self, form):
        protocol = form.protocol
        return "X" if protocol.mutSaturation else str(protocol.mutResidue.get())

    def getchainResidues(self, form):
        protocol = form.protocol
        structureHandler = emconv.AtomicStructHandler()
        structureHandler.read(protocol.inputAtomStruct.get().getFileName())
        structureHandler.getStructure()
        modelsLength, modelsFirstResidue = structureHandler.getModelsChains()
        chainResidues = {}

        for modelID, chains in modelsFirstResidue.items():
            for chainID, residues in chains.items():
                if chainID not in chainResidues:
                    chainResidues[chainID] = {}
                for residue in residues:
                    res_id = residue[0]
                    res_type = residue[1]
                    chainResidues[chainID][res_id] = res_type

        return chainResidues
    
    def getMutations(self, form):
        allRanPos = self.getPositions(form)
        aaTo = self.getaaTo(form)
        chainResidues = self.getchainResidues(form)
        mutations = []     
        
        for ranPos in allRanPos:
            ran = ranPos.split("-")
            for chain, residues_dict in chainResidues.items():
                for pos in range(int(ran[0]), int(ran[1]) + 1):
                    if pos in residues_dict:  
                        aaFrom = AA_THREE_TO_ONE[residues_dict[pos]]
                        mutation = '{}{}{}{}'.format(aaFrom, chain, pos, aaTo)
                        mutations.append(mutation)
        return mutations
    
    def show(self, form, *params):
        protocol = form.protocol
        mutations = self.getMutations(form)

        toMutateList = protocol.toMutateList.get()
        toMutateList += "\n" + "\n".join(mutations)
        form.setVar('toMutateList', toMutateList.strip())


class ClearMutationsFoldX(EmWizard):
  _targets = [(ProtocolFoldX, ['clearLabel'])]

  def show(self, form, *params):
    form.setVar('toMutateList', '')
