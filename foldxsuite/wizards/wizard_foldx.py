# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Natalia del Rey
# *              Judith Maestro Ciria
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

from foldxsuite.protocols import ProtocolDDGFoldX

from pwchem.wizards import SelectChainWizardQT, SelectResidueWizardQT, AddMutationsWizard, ClearMutationsWizard

SelectChainWizardQT().addTarget(protocol=ProtocolDDGFoldX,
                              targets=['mutChain'],
                              inputs=['inputAtomStruct'],
                              outputs=['mutChain'])

SelectChainWizardQT().addTarget(protocol=ProtocolDDGFoldX,
                              targets=['ROIChain'],
                              inputs=['inputAtomStruct'],
                              outputs=['ROIChain'])


class AddMutationsFoldX(AddMutationsWizard):
    _targets = [(ProtocolDDGFoldX, ['addMutation'])]


class ClearMutationsFoldX(ClearMutationsWizard):
  _targets = [(ProtocolDDGFoldX, ['clearLabel'])]
