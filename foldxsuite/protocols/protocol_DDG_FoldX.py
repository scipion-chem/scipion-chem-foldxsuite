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

"""
Wrapper around the FoldX method from https://foldxsuite.crg.es/documentation#manual
"""
import numpy as np
import os, re

from pyworkflow.constants import BETA
from pyworkflow.object import Object, Float, String
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message

from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct, SetOfStats
from pwchem.objects import SetOfStructROIs

import pwem.convert as emconv
from pwchem.utils.utils import cleanPDB

from foldxsuite import Plugin
from foldxsuite.constants import *

class ProtocolDDGFoldX(EMProtocol):
    """
    This protocol computes the change in free energy at the interface between two proteins
    when there is a mutation in one of the proteins. The result is returned standardized as 
    a z-score.
    """
    _label = 'DDG FoldX'
    _devStatus = BETA

    # -------------------------- DEFINE param functions ----------------------
    def _addMutationForm(self, form):
        form.addParam('multiPosition', params.BooleanParam, default=False,
                      label='Use a set of ROIs.',
                      help='Mutate and calculate the change in binding free energy (ΔΔG) '
                           'over a set of Regions Of Interest (ROIs).\nTo calculate ΔΔG '
                           'at specific positions, select "No" and directly specify the '
                           'mutations in "List of mutations".')
        
        form.addParam('ROIOrigin', params.EnumParam, default=0, condition='multiPosition',
                       label='Source of ROIs: ', choices=['Manual', 'SetOfStructROIs'],
                       help='Select the source of the regions of interest.')

        form.addParam('mutChain', params.StringParam, allowsNull=False, 
                      label='Chain to mutate', condition='ROIOrigin==0 and multiPosition',
                      help='Specify the protein chain to mutate.')
        
        form.addParam('RangPositions', params.StringParam, allowsNull=False,
                      label='Range of positions: ', condition='ROIOrigin==0 and multiPosition',
                      help='Specify the first and last index of each position range, separating '
                           'each range with a comma, i.e., "[FIRST_1]-[LAST_1], [FIRST_2]-[LAST_2]". '
                           'For example, "1-30, 50-70" will select for mutation all residues between '
                           'positions 1 and 30, and between positions 50 and 70 in the corresponding '
                           'chain.')
        
        form.addParam('inputStructROI', params.PointerParam, pointerClass="SetOfStructROIs",
                      label='Input structural ROI', condition='ROIOrigin==1 and multiPosition',
                      allowsNull=False, help='Select the source of the ROIs.') 

        form.addParam('mutSaturation', params.BooleanParam, default=True,
                       label='Saturation mutagenesis', condition='multiPosition',
                       help='Perform saturation mutagenesis, that is, replace each position '
                            'with each of the 20 protein-forming aminoacids (ACDEFGHIKLMNPQRSTVWY).')
        
        form.addParam('mutResidue', params.StringParam, allowsNull=False,
                      label="Residue to introduce", condition='multiPosition and not mutSaturation',
                      help='Define the substitute residue which will be introduced with its '
                           'one-letter code.\nFor the one-letter aminoacid code, see '
                           'https://foldxsuite.crg.eu/allowed-residues.')

        form.addParam('addMutation', params.LabelParam,
                      label='Add defined mutations', condition='multiPosition',
                      help='Add the defined mutations to the list of mutations below.')
        
        
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Structure templates')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                      label='Input atomic structure', allowsNull=False,
                      help='The atomic structure should have the protein-protein complex.')        
        
        group = form.addGroup('Define mutation')
        self._addMutationForm(group)
        group.addParam('toMutateList', params.TextParam, width=70,
                      default='', label='List of mutations:',
                      help='The syntax of a mutation is "[aaFrom][Chain][Position][aaTo]". For example, '
                           'the mutation "CA182Y", mutates position 182 of chain A that is a (C)ystein to '
                           'a t(Y)rosine.\nTo perform saturation mutagenesis (the amino acid is replaced '
                           'by each of the 20 protein-forming aminoacids), in the [aaTo] parameter specify '
                           '"X". For example, CA182X mutates Cys182 of the chain A to all protein-forming '
                           'aminoacids (ACDEFGHIKLMNPQRSTVWY).\nFor the one-letter aminoacid code, see '
                           'https://foldxsuite.crg.eu/allowed-residues.')       
        group.addParam('clearLabel', params.LabelParam,
                       label='Clear mutation list',
                       help='Clear mutations list')
        
        
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.computeDDG)
        self._insertFunctionStep(self.processResults)
        self._insertFunctionStep(self.calculateZScore)
        self._insertFunctionStep(self.createOutputStep)

    def computeDDG(self):
        fnPDB = "atomicStructure.pdb"
        cleanPDB(self.inputAtomStruct.get().getFileName(),fnPDB)
        
        fnMutL = []
        for i, line in enumerate(self.toMutateList.get().strip().split('\n')):
            pattern = re.compile(r'([A-Za-z]+)([A-Za-z]+)([^a-zA-Z]+)([A-Za-z]+)')
            match = re.match(pattern, line)
            if match:
                aaFrom, chain, position, aaTo = match.groups()
                mut = aaFrom + chain + position + "a"
            if mut not in fnMutL:
                fnMutL.append(mut)
        fnMut = ",".join(fnMutL)

        resultsDir = self._getExtraPath('Results_FoldX')
        if not os.path.exists(resultsDir):
            os.makedirs(resultsDir)

        args='--command=Pssm --pdb="%s" --positions="%s" --output-dir=%s'%(fnPDB, fnMut, resultsDir)
        Plugin.runFOLDX(self, args=args)
                
        os.remove(fnPDB)
    
    def processResults(self):
        pssmFile = os.path.join(self._getExtraPath('Results_FoldX'), 'PSSM_atomicStructure.txt')
        pssmProcess = self._getExtraPath('FOLDX_SM_FILE')
        
        outDdgSm = ""

        with open(pssmFile, "r") as foutput, open(pssmProcess, "w") as fddg:
            lines = foutput.readlines()
            residues = lines[0].strip().split()  
            
            for line in lines[1:]:
                parts = line.strip().split()
                mutationLabel = parts[0]
                energies = parts[1:]

                for i, residue in enumerate(residues):
                    mutation = f"{mutationLabel}{residue}"
                    energy = float(energies[i])
                    outDdgSm += f"{mutation}\t{energy}\n"

            outDdgSm = outDdgSm.rstrip()
            fddg.write(outDdgSm)

    def calculateZScore(self):
        pssmProcess = self._getExtraPath('FOLDX_SM_FILE')
        ddgUser = self._getExtraPath('FOLDX_ZSCORE_FILE')

        userMutations = self.toMutateList.get().strip().split('\n')

        with open(pssmProcess, "r") as f:
            muts = f.read().split("\n")
            mutDict = {line.split("\t")[0]: float(line.split("\t")[1]) for line in muts}

            # Calculating averages and standard deviations
            values = [mutDict[key] for key in mutDict]
            avg = np.mean(values)
            std = np.std(values)

            # Calculating Z-scores and consensus Z-scores
            allZscoresStr = "Mut\tddg\tzscore\n"  
            userZscoresStr = "Mut\tzscore\n"     
            for key in mutDict:
                ddg = mutDict[key]
                zscore = (ddg - avg) / std        
                mutDict[key] = zscore
                allZscoresStr += f"{key}\t{ddg}\t{zscore}\n"

                for userMut in userMutations:
                    if userMut.endswith("X"):
                        baseMut = userMut[:-1]  
                        if key.startswith(baseMut):
                            userZscoresStr += f"{key}\t{zscore}\n"
                    elif userMut == key:
                        userZscoresStr += f"{key}\t{zscore}\n"

        allZscoresStr = allZscoresStr.rstrip()  
        userZscoresStr = userZscoresStr.rstrip()  

        with open(pssmProcess, "w+") as fddg:
            fddg.write(allZscoresStr)

        with open(ddgUser, "w+") as fuser:
            fuser.write(userZscoresStr)

    def createOutputStep(self):
        foldxProcess = self._getExtraPath('FOLDX_SM_FILE')
        ddgUser = self._getExtraPath('FOLDX_ZSCORE_FILE')
        outputSet = SetOfStats.create(self.getPath())
        mutations = []
        zscoreMap = {}

        with open(ddgUser, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                if not line.strip():
                    continue
                cols = line.strip().split('\t')
                mutName = cols[0]
                mutations.append(mutName)
                if len(cols) > 2:
                    try:
                        zscoreMap[mutName] = float(cols[2])
                    except ValueError:
                        pass
        with open(foldxProcess, "r") as f:
            results = f.readlines()
        for line in results[1:]:
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            mutName = fields[0]
            if mutName not in mutations:
                continue
            item = Object()
            item.setObjLabel(label=mutName)
            item.mutation = String(mutName)
            try:
                item.ddg = Float(float(fields[1]))
            except (ValueError, IndexError):
                item.ddg = Float(0.0)
            if mutName in zscoreMap:
                item.zscore = Float(zscoreMap[mutName])
            else:
                try:
                    item.zscore = Float(float(fields[2]))
                except (ValueError, IndexError):
                    item.zscore = Float(0.0)
            outputSet.append(item)

        self._defineOutputs(outputStats=outputSet)
        self._defineTransformRelation(self.inputAtomStruct, outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []   

        structureHandler = emconv.AtomicStructHandler()
        structureHandler.read(self.inputAtomStruct.get().getFileName())
        structureHandler.getStructure()
        modelsLength, modelsFirstResidue = structureHandler.getModelsChains()
        
        validChains = set()
        chainResidues = {}

        for modelID, chains in modelsFirstResidue.items():
            for chainID, residues in chains.items():
                filteredResidues = [res for res in residues if res[1] != 'HOH']
                validChains.add(chainID)
                if chainID not in chainResidues:
                    chainResidues[chainID] = filteredResidues

        if not self.toMutateList.get().strip():
            errors.append('You have not added any mutation to the list. Do so using the "Add defined '
                          'mutations" wizard once you have defined it.')

        else:
            for i, line in enumerate(self.toMutateList.get().strip().split('\n')):
                pattern = re.compile(r'([A-Za-z]+)([A-Za-z]+)([^a-zA-Z]+)([A-Za-z]+)')
                match = re.match(pattern, line)

                if match:
                    aaFrom, chain, position, aaTo = match.groups()

                    if chain not in validChains:
                        errors.append(f'The chain "{chain}" of the mutation "{line}" is not present in the PDB file. '
                                        f'The PDB file contains the following chains: {", ".join(validChains)}.')
                    
                    elif not position.isdigit():
                        errors.append(f'The position of the mutation "{line}" must be an integer.')
                    
                    elif aaFrom not in AA_THREE_TO_ONE.values():
                        errors.append(f'The wild-type aminoacid of the mutation "{line}" does not ' 
                                        'exist or is not written with its one-letter code.')
                    
                    elif aaTo not in AA_THREE_TO_ONE.values():
                        errors.append(f'The mutant aminoacid of the mutation "{line}" does not ' 
                                        'exist or is not written with its one-letter code.')

                    else:   
                        position = int(position)             
                        residuesDict = {res[0]: res[1] for res in chainResidues[chain]}
                        if position not in residuesDict:
                            firstResidue = next(iter(residuesDict))
                            lastResidue = list(residuesDict)[-1]
                            errors.append(f'Position "{position}" in chain "{chain}" for mutation "{line}" is out of range. '
                                            f'The chain "{chain}" has positions from {firstResidue} to {lastResidue}.')
                        
                        elif AA_THREE_TO_ONE[residuesDict[position]] != aaFrom:
                            errors.append(f'The wild-type aminoacid "{aaFrom}" at position "{position}" in chain "{chain}" '
                                            f'for mutation "{line}" does not match the PDB file. The aminoacid at that position '
                                            f'is {residuesDict[position]} ({AA_THREE_TO_ONE[residuesDict[position]]}).')
                else:
                    errors.append(f'The mutation "{line}" does not have the 4 necessary parameters. ' 
                                   'Mutation format must be "[aaFrom][Chain][Position][aaTo]".')
        return errors

    def _summary(self):
        summary = []
        ddgFile = self._getExtraPath('FOLDX_ZSCORE_FILE')
        if os.path.exists(ddgFile):
            with open(ddgFile) as f:
              summary.append(f.read())    
        return summary
    

    def _methods(self):
        methods = []
        methods.append("Prediction of the binding free energy change (ΔΔG) for protein-protein interactions "
                       "due to a point mutation in an aminoacid using the program FoldX." 
                       "\nThe result is standardized as a z-score.")
        return methods
    
    def _citations(self):
        return ['schymkowitz2005']
    