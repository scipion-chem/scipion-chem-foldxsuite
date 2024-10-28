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

"""
Wrapper around the FoldX method from http://compbio.clemson.edu/saambe_webserver/
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

class ProtocolFoldX(EMProtocol):
    """
    This protocol computes the change in free energy at the interface between two proteins
    when there is a mutation in one of the proteins. The result is returned standardized as 
    a z-score.
    """
    _label = 'FoldX'
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

        results_dir = self._getExtraPath('Results_FoldX')
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)

        args='--command=Pssm --pdb="%s" --positions="%s" --output-dir=%s'%(fnPDB, fnMut, results_dir)
        Plugin.runFOLDX(self, args=args)
                
        os.remove(fnPDB)
    
    def processResults(self):
        pssm_file = os.path.join(self._getExtraPath('Results_FoldX'), 'PSSM_atomicStructure.txt')
        pssm_process = self._getExtraPath('FoldX_SM.tsv')
        
        out_ddg_sm = ""

        with open(pssm_file, "r") as foutput, open(pssm_process, "w") as fddg:
            lines = foutput.readlines()
            residues = lines[0].strip().split()  
            
            for line in lines[1:]:
                parts = line.strip().split()
                mutation_label = parts[0]
                energies = parts[1:]

                for i, residue in enumerate(residues):
                    mutation = f"{mutation_label}{residue}"
                    energy = float(energies[i])
                    out_ddg_sm += f"{mutation}\t{energy}\n"

            out_ddg_sm = out_ddg_sm.rstrip()
            fddg.write(out_ddg_sm)

    def calculateZScore(self):
        pssm_process = self._getExtraPath('FoldX_SM.tsv')
        ddg_user = self._getExtraPath('FoldX_zscore.tsv')

        user_mutations = self.toMutateList.get().strip().split('\n')

        with open(pssm_process, "r") as f:
            muts = f.read().split("\n")
            mut_dict = {line.split("\t")[0]: float(line.split("\t")[1]) for line in muts}

            # Calculating averages and standard deviations
            values = [mut_dict[key] for key in mut_dict]
            avg = np.mean(values)
            std = np.std(values)

            # Calculating Z-scores and consensus Z-scores
            all_zscores_str = "Mut\tddg\tzscore\n"  
            user_zscores_str = "Mut\tzscore\n"     
            for key in mut_dict:
                ddg = mut_dict[key]
                zscore = (ddg - avg) / std        
                mut_dict[key] = zscore
                all_zscores_str += f"{key}\t{ddg}\t{zscore}\n"

                for user_mut in user_mutations:
                    if user_mut.endswith("X"):
                        base_mut = user_mut[:-1]  
                        if key.startswith(base_mut):
                            user_zscores_str += f"{key}\t{zscore}\n"
                    elif user_mut == key:
                        user_zscores_str += f"{key}\t{zscore}\n"

        all_zscores_str = all_zscores_str.rstrip()  
        user_zscores_str = user_zscores_str.rstrip()  

        with open(pssm_process, "w+") as fddg:
            fddg.write(all_zscores_str)

        with open(ddg_user, "w+") as fuser:
            fuser.write(user_zscores_str)

    def createOutputStep(self):
        foldx_process = self._getExtraPath('FoldX_SM.tsv')  
        ddg_user = self._getExtraPath('FoldX_zscore.tsv')
      
        outputSet = SetOfStats.create(self.getPath())
        
        mutations = []
        with open(ddg_user, "r") as f:
            content = f.readlines()
            for line in content[1:]:  
                mut = line.split('\t')
                mutations.append(mut[0])

        with open(foldx_process, "r") as f:
            results = f.readlines()
        
        for line in results[1:]:                
            fields = line.strip().split("\t")
            
            if fields[0] in mutations:
                newItem = Object()
                newItem.setObjLabel(label=str(fields[0]))
                setattr(newItem, 'ddg', Float(fields[1]))   
                setattr(newItem, 'zscore', Float(fields[2])) 
                outputSet.append(newItem)

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
                filtered_residues = [res for res in residues if res[1] != 'HOH']
                validChains.add(chainID)
                if chainID not in chainResidues:
                    chainResidues[chainID] = filtered_residues

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
                        residues_dict = {res[0]: res[1] for res in chainResidues[chain]}
                        if position not in residues_dict.keys():
                            first_residue = list(residues_dict.keys())[0]
                            last_residue = list(residues_dict.keys())[-1]
                            errors.append(f'Position "{position}" in chain "{chain}" for mutation "{line}" is out of range. '
                                            f'The chain "{chain}" has positions from {first_residue} to {last_residue}.')
                        
                        elif AA_THREE_TO_ONE[residues_dict[position]] != aaFrom:
                            errors.append(f'The wild-type aminoacid "{aaFrom}" at position "{position}" in chain "{chain}" '
                                            f'for mutation "{line}" does not match the PDB file. The aminoacid at that position '
                                            f'is {residues_dict[position]} ({AA_THREE_TO_ONE[residues_dict[position]]}).')
                else:
                    errors.append(f'The mutation "{line}" does not have the 4 necessary parameters. ' 
                                   'Mutation format must be "[aaFrom][Chain][Position][aaTo]".')
        return errors

    def _summary(self):
        summary = []
        ddgFile = self._getExtraPath('FoldX_zscore.tsv')
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
    