import argparse
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from utils import *
import os


## it reads the output of the models, un-tokenises the predicted sequences and filters out unlikely metabolites
## -input_file: the csv file that has the input molecules (molecule ID and SMILES representations)
## -output_file: the filename where the processed predictions will be saved. It's a csv file. 
## predictions_file: the directory where the output of the models from the tranaslate_molecules script is saved
## predict_number: the number of predicted metabolite. It can be in [1,5,10,15]


def main(opt):
    input_file = opt.input_file
    output_file = opt.output_file
    predictions_file = opt.predictions_file
    predict_number = opt.predict_number
    figures_directory = 'Figures/'
    with open(predictions_file) as f_pred:
        pred_lines = [''.join(line.strip().split(' ')) for line in f_pred.readlines()]
    smiles2Name = {}
    smiles2metabolites = {}
    index = 0
    drug_lines = open(input_file).read().split('\n')
    if opt.visualise_molecules:
        if not os.path.exists(figures_directory):
            os.makedirs(figures_directory)

    for i in range(0, len(drug_lines) - 1):
        Name, smiles = drug_lines[i].split(',')
        if not check_smile(smiles):
            continue
        smiles = canonicalise_smile(smiles)
        smiles2Name[smiles] = Name
        predictions = set()
        for i in range(index, index + predict_number):
            predictions.add(pred_lines[i])
        index = index + predict_number
        processed, invalid, invalid_count,unrational_count = process_predictions(predictions, smiles, 0.25, True)
        if processed == []:
            smiles2metabolites[smiles] = ''
        else:
            smiles2metabolites[smiles] = processed
            drug = Chem.MolFromSmiles(smiles)
            preds = [Chem.MolFromSmiles(pred_smiles) for pred_smiles in processed]
            fig_dir = figures_directory + '/' + Name + '/'
            if not os.path.exists(fig_dir):
                os.makedirs(fig_dir)
            filename = fig_dir + Name + '.png'
            img = Draw.MolToFile(drug,filename,size=(500,500),wedgeBonds=False)
            prd_count = 1
            for prd in preds:
                filename = fig_dir + 'Metabolite' + str(prd_count) + '.png'
                img = Draw.MolToFile(prd, filename, size=(500, 500), wedgeBonds=False)
                prd_count = prd_count + 1

    table = ['Name', 'SMILES', 'predictions']
    for smiles in smiles2metabolites.keys():
        metabolites_str = ''
        Name = smiles2Name[smiles]
        metabolites = smiles2metabolites[smiles]
        for metabolite in metabolites:
            metabolites_str = metabolites_str + metabolite + ' '
        metabolites_str = metabolites_str[:-1]
        entry = [Name, smiles, metabolites_str]
        table = np.vstack((table, entry))

    with open(output_file, 'wb') as f:
        np.savetxt(f, table, fmt='%s', delimiter=',')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input_file', type=str, help='Input File')
    parser.add_argument('-output_file', type=str, help='Processed Predictions File')
    parser.add_argument('-predictions_file', type=str, help='Predictions Directory')
    parser.add_argument('-predict_number', type=int, help='The Number of Predict Metabolites')
    parser.add_argument('-visualise_molecules', type=bool, default=False, help='Visualise predicted metabolites')
    opt = parser.parse_args()
    main(opt)
