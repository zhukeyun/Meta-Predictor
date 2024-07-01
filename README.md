# MetaPredictor
MetaPredictor is a  end-to-end , prompt-based and transformer-base tool to predict human metabolites for small molecules. The methodology is described in detain in the paper _MetaPredictor: _in silico_ prediction of drug metabolites based on deep language models with prompt engineering_. The implementation of the Transformer model is based on the [OpenNMT toolkit](http://opennmt.net/OpenNMT-py/). The reaction data is processed based on the [RDKit](https://www.rdkit.org/) software.


## Installation
Create a conda environment:
```
conda create -n metapredictor python=3.8
source activate metapredictor
conda install rdkit -c rdkit
conda install pytorch=1.13.0
pip install -e .
```

## Predicting drug metabolites
### Download trained models
Step 1: Download the trained models from the folder /model.

### Prepare parent drug

Step 2: Prepare a csv file with the name or id of the parent compound in the first column and SMILES in the second column ( Example files can refer to input.csv ).Then canonicalise and tokenise SMILES of parent drug:
```
python prepare_input_file.py -input_file input.csv -output_file processed_data.txt
```
###  Translate

Step 3: Edit the bash file predict: define the source file (processed_data.txt ) and the beam size. The default beam size is 5 and the user can change it to get fewer or more predictions per molecule.  The user can also define the min and max length of the predicted SMILES sequences. Then, translate the input molecules into metabolites:
```
./translate.sh
```
### Get predictions

Step 4: Processing of prediction results and obtaining predictions:
```
python process_predictions.py -input_file input.csv -output_file predicted_metabolites.csv -beam_size ${beam} -visualise_predictions ${bool}
```
`input_file`  the name of the input csv file (same as the input file in step 2)  `output_file`  (optional) the name of the output csv file where the predictions will be saved. Default: predicted_metabolites.csv  `beam_size`  (optional) the beam size which has to be the same as the beam size used in step 3. Default: 5  `visualise_predictions`  (optional) if True then the predicted metabolites will be visualised using RDKit. Default: False.

This will generate a csv file with the predicted metabolites.
