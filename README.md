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
Step 1: Download the trained models from the folder ./model. And create two folders with paths ./model/SoM_identifier and ./model/metabolite_predictor to store the models separately.

### Prepare parent drug

Step 2: Prepare a csv file with the name or id of the parent compound in the first column and SMILES in the second column ( Example files can refer to input.csv ).Then canonicalise and tokenise SMILES of parent drug:
```
python prepare_input_file.py -input_file input.csv -output_file processed_data.txt
```
###  Translate and Get predictions

Step 3: Use the bash file to transalte and get predictions: you can choose different bash file according to the number of predicted metabolite needed. Each bash file has three parameters. You need to define the source file (processed_data.txt , same as the output file in step 2) for the first parameter. For the second parameter, You need to define the file path where the predictions will be saved. In the prediction files, you can get three files, two txt files store the predicted SoM and metabolite information, one csv stores parent compound names, SMILES and predicted metabolites. For the third parameter, you need to define the input file (same as the input file in step 2). If you want not to visualise predicted metabolites, you can change the visualise_predictions parameter in the last line of bash file. The user can also define the min and max length of the predicted SMILES sequences, number of predicted SoMs and number of predicted metabolites ( n_best ).
```
./predict-top1.sh  src_file  prediction_file_path  input_file
./predict-top1.sh  processed_data.txt  ./prediction  input.csv
```
