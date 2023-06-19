# CICERON
Classification bIoaCtive pEptides fRom micrObial fermeNtation

CICERON is a script for the functional classification of bioactive peptides specifically trained on BPs obtained from microbial fermentation. Starting from peptide sequences, nine binary classifiers assing a functional prediction to the bioactive peptide. 

## Installation
Create a conda environment with Python 3.9 and packages:

'''conda create -n ciceron python=3.9 pandas biopython scikit-learn numpy tensorflow'''

Enter the environment:

'''conda activate ciceron'''

Clone this git repository 

'''git clone https://github.com/BizzoTL/CICERON/'''

## Usage
Fasta protein files of interest must be present in the input folder.
From the cloned git folder, launch the following command:

"python3 model_prediction.py -i input_folder -o output_folder"

The results will be displayed in the output folder with the same name of the input file and the suffix chosen(facultative).
