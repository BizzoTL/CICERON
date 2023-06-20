# CICERON
Classification bIoaCtive pEptides fRom micrObial fermeNtation

CICERON is a script for the functional classification of bioactive peptides specifically trained on BPs obtained from microbial fermentation. Starting from peptide sequences, nine binary classifiers assign a functional prediction to the bioactive peptide. The functional classes that are predicted are the following: Antidiabetic, Antihypertensive, Antimicrobial, Antioxidant, Cardiovascular, Celiac disease, Immunomodulatory, Neuropeptides and Opiopid. For more information see the following paper: ...

## Installation
Create a conda environment with Python 3.9 and packages:

```conda create -n ciceron python=3.9 pandas biopython scikit-learn numpy tensorflow```

To enter the environment:

```conda activate ciceron```

Clone this git repository 

```git clone https://github.com/BizzoTL/CICERON/```

## Usage
Fasta protein files of interest must be present in the input folder.
From the cloned git folder, launch the following command:

```python3 model_prediction.py -i input_folder -o output_folder -s suffix```

The results will be displayed in the output folder with the same name as the input file and the suffix chosen(facultative).

Example:

Fasta file "Peptides.fasta" in the folder peptide_inputs:
```
>Peptide_1
ATGIQPP
>Peptide_2
MNIPYPYP
```

To run CICERON to predict the functional classification of these peptides use the following command line from the CICERON folder:

``` python3 model_prediction.py -i path/to/peptide_inputs -o path/to/predicted_peptides -s predicted```

This will return a file called "Peptides_predicted.csv" inside the folder "predicted_peptides".
