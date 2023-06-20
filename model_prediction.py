import pandas as pd
import numpy as np
import os
from tensorflow.keras.models import load_model
import tensorflow_addons as tfa
from Bio import SeqIO
from itertools import product
import pickle
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="CICERON, mircrobial fermentation peptide classifier", epilog = "https://github.com/BizzoTL/CICERON", formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input_folder", type=str, help="Folder containing the fasta file(s)",required=True)
    parser.add_argument("-o", "--output_folder", type=str, help="Output folder where to store the results",required=True)
    parser.add_argument("-s", "--suffix", default="_predicted_peptides", type=str, help="Suffix to add to the name of the input file to generate the final output")
    return parser.parse_args()

args = parse_args()
input_folder = args.input_folder
output_folder = args.output_folder
suffix = args.suffix

#amino acid position for sparse encoding
aa_pos = {"A":0,"C":1,"D":2,"E":3,"F":4,"G":5,"H":6,"I":7,"K":8,"L":9,"M":10,"N":11,"P":12,"Q":13,"R":14,"S":15,"T":16,"V":17,"W":18,"Y":19,"-":20}

#blosum62 substitution matrix
aa_blosum62 = {"A":np.asarray([4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2,-100]),
          "C":np.asarray([0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,-100]),
          "D":np.asarray([-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3,-100]),
          "E":np.asarray([-1,-4,-2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2,-100]),
          "F":np.asarray([-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3,-100]),
          "G":np.asarray([0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3,-100]),
          "H":np.asarray([-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2,-100]),
          "I":np.asarray([-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,-100]),
          "K":np.asarray([-1,-3,-1,-1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2,-100]),
          "L":np.asarray([-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,-100]),
          "M":np.asarray([-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,-100]),
          "N":np.asarray([-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2,-100]),
          "P":np.asarray([-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3,-100]),
          "Q":np.asarray([-1,-3,-0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,-100]),
          "R":np.asarray([-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2,-100]),
          "S":np.asarray([1,-1,0,0,-2,0,-1,-2,0,-2,-1,-1,-1,0,-1,4,1,-2,-3,-2,-100]),
          "T":np.asarray([0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2,100]),
          "V":np.asarray([0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,-100]),
          "W":np.asarray([-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2,-100]),
          "Y":np.asarray([-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7,-100]),
          "-":np.asarray([-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,100,-100,-100,0])}

#amino acid combination for 3-mer encoding
aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
combinations_list = ["".join(comb) for comb in product(aa, repeat=3)]

def threemers(peptide_list):  #threemers encoding function
    pept_vector = np.zeros((len(peptide_list),8000))
    counter_list = 0
    for peptide in peptide_list:
        counts = {aa_comb: 0 for aa_comb in combinations_list}
        for i in range(len(peptide)-2):
            window = peptide[i:i+3]
            if window in counts:
                counts[window] += 1
        pept_vector[counter_list]=list(counts.values())
        counter_list +=1
    return pept_vector

def denser_peptide(peptide_list): #dense encoding function
    pept_vector = np.zeros((len(peptide_list),100))
    counter_list = 0
    for peptide in peptide_list:
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i]=aa_pos[peptide[i]]
        counter_list +=1
    return pept_vector

def blosum_substitution(peptide_list): #blosum matrix encoding function
    pept_vector = np.zeros((len(peptide_list),100,21))
    counter = 0
    counter_list = 0
    for peptide in peptide_list:
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i]=aa_blosum62[peptide[i]]
            counter += 1
        for j in range(counter,100):
            pept_vector[counter_list,j]=aa_blosum62["-"]
        counter_list +=1
        counter = 0
    lista = np.zeros((len(peptide_list),2100))
    for i in range(len(pept_vector)):
        lista[i] = pept_vector[i].flatten()
    return lista

def sparse_peptide(peptide_list): #sparse peptide function
    pept_vector = np.zeros((len(peptide_list),100,21,1))
    counter_list = 0
    for peptide in peptide_list:
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i,aa_pos[peptide[i]]]=1
        for j in range(len(peptide),100):
            pept_vector[counter_list,j,aa_pos["-"]]=1
        counter_list +=1
    lista = np.zeros((len(peptide_list),2100))
    for i in range(len(pept_vector)):
        lista[i] = pept_vector[i].flatten()
    return lista

def vectorize_peptide(peptide_list): #sparse encoding for NN
    pept_vector = np.zeros((len(peptide_list),100,21,1))
    counter = 0
    counter_list = 0
    for peptide in peptide_list:
        if type(peptide) == float:
            continue
        if len(peptide) > 100:
            peptide = peptide[0:100]
        for i in range(len(peptide)):
            pept_vector[counter_list,i,aa_pos[peptide[i]]]=1
            counter += 1
        for j in range(counter,100):
            pept_vector[counter_list,j,aa_pos["-"]]=1
        counter_list +=1
        counter = 0
    counter = 0
    counter_list = 0
    if len(peptide_list) == 1:
        pept_vector = np.squeeze(pept_vector,axis=3)
    else:
        pept_vector = np.squeeze(pept_vector)
    return pept_vector

#load all the models(from sklearn and keras presaved models)
antidiabetic = pickle.load(open("Models/Antidiabetic.sav", 'rb'))
antihypertensive = pickle.load(open("Models/Antihypertensive.sav", 'rb'))
antimicrobial = load_model("Models/Antimicrobial.h5", compile = True)
antioxidant = pickle.load(open("Models/Antioxidant.sav", 'rb'))
cardiovascular = pickle.load(open("Models/Cardiovascular.sav", 'rb'))
celiac = pickle.load(open("Models/Celiac disease.sav", 'rb'))
immunomodulatory = pickle.load(open("Models/Immunomodulatory.sav", 'rb'))
neuropeptides = pickle.load(open("Models/Neuropeptides.sav", 'rb'))
opioid = load_model("Models/Opioid.h5", compile = True)

modelli = [antidiabetic,antihypertensive,antioxidant,cardiovascular,celiac,immunomodulatory,neuropeptides]
modelli_NN = [antimicrobial,opioid]
modelli_nomi = ["antidiabetic","antihypertensive","antioxidant","cardiovascular","celiac","immunomodulatory","neuropeptides"]
modelli_nomi_NN = ["antimicrobial", "opioid"]
for filename in os.listdir(input_folder+"/"):
        if filename.endswith('.fasta'):
            #prepare the dataframe to be used for storing the info
            input_file = input_folder+"/"+filename
            records = list(SeqIO.parse(input_file, "fasta"))
            data_dict = {"Header": [], "Sequence": []}
            for i, record in enumerate(records):
                data_dict["Header"].append(record.description)
                data_dict["Sequence"].append(str(record.seq))
            df = pd.DataFrame.from_dict(data_dict)
            df["Function"] = np.nan
            appaiati = []
            #search on the unfiltered database if there is a 100% correspondance
            full_db = pd.read_csv("Peptide_all.csv")
            for i in range(len(df["Sequence"])):
                correspondance = []
                for j in range(len(full_db["Sequence"])):
                    if df.iloc[i][1] == full_db.iloc[j][1]:
                        correspondance.append(full_db.iloc[j][0])
                        continue
                if len(correspondance) == 0:
                    correspondance.append("")
                appaiati.append(correspondance)
            df["Function"] = appaiati

            #encoding the sequences in different ways
            sparse_encoding = sparse_peptide(df["Sequence"])
            threemer_encoded = threemers(df["Sequence"])
            dense_encoding = denser_peptide(df["Sequence"])
            blosum_encoding = blosum_substitution(df["Sequence"])
            sparse_NN = vectorize_peptide(df["Sequence"])
            #predict using keras NN models
            for m in range(len(modelli_NN)):
                predizione = modelli_NN[m].predict(sparse_NN) 
                df[modelli_nomi_NN[m]] = predizione.tolist() 
            #predict with sklearn
            for m in range(len(modelli)):
                if modelli_nomi[m] in ["antioxidant","cardiovascular","celiac","immunomodulatory"]:
                    predizione = modelli[m].predict_proba(threemer_encoded)
                elif modelli_nomi[m] == "antidiabetic":
                    predizione = modelli[m].predict_proba(sparse_encoding)
                elif modelli_nomi[m] == "neuropeptides":
                    predizione = modelli[m].predict_proba(dense_encoding)
                elif modelli_nomi[m] == "antihypertensive":
                    predizione = modelli[m].predict_proba(blosum_encoding)
                df[modelli_nomi[m]] = predizione.tolist() 
            #save file with the same name as file in input+_predicted_peptides
            filename = filename.strip(".fasta")
            df.to_csv(output_folder+"/"+filename+"_"+suffix+".csv", index = False)
