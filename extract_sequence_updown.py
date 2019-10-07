#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 09:34:19 2019

@author: ls
"""
print('>==starting the sequence extracting==<')
import pandas as pd
import os
import argparse
import warnings
warnings.filterwarnings('ignore')
parser=argparse.ArgumentParser(description='the option for HAU genome for hirustum')                           
parser.add_argument('-gn','--gene_name',help='the name of the target genes',type=str)
parser.add_argument('-G','--genome_file',help='the genome file',type=str)
parser.add_argument('-gf','--gff_file',type=str,help='the gff file')
parser.add_argument('-rf','--result',type=str,help='the result dirctory')
parser.add_argument('-cl','--class_type',type=str,help='the type of sequence')
parser.add_argument('-up','--upstream',type=int,help='the number of upstream')
parser.add_argument('-down','--downstream',type=int,help='the number of downstream')
parser.add_argument('-o','--out',type=str,help='the format of output file')
args=parser.parse_args()
gene_name_inputfile=args.gene_name
genome_file_input=args.genome_file
gene_gff_input=args.gff_file
sequence_out=args.result
select_class=args.class_type
upstream=args.upstream
downstream=args.downstream
out=args.out
#look through the files input
if os.path.isfile(gene_name_inputfile):
    print('the gene name file exists, go on !')
    if os.access(gene_name_inputfile,os.R_OK):
        print('the gene name file is accessible')
    else:
        exit('gene name file is not compatible for read')
else:
    exit('you need to provide the gene names for extract')

if os.path.isfile(genome_file_input):
    print('the genome file exists, go on !')
    if os.access(genome_file_input,os.R_OK):
        print('the genome file is accessible')
    else:
        exit('the genome file is not compatible for read')
else:
    exit('you need to provide genome file')
if os.path.isfile(gene_gff_input):
    print('the gff file exists, go on!')
    if os.access(gene_gff_input,os.R_OK):
        print('gff file is accessible')
    else:
        exit('gff file is not compatible for read')
else:
    exit('you need to provide the gff file')
with open(gene_name_inputfile,'r') as file:
    gene_name_list=[]
    for line in file:
        line=line.strip()
        if line=='\n':
            continue
        else:            
            gene_name_list.append(line)
#create target_gene_gff
genome_gff=pd.read_table(gene_gff_input,header=None)  
gene_gff=genome_gff.loc[genome_gff[2]==select_class]
if gene_gff.shape[0]==0:
    exit('no select item in the gff')
gene_id=[]
for i in gene_gff[8]:
    i=i[3:18]
    gene_id.append(i)
gene_gff[8]=gene_id
target_gene_gff=gene_gff.loc[gene_gff[8].isin(gene_name_list)]
if target_gene_gff.shape[0]==0:
    exit('genes in gene name is not contained in gff file')
del genome_gff
del gene_gff
#seperate the genome.fa file
fa_name=[]
fa_seq=[]
seq=[]
with open(genome_file_input,'r') as file:
    for line in file:
        line=line.strip()
        if line=='\n':
            continue
        elif line[0]=='>':
            fa_name.append(line)
            s=''.join(fa_seq)
            seq.append(s)
            s=''
            fa_seq=[]
        else:
            fa_seq.append(line)
del seq[0]
seq.append(fa_seq)
fa_seq=[]
if len(seq)!=len(fa_name):
    exit('the genome file is not in correct format')
if os.path.isdir('seperate_chromosome')==True:
    pass
else:
    os.mkdir('seperate_chromosome')
    for i,j in zip(fa_name,seq):
    	with open('seperate_chromosome/%s.fa'%i[1:len(i)],'w') as file:
    	    j=''.join(j)
    	    file.write('>%s'%i+'\n')
    	    file.write(j)
#def get sequence function
def get_sequence(list_gene,s,up,down):
    if list_gene[0]-int(up)<=0:
        up=list_gene[0]
    if list_gene[1]+int(down)>=len(s):
        down=len(s)-list_gene[1]
    if list_gene[2]=='+':
        gene_seq=s[list_gene[0]-1:list_gene[1]]
        up_seq=s[list_gene[0]-int(up):list_gene[0]-1]
        down_seq=s[list_gene[1]:list_gene[1]+int(down)]
    else:
        table=str.maketrans('ATCG','TAGC')
        s=s.translate(table)
        s=s[::-1]
        gene_seq=s[list_gene[0]-1:list_gene[1]]
        up_seq=s[list_gene[0]-int(up):list_gene[0]-1]
        down_seq=s[list_gene[1]:list_gene[1]+int(down)]
    return gene_seq,up_seq,down_seq
genome_list=os.listdir('seperate_chromosome')
if os.path.isdir(sequence_out)==True:
    pass
else:  
    os.mkdir(sequence_out)
result_list=[]
for i in gene_name_list:
    with open('seperate_chromosome/%s.fa'%i[0:8],'r') as file:
        for line in file:
            line=line.strip()
            if line[0]=='>':
                continue
            else:
                s=line
    target_gff=target_gene_gff.loc[target_gene_gff[8]==i]
    list_input=[target_gff.iloc[0,3],target_gff.iloc[0,4],target_gff.iloc[0,6]]
    result=get_sequence(list_input,s,upstream,downstream)
    result_list.append(result)
if out=='merge':   
    with open(sequence_out+'/result.txt','w') as file:
        for i,j in zip(result_list,gene_name_list):
            file.write('gene%s: '%(j)+i[0]+'\n')
            file.write('up%s: '%(j)+i[1]+'\n')
            file.write('down%s:'%(j)+i[2]+'\n')
else:
    with open(sequence_out+'/gene.fa','w') as file:
        for i,j in zip(result_list,gene_name_list):
            file.write('>'+j+'\n')
            file.write(i[0]+'\n')
    with open(sequence_out+'/upstream.fa','w') as file:
        for i,j in zip(result_list,gene_name_list):
            file.write('>'+j+'\n')
            file.write(i[1]+'\n')
    with open(sequence_out+'/downstream.fa','w') as file:
        for i,j in zip(result_list,gene_name_list):
            file.write('>'+j+'\n')
            file.write(i[2]+'\n')
print('>===the sequences have been extracted')
        
