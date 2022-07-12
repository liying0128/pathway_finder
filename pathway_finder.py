from Bio import SeqIO
import os
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import swalign
from Bio.SeqRecord import SeqRecord
match=2
mismatch=-1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring)

def pathwayfinder(template_name,genome_name):
    allseq = [seq_record.seq for seq_record in SeqIO.parse(genome_name, "fasta")]
    allid= [seq_record.id for seq_record in SeqIO.parse(genome_name, "fasta")]
    templateseq=[seq_record.seq for seq_record in SeqIO.parse(template_name, "fasta")][0]
    templateid=[seq_record.id for seq_record in SeqIO.parse(template_name, "fasta")][0]
    tempscore=[]
    geneid=[]
    for k in range(len(allseq)):
        score=sw.align(templateseq,allseq[k]).score
        tempscore.append(score)
    for i in range(len(tempscore)):
        if tempscore[i]>len(templateseq)*0.1 and tempscore[i]>30:
            geneid.append(allid[i])
    plt.plot(tempscore,color='g',linewidth=4)
    plt.title(geneid)
    plt.ylim(10,int(len(templateseq)*2.1))
    plt.savefig(template_name+'.png',dpi=500)

genome_name='C:\\Users\\Administrator\\Desktop\\AA-seq\\Methylotenera.fasta'
pathway='C:\\Users\\Administrator\\Desktop\\AA-seq\\代谢途径-AA-seq\\'
os.chdir(pathway)
filelist=os.listdir()
for m in filelist[]:
    current_pathway=pathway+m+'\\'
    os.chdir(current_pathway)
    current_filelist=os.listdir()
    for n in current_filelist:
        template_name=n
        pathwayfinder(template_name, genome_name)
    