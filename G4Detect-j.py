#!/usr/bin/env python3
# -*- coding: utf-8 -*-:
  
  
import os
import pandas as pd
import re
import sys
import argparse
from pprint import pprint

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'G4RNAScreener Result Annotation')
    parser.add_argument('-PI','--PathofIndex',type = str)
    parser.add_argument('-PR','--PathofRawdata',type = str)
    parser.add_argument('-ToG4H','--Threshold_G4H',type = float)
    parser.add_argument('-W','--Window',type  = int)
    parser.add_argument('-S','--Step',type = int)
    return parser

def InputParam(arg):
    InputParam = {'ToG4H' : float(arg.Threshold_G4H), 'Window' : int(arg.Window),'Step' : int(arg.Step)}
    return InputParam
    
def MergeSequence(g4Seq, windowSeq, step, windowLength):
    if len(windowSeq) < windowLength:
		# if the sequence is shorter than the sequence of the windows
		# (generraly last windows above the thresolds)
        g4Seq += windowSeq[-(len(windowSeq)-(windowLength-step)):]
    else: #
        g4Seq += windowSeq[-step:] # take the stepsize added
    return g4Seq
    
def createG4ID(geneid, genename, startG4, endG4, strand):
    return str(geneid) + '|' + str(genename) + '|'+ str(startG4) + '|' + str(endG4) + '|' + str(strand)


def FindAndMergeG4(G4Detected,InputParam,PathofRawdata,PathofIndex):
    """ Annotation basic information of all sequence by index file created before.

    This function try to annotate all information to the raw result 
    returned by RNAG4 Screener.Including there names and IDs. 
    Especially, the genome position of gene containing pG4 window 
    sequence, and the relative position of each window sequence are 
    added for further analysis.
    """
    RawDataofScreener = pd.read_table(PathofRawdata , header = None, sep = '\t', dtype = 'str', engine = 'python')
    RawDataofScreener.columns = ['Line_number','Chromosome','Start_position','End_position','Sequence','Relative_Start_position','cGcC','G4H','G4NN']
    RawDataofScreener[['Start_position','End_position','Relative_Start_position']] = RawDataofScreener[['Start_position','End_position','Relative_Start_position']].astype(int)
    RawDataofScreener[['G4H']] = RawDataofScreener[['G4H']].astype(float)
    # import the result from G4RNAScreener
    IndexofGene = pd.read_table(PathofIndex , header = None, sep = '\t', dtype = 'str', engine = 'python')
    IndexofGene.columns = ['Chromosome','Start_position','End_position','Gene_name','grey','Strand','Gene_ID']
    IndexofGene[['Start_position','End_position']] = IndexofGene[['Start_position','End_position']].astype(int)
    # import the index file
    MergeResult = pd.merge(RawDataofScreener,IndexofGene,on=['Chromosome','Start_position','End_position'])
    # annotation by pandas merge function
    MergeResult.insert(6,'Relative_position_end',1)
    MergeResult['Relative_position_end'] = MergeResult['Sequence'].str.len() + MergeResult['Relative_Start_position'] -1
    # add windowsequence relative end position
    MergeResult.to_csv(PathofRawdata+"MergeIndex2Result.csv",sep = '\t')
    """
    next,to merge the overlap windows and output key dict  
    """
    ForwWindowPassed = False
    NowWindowPassed = False
    ReadLinetmp = MergeResult.to_dict('records')
    for ReadLine in ReadLinetmp:
    # find the window over shreshold
        if (ReadLine['G4H'] > InputParam['ToG4H']):
            NowWindowPassed = True
            if (ForwWindowPassed != NowWindowPassed):
                WindowGeneID = ReadLine['Gene_ID']
                G4Sequence = ReadLine['Sequence']
                ForwWindowPassed = NowWindowPassed
                startG4repo = ReadLine['Relative_Start_position']
                endG4repo = ReadLine['Relative_position_end']
                listedG4H = [ ReadLine['G4H'] ]
            else:
                G4Sequence = MergeSequence(G4Sequence, ReadLine['Sequence'], InputParam['Step'], InputParam['Window'])
                endG4repo = ReadLine['Relative_position_end']
                listedG4H.append(ReadLine['G4H'])
        if (ReadLine['G4H'] < InputParam['ToG4H'] or WindowGeneID != ReadLine['Gene_ID']):
            NowWindowPassed = False
            if (ForwWindowPassed != NowWindowPassed):
                meanG4H = sum(listedG4H)/len(listedG4H)
                ForwWindowPassed = NowWindowPassed
                G4ID = createG4ID(ReadLine['Gene_ID'],ReadLine['Gene_name'],startG4repo,endG4repo,ReadLine['Strand'])
                if G4ID not in G4Detected and ReadLine['Strand']:
                    G4Detected[G4ID] = [G4Sequence,str(meanG4H)]
    return G4Detected

def GetListedG4(allG4Detected):
    AllListedG4 = {}
    for descriptionG4 in allG4Detected:
        geneID = descriptionG4.split('|')[0]
        if geneID not in AllListedG4:
            ListedG4 = []
        else:
            ListedG4 = listeG4InGene.get(geneID)
        ListedG4.append(descriptionG4)
        AllListedG4[geneID] = ListedG4
    return AllListedG4
        
def main():
    G4Detected={}
    G4MergeandDetect = FindAndMergeG4(G4Detected,InputParam,PathofRawdata,PathofIndex)
    G4Infomation = pd.DataFrame.from_dict(G4MergeandDetect, orient = 'index')
    G4Infomation = G4Infomation.reset_index()
    tmp_colname = ['GeneInfo','g4Seq','G4H_score']
    #rename columns names
    G4Infomation.columns = tmp_colname
    # split gene infomation
    G4Infomation['Gene_ID'] = G4Infomation['GeneInfo'].map(lambda x:x.split('|')[0])
    G4Infomation['Gene_name'] = G4Infomation['GeneInfo'].map(lambda x:x.split('|')[1])
    G4Infomation['startG4position'] = G4Infomation['GeneInfo'].map(lambda x:x.split('|')[2])
    G4Infomation['endG4position'] = G4Infomation['GeneInfo'].map(lambda x:x.split('|')[3])
    G4Infomation['Strand'] = G4Infomation['GeneInfo'].map(lambda x:x.split('|')[4])
    G4Infomation = G4Infomation[['GeneInfo','Gene_ID','Gene_name','startG4position','endG4position','Strand','g4Seq','G4H_score']]
    G4Infomation.to_csv(PathofRawdata + "G4Infomation.csv",sep = '\t')
    # give G4 number per gene
    G4NumPerGene = G4Infomation['Gene_ID'].value_counts()
    G4NumPerGenelist = G4NumPerGene.to_frame().reset_index()
    tmp_colname2 = ['Gene_ID','G4Number']
    G4NumPerGenelist.columns = tmp_colname2
    #Annotate G4 number per gene list
    G4NumPerGeneAnnotmp = pd.merge(G4NumPerGenelist,G4Infomation,on = ['Gene_ID'])
    G4NumPerGeneAnno = G4NumPerGeneAnnotmp.drop_duplicates(subset = 'Gene_ID')
    G4NumPerGeneAnnoResult = G4NumPerGeneAnno[['Gene_ID','Gene_name','G4Number']]
    G4NumPerGeneAnnoResult.to_csv(PathofRawdata + "G4NumberPerGeneAnno.csv",sep = '\t')
    

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    PathofRawdata = arg.PathofRawdata
    PathofIndex = arg.PathofIndex
    InputParam = InputParam(arg)
    main()
