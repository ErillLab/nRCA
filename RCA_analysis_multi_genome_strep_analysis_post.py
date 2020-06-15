
import CodonUsagelib
import os
import sys
import time
from Bio import SeqUtils
import math


""" nRCA anaylsis
    Creates a normRelativeCodonAdaptationIndex object, generates an index from
    a reference set of genes (refset_file_name; FASTA format) and scores
    all the protein-coding genes for all provided accessions (one per line)
"""
# Copyright Ivan Erill (erill@umbc.edu)
# -*- coding: utf-8 -*-
from Bio import Entrez, SeqIO

def NCBIreadGBK(accession):
    """Gets an accession or GI number and downloads the GenBank record.
       Parses it with SeqIO and returns the sequence object.
    """
    net_handle = Entrez.efetch(db="nuccore",id=str(accession),
                                   rettype='gbwithparts', retmode="txt")
    gnome_record=SeqIO.read(net_handle, "genbank")
    net_handle.close()
    return gnome_record

def readGBK(filename):
    """Reads a GBK file and returns the object handle"""
    gnome_record = SeqIO.read(filename, "genbank")
    return gnome_record

def GBK_nRCA(genomeobject, cubindex, outfilename):
    """Takes in a GBK record object, a codon usage index and an output file.
       It goes through every CDS in the genome file, and scores its codon usage
       according to the provided index.
       
       Inputs:
       - genomeobject
           a parsed GBK genome object, as done by SeqIO.read (e.g. readGBK)
       - cubindex
           an instantiated CUB index
       - outfilename
           the output file name
    """
    
    #create interval set [0,0.05,0.10,0.15,0.20...] counter
    cnter = [0] * 20
    
    with open(outfilename,"a") as out_handle:
        ft_cnt=1
        #iterate features
        for feat in genomeobject.features:
            #get type features
            if feat.type in ['CDS', 'tRNA']:
                if feat.type == 'CDS':
                    nrca_val = nrca_index.nrca_for_gene(str(feat.location.extract(genomeobject).seq.upper()))
                    ps = int(math.floor(nrca_val*20))
                    cnter[ps]=cnter[ps] + 1
                else:
                    nrca_val = ''
                out_handle.write(genomeobject.description.split(',')[0] + ',')
                out_handle.write(genomeobject.id + ',')
                out_handle.write(str(SeqUtils.GC(genomeobject.seq)) + ',')
                out_handle.write(str(len(genomeobject.seq)) + ',')
                out_handle.write(feat.type + ',')
                if 'locus_tag' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['locus_tag'][0] + ',')
                elif 'note' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['note'][0] + ',')
                else:
                    out_handle.write('NO_LOCUS_TAG_' + str(ft_cnt) + ',')
                if 'protein_id' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['protein_id'][0] + ',')
                else:
                    out_handle.write('' + ',')
                out_handle.write(str(feat.location.nofuzzy_start) + '-')
                out_handle.write(str(feat.location.nofuzzy_end) + ' : ')
                out_handle.write(str(feat.location.strand))
                out_handle.write(',')
                out_handle.write(str(SeqUtils.GC(feat.location.extract(genomeobject).seq)) + ',')
                out_handle.write(str(len(feat.location.extract(genomeobject).seq)) + ',')
                out_handle.write(str(nrca_val))
                out_handle.write(',')
                if 'product' in feat.qualifiers:
                    out_handle.write(feat.qualifiers['product'][0])
                out_handle.write('\n')
            ft_cnt=ft_cnt+1
    
    #wite out histogram
    with open(outfilename.split('.')[0]+'_hist.csv',"a") as out_handle:
        out_handle.write(genomeobject.description.split(',')[0] + ',')
        out_handle.write(genomeobject.id + ',')
        out_handle.write(str(SeqUtils.GC(genomeobject.seq)) + ',')
        out_handle.write(str(len(genomeobject.seq)) + ',')
        for cnt in cnter:
            out_handle.write(str(cnt))
            out_handle.write(',')
        out_handle.write('\n')

    return(0)

################################################################################

Entrez.email = "EMAIL_HERE"
Entrez.apikey = "KEY_HERE"
	
infilename = 'streptophage_accessions_clean.txt'
outfilename = 'Streptophage.csv'
refset_file_name = 'S_griseus_ribosomal_CDS_NC_010572.fas'

#reset output file
f=open(outfilename,'w+')
f.write('Organism,GenomeID,GGC,GLen,Feat,LocusTag,ProtID,Location,GC,Len,nRCA,Product\n')
f.close()

f=open(outfilename.split('.')[0]+'_hist.csv','w+')
f.write('Organism,GenomeID,GGC,GLen,0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.0,0.85,0.90,0.95\n')
f.close()


# first make a nRCA object
nrca_index = CodonUsagelib.normRelativeCodonAdaptationIndex()
# now generate an index from a file
if os.path.exists(refset_file_name):
    nrca_index.generate_index(refset_file_name)
else:
    print("Cannot find the reference set file \n \
          Make sure it is in the appropriate folder")
    sys.exit()

with open(infilename,"r") as in_handle:
    lines=in_handle.readlines()
    for line in lines:
        gnomeID=line.strip('\n')
        genomefilename='genomes/' + gnomeID + '.gb'
        if not os.path.isfile(genomefilename):
            print "Downloading ", gnomeID
            mygenome = NCBIreadGBK(gnomeID)
            SeqIO.write(mygenome, genomefilename, "genbank")
            time.sleep(1)
        else:
            mygenome = readGBK(genomefilename)
        print "Processing ", gnomeID
        GBK_nRCA(mygenome, nrca_index, outfilename)

