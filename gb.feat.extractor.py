#!/usr/bin/env python

import os
import sys
import glob
from Bio import SeqIO
import argparse

def main():

    options=argparse.ArgumentParser(prog="gbFeatExtractor",
                usage=argparse.SUPPRESS,
                description='gbFeatExtractor: A tool for extracting sequence features from a GenBank annotation file',
                prefix_chars='-',
                add_help=True,
                epilog='Written by Chrispin Chaguza, St Jude Children\'s Research Hospital, 2026')

    options.add_argument('--in','-i',action='store',required=True,nargs=1,
                        metavar='inputfile',dest='inputfile',
                        help='Input sequence annotation file in GenBank format')
    options.add_argument('--out','-o',action='store',required=False,nargs=1,
                        metavar='outputfile',dest='outputfile',default='gb.features.tsv',
                        help='Output file containing extracted sequence features')
    options.add_argument('--seqfmt','-s',action='store',required=False,nargs=1,choices=["genbank","embl"],
                        metavar='seqformat',dest='seqformat',default="genbank",
                        help='Input sequence annotation format')
    options.add_argument('--outfmt','-m',action='store',required=False,nargs=1,
                        metavar='outputformat',dest='outputformat',default="tsv",
                        help='Output file containing extracted sequence features')
    options.add_argument('--feat','-f',action='store',required=False,nargs=1,
                        metavar='feattype',dest='feattype',default="cds",
                        help='Sequence feature type to extract')
    options.add_argument('--translate','-t',action='store_true',default=False,
                        dest='translate',help='Translate coding sequences')
    options.add_argument('--verbose','-v',action='store_true',default=False,
                        dest='verbose',help='Show progress')

    options=options.parse_args(args=None if sys.argv[2:] else ['--help'])


    cmdValues = {'inputFile': options.inputfile[0:][0],
                 'outputFile': options.outfilefile[0:][0] if isinstance(options.outputfile,list) else options.outputfile,
                 'seqformat': options.seqformat[0:][0] if isinstance(options.seqformat,list) else options.seqformat,
                 'outfmt': int(options.outputformat[0:][0]) if isinstance(options.outputformat,list) else options.outputformat,                 
                 'feattype': options.feattype[0:][0] if isinstance(options.feattype,list) else options.feattype,
                 'translate': options.translate,
                 'verbose': options.verbose}

    inputSeqFile=str(cmdValues["inputFile"]).lower()
    inputSeqFormat=str(cmdValues["seqformat"]).lower()
    outputFormat=str(cmdValues["outfmt"]).lower()
    seqFeatureType=str(cmdValues["feattype"]).lower()
    translateCodingSeq=cmdValues["translate"]
    verbose=cmdValues["verbose"]

    seqObj=SeqIO.read(inputSeqFile,inputSeqFormat)

    for seqFeat in seqObj.features:
        if seqFeat.type == "source":
            continue
        else:
            pass

        dnaSeq = []

        if seqFeat.location.strand == 1:
            for loc in seqFeat.location.parts:
                dnaSeq.append(str(seqObj.seq[seqFeat.location.start:seqFeat.location.end]))
        else:
            for loc in seqFeat.location.parts:
                dnaSeq.append(str(seqObj.seq[loc.start:loc.end].reverse_complement()))

        seqName = seqObj.id
        seqFeatType = seqFeat.type
        locusTagName = ""
        ProductName = ""
        GeneName = ""
        featSeqNuc = ''.join(dnaSeq)
        featSeqAA = ""
        featSeqLenNuc = len(featSeqNuc)
        featSeqLenAA = len(featSeqAA) 

        if translateCodingSeq==True and seqFeatType=="cds":
            featSeqAA = Seq.Seq(''.join(dnaSeq)).translate()
        else:
            featSeqAA = ""

        if 'locus_tag' in seqFeat.qualifiers.keys():
            locusTagName = seqFeat.qualifiers['locus_tag'][0]
        else:
            locusTagName = "Unknown"

        if 'gene' in seqFeat.qualifiers.keys():
            GeneName = seqFeat.qualifiers['gene'][0]
        else:
            GeneName = "Unknown"

        if 'product' in seqFeat.qualifiers.keys():
            ProductName = seqFeat.qualifiers['product'][0]
        else:
            ProductName = "Unknown"

        if outputFormat=="fasta":
            print(f">{seqName}___{locusTagName}___{GeneName} {ProductName}\n{featSeqNuc}")
        elif outputFormat=="csv":
            print(f"\"{seqName}\",\"{seqFeatType}\",\"{featSeqLenNuc}\",\"{locusTagName}\",\"{GeneName}\",\"{ProductName}\",{featSeqNuc}")
        elif outputFormat=="tsv":                        
            print(f"\"{seqName}\"\t\"{seqFeatType}\"\t\"{featSeqLenNuc}\"\t\"{locusTagName}\"\t\"{GeneName}\"\t\"{ProductName}\"\t{featSeqNuc}")
        else:
            pass

if __name__=="__main__":
    main()
