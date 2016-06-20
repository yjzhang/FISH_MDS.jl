#!/usr/bin/env python
import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description="""-h or --help for help page
This script is a utility to re-bin bedGraph format signal files,
such that they may be displayed in a 3D genome browser.  The re-binning
uses the bins from a Hi-C experiment and computes the weighted mean
of overlapping bedGraph signal.  In this case the weight is the number of base
pair overlap between the 'A' bedGraph feature and 'B' bin feature divided
by the total bases overlap of all 'A' features with a specific 'B' bin feature.
Therefore 'A' features with more overlap with the bin in 'B' will carry more
weight in resultant re-binned signal. The script requires bedtools and 
python module pybedtools to work properly.

Questions, Bug reporting: email steven.criscione@gmail.com

""",formatter_class=RawTextHelpFormatter)
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('input_signal', action= 'store', metavar='input_signal',
help="""Enter tab seperated bedGraph.  Example: signal.bedGraph or signal.bedGraph.gz
The bedGraph should have no header. The file optionally can be a url
see --url and optionally can be gzipped. The file must have 
chromosome, start, end, signal as the fields 1,2,3,4 with no additional fields.
example:
chr1	1    100 	106.5
chr1	300  400 	51.1
chr1	500  509 	53.2

For bigWig files the script can convert bigWig to bedGraph using the UCSC binary utility
found here: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/
the utility is bigWigToBedGraph this option is implemented by --bigWig 
to True.  The script searches for the directory ./kent_binaries/bigWigToBedGraph
within the current directory. If the input is bigWig it should have only four fields.
The file must have chromosome, start, end, signal with no additional fields.
""")
parser.add_argument('bins', action= 'store', metavar='bins',
help="""Enter tab seperated bed file of Hi-C bins.  Example: bins.bed
The bed should have chromosome, start, end fields as 1-3 and 
should not be gzipped. 
example:
chr1	1    1000
chr1	101  200
chr1	201  300
""")
parser.add_argument('resolution', action= 'store', metavar='resolution',
help="""Enter the resolution of the bins file in base pairs with no comma. 
example: 25000 """)
parser.add_argument('outname', action= 'store', metavar='outname',
help="""Enter the name of output re-binned bedGraph.  Example: rebin_signal.bedGraph
The output file is similar to a bedGraph, but it has a header. The columns are named:
chr,start,end,signal,log2.signal,rescaled.signal
Signal is the weighted average for the bin of input_signal.  log2.signal is the log2
transform of the signal.  Rescaled.signal is the signal rescaled to be between 0 and 1.
""")
parser.add_argument('--url', action= 'store', dest='url', metavar='url', default=False , type=bool,
help='Is the input bedGraph a url and not a local file? Default is  False change to True') 
parser.add_argument('--bigWig', action= 'store', dest='bigWig', metavar='bigWig', default=False , type=bool,
help='Is the input bigWig? Default is  False change to True') 
args = parser.parse_args()
####################################################################################################
import pybedtools
from pybedtools import BedTool
import subprocess
from subprocess import Popen, PIPE, STDOUT
import numpy as np
import os
import shlex
import scipy.stats

input_signal=args.input_signal
bins=args.bins
outname=args.outname
url=args.url
bigWig=args.bigWig
resolution=float(args.resolution)
####################################################################################################
def test_load(bigWig):
    try:
        subprocess.Popen("bedtools", stdout=PIPE, stderr=PIPE)
    except:
        raise Exception("Error: bedtools not in your $PATH")
    if bigWig==True:
	    try:
	        subprocess.Popen("./bigWigToBedGraph", stdout=PIPE, stderr=PIPE)
	    except:
	        raise Exception("Error: bigWigToBedGraph not in your $PATH")
test_load(bigWig)
####################################################################################################
def check_input():
	if url== True:
		if os.path.isfile(os.path.basename(input_signal) ) :
			raise Exception("Check the files within this folder, --url True was entered but file with same basename is already in the folder. Possibly remove and re-run")			
	if bigWig == False:
		if input_signal.endswith(".bedGraph"):
			pass
		elif input_signal.endswith(".bedgraph"):
			pass
		else:
			raise Exception("Expecting .bedGraph or .bedgraph suffix for input_signal file")
	if bigWig == True:
		if input_signal.endswith(".bigwig"):
			pass
		elif input_signal.endswith(".bigWig"):
			pass
		else:
			raise Exception("Expecting .bigwig or .bigWig suffix for input_signal file")
	if bins.endswith(".bed"):
		pass
	else:
		raise Exception("Expecting .bed suffix for bin file")    	
check_input()
####################################################################################################
def run_bash(command):
    p=subprocess.Popen(command,shell=True)
    p.wait()

def rebin_step1():
	if url == False:
		A=BedTool(input_signal)
		B=BedTool(bins)
		AB = A.intersect(B, wo=True)
		AB_inv=B.intersect(A, v=True)
		return AB,AB_inv
	elif url == True:
		to_download="'"+input_signal+"'"
		command= "wget "+to_download
		p = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
		stdout, stderr = p.communicate()
		ninput_signal=os.path.basename(input_signal) 
		if bigWig ==True:
			run_bash("./kent_binaries/bigWigToBedGraph "+ninput_signal+" "+ninput_signal.replace(".bigWig",".bedGraph").replace(".bigwig",".bedGraph"))
			run_bash("rm "+ninput_signal)
			ninput_signal=ninput_signal.replace(".bigWig",".bedGraph").replace(".bigwig",".bedGraph")	
		A=BedTool(ninput_signal)
		B=BedTool(bins)
		AB = A.intersect(B, wo=True)
		AB_inv=B.intersect(A, v=True)
		run_bash("rm "+ninput_signal)
		return AB,AB_inv

AB,AB_inv=rebin_step1()
####################################################################################################
def rebin_step2(AB,AB_inv):
	newbins={}
	newbinsigs={}
	newbins_inv={}
	for line in AB:
		line=str(line)
		line=line.strip("\n").split("\t")
		if len(line)!=8:
			raise Exception("bins should have 3 fields and bedGraph should have 4 fields")
		if not newbins.has_key(",".join(line[4:7])):
			newbins[",".join(line[4:7])]=[]
			newbinsigs[",".join(line[4:7])]=[]
		newbins[",".join(line[4:7])].append(float(line[7]))
		newbinsigs[",".join(line[4:7])].append(float(line[3]))
	for line in AB_inv:
		line=str(line)
		line=line.strip("\n").split("\t")
		if len(line)!=3:
			raise Exception("bins should have 3 fields")
		if not newbins_inv.has_key(",".join(line[0:3])):
			newbins_inv[",".join(line[0:3])]=[]
	return newbins,newbinsigs,newbins_inv

newbins,newbinsigs,newbins_inv=rebin_step2(AB,AB_inv)
####################################################################################################
def rebin_step3(newbins,newbinsigs,resolution,outname):
	rebin_array=[]
	for newbin in newbins.keys():
		weighted_avg=0
		for idx, weight in enumerate(newbins[newbin]):
			weighted_avg+=((weight/resolution)*(newbinsigs[newbin][idx]))
		rebin_array.append([newbin.split(",")+[weighted_avg]])
	return np.matrix(np.asarray(rebin_array))

rebin_array=rebin_step3(newbins,newbinsigs,resolution,outname)
nrebin_array= np.c_[rebin_array,(np.log2(np.asarray(rebin_array[:,3],dtype=float)))]

def zero_to_one(array,j):
	narray=np.asarray(array[:,j],dtype=float)
	rescaled=((narray-np.min(narray))/(np.max(narray)-np.min(narray)))
	return rescaled

def winsorising(array,j):
	narray=np.asarray(array[:,j],dtype=float)
	winsor=scipy.stats.mstats.winsorize(narray, limits=0.05)
	log2_winsor=np.log2(winsor)
	return winsor,log2_winsor

rescaled=zero_to_one(rebin_array,3)
nrebin_array= np.c_[nrebin_array,rescaled]
winsor,log2_winsor=winsorising(rebin_array,3)
nrebin_array=np.c_[nrebin_array,winsor]
nrebin_array=np.c_[nrebin_array,log2_winsor]

####################################################################################################
def rebin_step4(nrebin_array,outname,newbins_inv):
	with open(outname,"w") as fout:
		for i in nrebin_array.tolist():
			print >> fout, "\t".join(i)
		for key in newbins_inv.keys():
			print >> fout, "\t".join(key.split(","))+"\t0\t0\t0\t0\t0"
	fout.close()
	#run_bash("cat "+outname+" > check."+outname)
	run_bash("bedtools sort -i "+outname+" > fix."+outname)
	run_bash("mv  fix."+outname+" "+outname)
	with open("header."+outname,"w") as header:
		print >> header, "chr\tstart\tend\tsignal\tlog2.signal\trescaled.signal\ttrimmed.signal\tlog2trimmed.signal"
	header.close()
	run_bash("cat "+"header."+outname+" "+outname+" > fix."+outname)
	run_bash("rm "+"header."+outname)
	run_bash("mv  fix."+outname+" "+outname)
	
rebin_step4(nrebin_array,outname,newbins_inv)
print "Finished conversion to " +outname