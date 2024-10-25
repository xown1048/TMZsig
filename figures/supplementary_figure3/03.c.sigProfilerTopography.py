#####	SigProfilerTopography
#####	https://github.com/alexandrovlab/SigProfilerTopography
#####	https://osf.io/5unby/wiki/7.%20Output%20-%20Strand%20Asymmetry/

import sys
import os
import glob
import argparse
import pandas as pd
import re
import subprocess

'''
python 01.sigProfilerTopography.py \
-i /BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/sigProfilerTopography/part4/rev7.tmz \
-o /BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/sigProfilerTopography/part4/rev7.tmz/output
'''

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputdir',
					help='directory to input files',
					required=True,
					nargs='?',
					type=str,
					metavar='/path/to/input_directory',
					default='ExtractSample_input')

parser.add_argument('-o', '--outputdir',
					help='directory to store output files',
					required=True,
					nargs='?',
					type=str,
					metavar='/path/to/output_directory',
					default='ExtractSample_output')

args = parser.parse_args()
dir_input = os.path.realpath(args.inputdir)
dir_output = os.path.realpath(args.outputdir)

from SigProfilerTopography import Topography as topography

genome = 'GRCh38'
jobname = 'result'
numofSimulations = 10
replication_time_signal_file = '/BiO/Live/xown1048/.pyenv/versions/3.8.7/lib/python3.8/site-packages/SigProfilerTopography/lib/replication/wgEncodeUwRepliSeqK562WaveSignalRep1.GRCh38.wig'

def main_function():
	topography.runAnalyses(genome,
						dir_input,
						dir_output,
						jobname,
						numofSimulations,
						replication_time_signal_file=replication_time_signal_file,
						epigenomics=True,
						nucleosome=True,
						replication_time=True,
						strand_bias=True,
						processivity=True,
						discreet_mode=False,
						delete_unnecessary_files=False)

if __name__ == '__main__':
	main_function()


