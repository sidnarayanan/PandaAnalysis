#!/usr/bin/env python

from re import sub
from PandaCore.Tools.script import * 
from PandaCore.Tools.job_config import DataSample,convert_catalog

workdir = getenv('SUBMIT_WORKDIR')

args = parse('--infile', '--outfile', '--nfiles') 

fin = open(args.infile)
samples = convert_catalog(list(fin),as_dict=True)

fout = open(args.outfile,'w')
keys = sorted(samples)
if int(environ.get('SUBMIT_REVERSE', 0)):
	keys.reverse()
to_write = []
for k in keys:
	sample = samples[k]
	to_write += sample.get_config(args.nfiles,suffix='_%i')

counter=0
for c in to_write:
	fout.write(c%(counter,counter))
	counter += 1

logger.info('configBuilder.py','Submission will have %i jobs'%(counter))

fout.close()
