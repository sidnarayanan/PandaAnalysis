#!/usr/bin/env python

from glob import glob
from os import stat,getenv,system,path
from PandaCore.Tools.script import * 
from re import sub, match
from sys import argv
import argparse

args = parse(('--catalog',{'type':str, 'default':'/home/cmsprod/catalog/t2mit/pandaf/012'}),
             ('--user_catalog', STORE_TRUE),
             '--mc_catalog',
             '--data_catalog',
             ('--cfg',{'default':'common'}),
             '--outfile',
             ('--include', {'nargs':'+', 'default':None}),
             ('--exclude', {'nargs':'+', 'default':None}),
             ('--smartcache', STORE_TRUE),
             ('--force', STORE_TRUE),
             ('--max_files', {'type':int, 'default':-1}))

args.catalog = '/' + args.catalog.strip('/')
if not args.mc_catalog:
    args.mc_catalog = args.catalog
if not args.data_catalog:
    args.data_catalog = args.catalog
if not args.outfile:
    args.outfile = '/'.join(['/home', getenv('USER'), 'public_html'] + 
                            getenv('PANDA_CFG').replace('http://t3serv001.mit.edu/~','').split('/')[1:])

if args.cfg == 'leptonic':
    from PandaCore.Tools.process_leptonic import *
else:
    from PandaCore.Tools.process import *

user = getenv('USER')

class CatalogSample:
    def __init__(self,name,dtype,xsec):
        self.name = name
        self.dtype = dtype
        self.xsec = xsec
        self.files = []
    def add_file(self,f):
        self.files.append(f)
    def get_lines(self,smartcache_args=None,max_lines=-1):
        lines = []
        nickname = self.name+'_%i'
        for f in self.files:
            ds_ = f.split('/')[-2]
            f_ = f.split('/')[-1]
            book_ = '/'.join(args.catalog.split('/')[-2:])
            lines.append('{0:<25} {2:<10} {3:<15} {1}\n'.format(nickname,f,self.dtype,self.xsec)) 
            if smartcache_args is not None:
                smartcache_args.append(ds_)
        if max_lines > 0:
            lines = lines[:max_lines]
        return lines

def smartcache(arguments):
    book = args.catalog.strip('/').split('/')[-1]
    arguments = ' '.join(arguments)
    for a in [arguments]:
        cmd = ('python2.6 $(which dynamo-request) --panda %s --sample %s'%(book, a))
        # print cmd
        system(cmd)

def checkDS(nickname,include,exclude):
  included=False
  if include:
    for i in include:
      if i in nickname:
        included=True
        break
  else:
    included=True
  excluded=False
  if exclude:
    for e in exclude:
      if e in nickname:
        excluded=True
        break
  else:
    excluded=False
  return (included and not(excluded))

samples = {}

could_not_find = []

def cat(catalog, condition=lambda x : True): 
    global samples, could_not_find
    for d in sorted(glob(catalog+'/*')):
        dirname = d.split('/')[-1]
        if not condition(dirname):
            continue
        shortname = dirname.split('+')[0]
        try:
            properties = processes[shortname]
            found = True
        except KeyError:
            if args.force:
              found = False
              properties = (shortname,'MC',1)
            else:
              continue
        dtype = 'MC' if 'MINIAODSIM' in dirname else 'Data'
        if not checkDS(properties[0],args.include,args.exclude):
            continue
        if not found:
          could_not_find.append(shortname)
        if properties[0] not in samples:
            samples[properties[0]] = CatalogSample(*properties)
        logger.info(argv[0], 'Selecting %s'%properties[0])
        sample = samples[properties[0]]
        for rfpath in glob(d+'/RawFiles.*'):
            rawfile = open(rfpath)
            for line in rawfile:
                sample.add_file(line.split()[0])

cat(args.mc_catalog, lambda x : bool(match('.*SIM$', x)))
cat(args.data_catalog, lambda x : bool(match('.*AOD$', x)))
if args.user_catalog:
    cat(args.catalog.replace('cmsprod',user))
    cat(args.catalog.replace('cmsprod',user).replace('t2','t3'))

if len(could_not_find)>0:
    logger.warning(argv[0],"Could not properly catalog following datasets (force=%s)"%('True' if args.force else 'False'))
    for c in could_not_find:
        logger.warning(argv[0],'\t'+c)

cfg_file = open(args.outfile,'w')
lines = []
smartcache_args = [] if args.smartcache else None
for k in sorted(samples):
    sample = samples[k]
    lines += sample.get_lines(smartcache_args, args.max_files)
for iL in xrange(len(lines)):
    cfg_file.write(lines[iL]%(iL))

cfg_file.close()
logger.info(argv[0],'Cataloged %i files for %i datasets'%(len(lines),len(samples)))
logger.info(argv[0],'Output written to '+args.outfile)

if args.smartcache:
    smartcache_datasets = list(set(smartcache_args))
    logger.info(argv[0],'Making smartcache requests for files')
    smartcache(smartcache_datasets)
