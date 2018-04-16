#!/usr/bin/env python

from os import system,environ
from sys import exit,stdout, argv

which = argv[1]

deta = [1.,2.,2.4,2.8,3.2,3.6,4,4.4,4.8]
dphi = [3.14*x/10 for x in range(1,11)] 
mjj = range(100,1700,200)
pt1 = [80+20*x for x in range(5)]
pt2 = [40+10*x for x in range(9)]

datacards = {
        'cnc' : 'cnc_opt.txt',
        'mjj' : 'mjj_opt.txt',
        'deta' : 'deta_opt.txt',
        'dphi' : 'dphi_opt.txt',
        }

cuts = []
if which == 'cnc':
    cuts += ['deta>%.2f&&fabs(dphi)<%.2f&&mjj>%i'%(x,y,z) for x in deta for y in dphi for z in mjj]
elif which == 'mjj':
    cuts += ['deta>%.2f&&fabs(dphi)<%.2f'%(x,y) for x in deta for y in dphi]
elif which == 'deta':
    cuts += ['fabs(dphi)<%.2f&&mjj>%i'%(y,z) for y in dphi for z in mjj]
elif which == 'dphi':
    cuts += ['deta>%.2f&&mjj>%i'%(x,z) for x in deta for z in mjj]

#cuts += ['deta>%.2f&&mjj>%i'%(e,m) for e in deta for m in mjj]
#cuts += ['jot1Pt>%.2f&&jot2Pt>%.2f'%(x,y) for x in pt1 for y in pt2]
#cuts += ['deta>%.2f&&jot1Pt>%.2f'%(x,y) for x in deta for y in pt1]
#cuts += ['deta>%.2f&&jot2Pt>%.2f'%(x,y) for x in deta for y in pt2]
#cuts += ['mjj>%i&&jot1Pt>%.2f'%(x,y) for x in mjj for y in pt1]
#cuts += ['mjj>%i&&jot2Pt>%.2f'%(x,y) for x in mjj for y in pt2]
#cuts += ['deta>%.2f&&fabs(dphi)<%.2f'%(e,p) for e in deta for p in dphi]
#cuts += ['deta>%.2f&&fabs(minJetMetDPhi_withendcap)>%.2f'%(e,p) for e in deta for p in dphi]
#cuts += ['deta>%.2f&&mjj>%i'%(x,z) for x in deta for z in mjj]
#cuts += ['fabs(dphi)<%.2f&&mjj>%i'%(y,z) for y in dphi for z in mjj]
user = environ['USER']
outdir = environ['PANDA_VBFSCAN']+'/'+which+'/'
scramdir = environ['PANDA_FIT']

for cut in cuts:
    print "'%s' %s %s %s"%(cut,scramdir,outdir,datacards[which])
