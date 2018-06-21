d = {
    'JetHT_Run2016B'               : ['JetHT_Run2016B-03Feb2017_ver2-v2'],
    'JetHT_Run2016C'               : ['JetHT_Run2016C-03Feb2017-v1'],
    'JetHT_Run2016D'               : ['JetHT_Run2016D-03Feb2017-v1'],
    'JetHT_Run2016E'               : ['JetHT_Run2016E-03Feb2017-v1'],
    'JetHT_Run2016F'               : ['JetHT_Run2016F-03Feb2017-v1'],
    'JetHT_Run2016H'               : ['JetHT_Run2016H-03Feb2017_ver2-v1','JetHT_Run2016H-03Feb2017_ver3-v1'],

    'MET_Run2016B'               : ['MET_Run2016B-03Feb2017_ver2-v2'],
    'MET_Run2016C'               : ['MET_Run2016C-03Feb2017-v1'],
    'MET_Run2016D'               : ['MET_Run2016D-03Feb2017-v1'],
    'MET_Run2016E'               : ['MET_Run2016E-03Feb2017-v1'],
    'MET_Run2016F'               : ['MET_Run2016F-03Feb2017-v1'],
    'MET_Run2016G'               : ['MET_Run2016G-03Feb2017-v1'],
    'MET_Run2016H'               : ['MET_Run2016H-03Feb2017_ver3-v1'],
}

d['JetHT_AllEras'] = []
d['MET_AllEras'] = []
for e in 'BCDEFGH':
    if e == 'G':
        continue
    for pd in ['JetHT' , 'MET']:
        sample = '%s_Run2016%s'%(pd, e)
        if sample in d:
            d['%s_AllEras'%pd].extend(d[sample])
