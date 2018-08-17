d = {
    'JetHT_Run2016B'               : ['JetHT_Run2016B-03Feb2017_ver2-v2'],
    'JetHT_Run2016C'               : ['JetHT_Run2016C-03Feb2017-v1'],
    'JetHT_Run2016D'               : ['JetHT_Run2016D-03Feb2017-v1'],
    'JetHT_Run2016E'               : ['JetHT_Run2016E-03Feb2017-v1'],
    'JetHT_Run2016F'               : ['JetHT_Run2016F-03Feb2017-v1'],
    'JetHT_Run2016G'               : ['JetHT_Run2016G-03Feb2017-v1'],
    'JetHT_Run2016H'               : ['JetHT_Run2016H-03Feb2017_ver2-v1','JetHT_Run2016H-03Feb2017_ver3-v1'],

    'MET_Run2016B'               : ['MET_Run2016B-03Feb2017_ver2-v2'],
    'MET_Run2016C'               : ['MET_Run2016C-03Feb2017-v1'],
    'MET_Run2016D'               : ['MET_Run2016D-03Feb2017-v1'],
    'MET_Run2016E'               : ['MET_Run2016E-03Feb2017-v1'],
    'MET_Run2016F'               : ['MET_Run2016F-03Feb2017-v1'],
    'MET_Run2016G'               : ['MET_Run2016G-03Feb2017-v1'],
    'MET_Run2016H'               : ['MET_Run2016H-03Feb2017_ver3-v1'],

    'SingleMuon_Run2016B'               : ['SingleMuon_Run2016B-03Feb2017_ver2-v2'],
    'SingleMuon_Run2016C'               : ['SingleMuon_Run2016C-03Feb2017-v1'],
    'SingleMuon_Run2016D'               : ['SingleMuon_Run2016D-03Feb2017-v1'],
    'SingleMuon_Run2016E'               : ['SingleMuon_Run2016E-03Feb2017-v1'],
    'SingleMuon_Run2016F'               : ['SingleMuon_Run2016F-03Feb2017-v1'],
    'SingleMuon_Run2016G'               : ['SingleMuon_Run2016G-03Feb2017-v1'],
    'SingleMuon_Run2016H'               : ['SingleMuon_Run2016H-03Feb2017_ver2-v1', 'SingleMuon_Run2016H-03Feb2017_ver3-v1'],
}

pds = ['JetHT','MET','SingleMuon']
for pd in pds:
    d[pd+'_AllEras'] = []
for e in 'BCDEFGH':
    for pd in pds:
        sample = '%s_Run2016%s'%(pd, e)
        if sample in d:
            d['%s_AllEras'%pd].extend(d[sample])
for iov in ['BCD', 'EF', 'GH']:
    for pd in pds:
        sample = '%s_Run2016%s'%(pd, iov)
        d[sample] = []
        for e in iov:
            d[sample] += d['%s_Run2016%s'%(pd, e)]
