#!/usr/bin/env python

variables = []

formulae = [ (x.replace('ecfN','clf_ECFN'),'F', x) for x in 
    [   
            'ecfN_1_2_20/pow(ecfN_1_2_10,2.00)',
            'ecfN_1_3_40/ecfN_2_3_20',
            'ecfN_3_3_10/pow(ecfN_1_3_40,0.75)',
            'ecfN_3_3_10/pow(ecfN_2_3_20,0.75)',
            'ecfN_3_3_20/pow(ecfN_3_3_40,0.50)',
            'ecfN_1_4_20/pow(ecfN_1_3_10,2.00)',
            'ecfN_1_4_40/pow(ecfN_1_3_20,2.00)',
            'ecfN_2_4_05/pow(ecfN_1_3_05,2.00)',
            'ecfN_2_4_10/pow(ecfN_1_3_10,2.00)',
            'ecfN_2_4_10/pow(ecfN_2_3_05,2.00)',
            'ecfN_2_4_20/pow(ecfN_1_3_20,2.00)',
    ]
]
formulae.append(['clf_Tau32SD','F','tau32SD'])
formulae.append(['clf_HTTFRec','F','htt_frec'])

spectators = [
          ('gen_pt','F'),
          ('clf_MSD','F'),
              ]


