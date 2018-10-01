#!/usr/bin/env python

from re import sub 

variables = []

formulae = [ (x,'F',sub(r'_([0-9][0-9])', r'_\1[0]', x.replace('ecfN','fjECFN'))) for x in 
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
formulae.append(['tau32SD','F','fjTau32SD[0]'])
formulae.append(['htt_frec','F','fjHTTFRec[0]'])

spectators = [
          ('pt','F'),
          ('mSD','F'),
              ]



if __name__ == '__main__':
    from pprint import pprint
    pprint(variables)
    pprint(formulae)
    pprint(spectators)
