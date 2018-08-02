#!/usr/bin/env python

variables = []

formulae = [ (x,'F') for x in 
    [   
#            'ecfN_1_2_05/pow(ecfN_1_2_40,0.12)',
#            'ecfN_1_2_10/pow(ecfN_1_2_20,0.50)',
#            'ecfN_1_2_10/pow(ecfN_1_2_40,0.25)',
            'ecfN_1_2_20/pow(ecfN_1_2_10,2.00)',
#            'ecfN_1_2_20/pow(ecfN_1_2_40,0.50)',
#            'ecfN_1_2_40/pow(ecfN_1_2_20,2.00)',
#            'ecfN_1_3_20/pow(ecfN_1_2_05,4.00)',
#            'ecfN_1_3_20/ecfN_2_3_10',
#            'ecfN_1_3_20/pow(ecfN_3_3_05,1.33)',
#            'ecfN_1_3_40/pow(ecfN_1_2_20,2.00)',
#            'ecfN_1_3_40/pow(ecfN_2_3_10,2.00)',
            'ecfN_1_3_40/ecfN_2_3_20',
#            'ecfN_1_3_40/pow(ecfN_2_3_40,0.50)',
#            'ecfN_1_3_40/pow(ecfN_3_3_10,1.33)',
#            'ecfN_1_3_40/pow(ecfN_3_3_20,0.67)',
#            'ecfN_2_3_10/pow(ecfN_1_2_05,4.00)',
#            'ecfN_2_3_10/ecfN_1_3_20',
#            'ecfN_2_3_10/pow(ecfN_3_3_05,1.33)',
#            'ecfN_2_3_10/pow(ecfN_3_3_40,0.17)',
#            'ecfN_2_3_20/pow(ecfN_1_2_10,4.00)',
#            'ecfN_2_3_20/ecfN_1_3_40',
#            'ecfN_2_3_40/pow(ecfN_3_3_20,1.33)',
#            'ecfN_3_3_05/pow(ecfN_1_2_05,3.00)',
#            'ecfN_3_3_05/pow(ecfN_1_3_20,0.75)',
#            'ecfN_3_3_10/pow(ecfN_1_2_10,3.00)',
            'ecfN_3_3_10/pow(ecfN_1_3_40,0.75)',
            'ecfN_3_3_10/pow(ecfN_2_3_20,0.75)',
#            'ecfN_3_3_20/pow(ecfN_2_3_40,0.75)',
            'ecfN_3_3_20/pow(ecfN_3_3_40,0.50)',
#            'ecfN_3_3_40/pow(ecfN_1_2_40,3.00)',
            'ecfN_1_4_20/pow(ecfN_1_3_10,2.00)',
#            'ecfN_1_4_20/pow(ecfN_2_3_05,2.00)',
            'ecfN_1_4_40/pow(ecfN_1_3_20,2.00)',
            'ecfN_2_4_05/pow(ecfN_1_3_05,2.00)',
            'ecfN_2_4_10/pow(ecfN_1_3_10,2.00)',
            'ecfN_2_4_10/pow(ecfN_2_3_05,2.00)',
#            'ecfN_2_4_20/pow(ecfN_1_2_05,8.00)',
            'ecfN_2_4_20/pow(ecfN_1_3_20,2.00)',
#            'ecfN_2_4_20/pow(ecfN_3_3_05,2.67)',
#            'ecfN_2_4_20/pow(ecfN_2_4_40,0.50)',
#             'tau32SD',
    ]
]

spectators = [
#          ('tau32','F'),
#          ('tau32SD','F'),
          ('pt','F'),
          ('mSD','F'),
              ]


