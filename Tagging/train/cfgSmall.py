#!/usr/bin/env python

variables = []

formulae = [ (x,'F') for x in 
    [   
      'ecfN_1_4_10/pow(ecfN_1_3_05,2.00)',
      'ecfN_2_4_20/pow(ecfN_1_3_20,2.00)',
      'ecfN_2_4_10/pow(ecfN_1_3_10,2.00)',
      'ecfN_3_3_10/pow(ecfN_1_2_20,1.50)',
      'ecfN_2_4_10/ecfN_1_2_20',
      'ecfN_1_4_20/pow(ecfN_1_3_10,2.00)',
      'ecfN_2_4_05/pow(ecfN_1_3_05,2.00)',
      'ecfN_1_3_10/ecfN_2_3_05',
      'ecfN_3_3_10/pow(ecfN_3_3_20,0.50)',
      'ecfN_3_3_05/pow(ecfN_1_2_05,3.00)',
    ]
]

spectators = [
          ('pt','F'),
          ('mSD','F'),
              ]
