#!/usr/bin/env python

variables = [
            ]

formulae = []

spectators = [
          ('tau32','F'),
          ('tau32SD','F'),
          ('pt','F'),
          ('mSD','F'),
              ]

ecfnstring = 'ecfN_%i_%i_%.2i'
for N in [2,3,4]:
  for order in [1,2,3]:
    if N==4 and order==3:
      continue
    if N==2 and order!=1:
      continue
    for beta in [5,10,20]:
      variables.append((ecfnstring%(order,N,beta),'F'))

