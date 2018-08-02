#!/usr/bin/env python

class ECFN(object):
  def __init__(self,beta_,N_,i_):
    self.beta=beta_
    self.N=N_
    self.i=i_
    self.power=beta_*i_
  def __eq__(self,other):
    return ( (self.beta==other.beta) and (self.N==other.N) and (self.i==other.i) )

def makeECFString(e,power=1):
  base='ecfN_%i_%i_%.2i'%(e.i,e.N,int(e.beta*10))
  if power==1:
    return base
  else:
    return 'pow(%s,%.2f)'%(base,power)

ecfns = []
#for N in [2,3,4]:
for N in [3,4]:
  for order in [1,2,3]:
    if N==4 and order==3:
      continue
    if N==2 and order!=1:
      continue
    for beta in [.5,1.,2.,4.]:
      ecfns.append(ECFN(beta,N,order))

ratios = []
for numE in ecfns:
  numString=makeECFString(numE)
  for denE in ecfns:
#    if denE.N>numE.N or denE==numE:
    if denE.N>=numE.N or denE==numE:
      continue
    power = numE.power/denE.power
    r=numString+'/'
    if abs(power-1)<.01:
      r += makeECFString(denE)
    else:
      r += makeECFString(denE,power)
    ratios.append(r)

variables = [
            ]

spectators = [
#          ('tau32','F'),
#          ('tau32SD','F'),
          ('pt','F'),
          ('mSD','F'),
             ]

formulae = []

for r in ratios:
  if r=='ecfN_3_3_40/pow(ecfN_1_3_05,24.00)' or r=='ecfN_2_4_40/pow(ecfN_1_4_05,16.00)' or r=='ecfN_2_3_40/pow(ecfN_1_3_05,16.00)':
    continue
  formulae.append((r,'F'))
  if __name__ == "__main__":
    print r



