from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
    'mu':'(trigger&8)!=0',
}

metFilter='metFilter==1'

cuts = {
        'signal' : 'nJot>3 || (fjPt>400 && nJot>2)',
}

weights = {
        'signal' : '%f*normalizedWeight'
        }

