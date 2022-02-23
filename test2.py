#test2 
import numpy as np


def max_beta(Br):
    maxBta = [[],[]]
    i=0
    while i<len(Br):
        if len(maxBta[0]) < len(Br[i][0]):
            maxBta = Br[i]
        i=i+1
    return maxBta

Br = [   [[1],[]],[[],[]],[[],[1]]         ]

H = max_beta(Br)