#test2 
import numpy as np


def max_beta(Br):
    maxBta = [Br[0][0],Br[0][1]]
    i=0
    while i<len(Br):
        if len(maxBta[0]) < len(Br[i][0]):
            maxBta = Br[i]
        i=i+1
    return maxBta

Br = [[[], [3, 1, 0, 0, 0]], [[1,2,3,4], [6, 0, 0, 0, 0]], [[], [6, 0, 7, 0, 0]], [[], [6, 3, 0, 0, 0]]]

H = max_beta(Br)