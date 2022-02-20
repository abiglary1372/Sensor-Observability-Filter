# list of variable names
#
#

import pandas as pd
import numpy as np

G = []
Sr = []
A = []

def main():
    ObsvS = flt(Sr,G)
    if obs_check(ObsvS,A): print("observable schedule detected")
    
def flt(Vctrs,Bs):
    # main filtering algorithm
    GammaR = G
    So = []
    No = []
    z = []
    BM = []
    Br = []
    BetaR = []
    Bo = []
    m = 0
    l = 0
    Zeta = 0
    Cvr, BsCvrs = cover(Vctrs,Bs)
    if BsCvrs:
        while  GammaR != [] or Vctrs != []:
            l=0
            while l < Vctrs.size:
                Cvr[l]
                if BM == []:
                    BetaR = Cvr[l]
                    Br.append([BetaR,Vctrs[l]])
                else:
                    Cmn = conj(union(BM),Cvr[l])
                    BetaR = remove(Cvr[l],Cmn)
                    Br.append([BetaR,Vctrs[l]])
                l=l + 1      
             
            Bo = max_beta(Br)
            So.append(Bo[1])
            No.append(len(Bo[0]))
            BM.append(Bo[0])
            GammaR = remove(GammaR,Bo[0])
            Vctrs = remove(Vctrs,Bo[1])
            Zeta = Zeta + len(Bo[0]) 
            z.append(Zeta)
            m = m + 1
        return So

def cover(Vctrs,Bs):
    # this funvtion should check if  the sensor set is covering the basis first then 
    # outputs all the cover sets in a data structure for the res of the algorithm to use
    # if the sensor set os not covering a message should be returned stating this fact
    i=0
    Beta = []
    CvrsL = []
    InvBs = np.linalg.inv(np.array(Bs))
    l=0
    Chk=0
    BsCvrs = False
    while l < len(Vctrs):
        Alpha = np.dot(np.array(Vctrs[l]),InvBs)
        Alpha = Alpha.tolist()
        i=0
        while i<len(Alpha):
            if Alpha[i] != 0 :
                Beta.append(Bs[i])
            i=i+1
        CvrsL.append(Beta)
        l=l+1
    CvrU = union(CvrsL)
    l=0
    i=0
    while l < len(Bs):
        i=0
        while i < CvrU:
            if Bs[l] == CvrU[i]:
                Chk = Chk+1
            i=i+1
                
    if Chk == len(Bs):
        BsCvrs = True
    else:
        BsCvrs = False      
        
    return CvrsL, BsCvrs
        
def max_beta(Br):
    maxBta = [[],[]]
    i=0
    while i<len(Br):
        if len(maxBta[0]) < len(Br[i][0]):
            maxBta = Br[i]
        i=i+1
    return maxBta

# array/set operations 

def union(S):
    #recives a list with elemnts as lists and returns a list of non list elements 
    SU = []
    l=0
    i=0
    while i<len(S):
        l=0
        while l+i+1<len(S):
            if S[i] == S[l+i+1]:
                S.remove(S[l+i+1])
            else:
                l=l+1
        i=i+1
    
    l=0
    i=0
    while i<len(S):
        l=0
        while l<len(S[i]):
            SU.append(S[i][l])
            l=l+1
        i=i+1
    i=0
    l=0
    indx = []
    while i<len(SU):
        l=0
        while l+i+1<len(SU):
            if SU[i] == SU[l+i+1]:
                SU.remove(SU[l+i+1])
            else:
                l=l+1
        i=i+1
    return SU
    
def conj(Sa,Sb):
    Sab = []
    i=0
    l=0
    while i<len(Sa):
        l=0
        while l<len(Sb):
            if Sa[i]==Sb[l]:
                Sab.append(Sa[i])
                
            l=l+1
        i=i+1
    i=0
    l=0
    while i<len(Sab):
        l=0
        while l+i+1<len(Sab):
            if Sab[i] == Sab[l+i+1]:
                Sab.remove(Sab[l+i+1])
            else:
                l=l+1
        i=i+1
    return Sab

def remove(Sa, Sb): 
    ###**tested**###
    ################
    i=0
    l=0
    indx = []
    #detect
    while i<len(Sa):
        l=0
        while l<len(Sb):
            if Sa[i] == Sb[l]:
                indx.append(i)
            l=l+1
        i=i+1
  
    i=0
    ToRmv = []
    #save
    while i < len(indx):
        ToRmv.append(Sa[indx[i]])
        i=i+1
        
    #remove
    i=0
    while i < len(ToRmv):
        Sa.remove(ToRmv[i])
        i=i+1
    return Sa
    

# def obs_check(S,A):
    # at the end after filtering the sensor set this function will check the 
    # observability of the schedule using filtered set and returns a true or false value 
    # this function should be able to check any other random schedule as well and 
    # in general it should just check schedules and return true or flase
    
# def random_sch():
    # this function generate random schedules of requested time horizons using the 
    # sensor set or the result of the flt algorithm
  

if __name__ == "__main__":
    main()

## make the functions so that they just accept input and not use global varibles
## do not forget to test the unique ness and change of time step with small dimention systems and small time horizons 
## and discuss the numerical results 










