# list of variable names
#
#

import pandas as pd
import numpy as np

A = [[1,2,3,4,6],[0,4,5,3,7],[0,0,6,2,8],[0,0,0,3,0],[0,0,0,0,9]]
G = [[1,0,0,0,0],[4/3,1/3,-2/3,1,0],[2/3,1,0,0,0],[5/8,5/2,1,0,0],[83/30,61/15,8/3,0,1]]
Sr = [[3,1,0,0,0],[0,0,2,1,0],[0,0,0,0,3]] #covers
#Sr = [[3,1,0,0,0],[6,5,7,0,0],[6,4,7,0,0],[6,0,0,4,0],[6,3,0,0,5]] #covers
#Sr = [[3,1,0,0,0],[6,0,0,0,0],[6,0,7,0,0],[6,0,0,4,0],[6,3,0,0,0]] # does not cover
#Sr = [[1,1,1,1,1]]

def main():
   So, GammaR, VctrsL, BM2 = flt(Sr,G)
   return So, GammaR, VctrsL, BM2
   
    #if obs_check(ObsvS,A): print("observable schedule detected")
    
def flt(Vctrs,Bs):
    # main filtering algorithm
    GammaR = []
    VctrsL = []
    
    #copying
    i=0
    while i < len(Bs):
        GammaR.append(Bs[i])
        i=i+1
    #copying
    i=0
    while i < len(Vctrs):
        VctrsL.append(Vctrs[i])
        i=i+1  
    i=0
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
    #Cvr, BsCvrs = cover(VctrsL,GammaR)
    #if BsCvrs:
    while  (GammaR != [] and VctrsL != []):#######
        print("this is m",m)
        l=0
        print(VctrsL)
        while l < len(VctrsL):
            if BM == []:
                BetaR, BsCvrs = cover([VctrsL[l]],Bs)
                Br.append([BetaR[0],VctrsL[l]])
            else:
                Cvr, BsCvrs = cover([VctrsL[l]],Bs)
                Cmn = conj(union(BM),Cvr[0])
                print("this is Cmn\n",Cmn)
                print("this is Cvr[0]\n",Cvr[0])
                BetaR = remove(Cvr[0],Cmn)
                Br.append([BetaR,VctrsL[l]])     
                
            l=l + 1      
        print("this is Br\n",Br)
        #Bo = max_beta(Br)
        print('this is Bo\n', max_beta(Br))
        So.append(max_beta(Br)[1])
        No.append(len(max_beta(Br)[0]))
        BM.append(max_beta(Br)[0])
        GammaR = remove(GammaR,max_beta(Br)[0])
        print("*******************************")
        print("GammaR \n",GammaR)
        VctrsL = remove(VctrsL,[max_beta(Br)[1]])
        print("*******************************")
        print("VctrsL \n",VctrsL)
        Zeta = Zeta + len(max_beta(Br)[0]) 
        z.append(Zeta) 
        Bo.clear()
        Br = []
        m = m + 1
    return So, GammaR, VctrsL, BM
# remeber in the first iteration a sensor is removed from VctrsL but then reapears in the values in the next operation
# check why is that happening and if there is a problem with VctrsL

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
        Beta = []
        Alpha = np.dot(np.array(Vctrs[l]),InvBs)
        Alpha = Alpha.tolist()
        
        i=0
        while i<len(Alpha):
            #compensating numerical errors
            if abs(Alpha[i]) < 1e-15:
                Alpha[i]=0
            i=i+1

        
        
        i=0
        while i<len(Alpha):
            if Alpha[i] != 0 :
                Beta.append(Bs[i])
            i=i+1
        
        #removeing repeated elements
        i2=0
        l2=0
        while i2<len(Beta):
            l2=0
            while l2+i2+1<len(Beta):
                if Beta[i2] == Beta[l2+i2+1]:
                    Beta.remove(Beta[l2+i2+1])
                else:
                    l2=l2+1
            i2=i2+1
        ##########
        CvrsL.append(Beta)
        l=l+1
        
    #copying CvrsL :: because list are poointers and point to the value even
    #if the variable changes
    CvrsLU =[]
    i=0
    while i < len(CvrsL):
        CvrsLU.append(CvrsL[i])
        i=i+1
    ####
    CvrU = union(CvrsLU)
    l=0
    i=0
    while l < len(Bs):
        i=0
        while i < len(CvrU):
            if Bs[l] == CvrU[i]:
                Chk = Chk+1
            i=i+1
        l=l+1
                
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
    
    SCpy = []
    
    #copying
    i=0
    while i < len(S):
        SCpy.append(S[i])
        i=i+1
    #########
    
    SU = []
    l=0
    i=0
    while i<len(SCpy):
        l=0
        while l+i+1<len(SCpy):
            if SCpy[i] == SCpy[l+i+1]:
                SCpy.remove(SCpy[l+i+1])
            else:
                l=l+1
        i=i+1
    
    l=0
    i=0
    while i<len(SCpy):
        l=0
        while l<len(SCpy[i]):
            SU.append(SCpy[i][l])
            l=l+1
        i=i+1
    i=0
    l=0
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
    SaCpy = []
    SbCpy = []
    
    #copying
    i=0
    while i < len(Sa):
        SaCpy.append(Sa[i])
        i=i+1
    
    i=0
    while i < len(Sb):
        SbCpy.append(Sb[i])
        i=i+1
    #########

    i=0
    l=0
    indx = []
    #detect
    while i<len(SaCpy):
        l=0
        while l<len(SbCpy):
            if SaCpy[i] == SbCpy[l]:
                indx.append(i)
            l=l+1
        i=i+1
  
    i=0
    ToRmv = []
    #save
    while i < len(indx):
        ToRmv.append(SaCpy[indx[i]])
        i=i+1
        
    # remove repeated elements
    i=0
    l=0
    while i<len(ToRmv):
        l=0
        while l+i+1<len(ToRmv):
            if ToRmv[i] == ToRmv[l+i+1]:
                ToRmv.remove(ToRmv[l+i+1])
            else:
                l=l+1
        i=i+1
        
    #remove
    i=0
    while i < len(ToRmv):
        SaCpy.remove(ToRmv[i])
        i=i+1
    return SaCpy
    

# def obs_check(S,A):
    # at the end after filtering the sensor set this function will check the 
    # observability of the schedule using filtered set and returns a true or false value 
    # this function should be able to check any other random schedule as well and 
    # in general it should just check schedules and return true or flase
    
# def random_sch():
    # this function generate random schedules of requested time horizons using the 
    # sensor set or the result of the flt algorithm
  

if __name__ == "__main__":
   So, GammaR, VctrsL, BM2 = main()

## make the functions so that they just accept input and not use global varibles
## do not forget to test the unique ness and change of time step with small dimention systems and small time horizons 
## and discuss the numerical results 










