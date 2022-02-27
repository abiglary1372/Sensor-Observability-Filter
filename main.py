# what we want to do now is to designe the sensor set and choosing the coverset for each sensor to do so we eed to do numerical 
#sims and make suyre there are no numerical errors 

import math
from numpy import linalg
from numpy.linalg import matrix_rank
import numpy as np
#A here is the transpose of the system matrix
#Sr = [[3,1,0,0,0],[6,0,0,0,0],[6,0,7,0,0],[6,0,0,4,0],[6,3,0,0,0]] # does not cover
#Sr = [[0,1,0,0,0],[6,0,0,0,0],[0,0,2,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,10],[5,0,0,0,0]] # covers
#Sr = [[1,1,1,1,1]]
#Sr = [[3,1,0,0,0],[0,0,2,1,0],[0,0,0,0,3]] #covers
#Atrps = [[1,2,3,4,6],[0,4,5,3,7],[0,0,6,2,8],[0,0,0,3,0],[0,0,0,0,9]]
#Sr = [[3,1,0,0,0],[6,5,7,0,0],[6,4,7,0,0],[6,0,0,4,0],[6,3,0,0,5]] #covers
# Sr = [[1,0,1,0,0,4,0,0,0,0],
#        [1,0,0,0,0,0,1,0,0,0],
#        [0,0,1,0,0,0,3,0,0,0],
#        [0,1,0,4,0,0,0,0,0,0],
#        [0,0,0,0,5,0,0,2,0,0],
#        [1,0,0,0,0,6,0,0,0,0],
#        [0,0,0,0,0,0,7,0,0,0],
#        [0,2,0,0,0,0,0,8,0,0],
#        [0,0,0,2,0,0,1,0,0,0],
#        [1,0,17,0,0,1,0,0,0,2]]
# Sr = [[0.4201,0.0641,1.6619,0.6463,0.1268,0.6995,0,0.3664,0.0167,0]]
#G = [[1,0,0,0,0],[4/3,1/3,-2/3,1,0],[2/3,1,0,0,0],[8/5,5/2,1,0,0],[83/30,61/15,8/3,0,1]]
#G = [[1,0,0,0,0],[0.5547,0.8321,0,0,0],[0.5108,0.7982,0.3193,0,0],[0.7303,0.1826,-0.3651,0.5477,0],[0.4868,0.7155,0.4692,0,0.1759]]

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##example 1 (positive eigen values)
###########
# Atrps = [[1,0,0,0,6,0,12,21,22,0]
#           ,[0,4,0,0,7,0,14,0,16,0]
#           ,[0,0,6,2,0,11,0,45,67,89]
#           ,[0,0,0,3,7,10,45,12,0,0]
#           ,[0,0,0,0,9,0,77,0,99,0]
#           ,[0,0,0,0,0,11,0,21,0,0]
#           ,[0,0,0,0,0,0,12,0,0,23]
#           ,[0,0,0,0,0,0,0,21,22,23]
#           ,[0,0,0,0,0,0,0,0,22,23]
#           ,[0,0,0,0,0,0,0,0,0,23]]

# # the basis made of eigen vectors
# G=[[1,0,0,0,0,0,0,0,0,0],
#     [0,1,0,0,0,0,0,0,0,0],
#     [0,0,1,0,0,0,0,0,0,0],
#     [0,0,-0.5547,0.8321,0,0,0,0,0,0],
#     [0.3201,0.5976,0.3320,0.4980,0.4268,0,0,0,0,0],
#     [0,0,0.8602,0.3982,0,0.3186,0,0,0,0],
#     [0.3243,0.5202,0.1788,0.5364,0.5515,0,0.0215,0,0,0],
#     [0.1834,0,0.8359,0.3203,0,0.3669,0,0.1747,0,0],
#     [0.2100,0.0321,0.8309,0.3231,0.0634,0.3497,0,0.1832,0.0083,0],
#     [0.2318,0.0608,0.8217,0.3259,0.1245,0.3310,0.0015,0.1891,0.0164,0.0007]]

# # desining the sensor set coverage 
# Alpha = [[2,3,4,0,0,0,0,0,0,0],
#           [0,0,0,5,6,3,7,0,0,0],
#           [0,0,0,0,0,0,0,1,2,7],
#           [1,0,0,0,0,0,0,0,0,0],
#           [0,0,0,0,0,0,0,0,0,1],
#           [0,0,0,0,0,2,0,0,0,0]]
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##example 2 (one negative eigen value)
###########
Atrps=[[1,2,4,5,7],
        [0,11,4,8,7],
        [0,0,-8,0,3],
        [0,0,0,7,0],
        [0,0,0,0,6]]
G = [[1,0,0,0,0],
      [-68/171,-4/19,1,0,0],
      [33/35,-11/7,3/14,0,1],
      [1/6,-2,0,1,0],
      [1/5,1,0,0,0]]

# covering alpha
# Alpha =[[1,2,3,0,0],
#         [0,0,4,8,0],
#         [0,0,0,0,3],
#         [0,0,0,7,0],
#         [0,0,0,0,6],
#         [0,3,4,0,0]]

#not covering alpha

Alpha =[[0,2,0,0,0],
        [0,0,0,8,0],
        [0,0,0,3,3],
        [0,0,0,7,0],
        [0,0,0,0,6],
        [0,3,0,0,0]]

def main():
    
   Sr = snsr_fndr(G,Alpha)
   print(Sr)
   So, No, BsCvrsM, GammaRM,VctrsLM ,z = flt(Sr,G)
   rnk = obs_check(So,No,Atrps,G,z)
   if rnk == len(Atrps):
       print("observability matrix is full rank and the schedule exists")
   return So, Sr, No, BsCvrsM, GammaRM,VctrsLM ,z ,rnk
   
    
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
    Cvr = cover(VctrsL,GammaR)
    
    while  (GammaR != [] and VctrsL != []):
        l=0
        while l < len(VctrsL):
            if BM == []:
                BetaR = cover([VctrsL[l]],Bs)
                Br.append([BetaR[0],VctrsL[l]])
                # print("******************************")
                # print("Br first\n",Br)
            else:
                Cvr = cover([VctrsL[l]],Bs)
                Cmn = conj(union(BM),Cvr[0])
                BetaR = remove(Cvr[0],Cmn)
                Br.append([BetaR,VctrsL[l]])   
                # print("******************************")
                # print("Br\n",Br)
                
            l=l + 1      
        Bo = max_beta(Br)
        So.append(Bo[1])
        No.append(len(Bo[0]))
        BM.append(Bo[0])
        GammaR = remove(GammaR,Bo[0])
        VctrsL = remove(VctrsL,[Bo[1]])
        # print("******************************")
        # print('this is Bo \n',Bo )
        # print("******************************")
        # print('this is VctrsL\n',VctrsL )
        Zeta = Zeta + len(Bo[0]) 
        z.append(Zeta) 
        Br = []
        Bo = []
        m = m + 1
    
    if GammaR==[] or (GammaR==[] and VctrsL==[]):
        BsCvrsF = True
    elif VctrsL==[]:
        BsCvrsF = False
        
    return So, No,BsCvrsF, GammaR,VctrsL, z

def cover(Vctrs,Bs):
    # this funvtion should check if  the sensor set is covering the basis first then 
    # outputs all the cover sets in a data structure for the res of the algorithm to use
    # if the sensor set os not covering a message should be returned stating this fact
    i=0
    Beta = []
    CvrsL = []
    InvBs = np.linalg.inv(np.array(Bs))
    l=0
    while l < len(Vctrs):
        Beta = []
        Alpha = np.dot(np.array(Vctrs[l]),InvBs)
        Alpha = Alpha.tolist()

#compensating numerical errors
        i=0
        while i<len(Alpha):

            if abs(Alpha[i]) < 1e-10:
                Alpha[i]=0
            i=i+1
            
        #print("this is alpha\n",Alpha)
        
        
        
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
                
    return CvrsL
        
def max_beta(Br):
    maxBta = [Br[0][0],Br[0][1]]
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
    
def obs_check(S,No,A,Bs,z):
    
    l=0
    O=[]
    rnk = 0
    c=0
    while l < len(S):       
        while c<z[l]:
            o= np.dot(np.array(S[l]), np.linalg.matrix_power(np.transpose( np.array(A)), c) )
            o = o.tolist()
            O.append(o)
            c=c+1
        l=l+1
    rnk = matrix_rank(np.array(O))   
    #print(O)
    return rnk

def snsr_fndr(G,Alpha):
    
    #this function is used to design a sensor set for our numerical example
    Sv = []
    l=0
    while l < len(Alpha):
        s =np.dot(np.array(Alpha[l]),np.array(G))
        s=s.tolist()
        
        Sv.append(s)
        l=l+1
        
    return Sv
        
    
def schdl_rnk_chk(sgm,A,So):
    l=0
    O=[]
    while l<len(sgm):
        o= np.dot(np.array( So[sgm[l]] ), np.linalg.matrix_power(np.transpose( np.array(A)), l))
        o = o.tolist()
        O.append(o)
        l=l+1
        rnk = matrix_rank(np.array(O))
        
    return rnk
        
if __name__ == "__main__":
   So, Sr, No, BsCvrsO, GammaRO,VctrsLO ,z ,rnk = main()

# numerical results indicate that if sensor set dosnt cover the the basis there exist no observale schedule 
# numerical results also indicate that there are could multiplel possible sub sets for the observable schedule 
# numerical results so far suggest that it is possible to change the time steps 
# numerical results so far suggest that no matther how the sensors of the subset are used in the schedule the observabiliy is maintatined 












