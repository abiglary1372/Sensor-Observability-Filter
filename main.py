###############################################################################
#######################OBSERVABILITY OF A SENSOR SCHDULE#######################
###############################################################################
"""
TABLE OF PARAMETERS:
    
FUNCTION main():            A function where the correct sequence of function  
                            execution is managed.

1.FUNCTION flt():           Almost exactly the execution of the psudocode of 
                            the thesis document

2.FUNCTION cover():         finds the cover set of vector with respect to a 
                            basis

3.FUNCTION max_beta():      does exactly what function MAX() does in the psudo-
                            code of the thesis document, finding the paired 
                            element with fisrt elemnt of a pair that has the 
                            largest cardinality among other paired elements

4.FUNCTION union():         in this program at some point we will have a list of
                            sub cover sets (BM) in the form of lists where each
                            list contains basis elements in the form of lists
                            thus we have set of basis elements sets. this 
                            function recieves this set of sets and the takes 
                            union of the sets inside and outputs a list of 
                            basis elements(vectors)

5.FUNCTION conj():          This function finds the conjunction of two sets

6.FUNCTION remove():        This function accepts two Lists and removes the common 
                            elements of two lists and from the first one

7.FUNCTION obs_check():     After finding sensor set So this function builds 
                            the proposed schedule of the thesis and builds and 
                            claculates the  rank of the observability matrix
                           

8.FUNCTION snsr_fndr():     This Function is used to create a desired initial 
                            sensor set with a specific basis coverage  
"""

"""
TABLE OF VARIABLES:
    
    

"""

import math
from numpy import linalg
from numpy.linalg import matrix_rank
import numpy as np

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##example 3.3.1 
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

#covering alpha
Alpha =[[1,2,3,0,0],
        [0,0,4,8,0],
        [0,0,0,0,3],
        [0,0,0,7,0],
        [0,0,0,0,6],
        [0,3,4,0,0]]

#not covering alpha

# Alpha =[[0,2,0,0,0],
#         [0,0,0,8,0],
#         [0,0,0,3,3],
#         [0,0,0,7,0],
#         [0,0,0,0,6],
#         [0,3,0,0,0]]

## example 3.2

# Atrps=[[2,3,4],
#        [0,8,9],
#        [0,0,5]]

# G = [[1,0,0],
#      [-5/3,-3,1],
#      [0.5,1,0]]
#not covering
# Alpha =[[0,2,0],
#         [1,0,0],
#         [2,4,0]]
#covering
# Alpha =[[0,2,0],
#         [1,0,0],
#         [0,4,6]]

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

def main():
    
   Sr = snsr_fndr(G,Alpha)
   So, No, BsCvrsM, GammaRM,VctrsLM ,z = flt(Sr,G)
   rnk = obs_check(So,No,Atrps,G,z)
   if rnk == len(Atrps):
       print("********************************************************************")
       print("######  Observability matrix is full rank and the observable schedule exists!!  #####")
       print("********************************************************************")
       
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
            else:
                Cvr = cover([VctrsL[l]],Bs)
                Cmn = conj(union(BM),Cvr[0])
                BetaR = remove(Cvr[0],Cmn)
                Br.append([BetaR,VctrsL[l]])   
            l=l + 1      
            
        Bo = max_beta(Br)
        So.append(Bo[1])
        No.append(len(Bo[0]))
        BM.append(Bo[0])
        GammaR = remove(GammaR,Bo[0])
        VctrsL = remove(VctrsL,[Bo[1]])
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
    
    #copying the list values to a new variables
    i=0
    while i < len(S):
        SCpy.append(S[i])
        i=i+1
    #########
    
    SU = []
    l=0
    i=0
    # removing similar set elements
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
    
    #removing similar elements 
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












