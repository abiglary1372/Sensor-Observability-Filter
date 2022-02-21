# testing 
import numpy as np

A = [[1,2,3,4,6],[0,4,5,3,7],[0,0,6,2,8],[0,0,0,3,0],[0,0,0,0,9]]
G = [[1,0,0,0,0],[4/3,1/3,-2/3,1,0],[2/3,1,0,0,0],[5/8,5/2,1,0,0],[83/30,61/15,8/3,0,1]]
Sr = [[3,1,0,0,0],[6,5,7,0,0],[6,4,7,0,0],[6,0,0,4,0],[6,3,0,0,5]]
#Sr = [[1,0,0,0,0]]
#Sr = [[1,1,1,3,1]]

def main():
    R =cover(Sr,G)
    return  R
    
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
            if abs(Alpha[i])<1e-15:
                Alpha[i]=0
            i=i+1
        i=0
        while i<len(Alpha):
            if Alpha[i] != 0 :
                Beta.append(Bs[i])
            i=i+1
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
        
    return CvrsL


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

a =main()

#T=union(U)


# union not working 
# CoverU is suspicious
