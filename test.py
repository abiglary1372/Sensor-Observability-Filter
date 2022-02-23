import numpy as np

#Sr = [[3,1,0,0,0],[6,0,0,0,0],[6,0,7,0,0],[6,0,0,4,0],[6,3,0,0,0]]
Sr = [[3,1,0,0,0],[0,0,2,1,0],[0,0,0,0,3]]
Bo= [[[1, 0, 0, 0, 0], [1.3333333333333333, 0.3333333333333333, -0.6666666666666666, 1, 0], [0.6666666666666666, 1, 0, 0, 0], [0.625, 2.5, 1, 0, 0], [2.7666666666666666, 4.066666666666666, 2.6666666666666665, 0, 1]], [1, 1, 1, 1, 1]] 
Br= [[[[1, 0, 0, 0, 0], [1.3333333333333333, 0.3333333333333333, -0.6666666666666666, 1, 0], [0.6666666666666666, 1, 0, 0, 0], [0.625, 2.5, 1, 0, 0], [2.7666666666666666, 4.066666666666666, 2.6666666666666665, 0, 1]], [1, 1, 1, 1, 1]]]
BM= [[[1, 0, 0, 0, 0], [1.3333333333333333, 0.3333333333333333, -0.6666666666666666, 1, 0], [0.6666666666666666, 1, 0, 0, 0], [0.625, 2.5, 1, 0, 0], [2.7666666666666666, 4.066666666666666, 2.6666666666666665, 0, 1]]] 
GammaR = [[1,0,0,0,0],[4/3,1/3,-2/3,1,0],[2/3,1,0,0,0],[5/8,5/2,1,0,0],[83/30,61/15,8/3,0,1]]

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
        
CvrsL, BsCvrs = cover(Sr,GammaR)