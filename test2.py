#test2 
import numpy as np

G=[[1,0,0,0,0,0,0,0,0,0],
    [0,1,0,0,0,0,0,0,0,0],
    [0,0,1,0,0,0,0,0,0,0],
    [0,0,-0.5547,0.8321,0,0,0,0,0,0],
    [0.3201,0.5976,0.3320,0.4980,0.4268,0,0,0,0,0],
    [0,0,0.8602,0.3982,0,0.3186,0,0,0,0],
    [0.3243,0.5202,0.1788,0.5364,0.5515,0,0.0215,0,0,0],
    [0.1834,0,0.8359,0.3203,0,0.3669,0,0.1747,0,0],
    [0.2100,0.0321,0.8309,0.3231,0.0634,0.3497,0,0.1832,0.0083,0],
    [0.2318,0.0608,0.8217,0.3259,0.1245,0.3310,0.0015,0.1891,0.0164,0.0007]]
a1 = [0,0,0,0,0,0,0,0,2,0]

InvBs = np.linalg.inv(np.array(G))

Sv1 = np.dot(np.array(a1),np.array(G))

a1nw = np.dot(Sv1,InvBs)

i=0
while i<len(a1nw):
    #compensating numerical errors
    if abs(a1nw[i]) < 1e-10:
        a1nw[i]=0
    i=i+1
    
A=[[0,0,0,0,0,0,0,0,2,0],[0,1,0,0,0,3,0,0,2,0],[0,0,0,4,0,0,0,0,2,0]]

def snsr_fndr(G,Alpha):
    
    #this function is used to design a sensor set for our numerical example
    Sv = []
    l=0
    while l < len(Alpha):
        s =np.dot(np.array(Alpha[l]),np.array(G))
        s.tolist()
        
        Sv.append(s)
        l=l+1
        
    return Sv


Set = snsr_fndr(G,A)

