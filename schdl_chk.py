#schedual checker 
from numpy import linalg
from numpy.linalg import matrix_rank
import numpy as np
from random import randint

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##example 1
###########
Atrps = [[1,0,0,0,6,0,12,21,22,0]
          ,[0,4,0,0,7,0,14,0,16,0]
          ,[0,0,6,2,0,11,0,45,67,89]
          ,[0,0,0,3,7,10,45,12,0,0]
          ,[0,0,0,0,9,0,77,0,99,0]
          ,[0,0,0,0,0,11,0,21,0,0]
          ,[0,0,0,0,0,0,12,0,0,23]
          ,[0,0,0,0,0,0,0,21,22,23]
          ,[0,0,0,0,0,0,0,0,22,23]
          ,[0,0,0,0,0,0,0,0,0,23]]
So=[[4.1907, 7.227, 3.0507, 12.0979, 6.4213, 0.9558, 0.1505, 0.0, 0.0, 0.0], 
    [2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [2.226, 0.4898, 8.2496, 3.2478, 0.9983, 3.3833, 0.0105, 1.8648, 0.1314, 0.0049]]

#the schedule list sigma
sgm=[0,0,0,0,1,1,1,2,2,2]

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##example 2
###########
# Atrps=[[1,2,4,5,7],
#        [0,11,4,8,7],
#        [0,0,-8,0,3],
#        [0,0,0,7,0],
#        [0,0,0,0,6]]

# So=[[0.20467836257309946, -0.42105263157894735, 2.0, 0.0, 0.0], 
#     [5.104761904761904, -22.285714285714285, 0.8571428571428571, 8.0, 4.0], 
#     [0.6000000000000001, 3.0, 0.0, 0.0, 0.0]]
# sgm = [0,0,1,1,2]


def schdl_rnk_chk(sgm,Atrps,So):
    l=0
    O=[]
    while l<len(sgm):
        o= np.dot(np.array( So[sgm[l]] ), np.linalg.matrix_power(np.transpose( np.array(Atrps)), l))
        o = o.tolist()
        O.append(o)
        l=l+1
        rnk = matrix_rank(np.array(O))
        
    return rnk

# m is number of sensors in So and n is the schedule length

def rnd_schdl(m,n):
# this function generates random schedules such that all sensors in the set So are used 
    brk =True
    while brk ==True:
        sgm = []
        i=0
        while i<n:
            sgm.append(randint(0,m-1))
            i=i+1
        
        
    #checking to make sure all sensors in So are atleast 
    #for once present in the random schedule
        i=0
        t=0
        while i < m:
            if i in sgm :
                t=t+1
            i=i+1
            
        if t==m:
            brk =False
                
    return sgm
                


n=len(Atrps)
m=len(So)
sgm_rnd = rnd_schdl(m,n)

rnk = schdl_rnk_chk(sgm_rnd,Atrps,So)