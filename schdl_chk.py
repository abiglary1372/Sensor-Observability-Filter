#schedual checker 
from numpy import linalg
from numpy.linalg import matrix_rank
import numpy as np
from random import randint

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
##example 3.1
###########
# Atrps=[[1,2,4,5,7],
#         [0,11,4,8,7],
#         [0,0,-8,0,3],
#         [0,0,0,7,0],
#         [0,0,0,0,6]]


# So = [[3.0332, -5.1353, 2.6428, 0, 3], 
#   [5.1047, -22.2857, 0.8571, 8, 4], 
#   [0.6, 3, 0, 0, 0]]

# sgm = [0,0,0,1,2]
# sgm = [1,2,0,0,0]
# sgm = [0,1,0,0,2]
# sgm = [2,1,0,0,0]

# sgm = [0,2,0,1,1]


# So = [[1.1, -3.0, 0.0, 3.0, 0.0], 
#       [-0.7953216374269005, -0.42105263157894735, 2.0, 0.0, 0.0], 
#       [1.3333333333333333, -16.0, 0.0, 8.0, 0.0], 
#       [1.1666666666666665, -14.0, 0.0, 7.0, 0.0], 
#       [1.2000000000000002, 6.0, 0.0, 0.0, 0.0],
#       [-1.1929824561403508, -0.631578947368421, 3.0, 0.0, 0.0]]
# sgm = [0,0,1,5,6]


#example 3.2

So = [[-3.666666666666667, -6.0, 4.0], [1.0, 0.0, 0.0]]

SoIndx = [0,1]

Atrps=[[2,3,4],
       [0,8,9],
       [0,0,5]]


def main():
    all_schdl_tst()


def all_schdl_tst():
    
    for x1 in SoIndx:
        for x2 in SoIndx:
            for x3 in SoIndx:
                sgm = [x1,x2,x3]
                rnk = schdl_rnk_chk(sgm,Atrps,So)
                print(sgm,"\n")
                print(rnk,"\n")


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
                
def rnd_schdl_gnrl(m,n):
# this function generates random schedules 


    sgm = []
    i=0
    while i<n:
        sgm.append(randint(0,m-1))
        i=i+1
 
                
    return sgm

if __name__ == "__main__":
    main()

# n=len(Atrps)
# m=len(So)
# sgm_rnd = rnd_schdl_gnrl(m,n)

# rnk = schdl_rnk_chk(sgm,Atrps,So)