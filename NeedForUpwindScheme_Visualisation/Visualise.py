import numpy as np
import math
import matplotlib.pyplot as plt
def TDMA(A,b):														#solve Ax = b using TDM algorithm
    p=[]
    q=[]
    p.append(-1*A[0][1]/A[0][0])
    q.append(b[0]/A[0][0])
    x=np.zeros(len(b))
    for i in range(1,len(b)):
        if i!= len(b)-1:
            p.append(-1*A[i][i+1]/(A[i][i]+A[i][i-1]*p[i-1]))
        q.append((b[i]-A[i][i-1]*q[i-1])/(A[i][i]+A[i][i-1]*p[i-1]))
    x[len(b)-1]=q[len(b)-1]
    for i in range(len(b)-2,-1,-1):
        x[i]=p[i]*x[i+1]+q[i]
    return x 
def visulaise(peclet):
    length=peclet.size
    fig, axs = plt.subplots(length,figsize=(4,10))
    plt.ylim(-1,1)
    for iterator in range(length):
        n   =   100                                                 #No of nodes
        x   =   np.arange(0,n,1)                                    #value of x at every node
        pe  =   peclet[iterator]               
        delx=   0.1                                                 #space step
        Texact=np.zeros(n)                                          
        pec=pe*delx                                                 #cell peclet number
        aa=0.5+pec/4                                                #coefficients in TDMA formulation
        bb=0.5-pec/4                                                #coefficients in TDMA formulation
        b=np.zeros(n-2)                                             #b vector in Ax=b formulation
        b[n-3]=-1*bb                                            
        A=np.zeros((n-2,n-2),dtype='f')                             #coefficient matrix
        for i in range(n-2):
            A[i][i]=-1
            if i>0:
                A[i][i-1]=aa
            if i<n-3:
                A[i][i+1]=bb
        T=TDMA(A,b)                                                 #numerical solution
        T=np.insert(T,0,0)                                          #Adding boundary condition
        T=np.append(T,1)                                            #Adding boundary condition
        axs[iterator].plot(x,T,label='fds')
        try:
            for i in range(n):
                Texact[i]=(math.exp(pe*x[i]/10)-1)/(math.exp(pe)-1)
            for i in range(n):
                Texact[i]=Texact[i]/Texact[n-1]
            axs[iterator].plot(x,Texact,label='exact')
        except:
            plt.xlim(80,99)
            print("Peclet number being ",pe)
        axs[iterator].legend()
    plt.show()
peclet=np.array([0.01,10,20,50,100])
visulaise(peclet)
