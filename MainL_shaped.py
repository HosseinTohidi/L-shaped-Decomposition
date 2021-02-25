#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 00:17:27 2017

@author: hosseintohidi
"""
import sys
sys.path.insert(0, '/Users/hosseintohidi/Desktop/Lshape_code/')
from  knapsack2 import *
import numpy as np
import os
os.getcwd()
from gurobipy import *

import openpyxl
path='/Users/hosseintohidi/Desktop/Workbook1.xlsx'
wb=openpyxl.load_workbook(path)
wb.get_sheet_names()
sheet=wb.get_sheet_by_name('Sheet8')

# set parameters
m=10
n=2
k=4
kprime=k
omega=7

M=[i for i in range(m)]
N=[i for i in range(n)]
K=[i for i in range(k)]
Omega=[i for i in range(omega)]
a=[];b=[];c=[];q=[];w=[];teta=0;z=-10000
for i in M:
    a.append(sheet['B'+str(i+2)].value)
    c.append(sheet['A'+str(i+2)].value)

a=map(int,a)
c=map(int,c)
for j in N:
    b.append(sheet['C'+str(j+2)].value)
for s in K:
    w.append(sheet['D'+str(s+2)].value)
for kk in K:
    qtemp=[]
    for s in Omega:
        qtemp.append(sheet[chr(ord('E')+s)+str(kk+2)].value)
    
    q.append(map(int,qtemp))
b=map(int,b)
w=map(int,w)

model2s=Model()

def model2s_solve(k,sol):

        omegatemp=k
        model2s= Model('second-stage problem')

        xx=np.zeros((len(M),len(N)))
       
        for i in M:
            for j in N:
                
                st='X['+str(i)+','+str(j)+']'
                xx[i][j]=sol.get(st)

        y = model2s.addVars(K, N,Omega ,name = "Y",vtype = GRB.BINARY)
        #y = model2s.addVars(K, N,Omega ,name = "Y")

        obj2 = model2s.addVars(Omega ,name = "OBJ2")
        
        model2s.addConstrs((quicksum(a[i]*xx[i][j] for i in M)+quicksum(w[k]*y[k,j,omega] for k in K )<=b[j] for j in N for omega in Omega if omega==omegatemp) ,"C01")
        model2s.addConstrs((quicksum(y[k,j,omega] for j in N)<=1 for k in K for omega in Omega if omega==omegatemp) ,"C02")
        model2s.addConstrs((obj2[omega]==quicksum(q[k][omega]*y[k,j,omega] for j in N for k in K)for omega in Omega if omega==omegatemp) ,"C03")
        # Objective
        obj1=quicksum(c[i]*xx[i][j] for i in M for j in N)
        obj3=quicksum(1.0/len(Omega)*obj2[omega]for omega in Omega if omega==omegatemp)
        
        model2s.setObjective(obj3, GRB.MAXIMIZE) # maximize profit
        model2s.setParam( 'OutputFlag', False )

        model2s.update()        
        #model2s.write("model32.lp")
        
        model2s.optimize()
        
        sol2={}
        for v in model2s.getVars():
            sol2.update({v.Varname: v.X})
        
        return sol2
#solve the relaxed second stage problem to generate benders cuts
#create the T matrix
T=np.zeros([kprime+n,m*n])
nnn=0
for jj in N:
    for ii in M:
        T[jj,(ii*n)+jj]=a[ii]
    nnn+=1
#create H matrix
H=b+[1]*k    

def Benders_cut(k,sol,T,H):
        omegatemp=k
        modelB= Model('second-stage problem')

        xx=np.zeros((len(M),len(N)))
       
        for i in M:
            for j in N:
                
                st='X['+str(i)+','+str(j)+']'
                xx[i][j]=sol.get(st)

        #y = model2s.addVars(K, N,Omega ,name = "Y",vtype = GRB.BINARY)
        y = modelB.addVars(K, N,Omega ,name = "Y",ub=1,lb=0)

        obj2 = modelB.addVars(Omega ,name = "OBJ2")
        
        modelB.addConstrs((quicksum(a[i]*xx[i][j] for i in M)+quicksum(w[k]*y[k,j,omega] for k in K )<=b[j] for j in N for omega in Omega if omega==omegatemp) ,"C01")
        modelB.addConstrs((quicksum(y[k,j,omega] for j in N)<=1 for k in K for omega in Omega if omega==omegatemp) ,"C02")
        #modelB.addConstrs((obj2[omega]==quicksum(q[k][omega]*y[k,j,omega] for j in N for k in K)for omega in Omega if omega==omegatemp) ,"C03")
        # Objective
        obj1=quicksum(c[i]*xx[i][j] for i in M for j in N)
        obj3=quicksum(1.0/len(Omega)*quicksum(q[k][omega]*y[k,j,omega] for j in N for k in K)for omega in Omega if omega==omegatemp)
        
        modelB.setObjective(obj3, GRB.MAXIMIZE) # maximize profit
        modelB.setParam( 'OutputFlag', False )

        modelB.update()        
        modelB.write("model32.lp")
        
        modelB.optimize()
        dualB=modelB.getAttr('pi') # get the optimal dual variables' values

        return dualB,modelB.objVal
    

#Build model:        
model = Model('Two-stage stochastic multiple binary knapsack problem extensive form')

x = model.addVars(M, N, name = "X",vtype = GRB.BINARY)
y = model.addVars(K, N,Omega ,name = "Y",vtype = GRB.BINARY)
obj2 = model.addVars(Omega ,name = "OBJ2")
teta = model.addVar(name = "Teta")
model.addConstrs((quicksum(a[i]*x[i,j] for i in M)<=b[j] for j in N) ,"C1")
model.addConstrs((quicksum(x[i,j] for j in N)<=1 for i in M) ,"C2")

# Objective
obj1=quicksum(c[i]*x[i,j] for i in M for j in N)
obj3=quicksum(1.0/len(Omega)*obj2[omega]for omega in Omega)

#find L (here it is an upper bound)
L=findL(m,n,k,omega,c,a,b,w,q)
#L=1000
#model.addConstr((teta<=L) ,"C3")
#model.update()
check=True
BendersCheck=True
vv=0
nn=0
EE=[]
Ef=[]
objlist=[]
teta2=0
EQstar=0

while check:
   # print(nn,'EQstar',EQstar,'theta',teta2)
    nn=nn+1
    if vv==0:
        obj=obj1
        teta2=10000
    else:
        #teta = model.addVar(name = "Teta")
        obj=obj1+teta
        #model.update()
    model.setParam( 'OutputFlag', False )
    model.setObjective(obj, GRB.MAXIMIZE) # maximize profit
    model.update()
   # LogToConsole=0
    model.optimize()
    sol={}
    for v in model.getVars():
        sol.update({v.Varname: v.X})
    teta2=sol.get('Teta')    
    if vv==0:
        teta2=10000
    #solve model2s for k in K 
    sol3={}
    EQstar=0
    sxstar=[]
    xstarindex=[]
    ystarindex=[]

    sx=[]
    yindex=[]
    xindex=[]
    dual=[]
    objB=[]
    for k in Omega:
        sol3={}
        sol3=model2s_solve(k,sol)
        EQstar+=1.0/len(Omega)*sol3.get('OBJ2['+str(k)+']')
        #find benders cut
        dual.append(Benders_cut(k,sol,T,H)[0])
        objB.append(Benders_cut(k,sol,T,H)[1])
    e=0
    T=np.array(T)
    for k in Omega:
        for kkk in range(n+kprime):
            e+=1.0/len(Omega)*dual[k][kkk]*H[kkk]
    E=np.zeros([omega,m*n])
    for k in Omega:
        for iii in range(m*n):
            for kkk in range(n+kprime):
                E[k,iii]+=1.0/len(Omega)*dual[k][kkk]*T[kkk,iii]
                
    #check to see if benders cut can be generated:
    E=sum(E)   
    xx=np.zeros(m*n)
    nnn=0   
    for i in M:
        for j in N:
            st='X['+str(i)+','+str(j)+']'
            xx[nnn]=sol.get(st)
            nnn+=1
    ww=e-E.dot(xx)
    if teta2<= ww:
        BendersCheck=False
    
    E= np.reshape(E,[m,n])
   # print(E)
    print('BendersCheck:',BendersCheck)
        
    EE.append(EQstar)
    Ef.append(teta2)
    objlist.append(model.objVal)

    if abs(EQstar-teta2)<=1 or nn>50:
       check = False
    else:
        print(nn,'Obj=',round(model.objVal,2),'EQstar',round(EQstar,2),'theta',round(teta2,2))
        #add optimality cut
        for i in M:
            for j in N:
                st='X['+str(i)+','+str(j)+']'
               
                if sol.get(st)==1:
                    sxstar.append(st)
                    xstarindex.append(i)
                    ystarindex.append(j)

                else:
                    sx.append(st)
                    xindex.append(i)
                    yindex.append(j)
                    
        cutname='cut'+str(nn)  
        model.reset()         
        #model.addConstr((teta<=(EQstar-L)*(sum(x[i,j] for i,j in zip(xstarindex,ystarindex))-\
                                #sum(x[i,j] for i,j in zip(xindex,yindex))-len(sxstar))+EQstar),cutname)
        model.addConstr((teta<=(EQstar-L)*(sum(x[i,j] for i,j in zip(xstarindex,ystarindex))-sum(x[i,j] for i,j in zip(xindex,yindex))-len(xstarindex)+1)+L),cutname)
       
        #add benders cut
        
        bendcut='benderscut'+str(nn)
        #model.reset()
        if BendersCheck:
            model.addConstr(teta <=e+sum(E[i][j]*x[i,j] for i in M for j in N),bendcut)
        
        model.update()
        #model.write("test22.mps")
        model.write("test22.lp")   
    
    vv+=1
       
    
    
#print('the optimum solution is:',sol)
Ef[0]=2*Ef[1]


#for v in model.getVars():
 #   if v.X != 0:
 #      print("%s %f" % (v.Varname, v.X))
print('Objective function value=',model.objVal)
model.write("out.mst")
model.write("out1.sol")
model.write("iteration44.mps")
model.write("iteration45.lp")

# Lshaped Implementation

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.plot(range(nn),EE,label="EQstar")
plt.plot(range(nn),Ef,label="Theta")
plt.plot(np.arange(nn)[1:],objlist[1:],label="Obj. Function Value")

plt.legend(loc='upper right')
plt.show()

