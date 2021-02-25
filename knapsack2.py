#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 14:11:56 2017

@author: hosseintohidi
"""
from gurobipy import *

def findL(m, n, k, omega, c, a, b, w, q):
    M = [i for i in range(m)]
    N = [i for i in range(n)]
    K = [i for i in range(k)]
    Omega = [i for i in range(omega)]

    # Formulation in Extensive form
    modelL = Model('Two-stage stochastic multiple binary knapsack problem extensive form')
#    x = modelL.addVars(M, N, name = "X")
#    y = modelL.addVars(K, N,Omega ,name = "Y")
    x = modelL.addVars(M, N, name="X", vtype=GRB.BINARY)
    y = modelL.addVars(K, N, Omega, name="Y", vtype = GRB.BINARY)

    obj2 = modelL.addVars(Omega ,name = "OBJ2")



    modelL.addConstrs((quicksum(a[i]*x[i,j] for i in M)<=b[j] for j in N) ,"C1")

    modelL.addConstrs((quicksum(x[i,j] for j in N)<=1 for i in M) ,"C2")

    modelL.addConstrs((quicksum(a[i]*x[i,j] for i in M)+quicksum(w[k]*y[k,j,omega] for k in K)<=b[j] for j in N for omega in Omega) ,"C01")

    modelL.addConstrs((quicksum(y[k,j,omega] for j in N)<=1 for k in K for omega in Omega) ,"C02")

    modelL.addConstrs((obj2[omega]==quicksum(q[k][omega]*y[k,j,omega] for k in K for j in N)for omega in Omega) ,"C03")

#    modelL.addConstrs((x[i,j]>=0 for i in M for j in N), "R1")

#    modelL.addConstrs((x[i,j]<=1 for i in M for j in N), "R2")

#    modelL.addConstrs((y[k,j,omega]>=0 for k in K for j in N for omega in Omega), "R3")

#    modelL.addConstrs((y[k,j,omega]<=1 for k in K for j in N for omega in Omega), "R4")

    # Objective
    obj1=quicksum(c[i]*x[i,j] for i in M for j in N)
    obj3=quicksum(1.0/len(Omega)*obj2[omega]for omega in Omega)

    objj=obj3+obj1
    # run the model
    modelL.setObjective(obj3, GRB.MAXIMIZE) # maximize profit
    #optimizing the modelL and get variables as well as the SP model.
    modelL.update()
    modelL.optimize()
    modelL.write("outL.mst")
    modelL.write("outL.sol")
    modelL.write("modelL.lp")

    SOL={}

    for v in modelL.getVars():
        if v.X != 0:
            SOL.update({v.Varname: v.X})
            #print("%s %f" % (v.Varname, v.X))
    temp=0
    for ww in Omega:
        temp=temp+1.0/len(Omega)*SOL.get('OBJ2['+str(ww)+']')

    L=temp
    return L
