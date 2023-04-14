#!/usr/bin/env python3.6
"""
Created on Sat Mar 28 23:53:51 2020

@author: Alejandro Jodra
"""
import numpy as np
import math,os,error_handler

def unit_vector(v):
    mod = math.sqrt(prod_escalar(v,v)) ; v2 = [Point() for i in range(len(v))]
    try:    
        ndim = len(v)
        for i in range(ndim):
            v2[i].x = v[i].x/mod
            v2[i].y = v[i].y/mod
            v2[i].z = v[i].z/mod
    except TypeError:
        v2.x = v.x/mod
        v2.y = v.y/mod
        v2.z = v.z/mod
  
    return v2

def gs_cofficient(v1, v2):
    return np.dot(v2, v1) / np.dot(v1, v1)    

def multiply(cofficient, v):
    return list(map(lambda x : x * cofficient, v))

def proj(v1, v2): #Vector v2 proyectado en v1
    return np.array(multiply(gs_cofficient(v1, v2) , v1))

def Gram_Schmidt(X): #Vectores a procesar tienen que estar en las filas
    Y = []
    for i in range(len(X)):
        temp_vec = X[i]
        for inY in Y :
            proj_vec = proj(inY, X[i])
            #print "i =", i, ", projection vector =", proj_vec
            temp_vec = list(map(lambda x, y : x - y, temp_vec, proj_vec))
            #print "i =", i, ", temporary vector =", temp_vec
        Y.append(temp_vec)
    return np.array(Y) #Salida de los vectores en las filas

def norm_vect(X):
    try:
        Y = [X[i]/math.sqrt(X[i] @ X[i]) for i in range(len(X))]
    except TypeError:
        Y = X/math.sqrt(X @ X)
    return Y

def Second_grade_equation(a,b,c):
    raiz = (b**2)-(4*a*c)
    if raiz < 0.:
        Error = error_handler.NegativeRootError(raiz)
        return Error()
    x1 = (-b+math.sqrt((b**2)-(4*a*c)))/(2*a)
    x2 = (-b-math.sqrt((b**2)-(4*a*c)))/(2*a)

    return x1,x2

def kronecker(m,n):
    if m == n:
        kron = 1.
    else:
        kron = 0.

    return kron

def prod_escalar(V1,V2):
    x = 0.0
    y = 0.0
    z = 0.0   
    try:    
        ndim = len(V1)
        for i in range(ndim):
            x += V1[i].x * V2[i].x
            y += V1[i].y * V2[i].y
            z += V1[i].z * V2[i].z
        p = x + y + z    
    except TypeError:
        try:
            x = V1.x * V2.x
            y = V1.y * V2.y
            z = V1.z * V2.z
            p = x + y + z   
        except AttributeError:
            x = V1.x * V2[0]
            y = V1.y * V2[1]
            z = V1.z * V2[2]
            p = x + y + z 
    except:      
        print('Se ha producudo un error',)

    return p
