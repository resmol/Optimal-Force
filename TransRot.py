#!/usr/bin/env python3.6
"""
Created on Wed Feb  2 09:21:13 2022 in BRIHUEGA

@author: Alejandro Jodra
"""

import numpy as np
import math_utilities,math

class Point:
    """ Represent a point in a space of three dimensions.
    attributes: x,y,z."""

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return'%g, %g, %g' % (self.x, self.y, self.z)

def Inertia_matrix(atom):
    massdict = {1:1.007825,5:10.811,6:12.0011,7:14.003074,8:15.999,9:18.998,16:32.065,17:35.453,53:126.9044}
    numat = len(atom.cart)
    I_matrix = np.zeros((3,3), dtype = float)
#    L = np.zeros((3),dtype=float)

    I_matrix[0,0] =  sum([massdict[atom.atoms[n]] * (atom.cart[n].y**2 +  atom.cart[n].z**2) for n in range(numat)]) 
    I_matrix[0,1] = -sum([massdict[atom.atoms[n]]*atom.cart[n].x*atom.cart[n].y for n in range(numat)])
    I_matrix[0,2] = -sum([massdict[atom.atoms[n]]*atom.cart[n].x*atom.cart[n].z for n in range(numat)])

    I_matrix[1,0] = I_matrix[0,1]
    I_matrix[1,1] = sum([massdict[atom.atoms[n]] * (atom.cart[n].x**2 +  atom.cart[n].z**2) for n in range(numat)]) 
    I_matrix[1,2] = -sum([massdict[atom.atoms[n]]*atom.cart[n].y*atom.cart[n].z for n in range(numat)])

    I_matrix[2,0] = I_matrix[0,2]
    I_matrix[2,1] = I_matrix[1,2]
    I_matrix[2,2] = sum([massdict[atom.atoms[n]] * (atom.cart[n].x**2 +  atom.cart[n].y**2) for n in range(numat)]) 

    # Angular Momentum
#    lin_mom = [scalar_vector(massdict[atom[k].numat],atom[k].vel) for k in range(numat)]

#    lista = [cross(atom[n].cart,lin_mom[n]) for n in range(numat)]
#    L[0] = sum([lista[n][0] for n in range(numat)])
#    L[1] = sum([lista[n][1] for n in range(numat)])
#    L[2] = sum([lista[n][2] for n in range(numat)])

    return I_matrix 

def centro_masas(geom,matoms):
    
    massdict = {1:1.007825,5:10.811,6:12.0011,7:14.003074,8:15.999,9:18.998,16:32.065,17:35.453,53:126.9044}
    n = len(geom)
    Rc = Point()
    Msum = sum([massdict[matoms[i]] for i in range(n)])
    Rc.x = sum([geom[i].x*massdict[matoms[i]] for i in range(n)])
    Rc.y = sum([geom[i].y*massdict[matoms[i]] for i in range(n)])
    Rc.z = sum([geom[i].z*massdict[matoms[i]] for i in range(n)])
#    print(Rc.x/Msum,Rc.y/Msum,Rc.z/Msum)
    for i in range(n):
        geom[i].x = geom[i].x - Rc.x/Msum
        geom[i].y = geom[i].y - Rc.y/Msum
        geom[i].z = geom[i].z - Rc.z/Msum       
        
    return geom

def TransRot_coordinates(mol):

    numat = len(mol.cart)
    geom = centro_masas(mol.cart,mol.atoms)
#    geom = mol.cart
    RI = Inertia_matrix(mol)
    DI,VI = np.linalg.eig(RI)
    VI[:,0] = -VI[:,0] ; VI[:,1] = -VI[:,1]
#    VI = list(VI) ; DI = list(DI) 
#    VI[1] , VI[2] = VI[2] , VI[1]
#    DI[1] , DI[2] = DI[2] , DI[1]
#    VI = np.array(VI) ; DI = np.array(DI)    
    Q = np.zeros((6,3*numat), dtype = float)
    for i in range(numat):
        for j in range(3):
            l = 1 ; ll = 2
            Q[j,3*i+j] = 1./math.sqrt(numat)
            for k in range(3):
                Q[k+3,3*i+j] = (math_utilities.prod_escalar(geom[i],VI[l]) * VI[j,ll] - math_utilities.prod_escalar(geom[i],VI[ll]) * VI[j,l])
                l += 1 ; ll += 1
                if l > 2:
                    l = 0
                if ll > 2:
                    ll = 0

    for i in range(3,6):
        mod = math.sqrt(Q[i] @ Q[i])
        Q[i] = Q[i]/mod

    return Q

def TransRot_coordinates2(mol):
    numat = len(mol.cart)
#    print([[mol.cart[i].x,mol.cart[i].y,mol.cart[i].z] for i in range(numat)])
    geom = centro_masas(mol.cart,mol.atoms)
#    geom = mol.cart
#    print([[geom[i].x,geom[i].y,geom[i].z] for i in range(numat)])
    Q = np.zeros((6,3*numat), dtype = float)
    for i in range(numat):
        l = 0 ; ll = 0
        for j in range(3):
            Q[j,3*i+j] = 1./math.sqrt(numat)
            Q[3,3*i+j] = ll*geom[i].z + l*geom[i].y
            ll +=1
            if j == 2:
                ll -= 1 ; l = -1
            Q[4,3*i+j] = ll*geom[i].z + l*geom[i].x
            if j >= 1:
                l += 1
            ll = -ll
            Q[5,3*i+j] = ll*geom[i].y + l*geom[i].x

    for i in range(3,6):
        mod = math.sqrt(Q[i] @ Q[i])
        Q[i] = Q[i]/mod

    return Q
