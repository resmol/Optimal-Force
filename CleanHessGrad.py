#!/usr/bin/env python3.6
"""
Created on Wed Feb  2 09:21:13 2022 in BRIHUEGA

@author: Alejandro Jodra
"""

import numpy as np
import math_utilities,TransRot,math

class Point:
    """ Represent a point in a space of three dimensions.
    attributes: x,y,z."""

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return'%g, %g, %g' % (self.x, self.y, self.z)

class Molecule:
    """Contains the cartesian coordinates of the atoms
    atributes: atoms cart"""

    def __init__(self,atoms = [],all_atoms_cart = []):
        self.atoms = atoms
        self.cart = all_atoms_cart

def Hess_clean_proy(cord,atoms,Hess):
    mol = Molecule(atoms,[Point(cord[3*i]/1.88973,cord[3*i+1]/1.88973,cord[3*i+2]/1.88973) for i in range(len(atoms))])
    numat = len(mol.cart)
#    Q = TransRot.TransRot_coordinates(mol)
    Q = TransRot.TransRot_coordinates2(mol)
    Q = math_utilities.Gram_Schmidt(Q)
    for i in range(len(Q)):
        mod = math.sqrt(Q[i] @ Q[i])
        Q[i] = Q[i]/mod
    for i in range(len(Q)):
        coordinate = [Point(Q[i,3*j],Q[i,3*j+1],Q[i,3*j+2]) for j in range(numat)]
        file_jmol = 'Optimal_force/Trans_rot' + str(i) + '.xyz'
#        output_jmol(mol.atoms,mol.cart,coordinate,file_jmol)
#    A = np.outer(Q[0],Q[0]) + np.outer(Q[1],Q[1]) + np.outer(Q[2],Q[2]) + np.outer(Q[3],Q[3]) + np.outer(Q[4],Q[4]) + np.outer(Q[5],Q[5])
    A = np.transpose(Q) @ Q
    B = np.eye(3*numat,dtype=float) - A
    Hess = B @ Hess @ B
    u,v = np.linalg.eig(Hess)
    v = np.transpose(v)
    if isinstance(u[0],complex):
#        print('Cleaning the complex part')
        u = u.real
        v = v.real
    for i in range(3*numat-11,3*numat):
#        print(u[i])
        coordinate = [Point(v[i,3*j],v[i,3*j+1],v[i,3*j+2]) for j in range(numat)]
        file_jmol = 'Optimal_force/norm_mode' + str(i) + '.xyz'
#        output_jmol(mol.atoms,mol.cart,coordinate,file_jmol)
    for j in range(len(Q)):
        i = j+1
        v[3*numat-i] = Q[j]
        u[3*numat-i] = 0.1 #u[3*numat-i] son los autovalores para los modos normales de rotacion y vibracion
    Hess = np.transpose(v) @ np.diag(u) @ v
#    for i in v:
#        proy = A@i ; proy2 = B@i
#        print(math.sqrt(proy @ proy),math.sqrt(proy2 @ proy2),math.sqrt(i @ i))
    return Hess

def Vector_clean(cord,atoms,V): #Sirve para limpiar un vector fuera del plano XY
    mol = Molecule(atoms,[Point(cord[3*i]/1.88973,cord[3*i+1]/1.88973,cord[3*i+2]/1.88973) for i in range(len(atoms))])
    numat = len(mol.cart) ; numvect = len(V)
#    Q = TransRot.TransRot_coordinates(mol)
    Q = TransRot.TransRot_coordinates2(mol)
    Q = math_utilities.Gram_Schmidt(Q)
    for i in range(len(Q)):
        mod = math.sqrt(Q[i] @ Q[i])
        Q[i] = Q[i]/mod
    if len(V) != len(Q[0]):
        V = np.concatenate((np.zeros(len(Q[0])-len(V)),V),axis=0)
        i = 2 ; j = 3 ; k = 4 ; l = 3 ; m = 8 ; n = 3
#        i = 3-1 ; j = 2-1 ; k = 1-1 ; l = 2-1 ; m = 3-1 ; n = 1-1
        X = [[mol.cart[i].x-mol.cart[j].x,mol.cart[i].y-mol.cart[j].y,mol.cart[i].z-mol.cart[j].z],[mol.cart[k].x-mol.cart[l].x,mol.cart[k].y-mol.cart[l].y,mol.cart[k].z-mol.cart[l].z],[mol.cart[m].x-mol.cart[n].x,mol.cart[m].y-mol.cart[n].y,mol.cart[m].z-mol.cart[n].z]]
        Y = math_utilities.Gram_Schmidt(X)
        for i in range(len(Y)):
            Y[i]=Y[i]/math.sqrt(Y[i]@Y[i])

        B = np.zeros((3*numat,3*numat),dtype=float)
        for i in range(numat):
            B[3*i,3*i:3*i+3]=Y[0]
            B[3*i+1,3*i:3*i+3]=Y[1]
            B[3*i+2,3*i:3*i+3]=Y[2]

        for i in range(len(Q)):
            Q[i] = B@Q[i]

        B = np.zeros((3*numat,3*numat),dtype=float)
        for i in range(numat):
            B[i,3*i] = 1
            B[i+numat,3*i+1] = 1
            B[i+2*numat,3*i+2] = 1

        for i in range(len(Q)):
            Q[i] = B@Q[i]

    A = np.outer(Q[0],Q[0]) + np.outer(Q[1],Q[1]) + np.outer(Q[2],Q[2]) + np.outer(Q[3],Q[3]) + np.outer(Q[4],Q[4]) + np.outer(Q[5],Q[5])
#    A = np.outer(Q[0],Q[0]) + np.outer(Q[1],Q[1]) + np.outer(Q[2],Q[2])
    B = np.eye(3*numat,dtype=float) - A
    V = B @ V

    if 3*numat != numvect:
        return V[2*numat:]
    else:
        return V
