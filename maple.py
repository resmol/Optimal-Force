#!/usr/bin/env python3.6
"""
Created on Thu May 26 23:03:30 2022 in BRIHUEGA

@author: Alejandro Jodra
"""

import numpy as np

def output_maple(file,matrix):
    data = open(file,'w')
    data.write('%%MatrixMarket matrix coordinate real general\n')
    try:
        data.write('{:<3d} {:<3d} {:<7d}\n'.format(matrix.shape[0],matrix.shape[1],matrix.shape[0]*matrix.shape[1]))
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                data.write('{:<5d} {:<5d} {:17.9E}\n'.format(i+1,j+1,matrix[i,j]))
    except IndexError:
        data.write('{:<3d} 1   {:7d}\n'.format(matrix.shape[0],matrix.shape[0]))
        for i in range(matrix.shape[0]):
            data.write('{:<5d} 1   {:17.9E}\n'.format(i+1,matrix[i]))        
    data.close()

def MapleMatrix_data(file):
    data = open(file)
    lines = [line for line in data.read().split()]
    ndim = int(lines[5]) ; mdim = int(lines[6])
    Mat = np.empty((ndim,mdim), dtype = float)
    for i in range(ndim):
        for j in range(mdim):
            Mat[i,j] = lines[10+3*ndim*j+3*i]

    return Mat
