#!/usr/bin/env python3.6
"""
Created on Sat Mar 28 23:53:51 2020 in BRIHUEGA

@author: Alejandro Jodra
"""
import numpy as np
import math
import os

__metaclass__ = type

# Some Constants
delta = 1.e-5
au_time = 2.41888432e-2  # from femtoseconds to au
eh_ek = 3.157746e5       # relation Eh/Kb Hartrees/Boltzmann
au = 1822.88851543       # relation N/e
a0 = 0.529177249         # from bohr to angstrom
eh_kcal = 627.50469      # from Hartrees to Kcal/mol
kb = 3.1668152037e-6     # Boltzmann constant in au hartree/Kelvin
rad2deg = 180.0 / math.pi
deg2rad = math.pi / 180.0
angstrom2bohr = 1.88973

def atom_type(object):
    w = list()
    element = {1.:'H',5.:'B',6.:'C',7.:'N',8.:'O',9.:'F',16.:'S',17.:'Cl',35.:'Br',53.:'I'}
    for val in object:
        if val in element:
            w = w + list(element[val])
    return w

class Point:
    """ Represent a point in a space of three dimensions.
    attributes: x,y,z."""

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return'%g, %g, %g' % (self.x, self.y, self.z)

class Molecule_data:
    def __init__(self, gfile=""):
        self.path = gfile
        data = open(gfile)
        lines = [line for line in data.read().split()]
        data.close()
        try:
            n1 = lines.index('#p')
        except ValueError:       
            n1 = lines.index('#P')     
        massdict = {1:1.0079,5:10.811,6:12.0011,7:14.007,8:15.999,9:18.998,16:32.065,17:35.453,53:126.9044}
        calcOptions = dict(zip(lines[n1+1:n1+12],[True for i in lines[n1+1:n1+12]])) 
        self.calc = dict(zip(['opt','freq','force'],[False for i in ['optNo','freqNo','forceNo']]))      
        if 'opt' in calcOptions or 'Opt' in calcOptions or 'Opt=(TS,noeigentest,calcfc)' in calcOptions or 'opt=(TS,noeigentest,calcfc)' in calcOptions or 'opt=(addredundant,modredundant)' in calcOptions or 'Opt=(TS,noeigentest,calcAll)' in calcOptions or 'Opt=(TS,noeigentest,calcall)' in calcOptions or 'Opt=(TS,noeigentest,CalcAll)' in calcOptions:
            self.calc['opt'] = True
        if 'freq' in calcOptions or 'freq' in calcOptions or 'freq=(NoRaman)' in calcOptions:
            self.calc['freq'] = True      
        if 'opt(Maxcycle=1)' in calcOptions or 'force' in calcOptions or 'opt=ModRedundant' in calcOptions:
            self.calc['force'] = True
        
        n1 = lines.index('NAtoms=')
        self.numat=int(lines[n1+1])
           
        try:
            if self.calc['opt']:
                n1 = 0
                try: 
                    n1 = lines.index('Energy,',n1+1)
                    for i in range(1000):
                        try:
                            n1 = lines.index('Energy,',n1+1)
                            self.energy = float(lines[n1+3])
                        except ValueError:
                            break
                except ValueError:
                    for i in range(1000):
                        try:
                            n1 = lines.index('Done:',n1+1)
                            self.energy = float(lines[n1+3])
                        except ValueError:
                            break
            else:
                n1 = lines.index('Energy,')
                self.energy = float(lines[n1+3])
                n1 = lines.index('Done:')
                self.energy1 = float(lines[n1+3])     
                self.exc = self.energy1 - self.energy
        except ValueError:
            n1 = lines.index('Done:')
            self.energy = float(lines[n1+3])
            
        if self.calc['opt']:  
            try:      
                n2 = lines.index('Stationary')
                try:
                    n3 = lines.index('Standard',n2)
                    if lines[n3+1] == 'orientation:':
                        pass
                    else:
                        n3 = lines.index('Input',n2)
                except ValueError:
                    n3 = lines.index('Input',n2)
            except ValueError:
                n3 = lines.index('Input')
        else:
            try:
                n3 = lines.index('Input')
            except ValueError:
                n3 = lines.index('Standard')        
        geom_all_list = lines[n3+15:n3+15+(6*self.numat)] 
        self.cart = [Point(float(geom_all_list[i*6+3]),float(geom_all_list[i*6+4]),float(geom_all_list[i*6+5])) for i in range(self.numat)]
        self.cart_au = [Point(float(geom_all_list[i*6+3])*1.88973,float(geom_all_list[i*6+4])*1.88973,float(geom_all_list[i*6+5])*1.88973) for i in range(self.numat)]
        self.atoms = [int(geom_all_list[i*6+1]) for i in range(self.numat)]     
        self.AtMass = [massdict[self.atoms[i]] for i in range(self.numat)]    
        
    def logCI_data(self):
        data = open(self.path)  
        lines = [line for line in data.read().split()]
        n1 = lines.index('Multiplicity')
        n2 = lines.index('Variables:')
        conect_all = lines[n1+3:n2]
        numat = int((len(conect_all)+15)/8)
        conect_first = [(2,1),(3,int(conect_all[5]))]
        conect_bond = [(i+4,int(conect_all[9+8*i+1])) for i in range(numat-3)]
        conect_bond = conect_first + conect_bond
        conect_first = [(3,int(conect_all[5]),int(conect_all[7]))]
        conect_ang = [(i+4,int(conect_all[9+8*i+1]),int(conect_all[9+8*i+3])) for i in range(numat-3)]
        conect_ang = conect_first + conect_ang
        conect_dih = [(i+4,int(conect_all[9+8*i+1]),int(conect_all[9+8*i+3]),int(conect_all[9+8*i+5])) for i in range(numat-3)]
        self.bond = conect_bond
        self.ang = conect_ang
        self.dih = conect_dih
        self.conect = conect_bond + conect_ang + conect_dih                 
        try:
            n3 = lines.index('Parameters')
        except ValueError:
            n3 = lines.index('Standard')          
        geom_all_list = lines[n3+18:n3+18+(7*(3*numat-6))]
        self.q = [float(geom_all_list[i*7+2]) for i in range(3*numat-6)] 
        self.q[:len(conect_bond)] = [ i*angstrom2bohr for i in self.q[:len(conect_bond)]]
        self.q[len(conect_bond):] = [ i*deg2rad for i in self.q[len(conect_bond):]]
        
    def Read_Hess_int(self):
        data = open(self.path)
        lines = [line for line in data.read().split()] 
        if self.calc['freq']:    
            n4 = lines.index('Second')   
            n5 = lines.index('matrix:',n4)
            self.Hess_int = np.empty((3*self.numat-6,3*self.numat-6), dtype = float)
            m = int(round(((3*self.numat-6)/5-float(int((3*self.numat-6)/5)))*5))
            for i in range(int((3*self.numat-6)/5)):
                l = int(-0.5*i+0.5*i**2)           
                for j in range(3*self.numat-6-5*i):
                    n4 = j+1 ; n3 = n4
                    if j > 4:
                        n4 = 6 ; n3 = 5
                    for k in range(n3):
#                    print(i,j,l,k,n4)
                        self.Hess_int[j+5*i,i*5+k] = lines[n5+6+1*(j+1)+k+i*(6*(3*self.numat-6)-10+5)-5*6*l+j*5-(5*(n4-1)-(n4-1)**2-int(0.5*(n4-1))+int(0.5*(n4-1)**2))]    
            i = int((3*self.numat-6)/5) ; l = int(-0.5*i+0.5*i**2)
            for j in range(3*self.numat-6-5*i):
                n4 = j+1 ; n3 = n4
                if j > 4:
                    n4 = 6 ; n3 = 5
                for k in range(n3):            
                    self.Hess_int[j+5*i,i*5+k] = lines[n5+m+1+1*(j+1)+k+i*(6*(3*self.numat-6)-10+5)-5*6*l+j*5-(5*(n4-1)-(n4-1)**2-int(0.5*(n4-1))+int(0.5*(n4-1)**2))]
            for i in range(3*self.numat-6):
                for j in range(i):
                    self.Hess_int[j,i] = self.Hess_int[i,j] 
            n4 = lines.index('Eigenvectors')   
            n5 = lines.index('matrix:',n4)
            normMode_alllist = [] ; forceConst_alllist = []
            for i in range(int((3*self.numat-6)/5)):
                lista = lines[n5+13:n5+13+7*(3*self.numat-6)]
                lista2 = lines[n5+8:n5+8+5]
                normMode_alllist += lista
                forceConst_alllist += lista2
                n5 += 7*(3*self.numat-6) + 12   
            lista = lines[n5+2*m+3:n5+2*m+3+(m+2)*(3*self.numat-6)]
            lista2 = lines[n5+m+3:n5+m+3+m]
            normMode_alllist += lista
            forceConst_alllist += lista2
            self.normMode = [i for i in range(3*self.numat-6)] #Los modos normales estan en las filas
            self.forceConst = [i for i in range(3*self.numat-6)] 
            for i in range(3*self.numat-6):
                self.normMode[i] = [j for j in range(3*self.numat-6)]

            for i in range(int((3*self.numat-6)/5)):
                for k in range(5):
                    l = k + i*5
                    self.forceConst[l] = float(forceConst_alllist[l])
                    for j in range(3*self.numat-6):
                        self.normMode[l][j] = float(normMode_alllist[7*j+2+7*(3*self.numat-6)*i+k])
            i = int((3*self.numat-6)/5)
            for k in range(m):      
                l = k + i*5
                self.forceConst[l] = float(forceConst_alllist[l])
                for j in range(3*self.numat-6):
                    self.normMode[l][j] = float(normMode_alllist[(2+m)*j+2+7*(3*self.numat-6)*i+k])    
            self.normMode = np.array(self.normMode)        
#    if calc2 == 'freq' or calc3 == 'force':
 #   if calc3 == 'force': 
    def Read_force_int(self):
        data = open(self.path)
        lines = [line for line in data.read().split()]         
        n4 = lines.index('-DE/DX')
        self.grad_int = [-float(lines[n4+14+i*7]) for i in range(3*self.numat-6)]   
            
        data.close()  

    def Read_force(self):
        data = open(self.path)
        lines = [line for line in data.read().split()]
        if self.calc['freq'] or self.calc['force']:
#    if calc3 == 'force':
            try: 
                n4 = lines.index('restored')
                for i in range(1000):
                    try:
                        n4 = lines.index('restored',n4+1)
                    except ValueError:
                        break                    
            except ValueError:
                n4 = 0
                for i in range(1000):
                    try:    
                        n4 = lines.index('(Hartrees/Bohr)',n4+1)
                    except ValueError:
                        n4 = n4 - 9
                        break
            force_all_list = lines[n4+16:n4+16+5*self.numat]
            self.grad = [Point(-float(force_all_list[i*5+2]),-float(force_all_list[i*5+3]),-float(force_all_list[i*5+4])) for i in range(self.numat)]      
            
    def Read_Hess(self):
        data = open(self.path)
        lines = [line for line in data.read().split()]        
        if self.calc['freq']:    
            n4 = lines.index('coordinates:')   
            normMode_alllist = []
            forceConst_alllist = []
            for i in range(self.numat-2):
                lista = lines[n4+41:n4+41+11*self.numat]
                lista2 = lines[n4+21:n4+21+3]
                normMode_alllist += lista
                forceConst_alllist += lista2
                n4 += 11*self.numat + 40
            self.normMode = [i for i in range(3*self.numat-6)] #Los modos normales estan en las filas
        self.forceConst = [i for i in range(3*self.numat-6)] 
        for i in range(3*self.numat-6):
            self.normMode[i] = [j for j in range(self.numat)]

        for i in range(3*self.numat-6):
            for j in range(self.numat):
                self.normMode[i][j] = Point()
        
        for i in range(self.numat-2):
            for k in range(3):
                l = k + i*3
                self.forceConst[l] = float(forceConst_alllist[l])/15.57
                for j in range(self.numat):
                    self.normMode[l][j].x = float(normMode_alllist[11*j+2+11*self.numat*i+3*k])/math.sqrt(self.AtMass[j])
                    self.normMode[l][j].y = float(normMode_alllist[11*j+3+11*self.numat*i+3*k])/math.sqrt(self.AtMass[j])
                    self.normMode[l][j].z = float(normMode_alllist[11*j+4+11*self.numat*i+3*k])/math.sqrt(self.AtMass[j])
       
    def ExcitedStateGaussread(self,nStates=3):
        data = open(self.path)
        self.coeff = [i for i in range(nStates)] ; self.Energy = [i for i in range(nStates)] ; self.f = [i for i in range(nStates)] ;  self.Gap = [i for i in range(nStates)] ; self.Nat_exc = []
        lines = [line for line in data.read().split()]
        try:
            n1 = lines.index('#p')
        except ValueError:
            n1 = lines.index('#P')
        calc = lines[n1+1:n1+12]
        for i in range(10):
            if calc[i] == 'opt' or calc[i] == 'Opt' or calc[i] == 'Opt=(TS,noeigentest,calcfc)' or calc[i] == 'opt=(TS,noeigentest,calcfc)' or calc[i] == 'opt=(addredundant,modredundant)' or calc[i] == 'Opt=(TS,noeigentest,calcAll)' or calc[i] == 'Opt=(TS,noeigentest,calcall)' or calc[i] == 'Opt=(TS,noeigentest,CalcAll)':
                for j in range(1000):
                    try:
                        n1 = lines.index('strengths:',n1+1)
                        n0 = n1
                    except ValueError:
                        break
                break
            else:
                n0 = 0
        n1 = lines.index('strengths:',n0) ; n2 = lines.index('This',n1)
        num = 0
        if lines[n1+3] == 'symmetry':
            num = 7
        count = 0 ; num2 = [] ; count2 = 0
        for i in lines[n1+11+num:n2]:
            count += 1
            if count == 2:
                if i == '->':
                    num2.append(4)
                    par = 'YES'
                else:
                    num2.append(3)
                    par = 'NO'
            if count == 4 and par == 'YES':
                count = 0 ; count2 += 1
            if count == 3 and par == 'NO':
                count = 0 ; count2 +=1
        Orb = [lines[n1+11+num+i] for i in range(n2-n1-11-num)]
        self.coeff[0] = [2*float(Orb[sum(num2[:i])+num2[i]-1])**2 for i in range(count2)]
        self.Energy[0] = float(lines[n1+7+num]) ; self.Gap[0] = float(lines[n1+5+num]) ; aux = list(lines[n1+9+num]) ; self.f[0] = float(aux[2]+aux[3]+aux[4]+aux[5]+aux[6]+aux[7])
        for i in range(count2):
            if num2[i] == 3:
                #print(Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],coeff[0][i])
                self.Nat_exc.append([Orb[sum(num2[:i])],Orb[sum(num2[:i])+1]])
            if num2[i] == 4:
                #print(Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],Orb[sum(num2[:i])+2],coeff[0][i])
                self.Nat_exc.append([Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],Orb[sum(num2[:i])+2]])
        #print(' ')
        self.Nat_exc.append(' ')
        n3 = n2 + num
        for k in range(1,nStates):
            try:
                n100 = lines.index('nm',n3)
                n2 = lines.index('Excited',n3-1)-1 ; num3 = 0
                if k==1 and lines[n1+3] == 'symmetry':
                    num3 = 7
                try:
                    n100 = lines.index('nm',n3)+1
                    n3 = lines.index('nm',n100)-7
                    count = 0 ; num2 = [] ; count2 = 0
                    for i in lines[n2+11+num3:n3-num]:
                        count += 1
                        if count == 2:
                            if i == '->':
                                num2.append(4)
                                par = 'YES'
                            else:
                                num2.append(3)
                                par = 'NO'
                        if count == 4 and par == 'YES':
                            count = 0 ; count2 += 1
                        if count == 3 and par == 'NO':
                            count = 0 ; count2 +=1
                    Orb = [lines[n2+11+i+num3] for i in range(n3-n2-11-num-num3)]
                    self.coeff[k] = [2*float(Orb[sum(num2[:i])+num2[i]-1])**2 for i in range(count2)]
                    self.Energy[k] = float(lines[n2+7+num3]) ; self.Gap[k] = float(lines[n2+5+num3]) ; aux = list(lines[n2+9+num3]) ; self.f[k] = float(aux[2]+aux[3]+aux[4]+aux[5]+aux[6]+aux[7])
                    for i in range(count2):
                        if num2[i] == 3:
                            #print(Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],coeff[k][i])
                            self.Nat_exc.append([Orb[sum(num2[:i])],Orb[sum(num2[:i])+1]])
                        if num2[i] == 4:
                            #print(Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],Orb[sum(num2[:i])+2],coeff[k][i])
                            self.Nat_exc.append([Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],Orb[sum(num2[:i])+2]])
                    #print(' ')
                    self.Nat_exc.append(' ')
                except ValueError:
                    n3 = lines.index('SavETr:',n100)
                    count = 0 ; num2 = [] ; count2 = 0
                    for i in lines[n2+11:n3]:
                        count += 1
                        if count == 2:
                            if i == '->':
                                num2.append(4)
                                par = 'YES'
                            else:
                                num2.append(3)
                                par = 'NO'
                        if count == 4 and par == 'YES':
                            count = 0 ; count2 += 1
                        if count == 3 and par == 'NO':
                            count = 0 ; count2 +=1
                    Orb = [lines[n2+11+i] for i in range(n3-n2-11)]
                    self.coeff[k] = [2*float(Orb[sum(num2[:i])+num2[i]-1])**2 for i in range(count2)]
                    self.Energy[k] = float(lines[n2+7]) ; self.Gap[k] = float(lines[n2+5]) ; aux = list(lines[n2+9]) ; self.f[k] = float(aux[2]+aux[3]+aux[4]+aux[5]+aux[6]+aux[7])
                    for i in range(count2):
                        if num2[i] == 3:
                            #print(Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],self.coeff[k][i])
                            self.Nat_exc.append([Orb[sum(num2[:i])],Orb[sum(num2[:i])+1]])
                        if num2[i] == 4:
                            #print(Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],Orb[sum(num2[:i])+2],coeff[k][i])
                            self.Nat_exc.append([Orb[sum(num2[:i])],Orb[sum(num2[:i])+1],Orb[sum(num2[:i])+2]])
                    #print(' ')
                    self.Nat_exc.append(' ')
            except ValueError:
                break
        self.energy_st0 = self.energy-self.Gap[0]/27.2107
        
    def Homo_lumo_ener(self):    
        data = open(self.path)
        lines = [line for line in data.read().split()]
        n1 = lines.index('Current')
        n1 = lines.index('Alpha')
        n = int(lines[n1+5])
        self.HL_energies = np.array(lines[n1+6:n1+6+n],dtype=float)

    def get_vector(self,vector):
        " All the coordinates in one vector "
        vec_list = []
        for i in range(len(vector)):
            vec_list += [vector[i].x, vector[i].y, vector[i].z]

        return np.array(vec_list)

class Data_fchk:
    def __init__(self, gfile=""):
        self.cord = list()
        data = open(gfile)
        lines = [line for line in data.read().split()]

    #reading of the Cartesian Cordinates
        n1 = lines.index('Current')
#    print ('Number of atoms', int(int(lines[n1+5])/3))
        self.numat = int(int(lines[n1+5])/3)
        n2 = n1+6
        lista_cord = (lines[n2:n2+3*self.numat])
        self.cord = [float(w) for w in lista_cord]

        #reading atom type
        nN1 = lines.index('Nuclear')
        nN2 = nN1 + 5
        self.N_charges = [float(m) for m in lines[nN2:nN2+self.numat]]
        self.atomNumb = [int(float(m)) for m in lines[nN2:nN2+self.numat]]
        numbers = list(range(self.numat))
        self.typ = atom_type(self.N_charges)
        self.symb =[self.typ[m]+str(numbers[m]) for m in range(len(self.typ))]

        # Reading the Real atomic weights
        n5 = lines.index('Real')
        n6 = n5 + 6
        m_uma = [float(w)*au for w in lines[n6:n6+self.numat]]
        self.mass =  np.array(m_uma,dtype=float)

    # Reading the total energy
        n8 = lines.index('Total')
        n9 = n8 + 3
        self.energy = float(lines[n9])

    #reading of the Gradient in Cartesian Cordinates
        try:
            n3 = lines.index('Optimization')
            n3 = lines.index('Gradient')
            n3 = lines.index('Gradient',n3+1,len(lines))
            n4 = n3 + 4
            grad =(lines[n4:n4+(3*self.numat)])
            self.gx = np.array([float(x) for x in grad],dtype=float)
        except ValueError:
            n3 = lines.index('Gradient')
            n4 = n3 + 4
            grad =(lines[n4:n4+(3*self.numat)])
            self.gx = np.array([float(x) for x in grad],dtype=float)
    # Reading and Organized the Hessian cartesian Matrix
        try:
            n5 = lines.index('Constants')
            n6 = n5 + 4
            numberofdata = int(lines[n5+3])
            vector = (lines[n6:n6+numberofdata])
            ncart = 3*self.numat
            self.Hess = np.empty((ncart,ncart), dtype = float)
            for i in range(1,ncart+1):
                for j in range(1,i+1):
                    n = int((i*i - i)/2 + j)
                    self.Hess[i-1,j-1] =float(vector[n-1])
                    self.Hess[j-1,i-1] = self.Hess[i-1,j-1]
        except ValueError:
            data.close()
