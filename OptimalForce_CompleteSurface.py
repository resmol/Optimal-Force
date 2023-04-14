#!/usr/bin/env python3.6
"""
Created on Thu Jul  7 09:21:13 2022 in BRIHUEGA

@author: Alejandro Jodra
"""

###################################################################
## This is a module witch helps to calculate the optimal force    # 
##                                                                #
##                                                                #
###################################################################


import numpy as np
import os,math,Read_mol,maple,input_options,math_utilities,output_options,CleanHessGrad,sys
import random as rnd

__metaclass__= type

############### CLASS DEFINITION ###########################

class Point:
    """ Represent a point in a space of three dimensions.
    attributes: x,y,z."""

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return'%g, %g, %g' % (self.x, self.y, self.z)

def Extended_Hessian(H,g,gant,q):
    return H + np.outer(g-gant,g-gant)/((g-gant)@q) - np.outer(H@q,q@H)/(q@H@q)

def First_step(p,H,g,Fext,npaso,Inc):
    up = p/math.sqrt(p@p)
    HP = up @ (H @ H) @ up
    Hgp = (g @ H) @ up
    C = math_utilities.Second_grade_equation(HP,2*Hgp,g @ g - Fext**2)
    if isinstance(C,bool): return C
#    C = Fext/math.sqrt(HP)
    if Inc == 'Pos' or Inc == 'POS' or Inc == 'pos':
        q = C[0]*up 
    elif Inc == 'Neg' or Inc == 'NEG' or Inc == 'neg':
#        if npaso == 1:
        q = C[1]*up
#        else:
#            q = C1*up
#    print(C1,C2)

    return q

def Second_step(q,H1,H2,g1,g2):
    p = ((H1 @ H1) @ q) + H1 @ g1
    s = (g2 - g1 + (H2-H1) @ q)*1
    esc = math.acos(np.dot(s,p)/(math.sqrt(np.dot(s,s))*math.sqrt(np.dot(p,p))))*(180/math.pi)
    sp = s - math_utilities.proj(p,s)
#    print(math_utilities.gs_cofficient(gp/math.sqrt(gp@gp), g/math.sqrt(g@g)),math_utilities.gs_cofficient(p/math.sqrt(p@p), g/math.sqrt(g@g)))
    
    return sp, esc

def Third_step(p,q,H,g,Inc):
    if Inc == 'Pos' or Inc == 'POS' or Inc == 'pos':
        a = 0.0035
    elif Inc == 'Neg' or Inc == 'NEG' or Inc == 'neg':    
        a = -0.0035
    newq = a * p + q
    up = newq/math.sqrt(newq@newq)
    HP = up @ (H @ H) @ up
    gp = up @ g
    Hgp = (g @ H) @ up
    Fext = math.sqrt(q @ (H @ H) @ q + g @ g + 2*(g @ H @ q))
    C1,C2 = math_utilities.Second_grade_equation(HP,2*Hgp,g@g-Fext**2)
#    C1 = math.sqrt(q @ (H @ H) @ q)/(math.sqrt(HP)
    newq = C1*up

#    print(C1,C2)
    
    return newq

def Check_step(H,g1,g2):
    u = g2 - g1
    v = H@g1
    ang = math.acos(np.dot(u,v)/(math.sqrt(np.dot(u,u))*math.sqrt(np.dot(v,v))))*(180/math.pi)

    return ang

massdict = {'H':1,'B':5,'C':6,'N':7,'O':8,'F':9,'S':16,'Cl':17}

def Setup_variables(data1,data2,General_options):
    if int(General_options['APROXSTATE2']) == 2:
        if General_options['ANALYTICALSTATE2'] == False and int(General_options['PASO']) > 1:
            if General_options['EXTENDED HESSIAN 2'] == True:
                cord2, numat, st2_energy, st2_grad_cart = data2.cord, data2.numat, data2.energy, data2.gx
                H2ant = maple.MapleMatrix_data('Traj_Info/General/st2_hess_cart')
                g2ant = np.transpose(maple.MapleMatrix_data('Traj_Info/General/st2_grad_cart'))[0]
                disp = np.transpose(maple.MapleMatrix_data('Traj_Info/General/disp_cart'))[0]
                st2_hess_cart = Extended_Hessian(H2ant,st2_grad_cart,g2ant,disp)
            else:
                cord2, numat, st2_energy, st2_grad_cart, st2_hess_cart = data2.cord, data2.numat, data2.energy, data2.gx, data2.Hess
        else:
            cord2, numat, st2_energy, st2_grad_cart, st2_hess_cart = data2.cord, data2.numat, data2.energy, data2.gx, data2.Hess
    elif int(General_options['APROXSTATE2']) == 1:
        if General_options['ANALYTICALSTATE2'] == False and int(General_options['PASO']) > 1:
            cord2, numat, st2_energy, st2_grad_cart = data2.cord, data2.numat, data2.energy, data2.gx
            st2_hess_cart = np.zeros((3*numat,3*numat),dtype=float)
        else:
            cord2, numat, st2_energy, st2_grad_cart, st2_hess_cart = data2.cord, data2.numat, data2.energy, data2.gx, data2.Hess
            st2_hess_cart = np.zeros((3*numat,3*numat),dtype=float)

    if int(General_options['APROXSTATE1']) == 2:
        if General_options['ANALYTICALSTATE1'] == False and int(General_options['PASO']) > 1:
            if General_options['EXTENDED HESSIAN 1'] == True:
                cord1, numat, st1_energy, st1_grad_cart = data1.cord, data1.numat, data1.energy, data1.gx
                H1ant = maple.MapleMatrix_data('Traj_Info/General/st1_hess_cart')
                g1ant = np.transpose(maple.MapleMatrix_data('Traj_Info/General/st1_grad_cart'))[0]
                disp = np.transpose(maple.MapleMatrix_data('Traj_Info/General/disp_cart'))[0]
                st1_hess_cart = Extended_Hessian(H1ant,st1_grad_cart,g1ant,disp)
            else:
                cord1, numat, st1_energy, st1_grad_cart, st1_hess_cart = data1.cord, data1.numat, data1.energy, data1.gx, data1.Hess
        else:
            cord1, numat, st1_energy, st1_grad_cart, st1_hess_cart = data1.cord, data1.numat, data1.energy, data1.gx, data1.Hess
    elif int(General_options['APROXSTATE1']) == 1:
        if General_options['ANALYTICALSTATE1'] == False and int(General_options['PASO']) > 1:
            cord1, numat, st1_energy, st1_grad_cart = data1.cord, data1.numat, data1.energy, data1.gx
            st1_hess_cart = np.zeros((3*numat,3*numat),dtype=float)
        else:
            cord1, numat, st1_energy, st1_grad_cart, st1_hess_cart = data1.cord, data1.numat, data1.energy, data1.gx, data1.Hess
            st1_hess_cart = np.zeros((3*numat,3*numat),dtype=float)

    return cord2,st1_energy,st2_energy,st1_grad_cart,st2_grad_cart,st1_hess_cart,st2_hess_cart

def run():
    file = 'GENERAL.INPUT' ; General_options = input_options.read_options(ifile = file)
    file = 'GAUSSIAN.INPUT' ; Gaussian_options = input_options.read_options(ifile = file) 

    file_st1 = 'state0.fchk'
    file_st2 = 'state1.fchk'
    data_st1 = Read_mol.Data_fchk(file_st1)
    data_st2 = Read_mol.Data_fchk(file_st2)

    cord,st1_energy,st2_energy,st1_grad_cart,st2_grad_cart,st1_hess_cart,st2_hess_cart = Setup_variables(data_st1,data_st2,General_options)

    if General_options['CLEAN']:
        st1_hess_cart = CleanHessGrad.Hess_clean_proy(cord,data_st1.atomNumb,st1_hess_cart)
        st2_hess_cart = CleanHessGrad.Hess_clean_proy(cord,data_st1.atomNumb,st2_hess_cart)
        st1_grad_cart = CleanHessGrad.Vector_clean(cord,data_st1.atomNumb,st1_grad_cart)
        st2_grad_cart = CleanHessGrad.Vector_clean(cord,data_st1.atomNumb,st2_grad_cart)
    maple.output_maple('Traj_Info/General/st1_grad_cart',np.array(st1_grad_cart)); maple.output_maple('Traj_Info/General/st1_hess_cart',st1_hess_cart)
    maple.output_maple('Traj_Info/General/st2_grad_cart',np.array(st2_grad_cart)) ; maple.output_maple('Traj_Info/General/st2_hess_cart',st2_hess_cart)

    angle_gauss = Check_step(st1_hess_cart,st1_grad_cart,st2_grad_cart)

    Fextgauss = st1_grad_cart ; modgauss = math.sqrt(Fextgauss@Fextgauss) 

    ooptions = input_options.write_options(ifile = 'Output.Force')
    ooptions.write_option('npaso', General_options['PASO'])
    ooptions.write_option('angle', angle_gauss)
    ooptions.write_option('Gaussian Energy', (st2_energy-st1_energy)*627.5)
    ooptions.write_option('Gaussian Energy state0', (st1_energy))
    ooptions.write_option('Gaussian Energy state1', (st2_energy))
    ooptions.write_option('Gaussian modForce', modgauss*82.387)
    ooptions.write_option('Gaussian Force', Fextgauss)

    Fext = modgauss + (float(General_options['INCFORCE']))/82.387
    numat = len(data_st1.symb)
    coor = [Point(cord[3*j]/1.88973,cord[3*j+1]/1.88973,cord[3*j+2]/1.88973) for j in range(numat)]
#    print(npaso)

    if 'VECTOR' in General_options.opt_dict:
        v = General_options['VECTOR'].split(sep=',')
        v = np.array(v,dtype = 'float')

    if int(General_options['PASO']) == 1 or General_options['INC'] == 'NEG' or General_options['INC'] == 'neg' or General_options['INC'] == 'Neg':
        if 'VECTOR' in General_options.opt_dict:
            output_options.jmol(data_st1.typ,coor,v/math.sqrt(v@v),'Traj_Info/General/jmol_initial_geom')
            q = First_step(v,st1_hess_cart,st1_grad_cart,Fext,int(General_options['PASO']),General_options['INC'])
            if isinstance(q,bool):
                print('end program 2')
                output_options.print_header(ioptions = ooptions)
                sys.exit()            
        else:
            output_options.jmol(data_st1.typ,coor,st2_grad_cart/math.sqrt(st2_grad_cart@st2_grad_cart),'Traj_Info/General/jmol_initial_geom')
            q = First_step(st2_grad_cart,st1_hess_cart,st1_grad_cart,Fext,int(General_options['PASO']),General_options['INC'])
            if isinstance(q,bool):
                print('end program 2')
                output_options.print_header(ioptions = ooptions)
                sys.exit()

    else:
        disp = np.transpose(maple.MapleMatrix_data('Traj_Info/General/disp_cart'))[0]
        output_options.jmol(data_st1.typ,coor,disp/math.sqrt(disp@disp),'Traj_Info/General/jmol_initial_geom')
        q = First_step(disp,st1_hess_cart,st1_grad_cart,Fext,int(General_options['PASO']),General_options['INC'])
        if isinstance(q,bool):
            output_options.print_header(ioptions = ooptions)
            sys.exit()

#print(82.387*math.sqrt((q @ (st1_hess_cart @ st1_hess_cart) @ q) + (st1_grad_cart @ st1_grad_cart) + (2*(st1_grad_cart @ st1_hess_cart @ q))))
    for i in range(1000000):
        sp,ang = Second_step(q,st1_hess_cart,st2_hess_cart,st1_grad_cart,st2_grad_cart)
        mod = math.sqrt(sp@sp)
        if i == 0: print(mod)
#        print(mod,i)
        if mod <= 0.001 or i == 1000000-1:
            output_options.finish(i, 0.001, mod, ang)
            break
        q = Third_step(sp,q,st1_hess_cart,st1_grad_cart,General_options['INC'])
    
    maple.output_maple('Traj_Info/General/disp_cart',np.array(q))    
    cord += q 
    IncExc =  (st2_grad_cart - st1_grad_cart) @ q + 0.5*(q @ ((st2_hess_cart - st1_hess_cart) @ q))
#IncExc =  st2_grad_cart @ q 
    Fext = st1_hess_cart @ q + st1_grad_cart ; mod = math.sqrt(Fext@Fext)
    coor = [Point(cord[3*j]/1.88973,cord[3*j+1]/1.88973,cord[3*j+2]/1.88973) for j in range(numat)]

    ooptions.write_option('Predicted Energy', (st2_energy-st1_energy+IncExc)*627.5)
    ooptions.write_option('Predicted Energy state0', (st1_energy+st1_grad_cart @ q + 0.5*(q @ ((st1_hess_cart) @ q))))
    ooptions.write_option('Predicted Energy state1', (st2_energy+st2_grad_cart @ q + 0.5*(q @ ((st2_hess_cart) @ q))))
    ooptions.write_option('Predicted modForce', mod*82.387)
    ooptions.write_option('Predicted Force', Fext)
    output_options.print_header(ioptions = ooptions)

    output_options.gaussian_calc(Gaussian_options,data_st1.typ,coor)

    output_options.jmol(data_st1.typ,coor,st1_grad_cart/math.sqrt(st1_grad_cart@st1_grad_cart),'Traj_Info/General/jmol_Fopt_gauss')
    output_options.jmol(data_st1.typ,coor,Fext/math.sqrt(Fext@Fext),'Traj_Info/General/jmol_Fopt_predicted_NextStep')
    output_options.jmol(data_st1.typ,coor,q/math.sqrt(q@q),'Traj_Info/General/jmol_Fopt_disp')
    maple.output_maple('Traj_Info/General/Fpred_cart',np.array(Fext))
    maple.output_maple('Traj_Info/General/Fgauss_cart',np.array(st1_grad_cart))

run()
