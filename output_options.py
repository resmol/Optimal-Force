#!/usr/bin/env python3.6
"""
Created on Wed Feb  2 09:21:13 2022 in BRIHUEGA

@author: Alejandro Jodra
"""

from __future__ import print_function, division
import os,logging

width=80

def print_header(*args, **kwargs):
    print((ret_header(*args, **kwargs)))

def ret_header(title=None, ioptions=None, cfile=None, ver='1.0'):
    hstr  = width*'=' + '\n'

    if int(ioptions['npaso']) == 1:
        hstr += addlinec("BRIHUEGA FORCE %s"%ver)
        hstr += addlinec("Optimal force analysis")
        hstr += addlinec()
        hstr += addlinec("Author: Alejandro Jodra")
        hstr += addlinec("Contributions by: ")

    hstr += width*'-' + '\n'

    hstr += add_step(ioptions)

    hstr += addlinec()

    if not title==None:
        hstr += width*'-' + '\n'
        hstr += addlinec(title)

    hstr += width*'=' + '\n'

    return hstr

def addlinec(line=""):
    return "|" + line.center(width-2) + "|\n"

def addlinel(line="", lpad=5):
    return "|" + lpad*' ' + line.ljust(width-2-lpad) + "|\n"

def add_vectorXYZ(title, V):
    rstr = ''
    rstr += width*'-' + '\n'
    rstr += addlinec("%s"%title)
    rstr += addlinel('number     X          Y          Z')
    rstr += width*'-' + '\n' 
    try:
        for i in range(len(V)):
            rstr += addlinel("%2i    %10.7f %10.7f %10.7f"%(i+1,V[i].x,V[i].y,V[i].z))
    except:
        for i in range(int(len(V)/3)):
            rstr += addlinel("%2i    %10.7f %10.7f %10.7f"%(i+1,V[3*i],V[3*i+1],V[3*i+2]))
   
    return rstr

def add_step(ioptions):

    rstr = ''
    if 'npaso' in ioptions.opt_dict:
        rstr += addlinel("Step Number= %s"%ioptions['npaso'])
    if 'angle' in ioptions.opt_dict:
        rstr += addlinel("Vectors angle= %9.4f"%ioptions['angle'])
    if 'Gaussian Energy' in ioptions.opt_dict:
        rstr += addlinel("Gaussian Energy= %10.7f kcal/mol"%ioptions['Gaussian Energy'])
    if 'Gaussian Energy state0' in ioptions.opt_dict:
        rstr += addlinel("Gaussian Energy state0= %12.7f kcal/mol"%ioptions['Gaussian Energy state0'])
    if 'Gaussian Energy state1' in ioptions.opt_dict:
        rstr += addlinel("Gaussian Energy state1= %12.7f kcal/mol"%ioptions['Gaussian Energy state1'])
    if 'Predicted Energy' in ioptions.opt_dict:
        rstr += addlinel("Predicted Energy= %12.7f kcal/mol"%ioptions['Predicted Energy'])
    if 'Predicted Energy state0' in ioptions.opt_dict:
        rstr += addlinel("Predicted Energy state0= %12.7f kcal/mol"%ioptions['Predicted Energy state0'])
    if 'Predicted Energy state1' in ioptions.opt_dict:
        rstr += addlinel("Predicted Energy state1= %12.7f kcal/mol"%ioptions['Predicted Energy state1'])
    if 'Gaussian Force' in ioptions.opt_dict:
        rstr += add_vectorXYZ('Gaussian Force', ioptions['Gaussian Force'])
    if 'Gaussian modForce' in ioptions.opt_dict:
        rstr += addlinel("Gaussian modForce= %10.7f nN"%ioptions['Gaussian modForce'])
    if 'Predicted Force' in ioptions.opt_dict:
        rstr += add_vectorXYZ('Predicted Force', ioptions['Predicted Force'])
    if 'Predicted modForce' in ioptions.opt_dict:
        rstr += addlinel("Predicted modForce= %10.7f nN"%ioptions['Predicted modForce'])

    return rstr

def add_stoen(cfile):
    try:
        cfileb = os.path.basename(cfile)
    except TypeError:
        return ''

    rstr = ''
    if cfileb in ['analyze_tden.py', 'analyze_sden.py']:
        rstr += addlinec()
        rstr += addlinel("Transition density matrix analysis:", 3)
        rstr += addlinel("F. Plasser and H. Lischka")
        rstr += addlinel("J. Chem. Theory Comput. (2012), 8, 2777.")
        rstr += addlinec()
        rstr += addlinel("Transition and difference density matrix analysis:", 3)
        rstr += addlinel("F. Plasser, M. Wormit, A. Dreuw")
        rstr += addlinel("J. Chem. Phys. (2014), 141, 024106.")

    return rstr

def add_exciton(ioptions):
    try:
        prop_list = ioptions['prop_list']
    except TypeError:
        return ''

    rstr = ''

    if ('RMSeh' in prop_list) or ('dexc' in prop_list):
        rstr += addlinec()
        rstr += addlinel("Exciton analysis:", 3)
        rstr += addlinel("S. A. Baeppler, F. Plasser, M. Wormit, A. Dreuw")
        rstr += addlinel("Phys. Rev. A (2014), 90, 052521.")

    if 'RMSeh' in prop_list:
        rstr += addlinec()
        rstr += addlinel("Approximate RMSeh/dexc formula:", 3)
        rstr += addlinel("S. A. Mewes, J.-M. Mewes, A. Dreuw, F. Plasser")
        rstr += addlinel("Phys. Chem. Chem. Phys. (2016), 18, 2548.")

    if 'dH-E' in prop_list or 'Corr' in prop_list or 'sigH' in prop_list or 'sigE' in prop_list:
        rstr += addlinec()
        rstr += addlinel("Statistical analysis of excitations:", 3)
        rstr += addlinel("F. Plasser, B. Thomitzni, S. A. Baeppler et al.")
        rstr += addlinel("J. Comput. Chem. (2015), 36, 1609.")

    if ioptions['comp_dntos']:
        rstr += addlinec()
        rstr += addlinel("Conditional densities and DNTOs:", 3)
        rstr += addlinel("F. Plasser")
        rstr += addlinel("ChemPhotoChem (2019), DOI: 10.1002/cptc.201900014.")

    return rstr

def add_entanglement(ioptions):
    try:
        prop_list = ioptions['prop_list']
    except TypeError:
        return ''

    rstr = ''

    if ('S_HE' in prop_list) or ('Z_HE' in prop_list):
        rstr += addlinec()
        rstr += addlinel("Electron-hole entanglement:", 3)
        rstr += addlinel("F. Plasser")
        rstr += addlinel("J. Chem. Phys. (2016), 144, 194107.")

    return rstr

def add_cclib(ioptions):
    try:
        rtype = ioptions['rtype'].lower()
    except TypeError:
        return ''

    rstr = ''

    if rtype in ['cclib', 'gamess', 'orca']:
        rstr += addlinec()
        rstr += addlinel("cclib for structure parsing (http://cclib.github.io):", 3)
        rstr += addlinel("N. M. O'Boyle, A. L. Tenderholt, K. M. Langner")
        rstr += addlinel("J. Comput. Chem. (2008), 29, 839.")

    return rstr

def add_orbkit(ioptions):
    try:
        ok_use = ioptions['cube_orbitals'] or ioptions['comp_p_h_dens'] or ioptions['comp_rho0n']
    except TypeError:
        return ''

    rstr = ''
    if ok_use:
        rstr += addlinec()
        rstr += addlinel("orbkit for orbital/density plotting (http://orbkit.github.io):", 3)
        rstr += addlinel("G. Hermann, V. Pohl, J. C. Tremblay, B. Paulus, H.-C. Hege, A. Schild")
        rstr += addlinel("J. Comput. Chem. (2016), 37, 1511.")

    return rstr

def add_VIST(cfile):
    try:
        cfileb = os.path.basename(cfile)
    except TypeError:
        return ''

    rstr = ''
    if cfileb in ['plot_VIST.py']:
        rstr += addlinec()
        rstr += addlinel("Visualization of chemical shielding tensors (VIST):", 3)
        rstr += addlinel("F. Plasser, F. Gloecklhofer")
        rstr += addlinel("Eur. J. Org. Chem. (2021), DOI: 10.1002/ejoc.202100352.")

    return rstr

def gaussian_calc(ioptions,V1,V2):
    data = open('state0.com','w')
    data.write('%chk={}\n%mem={}\n%nproc={}\n#P {}\n\n{}\n\n{} {}\n'.format(ioptions['chk1'],ioptions['mem'],ioptions['nproc'],ioptions['calc1'],ioptions['comment'],ioptions['charge1'],ioptions['mult1']))
    for i in range(len(V2)):
        data.write('{:2s}     {:10.7f}    {:10.7f}    {:10.7f}\n'.format(V1[i],V2[i].x,V2[i].y,V2[i].z))
    data.write('  \n')
    data.close()

    data = open('state1.com','w')
    data.write('%chk={}\n%mem={}\n%nproc={}\n#P {}\n\n{}\n\n{} {}\n'.format(ioptions['chk2'],ioptions['mem'],ioptions['nproc'],ioptions['calc2'],ioptions['comment'],ioptions['charge2'],ioptions['mult2']))
    for i in range(len(V2)):
        data.write('{:2s}     {:10.7f}    {:10.7f}    {:10.7f}\n'.format(V1[i],V2[i].x,V2[i].y,V2[i].z))
    data.write('  \n')
    data.close()

def jmol(V1,V2,V3,file):
    data = open(file,'w')
    data.write('{}\n\n'.format(len(V1)))
    try:
        for i in range(len(V3)):
            data.write('{:2s}     {:10.7f}    {:10.7f}    {:10.7f}      {:10.7f}    {:10.7f}    {:10.7f}\n'.format(V1[i],V2[i].x,V2[i].y,V2[i].z,V3[i].x,V3[i].y,V3[i].z))
    except:
        for i in range(len(V1)):
            data.write('{:2s}     {:10.7f}    {:10.7f}    {:10.7f}      {:10.7f}    {:10.7f}    {:10.7f}\n'.format(V1[i],V2[i].x,V2[i].y,V2[i].z,V3[3*i],V3[3*i+1],V3[3*i+2]))
    data.close()

def finish(microiter, rmsdt, ndqt, angle):
    if ndqt > 1e-1:
        print("      newCartesian Iter: %i Failed to obtain coordinates (Threshold = %.6e |dQ| = %.6e, alpha = %.6e)\n" % (microiter, rmsdt, ndqt, angle))
    elif ndqt > 1e-3:
        print("      newCartesian Iter: %i Approximate coordinates obtained (Threshold = %.6e |dQ| = %.6e), alpha = %.6e\n" % (microiter, rmsdt, ndqt, angle))
    else:
        print("      newCartesian Iter: %i Cartesian coordinates obtained (Threshold = %.6e |dQ| = %.6e), alpha = %.6e\n" % (microiter, rmsdt, ndqt, angle))

