#!/usr/bin/env python
# coding: utf-8



import MDAnalysis as mda
import numpy as np


def cos_sim(a, b):
    """Takes 2 vectors a, b and returns the cosine similarity 
    """
    dot_product = np.dot(a, b) # x.y
    norm_a = np.linalg.norm(a) #|x|
    norm_b = np.linalg.norm(b) #|y|
    return dot_product / (norm_a * norm_b)

def dist_array(arr):
    """takes an array of positions (another array of length 2)  and 
    return an arryay composed of the norm of the positions and the positions
    eg; arr[N][3] is the input
    output will be Array[N][4] = [d, x, y, z]
    """
    ret_arr = np.zeros((len(arr),4))
    for i in range(len(arr)):
        ret_arr[i][0]=np.linalg.norm(arr[i])
        ret_arr[i][1]=arr[i][0]
        ret_arr[i][2]=arr[i][1]
        ret_arr[i][3]=arr[i][2]
    
    return ret_arr

def S_k(arr):
"""Calculate the translational tetrahedral order parameter given the array containing the distances of the four adjacent water oxygens to the the water (oxygen) of interest"""
    return 1-(1.0/3.0)*(((np.var(d_4, dtype=np.float32))/2*np.mean(d_4, dtype=np.float32))**2)


xtc_file = input("Enter the name of the .xtc file")
tpr_file = input("Enter the name of the .tpr file")


#r_min = 2.4000
r_max = 4.2000
#dr    = 0.6000 #it is user's choice
#r_thr =15.0000 #threshhold r value, also choice



import MDAnalysis.transformations as trans
"""The following lines are meant to center the protein in the box, it helps to prevent the boundary effects"""
u2 = mda.Universe(tpr_file, xtc_file)
protein2 = u2.select_atoms('protein')
not_protein2 = u2.select_atoms('not protein')
transforms = [trans.unwrap(protein2),
              trans.center_in_box(protein2, wrap=True),
              trans.wrap(not_protein2)]
u2.trajectory.add_transformations(*transforms)
"""code for centering ends here"""

with open('Tet_OP.txt', 'w') as f:
    t = 0.0
    for ts in u.trajectory:
        t = t + 1
        solv_shell_4A = u.atoms.select_atoms('name OW and around 4.2 protein')

        for hyd_water in solv_shell_4A:
            ne_water = u.select_atoms(f'name OW and around {r_max} index {hyd_water.id}')
            if(len(ne_water)) < 4: """If the number of waters within 4.2A is less than 4 (coordination number), skip the particular water molecule/step"""
                continue
            else:
                R_n = ne_water.positions
                R_0 = np.zeros(R_n.shape,dtype=np.float32)
                for i in range(len(R_n)):
                    for j in range(3):
                        R_0[i][j] = hyd_water.position[j]
                R_rel = R_n - R_0
                dist_R_rel = dist_array(R_rel)
                dist_nearest_4 = dist_R_rel[dist_R_rel[:,0].argsort()][0:4]
                R_4 = dist_nearest_4[:,1:]
                d_4 = dist_nearest_4[:,0:1]
                """orientational tetrahedral order parameter"""
                cos_sum = 0.0
                for i in range(3):
                    for j in range(i+1,4,1):
                        cos = cos_sim(R_4[i],R_4[j])
                        cos_sum = cos_sum + (cos + 1.0/3.0)*(cos + 1.0/3.0)
                q = 1.0 - (3.0/8.0)*cos_sum


                f.write(f'{t}\t{S_k(d_4)}\t{q}')
                f.write('\n')


