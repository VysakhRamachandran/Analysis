#!/usr/bin/env python
# coding: utf-8

# In[1]:


import MDAnalysis as mda
import numpy as np
import sys


# In[2]:


def cos_sim(a, b):
    """Takes 2 vectors a, b and returns the cosine similarity 
    """
    dot_product = np.dot(a, b) # x.y
    norm_a = np.linalg.norm(a) #|x|
    norm_b = np.linalg.norm(b) #|y|
    return dot_product / (norm_a * norm_b)


# In[3]:


def dist_array(arr):
    ret_arr = np.zeros((len(arr),4))
    for i in range(len(arr)):
        ret_arr[i][0]=np.linalg.norm(arr[i])
        ret_arr[i][1]=arr[i][0]
        ret_arr[i][2]=arr[i][1]
        ret_arr[i][3]=arr[i][2]
    
    return ret_arr


# In[4]:


u = mda.Universe('md_0_1.gro','md_0_1.xtc')


# In[5]:


solv_shell_4A = u.atoms.select_atoms('name OW and around 4.2 protein').positions


# In[11]:


r_min = 2.4000
r_max = 4.2000
dr    = 0.5000 #the user's choice
r_thr =15.0000 #threshhold r value, also choice
# q_array = np.zeros(len(solv_shell_4A))
# for w in range(len(solv_shell_4A)):
S_k=0
q_array = np.zeros(100) #for only 100 water molecules in this frame, ideally should be len(solv_shell_4A) --> all waters.
for w in range(100):
    coord_temp = np.zeros(3, dtype = np.float32)
    for i in range(3):
        coord_temp[i] = solv_shell_4A[w][i]
        x = coord_temp[0]
        y = coord_temp[1]
        z = coord_temp[2]
    r_1 = coord_temp
#locate the water molecule of interest
    select_string=f'name OW and (prop x=={x} and prop y=={y} and prop z=={z})'
    water_1=u.select_atoms(select_string)


    neighbour_list = u.select_atoms(f'name OW and around {r_max} {select_string}')
#if the number of neighbours within the cutoff of 4.2 Angstrom is less than 4,
#increment r_cut unitill a r_final
    while len(neighbour_list) < 4:     
        r = r_max + dr
        neighbour_list = u.select_atoms(f'name OW and around {r} {select_string}')
        if r > r_thr:
            print('Error in the coordinates')
            sys.exit('density too low for liquid water')
            break
    
    R_n = neighbour_list.positions
    R_0 = np.zeros(R_n.shape,dtype=np.float32)
    for i in range(len(R_n)):
        for j in range(3):
            R_0[i][j]=r_1[j]
    
    R_rel = R_n - R_0
    dist_R_rel = dist_array(R_rel)
    
    #sort out the nearest neighbour configurations.
    dist_nearest_4 = dist_R_rel[dist_R_rel[:,0].argsort()][0:4]
    d_4 = dist_nearest_4[:,0:1]
    #print(d_4) #remove print in actual code
    S_k=1-(1.0/3.0)*((np.var(d_4, dtype=np.float32))/(np.mean(d_4, dtype=np.float32))**2)   
    q_array[w] = S_k
       


# In[12]:


np.savetxt("qval_10.txt", q_array) 


# In[ ]:




