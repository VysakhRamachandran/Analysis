import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
u=mda.Universe("437d.gro")
rna_p=u.select_atoms('(nucleic or resname GTP) and (name O* or name P* or name N* or name C*) ')
print(rna_p)
n_rna_p = len(rna_p)
print(n_rna_p)
indexy=list(u.select_atoms('(nucleic or resname GTP) and (name O* or name P* or name N* or name C*) ').indices+1)
print(indexy)
final=np.zeros(shape= [n_rna_p,n_rna_p])
for ts in u.trajectory:
	rna_p=u.select_atoms('(nucleic or resname GTP) and (name O* or name P* or name N* or name C*) ')
	n_rna_p = len(rna_p)
##print('LID has {} residues and NMP has {} residues'.format(n_rna_p, n_k))

	dist_arr = distances.distance_array(rna_p.positions,rna_p.positions,box=u.dimensions)
	dist_arr[dist_arr<=3.9]=1
	dist_arr[dist_arr>3.9]=0
	final=final+dist_arr
#print(final)
#print(final.shape)
final=final*0.02

list1 = final.tolist()

for i, j in zip(list1,indexy):
    for k,z in zip(i,indexy):
        print(j,z,k)


#plt.imshow(final, cmap ='Greens')
#plt.show()
##list1 = final.tolist()

##for i, j in zip(list1,indexy):
##    for k,z in zip(i,indexy):
##        print(j,z,k)


#tick_interval=1
#fig, ax = plt.subplots()
#im = ax.imshow(dist_arr, cmap ='Greens')
#ax.set_yticks(np.arange(n_k)[::5])
#ax.set_xticks(np.arange(n_rna_p)[::20])
#ax.set_yticklabels([i for i in range(1,62)])
#ax.set_xticklabels(indexy, rotation='vertical')
#plt.ylabel('rna_o1p')
#plt.xlabel('k')
#plt.title('Distance between k and o2p at 3000')
#cbar = fig.colorbar(im)
#cbar.ax.set_ylabel('Distance (Angstrom)')
#plt.show()





