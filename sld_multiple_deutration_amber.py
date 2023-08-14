# Written by Mohd Ibrahim
import MDAnalysis as mda
import numpy as np
import argparse
import random
from tqdm import tqdm

parser = argparse.ArgumentParser(description='Calculate Neutron scattering length densities,at multiple\
                                 deutration levels')
parser.add_argument('-f', dest='xtc_file', type=str,default='rna.gro',help='xtc file')
parser.add_argument('-s', dest='tpr_file', type=str,default='rna.gro',help='tpr file')
parser.add_argument('-dz',dest='dz',type=float,default= 0.25,help='bin width in Angstroms')
parser.add_argument('-skip',dest='skip', type=int,default=1,help='skip every nth frame')
parser.add_argument('-deu','--list', nargs='+',dest='Deutration_level',help='Deuterium percentages, default is 0 8 100'
                    ,default=[0,100], type=float)

parser.add_argument('-b', dest='begin_time', type=float,default=0,help='Beginning time in ps')
parser.add_argument('-e', dest='end_time', type=float,default=0,help='End time in ps')

parser.add_argument('-o',dest='OutputFilename', type=str, help='outputfile')




        
def AssignScatteringLength(Atomname, deutration,component):

    deutration /= 100 ## Convert from percentage to fraction
    probability = random.random()

    if probability < deutration: HydrogenSL = 6.671

    else: HydrogenSL = -3.7390

    if component =='water':

        if Atomname == 'HW1' or Atomname=='HW2': return HydrogenSL

        if Atomname == 'OW':  return  5.803 
        if Atomname == 'NIO': return  3.63
        if Atomname == 'CXY': return  9.5770 

    if component =='all':

        if Atomname == 'HW1' or Atomname=='HW2': return HydrogenSL #6.671 :q!change to -3.7390 for normal water

        elif Atomname[0]=='O': return  5.803 

        elif Atomname[0]=='H': return -3.7390 

        elif Atomname[0]=='C' and Atomname != 'CXY': return 6.64

        elif Atomname[0]=='N' and Atomname != 'NIO':return 9.36

        elif Atomname[0]=='P': return 5.13

        elif Atomname == 'NIO': return  3.63

        elif Atomname == 'CXY': return  9.5770 

        else: print("Atom name '%s' not listed" %(Atom))
            

            
def WriteOutput(SLDBins, Deutration_level,OutputFilename):
     
    output_file = open(OutputFilename, 'w+')
    
    output_file.write("First column has distance along the z-direction in Angstroms, rest of the column\
 has Neutron \n SLDs at different Deuterium percentages depicted in each column. For each Deuterium level,\
SLD for \n full system [FS] and for just water and ions (JW) is evaluated. SLD units are (10^-5 A^-2)\n");

    output_file.write("\n%12s"%"Dist.[A]");

    for Col_name in Deutration_level:
        label = "Fs"+str(Col_name) + "%D2O"
        output_file.write("%12s"%label)
        label = "JW"+str(Col_name) + "%D2O"
        output_file.write("%12s"%label)

    output_file.write("\n\n");
    size = np.shape(SLDBins)

    for i in range(0,size[0]):
        for j in range(size[1]):output_file.write("%12.3f"%SLDBins[i,j])
        output_file.write("\n"); 
    output_file.close()
    
    return 0


#######Input paramters#################
args = parser.parse_args()
xtc_file = args.xtc_file
tpr_file = args.tpr_file
dz  = args.dz
OutputFilename = args.OutputFilename
Deutration_level = args.Deutration_level
####################################

u = mda.Universe(tpr_file, xtc_file)
N_atoms, Dimension, = len(u.atoms), u.dimensions,
N = int(Dimension[2]*1.3/dz) 
SLDBins    = np.zeros((N,2*len(Deutration_level)+1)) 


#SLDBins[:,0] = np.linspace(-Dimension[2]/2, Dimension[2]/2,N)
SLDBins[:,0] = np.linspace(0, Dimension[2]*1.3,N)
size = np.shape(SLDBins)

Slice_Volume = 0

###################################################

begin_time = args.begin_time
end_time   = args.end_time
skip       = args.skip
Timestep   = u.trajectory.ts.dt
begin    = round(begin_time/Timestep)
if end_time ==0: end_time = u.trajectory.totaltime
end   = round(end_time/Timestep)


for ts in tqdm(u.trajectory[begin:end:skip]):
  
    Atomname    =  u.atoms.names
    Resname     =  u.atoms.resnames
    Resids      =   u.atoms.resids
    Z_Position  =  u.atoms.positions[:,2]
    
    
    for atom_index in range(N_atoms):
        
        atomname    = Atomname[atom_index]
        resname     = Resname[atom_index]
        z_position  = Z_Position [atom_index] #+ Shift_z
        
        BinPosition = int(np.round(z_position/dz + dz/2))
        
        if BinPosition >= N: continue
            
        
            
      #  BinCounter[BinPosition] += 1
        deutration_index = 1
        
        for deutration in Deutration_level:

            
            SLDBins[BinPosition,deutration_index] +=   AssignScatteringLength(atomname, deutration,'all')
            
            if resname == 'SOL' or atomname=='CXY' or atomname=='NIO':
    
                SLDBins[BinPosition,deutration_index+1] +=   AssignScatteringLength(atomname, deutration,'water')
            
            deutration_index += 2
            
    Slice_Volume += u.dimensions[0]*u.dimensions[1]*dz
    

for i in  range(1,size[1]):  SLDBins[:,i] /= (Slice_Volume)            
            

WriteOutput(SLDBins, Deutration_level,OutputFilename)


