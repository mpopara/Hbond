# geometrical definition of Hbond employed here
# is inspired with S. Raschka paper doi: https://doi.org/10.1101/260612
# Here, I implemented three conditions to be simultaneously satisfied for Hbond to exist
# (distance H-A, angle D-H-A, angle H-A-PA; PA is acceptor-bound atom, aka pre-acceptor)
# One can add additional constraint, e.g. distance D-A)


from __future__ import division
import mdtraj
import numpy


def Hbond_occurence(trajectory_name, topology_name, donor, hydrogen, acceptor, pre_acceptor, outfilename_h):
    print(">>>> Reading %s"%(trajectory_name))
    traj = mdtraj.load(trajectory_name, top=topology_name)
    traj.unitcell_angles = None
    traj.unitcell_lengths = None
    traj.unitcell_vectors = None

    print(">>>> Calculating distances ...")
    proton_idx = traj.topology.select(hydrogen)
    acceptor_idx = traj.topology.select(acceptor)
    distances = mdtraj.compute_distances(traj, [(proton_idx[0], acceptor_idx[0])]) # in nm(!); native unit of mdtraj
    
    
    print(">>>> Calculating angles ...")
    D_coords = traj.xyz[:,traj.topology.select(donor),:]
    H_coords = traj.xyz[:,traj.topology.select(hydrogen),:]
    A_coords = traj.xyz[:,traj.topology.select(acceptor),:]
    PA_coords = traj.xyz[:,traj.topology.select(pre_acceptor),:]

    angles_D_H_A = numpy.zeros((len(traj.xyz),1))
    for i in range(len(traj.xyz)):
        d_HD = H_coords[i,:,:] - D_coords[i,:,:] 
        d_HA = H_coords[i,:,:] - A_coords[i,:,:] 
        e_HD = d_HD/numpy.sqrt(numpy.sum(d_HD * d_HD)) # unit vector
        e_HA = d_HA/numpy.sqrt(numpy.sum(d_HA * d_HA)) #unit vector
        angles_D_H_A[i] = numpy.degrees(numpy.arccos(numpy.sum(e_HD * e_HA)))
        
 
        
    angles_H_A_PA = numpy.zeros((len(traj.xyz),1))
    for i in range(len(traj.xyz)):
        d_AH = A_coords[i,:,:] - H_coords[i,:,:] #d1
        d_APA = A_coords[i,:,:] - PA_coords[i,:,:] #d2
        e_AH = d_AH/numpy.sqrt(numpy.sum(d_AH * d_AH)) #e1
        e_APA = d_APA/numpy.sqrt(numpy.sum(d_APA * d_APA)) #e2
        angles_H_A_PA[i] = numpy.degrees(numpy.arccos(numpy.sum(e_AH * e_APA)))
    
  

    print(">>>> Calculating h(t) ...")

    distances_c = distances<0.25  # distance H-A less then 2.5A
    angles_D_H_A_c = (angles_D_H_A>=120) & (angles_D_H_A<=180)
    angles_H_A_PA_c = (angles_H_A_PA>=90) & (angles_H_A_PA<=180)
    
    full_geometry = distances_c*angles_D_H_A_c*angles_H_A_PA_c #all three conditions must be satisfied simultaneously
    h = full_geometry *1 #convert true and fals to 1 and 0
                         # alternative way h = full_geometry.astype('uint8') 
                         # or h = full_geometry + 0

    
    fraction = ((numpy.count_nonzero(h))*100)/traj.n_frames
    print('Occurence is [%]:',fraction)     # total occurence of hydrogen bond 
    
    print(">>>> Writing %s"%(outfilename_h))
    outfile = open(outfilename_h,'w')
    outfile.write("##Frame\th(t)\n")
    for i in range(len(h)):
        outfile.write("%i\t%i\n"%(i+1,h[i][0]))
    outfile.close()            
     


def main():
    for traj_no in ["prod1"]:#,"prod2","prod3","prod4","prod5"]:
        trajectory_name = "trajectory_%s.dcd"%(traj_no)
        topology_name = "topology.pdb"
        outfilename_h = "15-111_Hbond_occurence_%s.dat"%(traj_no) #change the prefix on the file name
        

        donor = "resSeq 15 and name NZ"
        hydrogen = "resSeq 15 and name HZ1"
        acceptor = "resSeq 111 and name OE1"
        pre_acceptor = "resSeq 111 and name CD"
        Hbond_occurence(trajectory_name, topology_name, donor, hydrogen, acceptor, pre_acceptor, outfilename_h)
    
if __name__ == "__main__":
    main()

    