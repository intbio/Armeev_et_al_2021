import sys
sys.path.insert(0,'/home/armeev/projects/PyNAMod/')
import numpy as np
import pynamod
from pynamod.bp_step_geometry import get_params_for_single_step_debug
from pynamod.bp_step_geometry import rebuild_by_full_par_frame_numba,rebuild_by_full_par_frame
import MDAnalysis as mda
from MDAnalysis.analysis import align
#import io,requests

from tqdm.auto import tqdm,trange
import pandas as pd
from pynamod.visual_ngl import  show_ref_frames
from pynamod.non_DNA_geometry import get_obj_orientation_and_location,get_rotation_and_offset_ref_frame_to_obj
from scipy.spatial.distance import pdist, cdist


mda.Universe()
def gen_init_fiber_par_frame_blank(N_steps_in_NCP,N_NCPS,linker_length,BDNA_step=np.array([-0.032,-0.096,0.088,0.030,-15.120,-1.880,0.000, 0.436,3.353,-0.000 , 1.702, 35.666])):
    '''
    args:
    N_steps_in_NCP (int) amount of bp STEPS (N b.p. - 1) in ncp 
    N_NCPS (int) amount of NCPS on fiber
    linker_length (int or iterable of size N_NCPS + 1) linkers added between all NCPs and before and after first NCP
    BDNA_step = np.array([-0.032,-0.096,0.088,0.030,-15.120,-1.880,0.000, 0.436,3.353,-0.000 , 1.702, 35.666])
    Returns:
    np.array
    '''
    if not (linker_length is None):
        n=linker_length*(N_NCPS+1)+N_steps_in_NCP*N_NCPS        
        ncp_dyad_locs=np.arange(linker_length+N_steps_in_NCP//2,n,linker_length+N_steps_in_NCP,dtype=int)
        assert len(ncp_dyad_locs) == N_NCPS
        return (np.tile(BDNA_step,(n,1)),ncp_dyad_locs)

from numba import jit
@jit
def get_obj_orientation_and_location(ref_mat,ref_ori,rotation,offset):
    return(np.dot(ref_mat.T, rotation),np.dot(offset,ref_mat.T)+ref_ori)
    
def gen_fiber(ncp_loc,ncp_frame,all_steps,ncp_bp_frame,bp_num_map,
              R2=np.array([[-0.9194,  0.275,   0.2815], [-0.3555, -0.2733, -0.8939], [-0.1688, -0.9218,  0.349 ]]),
              of_vec=np.array([44.2273, -4.1654, -0.3646]),use_steps=[-70,70],get_dna_pos=False):
    all_steps=all_steps.copy()

    nucleosome_bp_steps=ncp_bp_frame[ncp_frame,bp_num_map[use_steps[0]+1]:bp_num_map[use_steps[1]+1]].reshape(-1,12)
    starts=ncp_loc+use_steps[0]+1
    mask=(np.repeat(np.arange(-use_steps[0]+use_steps[1]).reshape(1,-1),len(ncp_loc),axis=0)+starts.reshape(-1,1)).flatten()
    if (len(np.unique(mask)) != len(mask)) or (mask.max()>=all_steps.shape[0]) or (mask.min()<0):
        return(None)
        
    all_steps[mask]=nucleosome_bp_steps
    ref_frames=rebuild_by_full_par_frame_numba(all_steps)
    beads=np.zeros([len(ncp_loc),3])
    for j,i in enumerate(ncp_loc):
        o1 = ref_frames[i,3,:3]
        R1=ref_frames[i,:3,:3]
        mat,ori=get_obj_orientation_and_location(R1,o1,R2,of_vec)
        beads[j]=ori
    if get_dna_pos:
        #return only linker coordinates
        return(ref_frames,beads,ref_frames[np.in1d(np.arange(len(ref_frames)), mask)==False][:,3,:3])
        #return(ref_frames,beads,ref_frames[::5,3,:3])
   
    return(ref_frames,beads)

#def _gen_fiber

def calc_E(beads,Rm=60, eps=1):
    distances=pdist(np.array(beads).reshape((-1,3)))
    LJ=np.sum(eps*((Rm/distances)**12 - 2*(Rm/distances)**6))
    return(LJ)

def calc_E_w_dna(beads,dna_beads,Rm=60, eps=1,Rm_dna=10,eps_dna=0.5):
    #ncp_vs_ncp
    distances=pdist(np.array(beads).reshape((-1,3)))
    Enn=np.mean(eps*((Rm/distances)**12 - 2*(Rm/distances)**6))
    #linker_vs_linker
    distances_dna=pdist(np.array(dna_beads).reshape((-1,3)))
    Ell=np.mean(eps_dna*((Rm_dna/distances_dna)**12 - 2*(Rm_dna/distances_dna)**6))
    #ncp_vs_linker
    distances_dna_ncp=cdist(beads, dna_beads)
    Enl=np.mean(((eps_dna+Rm_dna)/2)*(((Rm+Rm_dna)/distances_dna_ncp)**12 - 2*((Rm+Rm_dna)/distances_dna_ncp)**6))
    
    
    LJ=Enn+Ell+Enl
    return(LJ)

def calc_E_w_dna_logistic(beads,dna_beads,Rm=60, eps=1,Rm_dna=10,eps_dna=0.5):
    #ncp_vs_ncp
    distances=pdist(np.array(beads).reshape((-1,3)))
    Enn=np.sum(eps*(1-1/(1+np.exp(-distances+Rm))))
  
    #linker_vs_linker
    distances_dna=pdist(np.array(dna_beads).reshape((-1,3)))
    Ell=np.sum(eps_dna*(1-1/(1+np.exp(-distances+Rm_dna))))
    
    distances_dna_ncp=cdist(beads, dna_beads)
    Enl=np.sum(((eps_dna+Rm_dna)/2)*(1-1/(1+np.exp(-distances+(Rm+Rm_dna)))))
    return(Enn+Ell+Enl)



from six import StringIO
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import align
def to_mda_traj(pos_traj,frame_traj,all_steps,ncp_bp_frames,bp_num_map,mute=True,alignTrj=False):
    #gen topology in pdb
    pdb_string=''
    bp_frames,octamers=gen_fiber(pos_traj[0],frame_traj[0],all_steps,ncp_bp_frames,bp_num_map)
    num=1
    pdb_string=pdb_string+'MODEL     1\n'
    for row in bp_frames[:,3,:3]/10:
        pdb_string=pdb_string+"ATOM  %5d  N   DNA A   1    %8.3f%8.3f%8.3f\n"%(num,row[0],row[1],row[2])
        num+=1
    for row in octamers/10:            
        pdb_string=pdb_string+"ATOM  %5d  O   NUC B   1    %8.3f%8.3f%8.3f\n"%(num,row[0],row[1],row[2])
        num+=1
    pdb_string=pdb_string+"ENDMDL\n"
    
    #gen trajectory np array
    coords=np.zeros((len(pos_traj),bp_frames.shape[0]+len(octamers),3))
    for frame,(positions,frames) in tqdm(enumerate(zip(pos_traj,frame_traj)),total=len(pos_traj),disable=mute):
         bp_frames,octamers=gen_fiber(positions,frames,all_steps,ncp_bp_frames,bp_num_map)
         coords[frame]=np.vstack([bp_frames[:,3,:3],octamers])
            
    u = mda.Universe(StringIO(pdb_string),coords/10.0,topology_format='pdb', format=MemoryReader, order='fac')
    
    if alignTrj:
        alignment = align.AlignTraj(u, u,in_memory=True)
        alignment.run()
    return(u)

def simulate_MC(ncp_bp_frames,bp_num_map,
               N_ncp=10,init_linker_length=16,ncp_bp_length=147,
               EN_initial=10000,linker_step_size=1,n_frames_md=1000,
               N_mc_evals=10000,N_results=None,MAX_mc_evals=np.inf,Rm=70, eps=1,Rm_dna=10,eps_dna=0.5,KT=1,E_type='LJ',mute=True):
    init_linker_length=init_linker_length+1
    linker_step_size=linker_step_size
    EN_initial=EN_initial
    pos_traj=[]
    frame_traj=[]
    EN_traj=[]

    all_steps,positions=gen_init_fiber_par_frame_blank(ncp_bp_length,N_ncp,init_linker_length)
    np.random.seed()
    #n_frames=len(full_data['1kx5_sym']['Frame'].unique())
    
    choice=np.random.choice(np.arange(n_frames_md),size=len(positions))   
    if N_results is None:
        pbar = trange(N_mc_evals,disable=mute)
        for i in pbar:
            np.random.seed()
            choice=np.random.choice(np.arange(n_frames_md),size=len(positions))   
            if linker_step_size==0:
                positions_new = positions
            else:
                pos_step=np.random.randint(-linker_step_size,linker_step_size+1, size = len(positions))
                positions_new = positions+pos_step

            result=gen_fiber(positions_new,choice,all_steps,ncp_bp_frames,bp_num_map,get_dna_pos=True)
            if not result is None:
                fiber,beads,beads_dna=result
                if E_type=='LJ':
                    EN_new=calc_E_w_dna(beads,beads_dna,Rm=Rm, eps=eps,Rm_dna=Rm_dna,eps_dna=eps_dna)
                elif E_type=='logistic':
                    EN_new=calc_E_w_dna_logistic(beads,beads_dna,Rm=Rm, eps=eps,Rm_dna=Rm_dna,eps_dna=eps_dna)

                Del_E = EN_new - EN_initial
                r = np.random.uniform(0,1)

                if Del_E < 0 or (not(np.isinf(np.exp(Del_E))) and (r  < np.exp(-Del_E/KT))):

                    positions=positions_new
                    EN_initial=EN_new
                    pos_traj.append(positions)
                    frame_traj.append(choice)
                    EN_traj.append(EN_new)
                    pbar.set_description(f"Found {len(EN_traj)}, E={EN_traj[-1]:.2f}")
    else:
        i=0
        j=0
        while (i < N_results) and (j<MAX_mc_evals):
            np.random.seed()
            choice=np.random.choice(np.arange(n_frames_md),size=len(positions))   
            if linker_step_size==0:
                positions_new = positions
            else:
                pos_step=np.random.randint(-linker_step_size,linker_step_size+1, size = len(positions))
                positions_new = positions+pos_step

            result=gen_fiber(positions_new,choice,all_steps,ncp_bp_frames,bp_num_map,get_dna_pos=True)
            j+=1
            if not result is None:
                fiber,beads,beads_dna=result
                if E_type=='LJ':
                    EN_new=calc_E_w_dna(beads,beads_dna,Rm=Rm, eps=eps,Rm_dna=Rm_dna,eps_dna=eps_dna)
                elif E_type=='logistic':
                    EN_new=calc_E_w_dna_logistic(beads,beads_dna,Rm=Rm, eps=eps,Rm_dna=Rm_dna,eps_dna=eps_dna)

                Del_E = EN_new - EN_initial
                r = np.random.uniform(0,1)

                if Del_E < 0 or (not(np.isinf(np.exp(Del_E))) and (r  < np.exp(-Del_E/KT))):

                    positions=positions_new
                    EN_initial=EN_new
                    pos_traj.append(positions)
                    frame_traj.append(choice)
                    EN_traj.append(EN_new)
                    i+=1
    return(pos_traj, frame_traj, EN_traj, all_steps)