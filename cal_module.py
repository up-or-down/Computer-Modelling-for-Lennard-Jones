import numpy as np
from particle3D import Particle3D
import lj_utils
import argparse
import math

def compute_separations(particles, box_size):
    '''
    Parameter:
        particles: list,Particles whose separations we want
        box_size:  float,length of one side of the box
        separations: 3d array,MIC vector separation of each pair of particles
    '''
    num = len(particles)
    separations = np.zeros((num,num,3))
    if not isinstance(particles,list):
        particles = [particles]
    if not isinstance(particles[0],Particle3D):
        print("Wrong type of input")
        return None
    half_L = box_size/2
    for i in range(num):
        for j in range(i):
#            if abs(particles[i].position[0]-particles[j].position[0])<2.5 and abs(particles[i].position[1]-particles[j].position[1])<2.5 and abs(particles[i].position[2]-particles[j].position[2])<2.5:
            separations[i,j,0] = (particles[i].position[0] - particles[j].position[0] + half_L)%box_size - half_L
            separations[i,j,1] = (particles[i].position[1] - particles[j].position[1] + half_L)%box_size - half_L
            separations[i,j,2] = (particles[i].position[2] - particles[j].position[2] + half_L)%box_size - half_L
            # else:
            #     separations[i,j,0] = 100
            #     separations[i,j,1] = 100
            #     separations[i,j,2] = 100
            separations[j,i] = -separations[i,j]
    return separations # 

def compute_forces_potential(separations):
    global flag
    shape = separations.shape
    forces = np.zeros((shape[0],3))
    potential = 0
    for i in range(shape[0]):
        for j in range(i):
            r2_ij = np.sum(separations[i,j]**2)
            if r2_ij < 6.25:
                forces[i] += 48*(r2_ij**(-7)-r2_ij**(-4)/2)*separations[i,j]
                if max(abs(forces[i]))>1000 and flag >0:
                    print("force of %d and %d out"%(i,j))
                    print(separations[i,j])
                    flag = 0
                forces[j] -= forces[i]
                potential += r2_ij**(-6)-r2_ij**(-3)

    potential *= 4
    return forces,potential

def generate_initial_conditions(N,rho,T):
    particles = []
    for i in range(N):
        particles.append(Particle3D("p_"+str(i),mass=1))
    box, full = lj_utils.set_initial_positions(rho,particles)
    if full: # 在固体时保持对称型
        lj_utils.set_initial_velocities(T,particles[:4],seed=100)
        for j in range(4,N):
            particles[j].velocity = particles[j%4].velocity # 每四个点为1组
    else:
        lj_utils.set_initial_velocities(T,particles,seed=100)
    return particles,box,full

def velocity_verlet_integrator(particles,forces,box_size,delta_t):
    natom = len(particles)
    delta_t = delta_t/10
    for i in range(natom):
        particles[i].position += particles[i].velocity*delta_t + forces[i]/particles[i].mass*delta_t*delta_t/2
        particles[i].position[0] = particles[i].position[0]%box_size
        particles[i].position[1] = particles[i].position[1]%box_size
        particles[i].position[2] = particles[i].position[2]%box_size
    separations = compute_separations(particles,box_size)
    new_forces, new_energy = compute_forces_potential(separations)
    for i in range(natom):
        particles[i].velocity += (forces[i]/particles[i].mass+new_forces[i]/particles[i].mass)*delta_t/2
    return new_forces, new_energy

# 3.1 Simulation code
def simulation(particles,box_size,time,delta_t,full):
    # compute all the separations between each pair of particles, using MIC,
    if not full: # gas
        natom = len(particles)
        separations = compute_separations(particles,box_size)
        # compute all the LJ forces at once from this array, using the cut-off radius(box_size)
        forces, potential = compute_forces_potential(separations)
        # update all the particle velocities at once, update all the particle positions at once, applying the PBC after they move
        new_force, new_energy = velocity_verlet_integrator(particles,forces,box_size,delta_t)
        positions = []
        velocities = []
        for i in range(natom):
            positions.append(particles[i].position)
            velocities.append(particles[i].velocity)
        positions = np.array(positions)
        velocities = np.array(velocities)
        kinetic = 0
        for i in range(natom):
            kinetic += particles[i].mass*np.sum(particles[i].velocity**2)/2
        total =  kinetic + new_energy
        # generate results as numpy arrays of the body positions and velocities, the time, and the energy across the simulation
        return positions, velocities, time+delta_t, total
    else: # situation for solid
        natom = len(particles)
        ndim = math.ceil((natom/4.0)**(1./3.)) # number of small box in each dimision
        box_size = box_size/ndim  # the smaller L of box in case of solid
        separations = compute_separations(particles[:4],box_size)  
        # compute all the LJ forces at once from this array, using the cut-off radius(box_size)
        forces, potential = compute_forces_potential(separations)
        # update all the particle velocities at once,
        new_force, new_energy = velocity_verlet_integrator(particles[:4],forces,box_size,delta_t)
        for i in range(ndim):
            for j in range(ndim):
                for k in range(ndim):
                    if i+j+k== 0:
                        continue
                    else:
                        index = 4*(i*ndim*ndim+j*ndim+k)
                        particles[index].position = particles[0].position + np.array([i,j,k])*box_size # 这是一个块内
                        particles[index+1].position = particles[1].position + np.array([i,j,k])*box_size
                        particles[index+2].position = particles[2].position + np.array([i,j,k])*box_size
                        particles[index+3].position = particles[3].position + np.array([i,j,k])*box_size
                        particles[index].velocity = particles[0].velocity 
                        particles[index+1].velocity = particles[1].velocity 
                        particles[index+2].velocity = particles[2].velocity 
                        particles[index+3].velocity = particles[3].velocity 
            # 此处计算separation，仍为每4个一组？根据对称性，每4个为一组
        positions = []
        velocities = []
        for i in range(natom):
            positions.append(particles[i].position)
            velocities.append(particles[i].velocity)
        positions = np.array(positions)
        velocities = np.array(velocities)
        kinetic = 0
        for i in range(natom):
            kinetic += particles[i].mass*np.sum(particles[i].velocity**2)/2
        total =  kinetic + new_energy
        # if np.round(time/delta_t) == 21 or np.round(time/delta_t) == 22 or np.round(time/delta_t) == 20 or np.round(time/delta_t) == 19:
        #     print("forces of"+str(np.round(time/delta_t))+"times\n",new_force)
        #     print("velocities of"+str(np.round(time/delta_t))+"times\n",velocities)

        #print(str(time/delta_t)+'  ',separations)
        #print(energy)
        # generate results as numpy arrays of the body positions and velocities, the time, and the energy across the simulation
        return positions, velocities, time+delta_t, total


# add a function to estimate these quantities as a function of time:
#   The kinetic, potential, and total energies of the system,
#   Particles’ mean square displacement (MSD),
# and save all of them to files for later analysis.
def estimate(particles,potential):
    # Make plots of the MSD for your solid and gas simulations, and try to find a temperature where the argon is a liquid
    natom = len(particles)
    

    return MSD, kinetic, potential, total
    pass

def run():
    parser = argparse.ArgumentParser(description="Full simulation")
    parser.add_argument("--dt", type=float, default=0.005, help="step size of simulation")
    parser.add_argument("--Nstep", type=int, default=2000, help="steps of simulation")
    parser.add_argument("--Num","-n", type=int, default=32, help="number of atoms in a box")
    parser.add_argument("--Temp", "-t",type=float, default=1.0, help="temprature of simulation")
    parser.add_argument("--rho", "-d",type=float, default=0.05, help="density of particles in reduced units")
    parser.add_argument("--output", "-o",type=str, default="./out.xyz", help="name of output files")
    args = parser.parse_args()
    time = 0
    
    # # basic test
    # particles, box_size = lj_utils.basic_test()
    # natom = 2
    
    # initial conditions
    particles, box_shape, full = generate_initial_conditions(args.Num,args.rho,args.Temp)
    natom = args.Num
    box_size = box_shape[0]  # get L of the model
    
    all_time = []
    all_energy = []
    with open(args.output,'w') as f:
        for i in range(args.Nstep):
            f.write(str(natom)+'\n')
            positions, velocities, time, energy = simulation(particles,box_size,time,args.dt,full) # do one simulation
            f.write("Point = "+str(i)+'\n')
            for j in range(natom):
                f.write(particles[j].name+' '
                        +str(particles[j].position[0])+' '
                        +str(particles[j].position[1])+' '
                        +str(particles[j].position[2])+'\n')
                #f.write(str(forces[j][0])+' '+str(forces[j][1])+' '+str(forces[j][2])+'\n')
            all_time.append(time)
            all_energy.append(energy)
        f.close()
    estimate(all_time,positions,velocities,)
    np.save("time.npy",all_time)
    np.save("energy.npy",all_energy)

if __name__ == '__main__':
    # default as for gas
    global flag
    flag = -1
    run()