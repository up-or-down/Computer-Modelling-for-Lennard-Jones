# Computer-Modelling-for-Lennard-Jones
## Files included:
cal_module.py: the python file for simulating LJ potential, including all of the functions needed for the experiment
particle3D.py: the python file defining the Particle3D class
li_utils.py, plot_xyz.py, test_basic_functions_lj.py : The provided python files 

solid.xyz: XYZ trajectory file for representing simulation run of solid material 
solid.gif: 3d display of solid.xyz
gas.xyz: XYZ trajectory file for representing simulation run of gas material 
gas.gif: 3d display of gas.xyz
out.xyz: XYZ trajectory file for basic test, showing two particles
out.gif: 3d display of out.xyz

## How to run my code:
input: parameters
output: time.npy, kinetic.npy, potential.npy, total.npy and MSD.npy  xxx.xyz(trajectory file), 
1、Solid Argon:Run the code with 32 particles, ρ = 1.0, T = 0.1, dt = 0.01 (in reduced units) for 1000 steps. 
for simulation, use the following command:
    
    python cal_module.py -n 32 -d 1.0 -t 0.1 --dt 0.01 -Nstep 1000 -o solid.xyz

in this case, input of the python file is the parameters, -n indicates number of particles, -d is the density, -t is the temperature, --dt and -Nstep is the simulation parameter, the time step and number of simulation steps.
-o is the name of output file, which is solid.xyz
![solid](https://github.com/up-or-down/Computer-Modelling-for-Lennard-Jones/blob/main/solid.gif)
2、Gaseous Argon:Run the code with 30 particles, ρ = 0.05, T = 1.0, dt = 0.005 (in reduced units) for 2000 steps. 
for simulation, use the following command:

    python cal_module.py -n 30 -d 0.05 -t 1.0 --dt 0.005 -Nstep 2000 -o gas.xyz

in this case, input of the python file is the parameters, -n indicates number of particles, -d is the density, -t is the temperature, --dt and -Nstep is the simulation parameter, the time step and number of simulation steps
-o is the name of output file, which is gas.xyz

3、for observable :
when run cal_module.py, other quantities such as kinetic, potential, total energies and MSD are saved as kinetic.npy, potential.npy, total.npy and MSD.npy separately.

4、for visualization
you can run plot_xyz.py to generate 3d animation of particles using the following command:
    
    python .\plot_xyz.py .\gas.xyz --3d -o gas.gif

in this case, the filename after python .\plot_xyz.py is the input file, --3d means plot in 3D, -o is the name of output file, which is gas.gif. In addition, you can also use -s to denote number of steps per frame of the animation
![gas](https://github.com/up-or-down/Computer-Modelling-for-Lennard-Jones/blob/main/gas.gif)
