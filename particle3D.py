import numpy as np

class Particle3D(object):
    def __init__(self,name,mass=1.0,position = [0,0,0],velocity=[0,0,0]):
        self.name = name
        self.mass = mass
        self.position = np.array(position)
        self.velocity = np.array(velocity)

    def __str__(self):
        return "Particle %s has mass %.2f, position is [%.2f,%.2f,%.2f], velocity is [%.2f,%.2f,%.2f]" %(
            self.name,self.mass,self.position[0],self.position[1],self.position[2],self.velocity[0],self.velocity[1],self.velocity[2])