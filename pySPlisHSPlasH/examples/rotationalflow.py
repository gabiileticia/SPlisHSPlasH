import sys
import os
os.chdir(r"C:\Users\gabri\Documents\thesis\SPlisHSPlasH")
sys.path.append(r"C:\Users\gabri\Documents\thesis\SPlisHSPlasH\build\lib\Release")  # where pysplishsplash.pyd lives

import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes
import numpy as np
import math
from scipy.spatial.transform import Rotation as R

def time_step_callback():  
    sim = sph.Simulation.getCurrent()
    # boundary = sim.getBoundaryModel(1)
    # animatedBody = boundary.getRigidBodyObject()
    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()

    a = 2
    fluid = sim.getFluidModel(0)
    n = fluid.numActiveParticles()
    if (t > 2.0 and t < 3.0):
        for i in range(n):
            pos = fluid.getPosition(i)  # numpy float32 array [x, y, z]
            x, y = pos[0], pos[2]
            
            fluid.setVelocity(i, [-y*a, 0, x*a])

def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create an imgui simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    base.setTimeStepCB(time_step_callback)

    # Get the scene and add objects
    scene = sph.Exec.SceneConfiguration.getCurrent().getScene()
    
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/cylinder.obj", translation=[0., 1.0, 0.], scale=[1., 2., 1.0], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))   
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin=[-0.7, 0.0, -0.7], boxMax = [0.7, 0.8, 0.7], mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    base.setValueFloat(base.DATA_EXPORT_FPS, 25.0)
    base.setValueString(base.PARTICLE_EXPORT_ATTRIBUTES,"m_is_vortex;m_lambdatwo;angular velocity;m_vorticity_init;m_angularAcceleration;density;m_direction_vort_dev;m_direction_velocity;velocity;m_masses;mass;vorticity_advected;vorticity_corrected_end;stream;vorticity_linear_field;vorticity_derivative;vorticity_dissipation;delta_velocity;m_vorticity_equation;m_a_adv;m_velocity_from_dfsph;m_velocity_corrected_end;m_position_from_dfsph;m_vorticity_rate_laplacian;m_vorticity_rate_gradient;m_gradV_x;m_gradV_y;m_gradV_z")
    base.activateExporter("VTK Exporter", True)
    base.setValueFloat(base.PAUSE_AT, 7.0) 
    
    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 2)
    sim.setVec3ValueReal(sim.GRAVITATION, [0,-9.81,0])
    
    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()

