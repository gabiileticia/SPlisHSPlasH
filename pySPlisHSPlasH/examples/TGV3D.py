import sys
import os
os.chdir(r"C:\Users\gabri\Documents\thesis\SPlisHSPlasH")
sys.path.append(r"C:\Users\gabri\Documents\thesis\SPlisHSPlasH\build\lib\Release")  # where pysplishsplash.pyd lives


# Set the working directory to the SPlisHSPlasH base path

import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes
import numpy as np
import math
from scipy.spatial.transform import Rotation as R

flag = 0

def time_step_callback():  
    global flag 
    sim = sph.Simulation.getCurrent()
    # boundary = sim.getBoundaryModel(0)
    # animatedBody = boundary.getRigidBodyObject()
    fluid = sim.getFluidModel(0)

    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()

    n = fluid.numActiveParticles()
    if (flag == 0):
        for i in range(n):
            pos = fluid.getPosition(i)  # numpy float32 array [x, y, z]
            x, y = pos[0], pos[2]

            tu = 10
            k = math.pi/2
            u = tu * np.sin(k*x) * np.cos(k*y)
            v = -tu * np.cos(k*x) * np.sin(k*y)
            w = 0.0

            vel = np.array([u, w, v], dtype=np.float32)
            fluid.setVelocity0(i, vel)
            fluid.setVelocity(i, vel)
        flag = 1

def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create an imgui simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    base.setTimeStepCB(time_step_callback)

    base.setValueFloat(base.DATA_EXPORT_FPS, 25.0)
    base.setValueString(base.PARTICLE_EXPORT_ATTRIBUTES, "angular velocity;m_vorticity_init;m_angularAcceleration;density;m_direction_vort_dev;m_direction_velocity;velocity;m_masses;mass;vorticity_advected;vorticity_corrected_end;stream;vorticity_linear_field;vorticity_derivative;vorticity_dissipation;delta_velocity;m_vorticity_equation;m_a_adv;m_velocity_from_dfsph;m_velocity_corrected_end;m_position_from_dfsph;m_vorticity_rate_laplacian;m_vorticity_rate_gradient;m_gradV_x;m_gradV_y;m_gradV_z")
    base.activateExporter("VTK Exporter", True)

    # Get the scene and add objects
    scene = sph.Exec.SceneConfiguration.getCurrent().getScene()
    scene.particleRadius = 0.025
    scene.sim2D = False

    # change camera position
    # base.setValueFloat(base.PAUSE_AT, 30.0)
    base.setVec3ValueReal(base.CAMERA_POSITION, [0,0,8])
    base.setVec3ValueReal(base.CAMERA_LOOKAT, [0,0,0])
    base.setValueInt(base.RENDER_WALLS, 1)
    base.setValueFloat(base.PAUSE_AT, 11.0) 

    # scene.materials.append(Scenes.MaterialData(id='Ghost'))
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 1.5, 0], scale=[4., 3., 4.], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=False, mapResolution=[25, 25, 25], isDynamic=False))

    # scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin = [0, 0, -0.32], boxMax = [6.527 , 6.527, 0.32], mode=0, initialVelocity=[0.0, 0.0, 0.0]))
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin = [-2, 0.0, -2], boxMax = [2 , 1.0, 2], mode=0, initialVelocity=[0.0, 0.0, 0.0]))
    # scene.fluidBlocks.append(Scenes.FluidBlock(id='Ghost', boxMin = [-7, -7, -0.32], boxMax = [-4 , -4, 0.32], mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 0)
    sim.setVec3ValueReal(sim.GRAVITATION, [0,-9.81,0])

    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()