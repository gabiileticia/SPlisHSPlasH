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

# vel_field = np.array([])

def time_step_callback():  
    sim = sph.Simulation.getCurrent()
    boundary = sim.getBoundaryModel(0)
    animatedBody = boundary.getRigidBodyObject()
    fluid = sim.getFluidModel(0)

    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()
    # print(t)
    if (t > 400.0):
        animatedBody.setAngularVelocity([0, 0, 0])
    else:
        animatedBody.setAngularVelocity([0, 0, 0.5])

    # vel = np.array(fluid.getFieldBuffer("velocity"), copy=False)
    # vorticity = np.array(fluid.getFieldBuffer("vorticity_linear_field"), copy=False)
    # for i in range(0,fluid.numActiveParticles()):
    #     vel_field[i,:] = vel[i,:]  

    animatedBody.animate()

def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create an imgui simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    base.setTimeStepCB(time_step_callback)

    base.setValueFloat(base.DATA_EXPORT_FPS, 1.0)
    base.setValueString(base.PARTICLE_EXPORT_ATTRIBUTES, "velocity;vorticity_linear_field;vorticity_current")
    base.activateExporter("VTK Exporter", True)

    # Get the scene and add objects
    scene = sph.Exec.SceneConfiguration.getCurrent().getScene()
    scene.particleRadius = 0.025
    scene.sim2D = True

    # change camera position
    # base.setValueFloat(base.PAUSE_AT, 30.0)
    base.setVec3ValueReal(base.CAMERA_POSITION, [0,0,8])
    base.setVec3ValueReal(base.CAMERA_LOOKAT, [0,0,0])
    base.setValueInt(base.RENDER_WALLS, 1)


    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="C:/Users/gabri/Documents/thesis/SPlisHSPlasH/data/models/sphere.obj", translation=[0., 0., 0.], scale=[0.5, 0.5, 0.5], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25], isAnimated=True))

    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin = [-0.32, -0.32, -0.32], boxMax = [0.32 , 0.32, 0.32], mode=0, initialVelocity=[0.5, 0.0, 0.0]))


    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 0)
    sim.setVec3ValueReal(sim.GRAVITATION, [0,0,0])

    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()


