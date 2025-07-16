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


def time_step_callback():  
    sim = sph.Simulation.getCurrent()
    boundary = sim.getBoundaryModel(1)
    animatedBody = boundary.getRigidBodyObject()

    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()
    #animatedBody.setAngularVelocity([0, 0, 10.0/180.0*math.pi])
    animatedBody.setVelocity([1., 0, 0])
    animatedBody.animate()

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
    scene.particleRadius = 0.025
    scene.sim2D = True

    # change camera position
    base.setValueFloat(base.PAUSE_AT, 30.0)
    base.setVec3ValueReal(base.CAMERA_POSITION, [0,0,8])
    base.setVec3ValueReal(base.CAMERA_LOOKAT, [0,0,0])
    base.setValueInt(base.RENDER_WALLS, 1)


    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 0.5, 0.], scale=[2., 3., 0.5], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25], isDynamic=False))
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[-14., 1.5, 0.], scale=[30., 1., 0.5], color=[0.1, 0.4, 0.5, 1.0], isAnimated=True, isWall=False, mapInvert=False))

    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin = [-1.0, -1.0, -0.3], boxMax = [1.0 , 1.0, 0.3], mode=0, initialVelocity=[0.0, 0.0, 0.0]))


    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 0)
    sim.setVec3ValueReal(sim.GRAVITATION, [0,0,0])

    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()

