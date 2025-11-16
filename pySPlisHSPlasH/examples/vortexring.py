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
    boundary = sim.getBoundaryModel(1)
    animatedBody = boundary.getRigidBodyObject()
    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()

    
    if (t > 3.0):
        animatedBody.setPosition([0, 3, 0])
        animatedBody.setVelocity([0, 0, 0])
        animatedBody.animate()
    elif (t > 2.0):
        animatedBody.setVelocity([1, 0, 0])
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
    
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 3.0, 0.], scale=[3., 6., 1.5], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))
    # scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 0.5, 0], scale=[2.0, 1.0, 2.0], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))

    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/cylinder.obj", translation=[-1.2, 0.5, 0], scale=[0.15, 0.05, 0.15], color=[0.5, 0.5, 0.5, 1.0], axis=[0,0,1], angle=math.pi/2.0, isAnimated=True, isWall=False, mapInvert=False, mapResolution=[25, 25, 25]))
    
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin=[-1.25, 0.0, -0.70], boxMax = [1.25, 0.7, 0.70], mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 2)
    
    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()

