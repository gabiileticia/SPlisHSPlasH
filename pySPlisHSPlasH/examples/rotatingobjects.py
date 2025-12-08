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

    boundary2 = sim.getBoundaryModel(2)
    animatedBody2 = boundary2.getRigidBodyObject()
    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()
    
    a = 3.0
    if (t < 0.7):
        animatedBody.setVelocity([ 0, -2.0, 0])
        animatedBody.animate()

        animatedBody2.setVelocity([ 0, -2.0, 0])
        animatedBody2.animate()
    elif (t < 10.0):
        pos = animatedBody.getPosition()
        x, y = pos[0], pos[2]
        animatedBody.setVelocity([-y*a, 0, x*a])
        animatedBody.animate()

        pos = animatedBody2.getPosition()
        x, y = pos[0], pos[2]
        animatedBody2.setVelocity([-y*a, 0, x*a])
        animatedBody2.animate()
    elif (t < 20.0):
        animatedBody.setVelocity([ 0, 5.0, 0])
        animatedBody.animate()

        animatedBody2.setVelocity([ 0, 5.0, 0])
        animatedBody2.animate()
    elif (t < 21.0):
        animatedBody.setVelocity([ 0, 0.0, 0])
        animatedBody.animate()

        animatedBody2.setVelocity([ 0, 0.0, 0])
        animatedBody2.animate()


def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create an imgui simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    base.setTimeStepCB(time_step_callback)

    base.setValueFloat(base.DATA_EXPORT_FPS, 25.0)
    base.setValueString(base.PARTICLE_EXPORT_ATTRIBUTES,"m_is_vortex;m_lambdatwo;angular velocity;m_vorticity_init;m_angularAcceleration;density;m_direction_vort_dev;m_direction_velocity;velocity;m_masses;mass;vorticity_advected;vorticity_corrected_end;stream;vorticity_linear_field;vorticity_derivative;vorticity_dissipation;delta_velocity;m_vorticity_equation;m_a_adv;m_velocity_from_dfsph;m_velocity_corrected_end;m_position_from_dfsph;m_vorticity_rate_laplacian;m_vorticity_rate_gradient;m_gradV_x;m_gradV_y;m_gradV_z")
    base.activateExporter("VTK Exporter", True)
    base.setValueFloat(base.STOP_AT, 50.0) 

    # Get the scene and add objects
    scene = sph.Exec.SceneConfiguration.getCurrent().getScene()
    scene.particleRadius = 0.01

    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 1.5, 0], scale=[3., 3., 3.0], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))
    # scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 0.5, 0], scale=[2.0, 1.0, 2.0], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))

    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/cylinder.obj", translation=[-1.1, 2.0, 0], scale=[0.15, 1.0, 0.15], color=[0.5, 0.5, 0.5, 1.0], axis=[0,0,1], isAnimated=True, isWall=False, mapInvert=False, mapResolution=[25, 25, 25]))
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/cylinder.obj", translation=[1.1, 2.0, 0], scale=[0.15, 1.0, 0.15], color=[0.5, 0.5, 0.5, 1.0], axis=[0,0,1], isAnimated=True, isWall=False, mapInvert=False, mapResolution=[25, 25, 25]))


    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin=[-1.45, 0.0, -1.45], boxMax = [1.45, 0.4, 1.45], mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 2)
    
    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()

