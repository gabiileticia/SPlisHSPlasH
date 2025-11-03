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
    if (flag < 5):
        for i in range(n):
            pos = fluid.getPosition(i)  # numpy float32 array [x, y, z]
            x, y = pos[0], pos[1]
            r = math.sqrt(x*x + y*y)
            tu = (2*math.pi)/(2*math.pi*r)
            w = 0.0
            r0 = 0.2
            if r <=  r0:
                tu = (2*math.pi*r)/(2*math.pi*r0*r0)
                u = - tu * y/ r
                v = tu * x / r
                vel = np.array([u, v, w], dtype=np.float32)
            else:
                tu = (2*math.pi)/(2*math.pi*r)
                u = - tu * y/ r
                v = tu * x / r
                vel = np.array([u, v, w], dtype=np.float32)
            fluid.setVelocity0(i, vel)
            fluid.setVelocity(i, vel)
        flag += 1

def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create an imgui simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    # base.setTimeStepCB(time_step_callback)

    base.setValueFloat(base.DATA_EXPORT_FPS, 10.0)
    base.setValueString(base.PARTICLE_EXPORT_ATTRIBUTES, "m_stream2d;m_stream2d_vol;m_stream3d;vorticity_current;density;m_direction_vort_dev;m_direction_velocity;velocity;mass;vorticity_advected;vorticity_corrected_end;stream;vorticity_linear_field;vorticity_derivative;vorticity_dissipation;delta_velocity;m_vorticity_equation;m_velocity_from_dfsph;m_velocity_corrected_end;m_position_from_dfsph;m_vorticity_rate_laplacian;m_vorticity_rate_gradient;m_gradV_x;m_gradV_y;m_gradV_z")
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
    base.setValueFloat(base.PAUSE_AT, 100.0) 

    # scene.materials.append(Scenes.MaterialData(id='Ghost'))
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 0, 0.], scale=[8., 8., 0.5], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=False, mapResolution=[25, 25, 25], isDynamic=False))

    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin = [-3.95, -3.95, -0.32], boxMax = [3.95 , 3.95, 0.32], mode=0, initialVelocity=[0.0, 0.0, 0.0]))
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin = [-10, -10, -0.32], boxMax = [-4.0 , -4, 0.32], mode=0, initialVelocity=[0.0, 0.0, 0.0]))
    # scene.fluidBlocks.append(Scenes.FluidBlock(id='Ghost', boxMin = [-7, -7, -0.32], boxMax = [-4 , -4, 0.32], mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 0)
    sim.setVec3ValueReal(sim.GRAVITATION, [0,0,0])

    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()