from __future__ import print_function
from __future__ import division
from past.utils import old_div

import numpy as np
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
from proteus import WaveTools as wt
from proteus.mbd.CouplingFSI import ProtChBody, ProtChSystem, ProtChMoorings, ProtChMesh
import py2gmsh

# dependencies for FSI
from proteus.mbd import CouplingFSI as fsi
import pychrono

import pycatenary


from proteus import MeshTools
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
from proteus.default_n import *
import ctypes
import os
import json
import numpy as np
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow




opts= Context.Options([
    ("T",40.0,"Final time for simulation"),
    ("dt_init",0.001,"initial time for simulation"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("nsave", 2, "number of timesteps saved per second"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),

    ("cfl",0.25,"Desired CFL restriction"),
    ("he",1.0,"he relative to Length of domain in x"),
    ("d_interface_he", 4.5, "refinement factore for a-w interface"),
    ("free_slip", True, "set boundaries to free slip"),
    ("Lz", 13.058, "Dimensions in Z direction [m]"),
    ("water_level", 6.096, "water level [m]"),
    ("wave_length", 15.02, "wavelength [m]"),
    ("chrono_dt", 0.0001, "water level [m]"),
    #Wave informations
    ("waves", True, "Enable waves [True/False]"),
    ("wind_speed", -7.72, "None or np.array([value, value, value]) for speed of wind"),
    ("wave_height", 1.14, "Height of waves, meters"),
    ("period", 3.75, "Wave period, seconds"),
    
    #Absorption Zone options
    ("absorption_zones", False, "Enable abs zones [True/False]"),
    ("abs_zone_length", 25.0, "Length of the absorption zones of the domain"),

    #Spring options
    ("springStiffness", 0.0 , "Spring constant, k: F=-kx"),
    ("deform_x", 0.0, "Initial deformation of the spring in the x-direction"),
    ("deform_y", 0.0, "Intial deformation of the springin the y-direction"),

    #Domain length options
    ("box_number", 3, "number of boxes"),
    ("width_multiplier", 3, "width of the domain is the width_multiplier*box_width"),
    ("Lz", 13.058, "Height of domain [m]"),

    #Mooring information
    ("moorings", False, "Option to enable moorings")

    ])

##############################################################################
####################### DIMENSIONS RELATED TO THE STL  #######################
####################### These are used for domain size #######################
####################### and STL translation input      #######################
##############################################################################

box_L = 21.2                 # Length of box, meters
box_W = 9                    # Width of box, meters
box_H = 1.5                  # Height of box, meters

box_shift = box_H/4.0       # Shift of box to get it at the right water level, meters
box_Spacing = 1.4            # Spacing between the boxes, meters




##############################################################################
#######################    MESH SIZING PARAMETERS      #######################
####################### These are used for smoothing   #######################
#######################          and meshing           #######################
##############################################################################

he = opts.he                            # background mesh size
d_interface_he = opts.d_interface_he    # refinement factor for a-w interface
d_min_mesh_size = he/d_interface_he     # smallest mesh size


##############################################################################
#######################                                #######################
#######################         PHYSICAL OPTIONS       #######################
#######################                                #######################
##############################################################################

# water density
rho_0 = 1025.
# water kinematic viscosity
nu_0 = 1.004e-6
# air density
rho_1 = 1.205
# air kinematic viscosity
nu_1 = 1.5e-5
# gravitational acceleration
g = np.array([0., 0., -9.81])


##############################################################################
#######################                                #######################
#######################            WAVE  OPTIONS       #######################
#######################                                #######################
##############################################################################

wave_period = opts.period
wave_height = opts.wave_height
water_level = opts.water_level
omega = 2.*3.14/wave_period
# wave options
if opts.waves is True:
    wave_direction = np.array([1., 0., 0.])
    wave_type = 'Fenton' 
    # number of Fourier coefficients
    Nfourier = 15
    wave = wt.MonochromaticWaves(period=wave_period,
                                 waveHeight=wave_height,
                                 mwl=water_level,
                                 depth=water_level,
                                 g=g,
                                 waveDir=wave_direction,
                                 waveType=wave_type,
                                 Nf=Nfourier)
    wavelength = wave.wavelength



##############################################################################
#######################                                #######################
#######################          DOMAIN  OPTIONS       #######################
#######################                                #######################
##############################################################################

## Setup the tank face metadata ##
domain = Domain.PiecewiseLinearComplexDomain()

Lx = wavelength*4.0+opts.box_number*box_L+(opts.box_number-1)*box_Spacing
Ly = box_W*opts.width_multiplier
tank = st.Tank3D(domain, dim=[Lx,Ly,opts.Lz])


tank.setSponge(x_n=float(wavelength),x_p=float(wavelength*2.0))
boundaryTags = {'y-' : 1,
                'x+' : 2,
                'y+' : 3,
                'x-' : 4,
                'z+' : 5,
                'z-' : 6,
                'sponge' : 7,
               }






### Barge Setup ###
# Specify the input CAD files
s_geom_filename_section = 'box.stl'
## Set the tank physical properties ##
free_x = np.array([1., 1., 1.])  # Translational DOFs
free_r = np.array([1., 1., 1.])  # Rotational DOFs
### Define the setup properties ###
chrono_dt = opts.chrono_dt
sampleRateChrono = 0.0
# Set the domain dictionary
c_region_dictionary = {'water': 0, 'box': 1}

roll_axis = (1., 0., 0.)
pitch_axis = (0., 1., 0.)
heave_axis = (0., 0., 1.)


for i in range(opts.box_number):
    # Open the STL file
    caisson1 = st.ShapeSTL(domain, s_geom_filename_section)
    # Set the regions in the STL object
    caisson1.setRegions([[0.0, 0.0, 0.0]], [c_region_dictionary['box']])  # This is targeted to the initial COG
    caisson1.setHoles([[0.0, 0.0, 0.0]])
    caisson1.holes_ind = [0]
    #Set the barycenter for the object
    caisson1.setBarycenter([0.0,0.0,0.0])
    # Translate the section to be internal to the tank
    caisson1.translate([2*wavelength+box_L/2.0+i*box_L+(i-1)*box_Spacing,box_W*opts.width_multiplier/2.0 , water_level+box_shift])
    # Set the section into the tank as a child object
    tank.setChildShape(caisson1)
    # Set boundary condition for pier section
    for key, bc in caisson1.BC.items():
        bc.setNoSlip()



#   ____ _
#  / ___| |__  _ __ ___  _ __   ___
# | |   | '_ \| '__/ _ \| '_ \ / _ \
# | |___| | | | | | (_) | | | | (_) |
#  \____|_| |_|_|  \___/|_| |_|\___/
# Chrono

# SYSTEM

## Set the chrono properties ##
o_chrono_system = ProtChSystem(sampleRate=sampleRateChrono)  # @@ Accepted the default inputs during the initial conversion
o_chrono_system.ChSystem.Set_G_acc(pychrono.ChVectorD(0., 0., -9.81))  # [m/s^2]
refinement_grading = 1.2
# Setup the Chrono solution properties
o_chrono_system.setTimeStep(chrono_dt)
o_chrono_system.build_kdtree = True
chsystem = o_chrono_system.getChronoObject()
#solver = pychrono.ChSolverMINRES()
#solver = solver = pychrono.ChSolverBB()
solver = pychrono.ChSolverSparseLU()
#solver.SetVerbose(True)
chsystem.SetSolver(solver)
o_chrono_system.ChSystem.SetTimestepperType(pychrono.ChTimestepper.Type_EULER_IMPLICIT_LINEARIZED)


d_section_single_mass = box_H*box_W*box_L*rho_0/4.0     # [kg]
d_section_single_Ixx = 508787.601      # [kg/m^2]
d_section_single_Iyy = 2760531.792       # [kg/m^2]
d_section_single_Izz = 3241817.361      # [kg/m^2]



for i in range(opts.box_number):
    #-------------------------BODY1-----------------------
    o_chrono_body_section = ProtChBody(o_chrono_system)
    o_chrono_body_section.attachShape(domain.shape_list[i+1])
    o_chrono_body_section.setConstraints(free_x=free_x, free_r=free_r)
    # Set the body physical properties
    o_chrono_body_section.ChBody.SetMass(d_section_single_mass)
    o_inertia_section = pychrono.ChVectorD(d_section_single_Ixx, d_section_single_Iyy, d_section_single_Izz)
    o_chrono_body_section.ChBody.SetInertiaXX(o_inertia_section)
    o_chrono_body_section.setRecordValues(all_values=True)

# Fix the first section within the domain
o_chrono_system.subcomponents[opts.box_number-1].ChBody.SetBodyFixed(True)  # Fix first section



CH_C_PI = np.pi
# In most CADs the Y axis is horizontal, but we want it vertical.
# So define a root transformation for rotating all the imported objects.
rotation1 = pychrono.ChQuaternionD()
rotation1.Q_from_AngAxis(-CH_C_PI / 2, pychrono.ChVectorD(1, 0, 0))  # 1: rotate 90° on X axis

for i in range(opts.box_number-1):
    # Define the cylindrical link lock joint
    o_link_lock_section = pychrono.ChLinkLockRevolute()

    # Define the two points give the rotation axis. This is at the pin location half way between the sections
    # Positioning is automatically calculated to the joint.
    o_first_point = pychrono.ChCoordsysD(pychrono.ChVectorD(2*wavelength+box_L/2.0+i*box_L+i*(box_Spacing+box_L)/2.0,
                                                            box_W*opts.width_multiplier/2.0,
                                                            water_level+box_shift-box_H+1.75), rotation1)                                                                                    # Shift y points to be edges of section, and z to be displaced from new top of section.

    ## Initialize the rotational constraint. This is intended to be around the
    o_link_lock_section.Initialize(o_chrono_system.subcomponents[i].ChBody,
                                o_chrono_system.subcomponents[i+1].ChBody,
                                o_first_point)
#
    #mm = pychrono.ChVectorD(0,0,0)
    #mr = pychrono.ChMatrix33D()
    #mr.SetMatr([[1,0,0],[0,0,-1],[0,1,0]])
    #markTransf = pychrono.ChFrameMovingD(mm,mr)
#
    #o_link_lock_section.GetMarker1().ConcatenatePreTransformation(markTransf)
    #o_link_lock_section.GetMarker2().ConcatenatePreTransformation(markTransf)

    # Add the joint to the system
    o_chrono_system.ChSystem.Add(o_link_lock_section)

if opts.moorings==True:
    #-------------MOORINGS------------------------------------------------------------------------------------------------------------------#
    # variables
    # length
    L = 12.# m
    # submerged weight
    w = 0.0778  # kg/m
    # equivalent diameter (chain -> cylinder)
    d = 32.34e-3 # m
    # unstretched cross-sectional area
    A0 = (np.pi*d**2/4)
    # density
    dens = w/A0+rho_0
    # number of elements for cable
    nb_elems = 20
    # Young's modulus
    #E = (753.6e6)/50**3/A0
    E = 5.44e10

    # fairleads coordinates
    #fairlead_center = np.array([opts.Lx/2, opts.Ly/2, water_level - 0.045 + caisson_dim[2]/2])
    #fairlead1 = fairlead_center+np.array([0., (caisson_dim[1] + caisson_dim[1]/2 + d_section_side), 0.])
    #fairlead2 = fairlead_center+np.array([0., -(caisson_dim[1] + caisson_dim[1]/2 + d_section_side), 0.])
    # anchors coordinates
    anchor1 = np.array([opts.Lx/2, fairlead1[1]+1., 0.])
    anchor2 = np.array([opts.Lx/2, fairlead2[1]-1., 0.])

    # quasi-statics for finding shape of cable
    from pycatenary.cable import MooringLine
    # create lines
    EA = E*A0
    cat1 = MooringLine(L=L,
                    w=w*9.81,
                    EA=EA,
                    anchor=anchor1,
                    fairlead=fairlead1,
                    nd=2,
                    floor=True)
    cat2 = MooringLine(L=L,
                    w=w*9.81,
                    EA=EA,
                    anchor=anchor2,
                    fairlead=fairlead2,
                    nd=2,
                    floor=True)
    cat1.computeSolution()
    cat2.computeSolution()

    # ANCHOR1
    # arbitrary body fixed in space
    body1 = fsi.ProtChBody(o_chrono_system)
    body1.barycenter0 = np.zeros(3)
    # fix anchor in space
    body1.ChBody.SetBodyFixed(True)

    # ANCHOR2
    # arbitrary body fixed in space
    body2 = fsi.ProtChBody(o_chrono_system)
    body2.barycenter0 = np.zeros(3)
    # fix anchor in space
    body2.ChBody.SetBodyFixed(True)

    # MESH
    # initialize mesh that will be used for cables
    #mesh = fsi.ProtChMesh(o_chrono_system)

    # FEM CABLES
    # moorings line 1
    #m1 = fsi.ProtChMoorings(system=o_chrono_system,
    #                        mesh=mesh,
    #                        length=np.array([L]),
    #                        nb_elems=np.array([nb_elems]),
    #                        d=np.array([d]),
    #                        rho=np.array([dens]),
    #                        E=np.array([E]))
    #m1.setName(b'mooring1')
    # send position functions from catenary to FEM cable
    #m1.setNodesPositionFunction(cat1.s2xyz, cat1.ds2xyz)
    # sets node positions of the cable
    #m1.setNodesPosition()
    # build cable
    #m1.buildNodes()
    # apply external forces
    #m1.setApplyDrag(True)
    #m1.setApplyBuoyancy(True)
    #m1.setApplyAddedMass(True)
    # set fluid density at cable nodes
    #m1.setFluidDensityAtNodes(np.array([rho_0 for i in range(m1.nodes_nb)]))
    # sets drag coefficients
    #m1.setDragCoefficients(tangential=1.15, normal=0.213, segment_nb=0)
    # sets added mass coefficients
    #m1.setAddedMassCoefficients(tangential=0.269, normal=0.865, segment_nb=0)
    # small Iyy for bending
    #m1.setIyy(0., 0)
    # attach back node of cable to body
    #m1.attachBackNodeToBody(o_chrono_system.subcomponents[3]) #**********************************************
    # attach front node to anchor
    #m1.attachFrontNodeToBody(body1)

    # mooring line 2
    #m2 = fsi.ProtChMoorings(system=o_chrono_system,
    #                        mesh=mesh,
    #                        length=np.array([L]),
    #                        nb_elems=np.array([nb_elems]),
    #                        d=np.array([d]),
    #                        rho=np.array([dens]),
    #                        E=np.array([E]))
    #m2.setName(b'mooring2')
    # send position functions from catenary to FEM cable
    #m2.setNodesPositionFunction(cat2.s2xyz, cat2.ds2xyz)
    # sets node positions of the cable
    #m2.setNodesPosition()
    # build cable
    #m2.buildNodes()
    # apply external forces
    #m2.setApplyDrag(True)
    #m2.setApplyBuoyancy(True)
    #m2.setApplyAddedMass(True)
    # set fluid density at cable nodes
    #m2.setFluidDensityAtNodes(np.array([rho_0 for i in range(m2.nodes_nb)]))
    # sets drag coefficients
    #m2.setDragCoefficients(tangential=1.15, normal=0.213, segment_nb=0)
    # sets added mass coefficients
    #m2.setAddedMassCoefficients(tangential=0.269, normal=0.865, segment_nb=0)
    # small Iyy for bending
    #m2.setIyy(0., 0)
    # attach back node of cable to body
    #m2.attachBackNodeToBody(o_chrono_system.subcomponents[4]) #**********************************************
    # attach front node to anchor
    #m2.attachFrontNodeToBody(body2)

    # SEABED
    # create a box
    #seabed = pychrono.ChBodyEasyBox(100., 0.2, 1., 1000, True)
    # move box
    #seabed.SetPos(pychrono.ChVectorD(0., -0.1-d*2, 0.))
    # fix boxed in space
    #seabed.SetBodyFixed(True)
    # add box to system
    #o_chrono_system.ChSystem.Add(seabed)

    # CONTACT MATERIAL
    # define contact material for collision detection
    material = pychrono.ChMaterialSurfaceSMC()
    material.SetKn(3e6)  # normal stiffness
    material.SetGn(1.)  # normal damping coefficient
    material.SetFriction(0.3)
    material.SetRestitution(0.2)
    material.SetAdhesion(0)

    # add material to objects
    #seabed.SetMaterialSurface(material)
    #m1.setContactMaterial(material)
    #m2.setContactMaterial(material)


#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

# TANK
smoothing = 1.5*he/d_interface_he
dragAlpha = 5*omega/(1e-6)

tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=2, wind_speed=np.array([opts.wind_speed, 0., 0.]))  # Setting wavetools boundary

if opts.waves is True:
    tank.setGenerationZones(x_n=True,
                        waves=wave,
                        smoothing=smoothing,
                        dragAlpha=dragAlpha)

    tank.setAbsorptionZones(x_p=True,
                        dragAlpha=dragAlpha)

# open top
tank.BC['z+'].setAtmosphere()
if opts.free_slip:
    tank.BC['z-'].setFreeSlip()
    tank.BC['y-'].setFreeSlip()
    tank.BC['y+'].setFreeSlip()
    tank.BC['x-'].setFreeSlip()
    tank.BC['x+'].setFreeSlip()
else:  # no slip
    tank.BC['z-'].setNoSlip()
    tank.BC['x-'].setNoSlip()
    tank.BC['y-'].setNoSlip()
    tank.BC['y+'].setNoSlip()
    tank.BC['x+'].setNoSlip()
for bc in tank.BC_list:
    bc.setFixedNodes()

tank.BC['sponge'].setNonMaterial()
tank.BC['sponge'].setFixedNodes()
tank.BC['y+'].setFixedNodes()#setTank() #Allows nodes to slip freely
tank.BC['y-'].setFixedNodes()#setTank() #Allows nodes to slip freely
tank.BC['z+'].setFixedNodes()#setTank() #Allows nodes to slip freely
tank.BC['z-'].setFixedNodes()#setTank() #Allows nodes to slip freely

if opts.absorption_zones is True:
    tank.setAbsorptionZones(x_n=True,
                        dragAlpha=dragAlpha)
    tank.setAbsorptionZones(x_p=True,
                        dragAlpha=dragAlpha)
    tank.setAbsorptionZones(y_n=True,
                        dragAlpha=dragAlpha)
    tank.setAbsorptionZones(y_p=True,
                        dragAlpha=dragAlpha)

#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
nd = domain.nd

class P_IC:
    def uOfXT(self, x, t):
        p_L = 0.0
        phi_L = tank.dim[nd-1] - water_level
        phi = x[nd-1] - water_level
        p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
        return p
class U_IC:
    def uOfXT(self, x, t):
        return 0.0
class V_IC:
    def uOfXT(self, x, t):
        return 0.0
class W_IC:
    def uOfXT(self, x, t):
        return 0.0
class VF_IC:
    def uOfXT(self, x, t):
        return smoothedHeaviside(smoothing, x[nd-1]-water_level)
class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - water_level

# instanciating the classes for *_p.py files
initialConditions = {'pressure': P_IC(),
                     'vel_u': U_IC(),
                     'vel_v': V_IC(),
                     'vel_w': W_IC(),
                     'vof': VF_IC(),
                     'ncls': PHI_IC(),
                     'rdls': PHI_IC()}

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|

### Setup the mesh of the system ###
# Set the mesh options
mesh_fileprefix = 'mesh_test'
domain.use_gmsh = True
domain.MeshOptions.he = opts.he
domain.MeshOptions.setTriangleOptions()
domain.MeshOptions.genMesh = False
domain.MeshOptions.use_gmsh = domain.use_gmsh
domain.MeshOptions.setOutputFiles(name=mesh_fileprefix)
domain.MeshOptions.setElementSize(he=opts.he)
domain.MeshOptions.setTriangleOptions()
domain.geofile = mesh_fileprefix
# Assemble the domain
st.assembleDomain(domain)
## Create the mesh ##
# Create the mesh
mesh = py2gmsh.geometry2mesh(domain)

# create field for air water interface at
f1 = py2gmsh.Field.Box(mesh=mesh)
f1.VIn = opts.he/opts.d_interface_he
f1.VOut = opts.he
f1.XMax = Lx+wavelength*10.0
f1.XMin = 0.0-wavelength*10.0
f1.YMax = Ly+opts.abs_zone_length
f1.YMin = 0.0-opts.abs_zone_length
f1.ZMax = water_level+opts.he/opts.d_interface_he*1.75
f1.ZMin = water_level-opts.he/opts.d_interface_he*1.75
f1.Thickness = opts.Lz/2.

fmin = py2gmsh.Field.Min(mesh=mesh)
fmin.FieldsList = [f1]

#create a counter to store the number of fields
i_gmshFieldCounter = 2

# set the background field as minimum field
mesh.setBackgroundField(fmin)
# set max element size
mesh.Options.Mesh.CharacteristicLengthMax = opts.he
mesh.Coherence = True
#mesh.Options.Mesh.OptimizeNetgen = 1
max_flag = max(domain.vertexFlags + domain.segmentFlags + domain.facetFlags)
mesh.writeGeo(mesh_fileprefix + '.geo')

flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
for s in o_chrono_system.subcomponents:
    if type(s) is ProtChBody:
        for i in s.boundaryFlags:
            flags_rigidbody[i] = 1

#####################################################
###  _   _                           _            ###
### | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___   ###
### |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|  ###
### | |\  | |_| | | | | | |  __/ |  | | (__\__ \  ###
### |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/  ###
#####################################################


   
# todo: re-enable everything after this
addedMass = True
movingDomain = True
checkMass = False
applyCorrection = True
applyRedistancing = True
freezeLevelSet = True
# ----------------------------------------------------
# Time stepping and velocity
# ----------------------------------------------------
weak_bc_penalty_constant = 100.0 #what is this doing??? Maybe I should delete this... same with below -Daniel B
nDTout = int(opts.T * opts.nsave)
timeIntegration = opts.timeIntegration
if nDTout > 0:
    dt_out = (opts.T - opts.dt_init) / nDTout
else:
    dt_out = 0
runCFL = opts.cfl
dt_fixed = opts.dt_fixed



# ----------------------------------------------------
#  Discretization -- input options
useOldPETSc = False
useSuperlu = not True
spaceOrder = 1
useHex = False
useRBLES = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False

ns_closure = 2
# 1 -- K-Epsilon
# 2 -- K-Omega, 1998
# 3 -- K-Omega, 1988
# Input checks
if spaceOrder not in [1, 2]:
    print("INVALID: spaceOrder" + spaceOrder)
    sys.exit()
if useRBLES not in [0.0, 1.0]:
    print("INVALID: useRBLES" + useRBLES)
    sys.exit()
if useMetrics not in [0.0, 1.0]:
    print("INVALID: useMetrics")
    sys.exit()
#  Discretization
nd = domain.nd
if spaceOrder == 1:
    hFactor = 1.0
    if useHex:
        basis = C0_AffineLinearOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 3)
    else:
        basis = C0_AffineLinearOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 3)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 3)
        # elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor = 0.5
    if useHex:
        basis = C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd - 1, 4)
    else:
        basis = C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd, 4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd - 1, 4)
sc = 0.5  # default: 0.5. Test: 0.25
sc_beta = 1.5  # default: 1.5. Test: 1.
epsFact_consrv_diffusion = 1.  # default: 1.0. Test: 0.1
ns_forceStrongDirichlet = False
backgroundDiffusionFactor = 0.01
sc = 0.5
sc_beta = 1.5
if useMetrics:
    ns_shockCapturingFactor = sc
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = sc
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = sc_beta
    vof_shockCapturingFactor = sc
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = sc_beta
    rd_shockCapturingFactor = sc
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.
    redist_Newton = True  # False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = True  # False
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = True  # False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = sc_beta
else:
    ns_shockCapturingFactor = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref = 1.0
    ls_sc_beta = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = 1.0
    rd_shockCapturingFactor = 0.9
    rd_lag_shockCapturing = False
    epsFact_density = 1.5
    epsFact_viscosity = epsFact_curvature = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False  # True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True  # False
    kappa_sc_uref = 1.0
    kappa_sc_beta = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True  # False
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = 1.0
tolfac = 0.0005
ns_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size ** 2)
vof_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size ** 2)
ls_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size ** 2)
mcorr_nl_atol_res = max(1.0e-6, 0.1 * tolfac * d_min_mesh_size ** 2)
rd_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size)
kappa_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size ** 2)
dissipation_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size ** 2)
mesh_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size ** 2)
am_nl_atol_res = max(1.0e-6, tolfac * d_min_mesh_size ** 2)
# turbulence
#nsclosure: 1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega
useRANS =0
ns_closure=2
if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4
def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = opts.Lz - water_level
    phi = x[nd - 1] - water_level
    return p_L - g[nd - 1] * (rho_0 * (phi_L - phi) + (rho_1 - rho_0) * (
                smoothedHeaviside_integral(epsFact_consrv_heaviside * d_min_mesh_size, phi_L)
                - smoothedHeaviside_integral(epsFact_consrv_heaviside * d_min_mesh_size, phi)))
