import numpy as np
from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
import proteus.TwoPhaseFlow.TwoPhaseFlowProblem as TpFlow
import proteus.TwoPhaseFlow.utils.Parameters as Parameters
from proteus import WaveTools as wt
import math

# dependencies for FSI
from proteus.mbd import CouplingFSI as fsi
import pychrono

# general options

opts= Context.Options([
    ("T",30.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.5,"Desired CFL restriction"),
    ("he",0.01,"he relative to Length of domain in x"),
    ("fr",1.0,"Forcing frequency ratio"),
    ("fnx",0.6104,"Natural frequency of sway motion"),
    ("fny",0.6104,"Natural frequency of heave motion"),
    ("fnz",0.6104,"Natural frequency of roll motion"),
    ("TLD_type",True,"True if the TLD is on"),
    ("mooring",False,"True if the mooring lines are attached"),
    ("ic_angle",0.,"Initial pitch angle of the floating platform (deg)"),
    ])

T = opts.T
sampleRate = opts.dt_output
dt_init = 0.001
cfl = opts.cfl
he = opts.he
a=1
# for ALE formulation
movingDomain = True
# for added mass stabilization
addedMass = True

# physical options
# water density
rho_0 = 998.2
# water kinematic viscosity
nu_0 = 1.004e-6
# air density
rho_1 = 1.205
# air kinematic viscosity
nu_1 = 1.5e-5
# gravitational acceleration
g = np.array([0., -9.81, 0.])

# body options
fixed = False

# forcing frequency
fr = opts.fr
fn = opts.fnz
fc = fr*fn

# Main structure dimensions
body_w = 1.
body_h = 0.5

# TLD parameters
tld_w = 0.62 # tank length
tld_h = 0.2  # tank height
ht = 0.049   # water depth of TLD
tld_t = 0.05
ic_angle = (opts.ic_angle/180.)*math.pi

# wave channel
water_level = 1.
water_length = 8.

# Regular wave parameters
wave_period = 1/(fc)
wave_height = 0.05
wave_direction = np.array([1., 0., 0.])
wave_type = 'Fenton'  #'Linear'
# number of Fourier coefficients
Nf = 8
wave = wt.MonochromaticWaves(period=wave_period,
                             waveHeight=wave_height,
                             mwl=water_level,
                             depth=water_level,
                             g=g,
                             waveDir=wave_direction,
                             waveType=wave_type,
                             Nf=8)
wavelength = wave.wavelength

#  ____                        _
# |  _ \  ___  _ __ ___   __ _(_)_ __
# | | | |/ _ \| '_ ` _ \ / _` | | '_ \
# | |_| | (_) | | | | | | (_| | | | | |
# |____/ \___/|_| |_| |_|\__,_|_|_| |_|
# Domain
# All geometrical options go here (but not mesh options)

domain = Domain.PlanarStraightLineGraphDomain()

# ----- SHAPES ----- #

# TANK (=wave channel)
tank = st.Tank2D(domain, dim=(water_length, 2.*water_level))
#tank = st.Tank2D(domain, dim=(2.*wavelength, 2.*water_level))

# SPONGE LAYERS
# generation zone: 1 wavelength
# absorption zone: 2 wavelengths
#tank.setSponge(x_n=2.)
tank.setSponge(x_n=wavelength, x_p=wavelength)

# FLOATING STRUCTURE (main structure)
caisson = st.Rectangle(domain, dim=(body_w, body_h), coords=(0., 0.))
# set barycenter in middle of caisson
caisson.setBarycenter([0., 0.])
# caisson is considered a hole in the mesh
caisson.setHoles([[0., 0.]])
# 2 following lines only for py2gmsh
caisson.holes_ind = np.array([0])
tank.setChildShape(caisson, 0)
# translate caisson to middle of the tank
caisson.translate(np.array([0.5*water_length, water_level]))
#caisson.translate(np.array([wavelength, water_level]))

# PMTLD (water tank on the main structure)
if opts.TLD_type:
    tld_dims = (tld_w+2.*tld_t, tld_h+tld_t)
    # define vertices
    tld_vertices = np.array([
        [0., 0.],
        [tld_dims[0], 0.],
        [tld_dims[0], tld_dims[1]],
        [tld_dims[0]-tld_t, tld_dims[1]],
        [tld_dims[0]-tld_t, tld_t],
        [tld_t, tld_t],
        [tld_t, tld_dims[1]],
        [0., tld_dims[1]],
    ])
    # give flags to vertices (1 flag per vertex, here all have the same flag)
    tld_vertexFlags = np.array([1 for ii in range(len(tld_vertices))])
    # define segments
    tld_segments = np.array([[ii-1, ii] for ii in range(1, len(tld_vertices))])
    # add last segment
    tld_segments = np.append(tld_segments, [[len(tld_vertices)-1, 0]], axis=0)
    # give flags to segments (1 flag per segment, here all have the same flag)
    tld_segmentFlags = np.array([1 for ii in range(len(tld_segments))])
    # define regions inside the body
    tld_regions = np.array([[0.5*tld_dims[0],0.5*tld_t]]) 
    tld_regionFlags = np.array([2])
    # define holes inside the body
    tld_holes = np.array([[0.5*tld_dims[0],0.5*tld_t]])
    tld_boundaryTags = {'wall': 1}
    # barycenter
    tld_barycenter = np.array([0.5*tld_dims[0], 0.5*tld_t, 0.])
    pmtld = st.CustomShape(
        domain=domain,
        vertices=tld_vertices,
        vertexFlags=tld_vertexFlags,
        segments=tld_segments,
        segmentFlags=tld_segmentFlags,
        regions=tld_regions,
        regionFlags=tld_regionFlags,
        holes=tld_holes,
        boundaryTags=tld_boundaryTags,
        barycenter=tld_barycenter,
    )
    # translation and rotation of PMTLD
    pmtld.translate(np.array([0.5*tank.dim[0]-0.5*tld_dims[0], water_level+0.5*body_h]))
    pmtld.rotate(rot = ic_angle)

#   ____ _
#  / ___| |__  _ __ ___  _ __   ___
# | |   | '_ \| '__/ _ \| '_ \ / _ \
# | |___| | | | | | (_) | | | | (_) |
#  \____|_| |_|_|  \___/|_| |_|\___/
# Chrono

# SYSTEM

# create system
system = fsi.ProtChSystem()
# access chrono object
chsystem = system.getChronoObject()
# communicate gravity to system
# can also be set with:
# system.ChSystem.Set_G_acc(pychrono.ChVectorD(g[0], g[1], g[2]))
system.setGravitationalAcceleration(g)
# set maximum time step for system
system.setTimeStep(1e-4)

solver = pychrono.ChSolverMINRES()
chsystem.SetSolver(solver)


# BODY-1: MAIN FLOATING STRUCTURE

# create floating body
body = fsi.ProtChBody(system=system)
# give it a name
body.setName(b'my_body')
# attach shape: this automatically adds a body at the barycenter of the caisson shape
body.attachShape(caisson)
# set 2D width (for force calculation)
body.setWidth2D(1.)
# access chrono object
chbody = body.getChronoObject()
# impose constraints
chbody.SetBodyFixed(fixed)
free_x = np.array([1., 1., 0.]) # translational
free_r = np.array([0., 0., 1.]) # rotational
body.setConstraints(free_x=free_x, free_r=free_r)
# access pychrono ChBody
# set mass
# can also be set with:
# body.ChBody.SetMass(14.5)

# set main structure density, mass, and moment of inertia
thob = 500.
mb = thob*body_w*body_h
body.setMass(mb)

mbi = mb*(body_w*body_w+body_h*body_h)/12.
body.setInertiaXX(np.array([1., 1., mbi]))

# set inertia
# can also be set with:
# body.ChBody.setInertiaXX(pychrono.ChVectorD(1., 1., 0.35))
# body.setInertiaXX(np.array([1., 1., 0.35*body.getMass()/14.5]))

# record values
body.setRecordValues(all_values=True)



# BODY-2: TLD

if opts.TLD_type:
    body = fsi.ProtChBody(system=system)
    # give it a name
    body.setName(b'my_tld')
# attach shape: this automatically adds a body at the barycenter of the caisson shape
    body.attachShape(pmtld)
# set 2D width (for force calculation)
    body.setWidth2D(1.)
# access chrono object
    chbody = body.getChronoObject()
# impose constraints
    chbody.SetBodyFixed(fixed)
    free_x = np.array([1., 1., 0.]) # translational
    free_r = np.array([0., 0., 1.]) # rotational
    body.setConstraints(free_x=free_x, free_r=free_r)
# access pychrono ChBody
# set mass
# can also be set with:
# body.ChBody.SetMass(14.5)
# set main structure density, mass, and moment of inertia
    thob = 300.
    m1 = thob*tld_t*(tld_w-2.*tld_t)
    m2 = thob*tld_t*tld_h
    mb = m1+m2
    body.setMass(mb)

    mbi1 = thob*tld_h*tld_w*(4.*tld_h*tld_h+tld_w*tld_w)/12.
    mbi2 = thob*(tld_h-tld_t)*(tld_w-2.*tld_t)*(tld_h*tld_h+tld_w*tld_w)/12.+tho*(tld_h-tld_t)*(tld_w-2.*tld_t)*(0.5*tld)**2
    mbi = mbi1 - mbi2
    body.setInertiaXX(np.array([1., 1., mbi]))


#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

# CAISSON
# set no-slip conditions on caisson
for tag, bc in caisson.BC.items():
    bc.setNoSlip()

# PMTLD
# set no-slip conditions on PMTLD
if opts.TLD_type:
    for tag, bc in pmtld.BC.items():
        bc.setNoSlip()

# TANK
# atmosphere on top
tank.BC['y+'].setAtmosphere()
# free slip on bottom
tank.BC['y-'].setFreeSlip()
# free slip on the right
tank.BC['x+'].setFreeSlip()
# non material boundaries for sponge interface
#tank.BC['x-'].setFreeSlip()
tank.BC['sponge'].setNonMaterial()

# fix in space nodes on the boundaries of the tank
for tag, bc in tank.BC.items():
    bc.setFixedNodes()

# WAVE AND RELAXATION ZONES

smoothing = he*1.5
dragAlpha = 5*2*np.pi/wave_period/(1.004e-6)
tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave,
                                               smoothing=smoothing,
                                               vert_axis=1)
tank.setGenerationZones(x_n=True,
                        waves=wave,
                        smoothing=smoothing,
                        dragAlpha=dragAlpha)
tank.setAbsorptionZones(x_p=True,
                        dragAlpha=dragAlpha)

#  ___       _ _   _       _    ____                _ _ _   _
# |_ _|_ __ (_) |_(_) __ _| |  / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
#  | || '_ \| | __| |/ _` | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
#  | || | | | | |_| | (_| | | | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |___|_| |_|_|\__|_|\__,_|_|  \____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
# Initial Conditions

from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral
smoothing = 1.5*he
nd = domain.nd

x_limit = (0.5*tank.dim[0]-0.5*tld_w, 0.5*tank.dim[0]+0.5*tld_w)
y_limit = (water_level, water_level+tld_h)


class P_IC:
    def uOfXT(self, x, t):
            p_L = 0.0
            phi_L = tank.dim[nd-1] - water_level
            phi = x[nd-1] - water_level
            p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
            return p

class P_IC_TLD:
    def uOfXT(self, x, t):
        if x[0]<=x_limit[0]:
            p_L = 0.0
            phi_L = tank.dim[nd-1] - water_level
            phi = x[nd-1] - water_level
            p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
            return p
        elif x[0]>=x_limit[1]:
            p_L = 0.0
            phi_L = tank.dim[nd-1] - water_level
            phi = x[nd-1] - water_level
            p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
            return p
        else:
            p_L = 0.0
            phi_L = tank.dim[nd-1] - (water_level+tld_h)
            phi = x[nd-1] - (water_level+tld_h)
            p = p_L -g[nd-1]*(rho_0*(phi_L - phi)
                          +(rho_1 -rho_0)*(smoothedHeaviside_integral(smoothing,phi_L)
                                                -smoothedHeaviside_integral(smoothing,phi)))
            return p


class zero:
    def uOfXT(self, x, t):
        return 0.0

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

class VF_IC_TLD:
    def uOfXT(self, x, t):
        if x[0]<=x_limit[0]:
            return smoothedHeaviside(smoothing, x[nd-1]-water_level)
        elif x[0]>=x_limit[1]:
            return smoothedHeaviside(smoothing, x[nd-1]-water_level)
        else:
            return smoothedHeaviside(smoothing, x[nd-1]-(water_level+tld_h)) # fill water in TLD


class PHI_IC:
    def uOfXT(self, x, t):
        return x[nd-1] - water_level

class PHI_IC_TLD:
    def uOfXT(self, x, t):
        if x[0]<=x_limit[0]:
            return x[nd-1] - water_level
        elif x[0]>=x_limit[1]:
            return x[nd-1] - water_level
        else:
            return x[nd-1] - (water_level+tld_h) # fill water in TLD

# instanciating the classes for *_p.py files
if opts.TLD_type:
    initialConditions = {'pressure': P_IC_TLD(),
                         'vel_u': U_IC(),
                         'vel_v': V_IC(),
                         'vel_w': W_IC(),
                         'vof': VF_IC(),
                         'ncls': PHI_IC(),
                         'rdls': PHI_IC()}
else:
    initialConditions = {'pressure': P_IC(),
                         'vel_u': U_IC(),
                         'vel_v': V_IC(),
                         'vel_w': W_IC(),
                         'vof': VF_IC_TLD(),
                         'ncls': PHI_IC_TLD(),
                         'rdls': PHI_IC_TLD()}                         
                         
                         

#  __  __           _        ___        _   _
# |  \/  | ___  ___| |__    / _ \ _ __ | |_(_) ___  _ __  ___
# | |\/| |/ _ \/ __| '_ \  | | | | '_ \| __| |/ _ \| '_ \/ __|
# | |  | |  __/\__ \ | | | | |_| | |_) | |_| | (_) | | | \__ \
# |_|  |_|\___||___/_| |_|  \___/| .__/ \__|_|\___/|_| |_|___/
#                                |_|


domain.MeshOptions.genMesh = True
domain.MeshOptions.he = he
mesh_fileprefix = 'mesh'
domain.MeshOptions.setOutputFiles(mesh_fileprefix)

st.assembleDomain(domain)


#  _   _                           _
# | \ | |_   _ _ __ ___   ___ _ __(_) ___ ___
# |  \| | | | | '_ ` _ \ / _ \ '__| |/ __/ __|
# | |\  | |_| | | | | | |  __/ |  | | (__\__ \
# |_| \_|\__,_|_| |_| |_|\___|_|  |_|\___|___/
# Numerics

myTpFlowProblem = TpFlow.TwoPhaseFlowProblem()

myTpFlowProblem.outputStepping.final_time = T
myTpFlowProblem.outputStepping.dt_output=sampleRate
myTpFlowProblem.outputStepping.dt_init=dt_init

myTpFlowProblem.domain = domain

myTpFlowProblem.SystemNumerics.cfl = cfl
myTpFlowProblem.SystemNumerics.useSuperlu=False
# Necessary for moving domains
myTpFlowProblem.SystemPhysics.movingDomain = movingDomain

params = myTpFlowProblem.SystemPhysics

# PHYSICAL PARAMETERS
params['rho_0'] = rho_0  # water
params['rho_1'] = rho_1 # air
params['nu_0'] = nu_0  # water
params['nu_1'] = nu_1  # air
params['gravity'] = g
params['surf_tension_coeff'] = 0.0

params.addModel(Parameters.ParametersModelMoveMeshElastic,'move')
params.useDefaultModels()
# added mass estimation
if addedMass is True:
    params.addModel(Parameters.ParametersModelAddedMass,'addedMass')

m = myTpFlowProblem.SystemPhysics.modelDict

m['move'].p.initialConditions['hx'] = zero()
m['move'].p.initialConditions['hy'] = zero()
m['flow'].p.initialConditions['p'] = zero()
m['flow'].p.initialConditions['u'] = zero()
m['flow'].p.initialConditions['v'] = zero()
m['vof'].p.initialConditions['vof'] = VF_IC()
m['ncls'].p.initialConditions['phi'] = PHI_IC()
m['rdls'].p.initialConditions['phid'] = PHI_IC()
m['mcorr'].p.initialConditions['phiCorr'] = zero()
m['addedMass'].p.initialConditions['addedMass'] = zero()

# ADD RELAXATION ZONES TO AUXILIARY VARIABLES
m['flow'].auxiliaryVariables += domain.auxiliaryVariables['twp']
# ADD SYSTEM TO AUXILIARY VARIABLES
m['flow'].auxiliaryVariables += [system]
m['flow'].p.coefficients.NONCONSERVATIVE_FORM=0
if addedMass is True:
    # passed in added_mass_p.py coefficients
    m['addedMass'].auxiliaryVariables += [system.ProtChAddedMass]
    max_flag = 0
    max_flag = max(domain.vertexFlags)
    max_flag = max(domain.segmentFlags+[max_flag])
    max_flag = max(domain.facetFlags+[max_flag])
    flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
    for s in system.subcomponents:
        if type(s) is fsi.ProtChBody:
            for i in s.boundaryFlags:
                flags_rigidbody[i] = 1
    m['addedMass'].p.coefficients.flags_rigidbody = flags_rigidbody

