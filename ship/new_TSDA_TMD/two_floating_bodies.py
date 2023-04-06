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

#class spring_function(pychrono.ForceFunctorP):

#    def __call__(self, d_time, d_rest_length, d_length, d_vel, o_spring_object):
#        """
#        Overrides the spring calculation class with a custom calculation function. This current sets any value over the
#        spring maximum as the maximum to mnimize the the affect of the broken flexor and the simulation.

#        Parameters
#        ----------
#        d_time: float
#            Seconds
#        d_rest_length: float
#            meters
#        d_length: float
#            metters
#        d_vel: float
#            meeers per second
#        o_spring_object: objecti
#            Sping chrono object

#        Returns
#        -------
#        d_force: float
#            netwons

#        """

#        d_spring_constant = 4000.0 # [N/m] Spring force after initial tolerance

#        d_force = -d_spring_constant * (d_length-d_rest_length)

#        return d_force

opts= Context.Options([
    ("final_time",60.0,"Final time for simulation"),
    ("dt_output",0.1,"Time interval to output solution"),
    ("cfl",0.5,"Desired CFL restriction"),
    ("he",0.05,"he relative to Length of domain in x"),
    ("fr",1.0,"Forcing frequency ratio"),
    ("fnx",0.6104,"Natural frequency of sway motion of the main structure"),
    ("fny",0.6104,"Natural frequency of heave motion of the main structure"),
    ("fnz",0.6104,"Natural frequency of roll motion of the main structure"),
    ("mooring",False,"True if the mooring lines are attached"),
    ("ic_angle",0.,"Initial pitch angle of the floating platform (deg)"),
    ])

# general options
# sim time
T = opts.final_time
# initial step
dt_init = 0.001
# CFL value
cfl = opts.cfl
# mesh size
he = opts.he
# rate at which values are recorded
sampleRate = 0.05

# physical configurations
mooring = opts.mooring

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

# forcing frequency (aim at roll motion)
fr = opts.fr
fn = opts.fnz
fc = fr*fn

# wave channel
water_level = 10.
water_length = 30.

# wave options
wave_period = 1./fc
wave_height = 0.1
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

# Space between TLD and main structure
spacing = 1.


# Main structure dimensions
body_w1 = 8.
body_w2 = 1.
body_h1 = 0.5
body_h2 = 1.5

#  --------w1--------
#  |                |
#  |                h1
#  \               /
#   \             / h2
#    \-----w2----/

thob = 400.
mb1 = thob*(body_w1*body_h1+0.5*(body_w1+body_w2)*body_h2)
by = rho_0*0.5*(body_w1+body_w2)*body_h2

if mb1 < by:
    w_temp = (2.*mb1/rho_0/body_h2*(body_w1-body_w2)+body_w2**2)**0.5
    yst = 2.*mb1/rho_0/(w_temp+body_w2)
else:
    yst = (mb1-by)/rho_0/body_w1+body_h2

ic_angle = (opts.ic_angle/180.)*math.pi


# TMD dimension, assume a rectangular box
tld_lx = body_w1
tld_ly = 0.1*tld_lx
tld_tho = 100.

# Design vertical spring based on TMD theory
mb2 = tld_lx*tld_ly*tld_tho
rm = mb2/mb1
ft = 1./(1.+rm)
xi_opt = (3.*rm/8./(1.+rm))**0.5
keq = mb2*(2.*np.pi*fny*ft)**2
ceq = 2.*mb2*(2.*np.pi*fny*ft)
cosa = spacing**2/(spacing**2+(0.5*tld_lx+0.5*body_w1)**2) # square of cosine angle of spring and dashpot
ki = keq/2./cosa
ci = ceq/2./cosa

# TANK
#tank = st.Tank2D(domain, dim=(2*wavelength, 2*water_level))
tank = st.Tank2D(domain, dim=(tank_length, 2.*water_level))

# SPONGE LAYERS
# generation zone: 1 wavelength
# absorption zone: 2 wavelengths
tank.setSponge(x_n=wavelength, x_p=wavelength)

# MAIN STRUCTURE
#caisson1 = st.Rectangle(domain, dim=(fb_lx, fb_ly), coords=(0., 0.))

w05 = 0.5*(body_w1-body_w2)
vertices = np.array([
    [w05,         0.], # 0
    [w05+body_w2, 0.], # 1
    [body_w1,     body_h2], # 2
    [body_w1,     body_h1+body_h2], # 3
    [0.,          body_h1+body_h2], # 4
    [0.,          body_h2], # 5
])

# give flags to vertices (1 flag per vertex, here all have the same flag)
vertexFlags = np.array([1 for ii in range(len(vertices))])
# define segments
segments = np.array([[ii-1, ii] for ii in range(1, len(vertices))])
# add last segment
segments = np.append(segments, [[len(vertices)-1, 0]], axis=0)
# give flags to segments (1 flag per segment, here all have the same flag)
segmentFlags = np.array([1 for ii in range(len(segments))])
# define regions inside the body
regions = np.array([[0.5*body_w1, body_h2]])
regionFlags = np.array([1])
# define holes inside the body
holes = np.array([[0.5*body_w1, body_h2]])
regionFlags = np.array([1])
boundaryTags = {'wall': 1}
# barycenter
arm = body_w1*body_h1*(body_h2+0.5*body_h1)+0.5*(body_w1+body_w2)*body_h2*(body_w2+2.*body_w1)*body_h2/(body_w1+body_w2)/3.
area = body_w1*body_h1+0.5*(body_w1+body_w2)*body_h2
gy = arm/area
barycenter = np.array([0.5*body_w1, gy, 0.])

caisson1 = st.CustomShape(
    domain=domain,
    vertices=vertices,
    vertexFlags=vertexFlags,
    segments=segments,
    segmentFlags=segmentFlags,
    regions=regions,
    regionFlags=regionFlags,
    holes=holes,
    boundaryTags=boundaryTags,
    barycenter=barycenter,
)

# set barycenter in middle of caisson
#caisson.setBarycenter([0., 0.])
# caisson is considered a hole in the mesh
#caisson.setHoles([[0., 0.]])

# 2 following lines only for py2gmsh
caisson1.holes_ind = np.array([0])
tank.setChildShape(caisson1, 0)
# translate caisson to middle of the tank
caisson1.translate(np.array([0.5*water_length-0.5*body_w1, water_level-yst])) # initial motion is getting down
caisson1.rotate(rot = ic_angle)


#DAMPER
caisson2 = st.Rectangle(domain, dim=(tld_lx, tld_ly), coords=(0., 0.))
# set barycenter in middle of caisson
caisson2.setBarycenter([0., 0.])
# caisson is considered a hole in the mesh
caisson2.setHoles([[0., 0.]])
# 2 following lines only for py2gmsh
caisson2.holes_ind = np.array([0])
tank.setChildShape(caisson2, 0)
# translate caisson to middle of the tank
caisson2.translate(np.array([0.5*tank_length, water_level+body_h1+body_h2-yst+spacing+0.5*tld_ly]))



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

# BODY

# create floating body
body = fsi.ProtChBody(system=system)
# give it a name
body.setName(b'my_body1')
# attach shape: this automatically adds a body at the barycenter of the caisson shape
body.attachShape(caisson1)
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
#mb1 = fb_lx*fb_ly*fb_tho
body.setMass(mb1)
# set inertia
# can also be set with:
# body.ChBody.setInertiaXX(pychrono.ChVectorD(1., 1., 0.35))
ib1 = 0.8*mb*(body_w1**2+(body_h1+body_h2)**2)/12. # very rough estimation
body.setInertiaXX(np.array([1., 1., ib1]))
# record values
body.setRecordValues(all_values=True)

# create attached body
body = fsi.ProtChBody(system=system)
# give it a name
body.setName(b'my_body2')
# attach shape: this automatically adds a body at the barycenter of the caisson shape
body.attachShape(caisson2)
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
#mb2 = tld_lx*tld_ly*tld_tho
body.setMass(mb2)
# set inertia
# can also be set with:
# body.ChBody.setInertiaXX(pychrono.ChVectorD(1., 1., 0.35))
ib2 = mb2*(tld_lx**2+tld_ly**2)/12.
body.setInertiaXX(np.array([1., 1., ib2]))
# record values
body.setRecordValues(all_values=True)

# SPRING 1 (left tilted)
TSDA1 = pychrono.ChLinkTSDA()
body1_point = pychrono.ChVectorD(0.5*tank_length-0.5*body_w1,water_level+body_h1+body_h2-yst,0.0)
body2_point = pychrono.ChVectorD(0.5*tank_length+0.5*tld_lx,water_level+body_h1+body_h2-yst+spacing,0.0)

TSDA1.Initialize(system.subcomponents[0].ChBody,
                                    system.subcomponents[1].ChBody,
                                    False, body1_point, body2_point, auto_rest_length=True)
#o_spring_force_functor = spring_function()
#TSDA.RegisterForceFunctor(o_spring_force_functor)
TSDA1.SetSpringCoefficient(ki)
TSDA1.SetDampingCoefficient(ci)
system.ChSystem.Add(TSDA1)


# SPRING 2 (middle vertical)
#TSDA2 = pychrono.ChLinkTSDA()
#body1_point = pychrono.ChVectorD(0.5*tank_length,water_level+body_h1+body_h2-yst,0.0)
#body2_point = pychrono.ChVectorD(0.5*tank_length,water_level+body_h1+body_h2-yst+spacing,0.0)

#TSDA2.Initialize(system.subcomponents[0].ChBody,
#                                    system.subcomponents[1].ChBody,
#                                    False, body1_point, body2_point, auto_rest_length=True)
##o_spring_force_functor = spring_function()
##TSDA.RegisterForceFunctor(o_spring_force_functor)
#TSDA2.SetSpringCoefficient(10000000.)
#TSDA2.SetDampingCoefficient(40.)
#system.ChSystem.Add(TSDA2)

# SPRING 3 (right tilted)
TSDA3 = pychrono.ChLinkTSDA()
body1_point = pychrono.ChVectorD(0.5*tank_length+0.5*body_w1,water_level+body_h1+body_h2-yst,0.0)
body2_point = pychrono.ChVectorD(0.5*tank_length-0.5*tld_lx,water_level+body_h1+body_h2-yst+spacing,0.0)

TSDA3.Initialize(system.subcomponents[0].ChBody,
                                    system.subcomponents[1].ChBody,
                                    False, body1_point, body2_point, auto_rest_length=True)
#o_spring_force_functor = spring_function()
#TSDA.RegisterForceFunctor(o_spring_force_functor)
TSDA3.SetSpringCoefficient(ki)
TSDA3.SetDampingCoefficient(ci)
system.ChSystem.Add(TSDA3)


#  ____                        _                   ____                _ _ _   _
# | __ )  ___  _   _ _ __   __| | __ _ _ __ _   _ / ___|___  _ __   __| (_) |_(_) ___  _ __  ___
# |  _ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | | |   / _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
# | |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | |__| (_) | | | | (_| | | |_| | (_) | | | \__ \
# |____/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |\____\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
#                                           |___/
# Boundary Conditions

# CAISSON

# set no-slip conditions on caisson
for tag, bc in caisson1.BC.items():
    bc.setNoSlip()
for tag, bc in caisson2.BC.items():
    bc.setNoSlip()

# TANK

# atmosphere on top
tank.BC['y+'].setAtmosphere()
# free slip on bottom
tank.BC['y-'].setFreeSlip()
# free slip on the right
tank.BC['x+'].setFreeSlip()
# non material boundaries for sponge interface
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

class P_IC:
    def uOfXT(self, x, t):
        p_L = 0.0
        phi_L = tank.dim[nd-1] - water_level
        phi = x[nd-1] - water_level
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


domain.MeshOptions.genMesh = True
domain.MeshOptions.he = he
mesh_fileprefix = 'mesh'
domain.MeshOptions.setOutputFiles(mesh_fileprefix)
#
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
myTpFlowProblem.SystemNumerics.useSuperlu=True #False


#outputStepping = TpFlow.OutputStepping(
#    final_time=T,
#    dt_init=dt_init,
#    # cfl=opts.cfl,
#    dt_output=sampleRate,
#    nDTout=None,
#    dt_fixed=None,
#)

#myTpFlowProblem = TpFlow.TwoPhaseFlowProblem(
#    ns_model=None,
#    ls_model=None,
#    nd=domain.nd,
#    cfl=cfl,
#    outputStepping=outputStepping,
#    structured=False,
#    he=he,
#    domain=domain,
#    initialConditions=initialConditions,
#    useSuperlu=True
#)


# Necessary for moving domains
#myTpFlowProblem.movingDomain = movingDomain
#params = myTpFlowProblem.Parameters

myTpFlowProblem.SystemPhysics.movingDomain = movingDomain
params = myTpFlowProblem.SystemPhysics

# PHYSICAL PARAMETERS
params['rho_0'] = rho_0  # water
params['rho_1'] = rho_1 # air
params['nu_0'] = nu_0  # water
params['nu_1'] = nu_1  # air
params['gravity'] = g
params['surf_tension_coeff'] = 0.0

#params.physical.densityA = rho_0  # water
#params.physical.densityB = rho_1  # air
#params.physical.kinematicViscosityA = nu_0  # water
#params.physical.kinematicViscosityB = nu_1  # air
#params.physical.gravity = g
#params.physical.surf_tension_coeff = 0.

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


#m = myTpFlowProblem.Parameters.Models

# MODEL PARAMETERS
#ind = -1
# first model is mesh motion (if any)
#if movingDomain:
#    m.moveMeshElastic.index = ind+1
#    ind += 1
# navier-stokes
#m.rans2p.index = ind+1
#ind += 1
# volume of fluid
#m.vof.index = ind+1
#ind += 1
# level set
#m.ncls.index = ind+1
#ind += 1
# redistancing
#m.rdls.index = ind+1
#ind += 1
# mass correction
#m.mcorr.index = ind+1
#ind += 1
# added mass estimation
#if addedMass is True:
#    m.addedMass.index = ind+1
#    ind += 1

# ADD RELAXATION ZONES TO AUXILIARY VARIABLES
#m.rans2p.auxiliaryVariables += domain.auxiliaryVariables['twp']
# ADD SYSTEM TO AUXILIARY VARIABLES
#m.rans2p.auxiliaryVariables += [system]
#m.rans2p.p.coefficients.NONCONSERVATIVE_FORM=0
#if addedMass is True:
#    # passed in added_mass_p.py coefficients
#    m.addedMass.auxiliaryVariables += [system.ProtChAddedMass]
#    max_flag = 0
#    max_flag = max(domain.vertexFlags)
#    max_flag = max(domain.segmentFlags+[max_flag])
#    max_flag = max(domain.facetFlags+[max_flag])
#    flags_rigidbody = np.zeros(max_flag+1, dtype='int32')
#    for s in system.subcomponents:
#        if type(s) is fsi.ProtChBody:
#            for i in s.boundaryFlags:
#                flags_rigidbody[i] = 1
#    m.addedMass.p.coefficients.flags_rigidbody = flags_rigidbody
