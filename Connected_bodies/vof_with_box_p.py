from builtins import object
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.mprans import VOF
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.opts.T

LevelModelType = VOF.LevelModel
if ct.useOnlyVF:
    RD_model = None
    LS_model = None
else:
    RD_model = 3
    LS_model = 2

coefficients = VOF.Coefficients(LS_model=int(ct.movingDomain)+LS_model,
                                V_model=int(ct.movingDomain)+0,
                                RD_model=int(ct.movingDomain)+RD_model,
                                ME_model=int(ct.movingDomain)+1,
                                checkMass=True,
                                useMetrics=ct.useMetrics,
                                epsFact=ct.epsFact_vof,
                                sc_uref=ct.vof_sc_uref,
                                sc_beta=ct.vof_sc_beta,
                                movingDomain=ct.movingDomain)

dirichletConditions = {0: lambda x, flag: domain.bc[flag].vof_dirichlet.init_cython()}

advectiveFluxBoundaryConditions = {0: lambda x, flag: domain.bc[flag].vof_advective.init_cython()}

diffusiveFluxBoundaryConditions = {0: {}}

class VF_IC(object):
    def uOfXT(self, x, t):
        for i in range(len(ct.box_start)):
            if x[0]>ct.box_start[i][0] and x[0]<ct.box_end[i][0] and x[1]>ct.box_start[i][1] and x[1]<ct.box_end[i][1] and x[2]>ct.box_start[i][2] and x[2]<ct.box_end[i][2]:
                return smoothedHeaviside(ct.epsFact_consrv_heaviside*ct.opts.he/2.0/ct.opts.d_interface_he,x[2]-ct.box_water_level)
            else:
                return smoothedHeaviside(ct.epsFact_consrv_heaviside*ct.opts.he/2.0/ct.opts.d_interface_he,x[2]-ct.water_level)   

initialConditions = {0: VF_IC()}
