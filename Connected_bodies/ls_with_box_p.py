from builtins import object
from proteus.default_p import *
from proteus.mprans import NCLS
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.opts.T

LevelModelType = NCLS.LevelModel

coefficients = NCLS.Coefficients(V_model=int(ct.movingDomain)+0,
                                 RD_model=int(ct.movingDomain)+3,
                                 ME_model=int(ct.movingDomain)+2,
                                 checkMass=False,
                                 useMetrics=ct.useMetrics,
                                 epsFact=ct.epsFact_consrv_heaviside,
                                 sc_uref=ct.ls_sc_uref,
                                 sc_beta=ct.ls_sc_beta,
                                 movingDomain=ct.movingDomain)

dirichletConditions = {0: lambda x, flag: None}

advectiveFluxBoundaryConditions = {}

diffusiveFluxBoundaryConditions = {0: {}}

class PHI_IC(object):
    def uOfXT(self, x, t):
        for i in range(len(ct.box_start)):
            if x[0]>ct.box_start[i][0] and x[0]<ct.box_end[i][0] and x[1]>ct.box_start[i][1] and x[1]<ct.box_end[i][1] and x[2]>ct.box_start[i][2] and x[2]<ct.box_end[i][2]:
                return x[2] - ct.box_water_level     
            else:
                return x[2] - ct.water_level   


initialConditions = {0: PHI_IC()}
