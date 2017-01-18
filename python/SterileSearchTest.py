import SterileSearchPy as ssp
import numpy as np

# constructing a data path object

dp=ssp.DataPaths

# constructing a steering object

steer=ssp.SteeringParams

# constructing a sterile neutrino parameter object

nus=ssp.SterileNuParams

# constructing a sterilizer object

guy=ssp.Sterilizer(dp,steer,nus)

