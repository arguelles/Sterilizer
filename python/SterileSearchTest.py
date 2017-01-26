import SterileSearchPy as ssp
import numpy as np

# constructing a data path object

dp=ssp.DataPaths()
dp.squids_files_path="/home/carguelles/work/Sterilizer/conventional_fluxes/"
dp.prompt_squids_files_path="/home/carguelles/work/Sterilizer/prompt_fluxes/"

steer=ssp.SteeringParams()
nus=ssp.SterileNuParams()
steer.ReadCompact=False

guy=ssp.Sterilizer(dp,steer,nus)

