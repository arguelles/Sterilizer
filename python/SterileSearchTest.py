import SterileSearchPy as ssp
import numpy as np

# constructing a data path object

dp=ssp.DataPaths()

steer=ssp.SteeringParams()
nus=ssp.SterileNuParams()

dp.squids_files_path="/data/user/bjones/Sterilizer/Sterilizer/test_data/"
steer.ReadCompact=False

guy=ssp.Sterilizer(dp,steer,nus)


