import numpy as np
import subprocess
import emcee
import time
import sys
import tqdm
import SterileSearchPy as ssp

dp=ssp.DataPaths()
dp.squids_files_path="/home/carguelles/work/Sterilizer/conventional_fluxes/"
dp.prompt_squids_files_path="/home/carguelles/work/Sterilizer/prompt_fluxes/"

steer=ssp.SteeringParams()
steer.ReadCompact=False

nus=ssp.SterileNuParams()

sterilizer=ssp.Sterilizer(dp,steer,nus)

#calculate likelihood from c++
def llhCPP(theta):
    nuisance = ssp.Nuisance()
    #normalization, astroflux, promptflux, crslope, domeff, pik, nunubarr, zc = theta
    normalization, crslope, domeff, pik, nunubarr, zc = theta

    nuisance.normalization=normalization
    nuisance.astroFlux=0.
    nuisance.promptFlux=0.
    nuisance.crSlope=crslope
    nuisance.domEfficiency=domeff/0.9-1.
    nuisance.piKRatio=pik
    nuisance.nuNubarRatio=nunubarr
    nuisance.zenithCorrection=zc

    return -sterilizer.EvalLLH(nuisance)

def lnprior(theta):
    return 0.0

def lnprob(theta):
	lp = lnprior(theta)
	if not np.isfinite(lp):
		return -np.inf
	return lp + llhCPP(theta)

## MCMC business

tt = time.time()
print("Initializing walkers")
ndim = 6
nwalkers = 50

p0_base = [1.,0.,1.,1.,1., 0.]
p0_std = [0.3, 0.05, 0.1, 0.1, 0.1, 0.1]
p0 = np.random.normal(p0_base, p0_std, size=[nwalkers, ndim])

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=10)

print("Running burn-in")

pos, prob, state = sampler.run_mcmc(p0, 500)
sampler.reset()

nsteps = 1000

# sampler.run_mcmc(pos,500) #regular run
for _ in tqdm.tqdm(sampler.sample(pos, iterations=nsteps), total=nsteps):
    pass
print("Time elapsed", time.time()-tt)

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

np.savetxt("sterile_chain.dat",samples)
import corner
#fig = corner.corner(samples, labels=["$log(ReC_{\mu\tau})$", "$log(ImagC_{\mu\tau})$", "$log(C_{\mu\mu})$"])
fig = corner.corner(samples)
fig.savefig("./triangle_full.png")

np.savetxt("sterile_chain_fc.dat",sampler.flatchain)
