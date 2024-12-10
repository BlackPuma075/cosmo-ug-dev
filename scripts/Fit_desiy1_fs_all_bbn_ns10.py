#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib.pyplot as plt
from desilike.theories.galaxy_clustering import FOLPSAXTracerPowerSpectrumMultipoles, DirectPowerSpectrumTemplate
from desilike.observables.galaxy_clustering import TracerPowerSpectrumMultipolesObservable 
from desilike.observables import ObservableCovariance
from desilike.emulators import EmulatedCalculator, Emulator, TaylorEmulatorEngine
from desilike.likelihoods import ObservablesGaussianLikelihood, SumLikelihood
from desilike.theories import Cosmoprimo
from cosmoprimo.fiducial import DESI
from desilike import setup_logging

######### Settings ######

#model: LCDM or HS
model = 'LCDM'

#Biasing and EFT parametrization: 'physical' or 'default' (non-physical)
prior_basis = 'physical' #Prior to be used 

k_max = 0.20

#width for EFT and SN paramss
width_EFT = 12.5
width_SN0 = 2.0
width_SN2 = 5.0

sampler = 'cobaya'


# List of tracers
tracers = ['BGS', 'LRG1', 'LRG2', 'LRG3', 'ELG', 'QSO']  # Add more tracers as needed
chain_name = f'Chains/ug_desiy1_rept_fs-all_schoneberg2024-bbn_planck2018-ns10_physprior_LCDM_'f'kmax{str(k_max).replace(".", "p")}' #fn to save the chain from the sampler

# Define file paths for each tracer (in the same order as the tracers list above)
tracer_params = {
    0: {
        'data_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/corrected/pkpoles_corrected_BGS_BRIGHT-21.5_GCcomb_z0.1-0.4_thetacut0.05.npy',
        'wmatrix_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/wmatrix_smooth_BGS_BRIGHT-21.5_GCcomb_z0.1-0.4_thetacut0.05.npy',
        'covariance_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/cov_2pt/ezmock/v1/covariance_power_BGS_BRIGHT-21.5_GCcomb_z0.1-0.4_default_FKP_lin_thetacut0.05.npy'
    },
    1: {
        'data_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/corrected/pkpoles_corrected_LRG_GCcomb_z0.4-0.6_thetacut0.05.npy',
        'wmatrix_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/wmatrix_smooth_LRG_GCcomb_z0.4-0.6_thetacut0.05.npy',
        'covariance_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/cov_2pt/ezmock/v1/covariance_power_LRG_GCcomb_z0.4-0.6_default_FKP_lin_thetacut0.05.npy'
    },
    2: {
        'data_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/corrected/pkpoles_corrected_LRG_GCcomb_z0.6-0.8_thetacut0.05.npy',
        'wmatrix_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/wmatrix_smooth_LRG_GCcomb_z0.6-0.8_thetacut0.05.npy',
        'covariance_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/cov_2pt/ezmock/v1/covariance_power_LRG_GCcomb_z0.6-0.8_default_FKP_lin_thetacut0.05.npy'
    },
    3: {
        'data_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/corrected/pkpoles_corrected_LRG_GCcomb_z0.8-1.1_thetacut0.05.npy',
        'wmatrix_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/wmatrix_smooth_LRG_GCcomb_z0.8-1.1_thetacut0.05.npy',
        'covariance_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/cov_2pt/ezmock/v1/covariance_power_LRG_GCcomb_z0.8-1.1_default_FKP_lin_thetacut0.05.npy'
    },
    4: {
       'data_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/corrected/pkpoles_corrected_ELG_LOPnotqso_GCcomb_z1.1-1.6_thetacut0.05.npy',
        'wmatrix_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/wmatrix_smooth_ELG_LOPnotqso_GCcomb_z1.1-1.6_thetacut0.05.npy',
        'covariance_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/cov_2pt/ezmock/v1/covariance_power_ELG_LOPnotqso_GCcomb_z1.1-1.6_default_FKP_lin_thetacut0.05.npy'
    },
    5: {
        'data_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/corrected/pkpoles_corrected_QSO_GCcomb_z0.8-2.1_thetacut0.05.npy',
        'wmatrix_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/baseline_2pt/pk/wmatrix_smooth_QSO_GCcomb_z0.8-2.1_thetacut0.05.npy',
        'covariance_fn': '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/unblinded/desipipe/cov_2pt/ezmock/v1/covariance_power_QSO_GCcomb_z0.8-2.1_default_FKP_lin_thetacut0.05.npy'
    }
}    



#Define a cosmology to get sigma_8, Omega_m and fR0
cosmo = Cosmoprimo(engine='class')
cosmo.init.params['H0'] = dict(derived=True)
cosmo.init.params['Omega_m'] = dict(derived=True)
cosmo.init.params['sigma8_m'] = dict(derived=True)
fiducial = DESI() #fiducial cosmology


#Define tracer types and their corresponding redshifts
tracer_redshifts = {
    'LRG1': 0.510, 'LRG1': 0.510,
    'LRG2': 0.706, 'LRG2': 0.706,
    'LRG3': 0.930, 'LRG3': 0.930,
    'QSO': 1.491, 'QSO': 1.491,
    'ELG': 1.317, 'ELG': 1.317,
    'BGS': 0.295, 'BGS': 0.295
}

#Initialize an empty list to store the theory objects
theories = []

#Iterate over each tracer and create the corresponding theory object
for tracer in tracers:
    if tracer in tracer_redshifts:
        z = tracer_redshifts[tracer]
    else:
        print(f'Invalid tracer: {tracer}. Skipping.') 
        continue

    #Create the template and theory objects
    template = DirectPowerSpectrumTemplate(fiducial = fiducial, z=z) #cosmology and fiducial cosmology defined above
    
    theory = FOLPSAXTracerPowerSpectrumMultipoles(template=template, prior_basis = prior_basis) #Add the prior_basis='physical' argument to use physically motivated priors

    #Update cosmo priors
    for param in ['n_s', 'h','omega_cdm', 'omega_b', 'logA', 'tau_reio']:
        template.params[param].update(fixed = False)
        if param == 'tau_reio':
            template.params[param].update(fixed = True)
        if param == 'n_s':
            #cosmo.params[param].update(fixed = True)
            template.params[param].update(prior={'dist': 'norm', 'loc': 0.9649, 'scale': 0.042})
        if param == 'omega_b':
            template.params[param].update(prior={'dist': 'norm', 'loc': 0.02218, 'scale': 0.00055})
        if param == 'h':
            template.params[param].update(prior = {'dist':'uniform','limits': [0.5,0.9]})
        if param == 'omega_cdm':
            template.params[param].update(prior = {'dist':'uniform','limits': [0.05, 0.2]})
        if param == 'logA':
            template.params[param].update(prior = {'dist':'uniform','limits': [2.0, 4.0]})
    
    #Update bias and EFT priors
    #theory.params['bs'].update(fixed=True)
    if prior_basis == 'physical':
        theory.params['b1p'].update(prior = {'dist':'uniform','limits': [1e-5, 10]})
        theory.params['b2p'].update(prior = {'dist':'uniform','limits': [-50, 50]})
        theory.params['bsp'].update(prior = {'dist': 'norm', 'loc': 0, 'scale': 20})
        theory.params['alpha0p'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_EFT})
        theory.params['alpha2p'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_EFT})
        theory.params['alpha4p'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_EFT})
        theory.params['sn0p'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_SN0})
        theory.params['sn2p'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_SN2})
    else:
        theory.params['b1'].update(prior = {'dist':'uniform','limits': [1e-5, 10]})
        theory.params['b2'].update(prior = {'dist':'uniform','limits': [-50, 50]})
        theory.params['bs'].update(prior = {'dist': 'norm', 'loc': 0, 'scale': 20})
        theory.params['alpha0'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_EFT})
        theory.params['alpha2'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_EFT})
        theory.params['alpha4'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_EFT})
        theory.params['sn0'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_SN0})
        theory.params['sn2'].update(prior={'dist': 'norm', 'loc': 0, 'scale': width_SN2})
    
    #Append the theory object to the list
    theories.append(theory)

#Print used priors
for param in theory.all_params:
    print(param,':',theory.all_params[param].prior)


# In[ ]:


#Define a function to create an observable
def create_observable(data_fn, wmatrix_fn, covariance_fn, theory, tracer_index):
    #Load and process covariance
    covariance = ObservableCovariance.load(covariance_fn)
    covariance = covariance.select(xlim=(0.02, k_max), projs=[0, 2])
    
    #Create and return the observable
    return TracerPowerSpectrumMultipolesObservable(
        data=data_fn,
        covariance=covariance,
        klim={ell: [0.02, k_max, 0.005] for ell in [0, 2]},
        theory=theories[tracer_index],
        wmatrix=wmatrix_fn,
        kin=np.arange(0.001, 0.35, 0.001),
    )

#Create observables for each tracer
observables = [create_observable(params['data_fn'], params['wmatrix_fn'], params['covariance_fn'],theories, i) 
                for i, params in tracer_params.items()]


#Create an emulated theory for each tracer
for i in range(len(theories)):
    theories[i] = observables[i].wmatrix.theory
    emulator = Emulator(theories[i].pt, engine=TaylorEmulatorEngine(method = 'finite', order = 3))
    #emulator.save('Emulator/FOLPSAX_mf_Taylor_o4_LRG1')
    emulator.set_samples()
    emulator.fit()
    
    theories[i].init.update(pt = emulator.to_calculator())

print('All theories have been emulated succesfully')

#Analytic marginalization over eft and nuisance parameters
for i in range(len(theories)): 
    if prior_basis == 'physical':
        params_list = ['alpha0p', 'alpha2p', 'alpha4p', 'sn0p', 'sn2p']
    else:
        params_list = ['alpha0', 'alpha2', 'alpha4', 'sn0', 'sn2']

    for param in params_list:    
        theories[i].params[param].update(derived = '.marg')

        
#Rename the eft and nuisance parameters to get a parameter for each tracer (i.e. QSO_alpha0, QSO_alpha2, BGS_alpha0,...)        
for i in range(len(theories)):    
    for param in theories[i].init.params:
        # Update latex just to have better labels
        param.update(namespace='{}'.format(tracers[i])) 
        

#Create a likelihood per theory object
setup_logging()
Likelihoods = []
for i in range(len(theories)):
    Likelihoods.append(ObservablesGaussianLikelihood(observables = [observables[i]]))


#Sum the likelihoods and initialize
likelihood = SumLikelihood(likelihoods = (Likelihoods))

likelihood()


#Run the sampler and save the chain
from desilike.samplers import EmceeSampler, MCMCSampler

if sampler == 'cobaya':
    sampler = MCMCSampler(likelihood ,save_fn = chain_name)
    sampler.run(check={'max_eigen_gr': 0.03})
    
else:
    sampler = EmceeSampler(likelihood ,save_fn = chain_name)
    sampler.run(check={'max_eigen_gr': 0.3})

