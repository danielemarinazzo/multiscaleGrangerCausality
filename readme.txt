msGC - MATLAB toolbox for the computation of multiscale Granger Causality
based on Faes et al. 2017 [1]

Main functions:
- mscg.m: computes Granger causality at a given scale from assigned VAR parameters
- iss_PV: computes partial variances for a state space model from innovations form SS parameters
- iss_ss2iss.m: computes innovations form SS parameters from SS parameters
- iss_varma2iss.m: computes innovations form SS parameters from VARMA parameters
- iss_ds.m: computes innovations form SS parameters for a downsampled SS model

External functions
- egc_gcMVAR.m (LinReg, linReg_Ftest, SetLag, buildvectors): functions to compute standard time domain Granger Causality
- eMVAR_idMVAR.m (choldiag, InstModelFilter, mos_idMVAR, MVARfilter): functions to identify and simulate a MVAR process
- surriaaft.m (surrshuf): functions to generate IAAFT surrogates

Scripts
- test_ssGC.m: computation of multiscale GC with state-space models for simple simulated VAR processes
- PRE2017_simulationA.m: realizes simulation A of [1]
- PRE2017_simulationB.m: realizes simulation b of [1]
- PRE2017_application_modern.m: realizes application to modern climate data performed in [1]
- PRE2017_application_paleo.m: realizes application to paleolithic climate data performed in [1]

Data:
- data_modern.mat (data_modern_L1-filtered): climatological time series analyzed in [1] (modern climate data)
- data_800.mat: climatological time series analyzed in [1] (paleolithic climate data)

[1] L Faes, S Stramaglia, G Nollo, D Marinazzo. Multiscale Granger causality. Physical Review E. 2017