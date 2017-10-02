# multiscaleGrangerCausality
In the study of complex physical and biological systems represented by multivariate stochastic processes, an issue of great relevance is the description of the system dynamics spanning multiple temporal scales. While methods to assess the dynamic complexity of individual processes at different time scales are well-established, multiscale analysis of directed interactions has never been formalized theoretically, and empirical evaluations are complicated by practical issues such as filtering and downsampling.
We extend the very popular measure of Granger causality (GC), a prominent tool for assessing directed lagged interactions between joint processes, to quantify information transfer across multiple time scales. We show that the multiscale processing of a vector autoregressive (AR) process introduces a moving average (MA) component, and describe how to represent the resulting ARMA process using state space (SS) models and to combine the SS model parameters for computing exact GC values at arbitrarily large time scales. We exploit the theoretical formulation to identify peculiar features of multiscale GC in basic AR processes, and demonstrate with numerical simulations the much larger estimation accuracy of the SS approach compared with pure AR modeling of ltered and downsampled data. The improved computational reliability is exploited to disclose meaningful multiscale patterns of information transfer between global temperature and carbon dioxide concentration time series, both in paleoclimate and in recent years. 
The msGC Matlab toolbox reroduces algorithms, simulations and real data analysis reported in [1].

Main functions
--------------------
- mscg.m: computes Granger causality at a given scale from assigned VAR parameters
- iss_PV: computes partial variances for a state space model from innovations form SS parameters
- iss_ss2iss.m: computes innovations form SS parameters from SS parameters
- iss_varma2iss.m: computes innovations form SS parameters from VARMA parameters
- iss_ds.m: computes innovations form SS parameters for a downsampled SS model
External functions
--------------------
- egc_gcMVAR.m (LinReg, linReg_Ftest, SetLag, buildvectors): functions to compute standard time domain Granger Causality
- eMVAR_idMVAR.m (choldiag, InstModelFilter, mos_idMVAR, MVARfilter): functions to identify and simulate a MVAR process
- surriaaft.m (surrshuf): functions to generate IAAFT surrogates
Scripts
--------------------
- test_ssGC.m: computation of multiscale GC with state-space models for simple simulated VAR processes
- PRE2017_simulationA.m: realizes simulation A of [1]
- PRE2017_simulationB.m: realizes simulation b of [1]
- PRE2017_application_modern.m: realizes application to modern climate data performed in [1]
- PRE2017_application_paleo.m: realizes application to paleolithic climate data performed in [1]
Data
--------------------
- data_modern.mat (data_modern_L1-filtered): climatological time series analyzed in [1] (modern climate data)
- data_800.mat: climatological time series analyzed in [1] (paleolithic climate data)
NOTE: the iss_ds and ss2iss functions are taken from the State-Space Granger Causality Matlab Toolbox  - http://users.sussex.ac.uk/~lionelb/downloads/ssgc.zip

[1] L Faes, L Faes, S Stramaglia, G Nollo, D Marinazzo. Multiscale Granger causality. Physical Review E. 2017. https://arxiv.org/abs/1703.08487
