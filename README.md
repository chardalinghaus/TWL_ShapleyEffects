# Assessing coastal flooding uncertainties using global sensitivity analysis
Use of a variance-based global sensitivity analysis called ”Shapley effects” to quantify the contributing uncertainty in total water levels by mean sea level, astronomical tide, storm surge and runup/wave dynamics, both historically and for future projections up to 2099. The study areas are Ōrewa and Muriwai beaches, located on the East and West coasts of New Zealand, characterized by entirely different processes and conditions.

Sea-level rise scenarios, astronomical tides, storm surge and wave data used in this study are publicly available. Sea-level rise scenarios, including vertical land motion, are provided by the [New Zealand SeaRise Project](https://www.searise.nz/). Astronomical tides are provided by [NIWA](https://tides.niwa.co.nz/). Storm surge and wave data can be visualized and downloaded on an interactive web dashboard from [Coastal and Ocean Collective](https://coastalhub.science/data). 

**SLR_WeibullDistribution.py**: After downloaded, we fitted a Weibull distribution for each sea level rise scenario, generating probability distributions. Using linear interpolation, we generated a comprehensive projected uniform sea level rise distribution with 5000 values per year. 
**AT_TidePredictions.py**: Extraction of tidal constituents, reconstruction and prediction of future tides.
**MonteCarloSimulations.py**: After preparing the data, we sampled and estimated 100,000 values of total water level per year.

The implementation of Shapley Effects (**Shapley_Effects.R**) was realized by applying the sobolshap_knn() ([Broto et al., 2020](https://epubs.siam.org/doi/10.1137/18M1234631); [Song et al., 2016](https://epubs.siam.org/doi/10.1137/15M1048070)) algorithm implemented in the [sensitivity](https://cran.r-project.org/web/packages/sensitivity/index.html) R package. Finally, the results can be visualized with **Shapley_Effects_Plot.py**.

Note that the scripts here are for near-future scenarios but can easily be modified for other periods.
