# HMDS Modeling

Implements and analyzes the HMDS model.

``modeling.R`` - fits the HMDS model in Rstan (``HMDS.stan``) and saves the results for the Beethoven symphony data (``Beethoven-hellinger.rda``)
for the tempo, dynamics and spectral flatness/timbre metrics

``model-checks.R`` - performs MCMC approximation and goodness-of-fit checks for the HMDS results

``results-analysis.R`` - analyzes model fit and reproduces results presented in the paper


## Data

``Beethoven-hellinger.rda`` contains pairwise distance matrices (calculated using the Hellinger distance) for various musical features, between 10 orchestras across all 9 Beethoven symphonies (with movements separated as different pieces)

The array is 10 x 10 x 37 x 23, where 10 is the number of orchestras, 37 is the number of pieces and 23 is the total number of musical metrics considered
All dimensions of the array are labeled, but the 23 musical metrics are (in order):

1. ``tempo`` - tempo metric
2. ``volume`` - dynamics metric
3. ``SF`` - spectral flatness or timbre metric, over the entire spectrum
4. ``spec1`` - ``spec10`` - dynamics metric for a subset/group of frequency ranges, roughly dynamics balance across the spectrum
5. ``SF1`` - ``SF10`` - spectral flatness metric for a subset/group of frequency ranges
