# MOSCATEL
Rough description:
Command-line implementation of basic MuSCAT photometry pipeline and transit lightcurve analysis

## Progress
2017/03/21: basic analysis during Okayama observation
2017/04/09: basic high-level scripting

## Sample run
### Part 1 Photometry
e.g. python moscatel --band_idx=2 --skip_every=10
### Part 2 Lightcurve analysis
e.g. python analysis --target=b --ref=a star=1

See also other plotting helper functions in /scripts/utils.py.

TO DO: 
1. define centroids (tuples) either in config.dat or as user input
2. low-level control on parameters: e.g. radius of aperture, annuli etc.
3. incorporate advanced analysis e.g. mcmc and GP modeling
4. upgrade to classes
5. fix setup.py
