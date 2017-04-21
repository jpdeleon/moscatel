# MOSCATEL
Rough description:
Command-line implementation of basic MuSCAT photometry pipeline and transit lightcurve analysis

## Progress
2017/03/21: basic analysis during Okayama observation
2017/04/09: basic high-level scripting
2017/04/20: created setup.py

## Sample run
### Part 1 Photometry

* e.g. $ python moscatel-phot --band_idx=2 --skip_every=5

### Part 2 Lightcurve analysis
* e.g. $ python moscatel-analysis --target=b --ref=a star=0

See also other plotting helper functions in /moscatel/utils.py.

TO DO: 
1. define input/output directories, centroids (tuples), etc in config.yaml
2. implement low-level control: e.g. radius of aperture, annuli etc.
3. incorporate advanced analysis e.g. mcmc and GP modeling
4. upgrade to classes
