# Boltzmann Codes
Boltzmann codes are used for calculate cosmological quantities given a set of parameters. [CAMB](https://camb.info/) and 
[CLASS](https://lesgourg.github.io/class_public/class.html) are boltzmann codes that compute the evolution of linear perturbations,
Cosmic Microwave Background (CMB), lensing, transfer and correlation functions, power spectra, etc.
## Usage
You need the python wrapper to be able to import camb and class. boltzmann_code.py contains functions that calculates different
cosmological quantities. Every function works depending of the boltzmann code. The main function its 'run', that initialize
CAMB or CLASS setting the cosmology and power spectra. Every other function get 'run' as a parameter to compute the other quantities.
