# Boltzmann Codes
Boltzmann codes are used for calculate cosmological quantities given a set of parameters. [CAMB](https://camb.info/) and 
[CLASS](https://lesgourg.github.io/class_public/class.html) are boltzmann codes that compute the evolution of linear perturbations,
Cosmic Microwave Background (CMB), lensing, transfer and correlation functions, power spectra, etc.
## Usage
You need the python wrapper to be able to import camb and class. boltzmann_code.py contains functions that calculates different
cosmological quantities. Every function works depending of the boltzmann code. The main function its 'run', this function initialize
CAMB or CLASS setting the cosmology and power spectra. Every other function get 'run' as a parameter to compute the other quantities.
```python
#EXAMPLE
import boltzmann_codes

zmin, zmax, step = 0, 3, 0.05
redshifts = np.arange(zmin, zmax, step)
code = 'camb'
hubble_constant, omega_baryon, omega_dark_matter = 70, 0.04, 0.24
camb_code = boltzmann_codes.run(redshifts, code, hubble_constant, omega_baryon, omega_dark_matter)
hubble_parameter = boltzmann_codes.HubbleParameter(redshifts, camb_code, code)
```
cosmology.py it is an example of the usage of boltzmann_codes.py.
## License
[MIT](https://choosealicense.com/licenses/mit/)
