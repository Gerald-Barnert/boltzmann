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

zz = np.arange(zmin=0, zmax=3, step=0.05)
camb_code = boltzmann_codes.run(redshift=zz, boltzamnn_code="camb", hubble_constant=70, 
                                omega_baryion=0.04, omega_dark_matter=0.24)
hubbleparameter = boltzmann_codes.HubbleParameter(zz, camb_code, "camb")
```
cosmology.py it is an example of the usage of boltzmann_codes.py.
## License
[MIT](https://choosealicense.com/licenses/mit/)
