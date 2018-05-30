import astropy.units as q
import pickle
import mcmc_fit.mcmc_fit as mc
import numpy as np
import time
start_time = time.time()


# Read in SED of 1256-0224
SED_path = '../Atmospheres_paper/Data/Gaia1256-0224 (L3.5sd) SED.txt'  # Move up to Modules and into folder
w, f, e = np.loadtxt(SED_path, delimiter=" ", unpack=True)

# Read in only NIR spectrum
wavelength, flux, err = np.loadtxt('SED1256_nir_Fire.txt', delimiter="\t", unpack=True, skiprows=1)

# ----------------------- If you have a pickle file of the model grid (it must have no units!) --------------------
# Load the pickle file of the model grid
# pickle_path = '/Users/eileengonzales/Dropbox/BDNYC/BDNYCdb_copy/SXD_r2000_Marley.pkl'
# file = open(pickle_path, 'rb')  # rb means read binary
# models = pickle.load(file)
# file.close()
#
# # Make the wavelength array the proper length
# models['wavelength'] = [models['wavelength']] * len(models['flux'])
#
# # This makes the model grid
# mg = mc.make_model_db('btsettl', 'model_atmosphere_db', model_grid=models, grid_data='spec',
#                       param_lims=[('teff', 2000, 2500, 50), ('logg', 3.5, 5.5, 0.5)], fill_holes=False, bands=[],
#                       rebin_models=w, use_pandas=False)

# --------------------------- If you are getting grid from the Model Database ---------------------------------------
# Testing to run from calling the database and not a pickle file
model_atmosphere_db = '/Users/eileengonzales/Dropbox/BDNYC/BDNYCdb_copy/model_atmospheres.db'
# Get the BtSettl grid
mg = mc.make_model_db('bt_settl_2013', model_atmosphere_db, grid_data='spec',
                      param_lims=[('teff', 2000, 2500, 50), ('logg', 3.5, 5.5, 0.5)], fill_holes=False, bands=[],
                      rebin_models=wavelength, use_pandas=False)

# Get Marley & Saumon. f_sed is only 2 for temps above 1600 K
mg2 = mc.make_model_db('marley_saumon', model_atmosphere_db, grid_data='spec',
                       param_lims=[('teff', 2000, 2400, 50), ('logg', 3.5, 5.5, 0.5), ('f_sed', 0, 2, 1)],
                       fill_holes=False, bands=[], rebin_models=w, use_pandas=False)

# Adding units to the arrays for either way grid is generated
w = w * q.um
f = f * q.erg/q.AA/q.cm**2/q.s
e = e * q.erg/q.AA/q.cm**2/q.s

wavelength = wavelength * q.um
flux = flux * q.erg/q.AA/q.cm**2/q.s
err = err * q.erg/q.AA/q.cm**2/q.s

# This does the actual fitting to the spectra
# def fit_spectrum(raw_spectrum, model_grid, model_grid_name, shortname, walkers, steps, mask=[], db='', extents=None,
#                 object_name='Test', log=False, plot=True, prnt=True, generate=True, outfile=None):
# See line 212 fo mcmc_fit for more details

# test Run to see if the code is working properly
# bdsamp = mc.fit_spectrum([w, f, e], mg, 'bt_settl_2013', '1256-0224', 1, 5, mask=[], db='', extents=None,
#                          object_name='J1256', log=False, plot=True, prnt=True, generate=True, outfile=None)


# Full Run on J1256 trying 3 walkers and 100 steps to get an idea of how long the run takes.
# bdsamp = mc.fit_spectrum([w, f, e], mg, 'bt_settl_2013', '1256-0224', 2, 5, mask=[], db='', extents=None,
#                          run_name='J1256_test', log=False, plot=True, prnt=True, generate=True, outfile=None)

bdsamp_nir = mc.fit_spectrum([wavelength, flux, err], mg, 'bt_settl_2013', '1256-0224', 2, 5, mask=[], db='', range=None,
                         run_name='J1256_nironly', log=False, plot=True, prnt=True, generate=True, outfile=None)

# bdsamp2 = mc.fit_spectrum([w, f, e], mg2, 'marley_saumon', '1256-0224', 5, 100, mask=[], db='', range=None,
#                           run_name='J1256_marley', log=False, plot=True, prnt=True, generate=True, outfile=None)

# determine time it takes code to run 2, walkers 5 steps 11.07 minutes. 3 walkers, 100 steps -7.8 hours
# nir only 2,5 = 7.8 minutes
print "time elapsed: {:.2f}s".format(time.time() - start_time)
