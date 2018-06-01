import astropy.units as q
import pickle
import mcmc_fit.mcmc_fit as mc
import numpy as np
import time
start_time = time.time()
#from synth_fit.Calc_chisquare_modifiedEG import *
from synth_fit.calc_chisq import *


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
                      param_lims=[('teff', 1200, 3000, 50), ('logg', 3.5, 5.5, 0.5)], fill_holes=False, bands=[],
                      rebin_models=w, use_pandas=False)

# Get Marley & Saumon. f_sed is only 2 for temps above 1600 K
mg2 = mc.make_model_db('marley_saumon', model_atmosphere_db, grid_data='spec',
                       param_lims=[('teff', 1600, 2400, 50), ('logg', 3.5, 5.5, 0.5), ('f_sed', 0, 2, 1)],
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

# --------------- Probably won't run the MCMC and will just do a Chi-sqaured analysis for the report. ----------------
# Full Run on J1256 trying 3 walkers and 100 steps to get an idea of how long the run takes.
bdsamp = mc.fit_spectrum([w, f, e], mg, 'bt_settl_2013', '1256-0224chi', 2, 5, mask=[], db='', range=None,
                         run_name='J1256_testchi', log=False, plot=True, prnt=True, generate=True, outfile=None)

bdsamp_nir = mc.fit_spectrum([wavelength, flux, err], mg, 'bt_settl_2013', '1256-0224', 2, 5, mask=[], db='', range=None,
                         run_name='J1256_nironly', log=False, plot=True, prnt=True, generate=True, outfile=None)

bdsamp2 = mc.fit_spectrum([w, f, e], mg2, 'marley_saumon', '1256-0224test', 5, 100, mask=[], db='', range=None,
                          run_name='J1256_test', log=False, plot=True, prnt=True, generate=True, outfile=None)

# determine time it takes code to run 2, walkers 5 steps 11.07 minutes. 3 walkers, 100 steps -7.8 hours
# nir only 2,5 = 7.8 minutes
print "time elapsed: {:.2f}s".format(time.time() - start_time)

# Run the Chi-squared individually
test_all(w,f,e,mg2,['teff','logg','f_sed','k_zz'],smooth=False, resolution=None, shortname='J1256-0224_tf')
#Output:  (array([  2.40000000e+03,   4.00000000e+00,   2.00000000e+00, 0.00000000e+00]), 0.17684231505573014)

test_all(w,f,e,mg,['teff','logg'],smooth=False, resolution=None, shortname='J1256-0224_tfbt')
# Output: (array([ 1650.,     4.]), 0.014132433189028338)


# -------- Function needed to plot the model ----------------

def rebin_spec(spec, wavnew, waveunits='um'):
    from pysynphot import spectrum, observation
    # Gives same error answer: Err = np.array([np.sqrt(sum(spec[2].value[idx_include(wavnew,[((wavnew[0] if n==0 else
    #  wavnew[n-1]+wavnew[n])/2,wavnew[-1] if n==len(wavnew) else (wavnew[n]+wavnew[n+1])/2)])]**2)) for n in range(
    # len(wavnew)-1)])*spec[2].unit if spec[2] is not '' else ''
    if len(spec) == 2:
        spec += ['']
    try:
        Flx, Err, filt = spectrum.ArraySourceSpectrum(wave=spec[0].value,
                                                      flux=spec[1].value), spectrum.ArraySourceSpectrum(
            wave=spec[0].value, flux=spec[2].value) if spec[2] else '', spectrum.ArraySpectralElement(spec[0].value,
                                                                                                      np.ones(
                                                                                                          len(spec[0])),
                                                                                                      waveunits=waveunits)
    except:
        spec, wavnew = [i * q.Unit('') for i in spec], wavnew * q.Unit('')
        Flx, Err, filt = spectrum.ArraySourceSpectrum(wave=spec[0].value,
                                                      flux=spec[1].value), spectrum.ArraySourceSpectrum(
            wave=spec[0].value, flux=spec[2].value) if spec[2] else '', spectrum.ArraySpectralElement(spec[0].value,
                                                                                                      np.ones(
                                                                                                          len(spec[0])),
                                                                                                      waveunits=waveunits)
    return [wavnew, observation.Observation(Flx, filt, binset=wavnew.value, force='taper').binflux * spec[1].unit,
            observation.Observation(Err, filt, binset=wavnew.value, force='taper').binflux * spec[2].unit if spec[
                2] else np.ones(len(wavnew)) * spec[1].unit]

# -------------------------------------


# Pull out the models in order to plot them up.
ma_db = astrodb.Database(model_atmosphere_db)
ms_nocloud = ma_db.query("select wavelength, flux from marley_saumon where logg=4.0 and teff=2400 and f_sed=2 and k_zz=0", fmt='dict')
bt_cloud = ma_db.query("select wavelength, flux from bt_settl_2013 where logg=4.0 and teff=1650", fmt='dict')


# pull out the wavelenght and flux
flux_nc=ms_nocloud[0]['flux']
wave_nc=ms_nocloud[0]['wavelength']
flux_cl=bt_cloud[0]['flux']
wave_cl=bt_cloud[0]['wavelength']

# Normalize to max flux and smooth to same resolution
norm_flnc=flux_nc/max(flux_nc)
norm_flcl=flux_cl/max(flux_cl)

# create an uncertainty array for the model
unc = np.ones(len(norm_flnc))
unc2 = np.ones(len(norm_flcl))

spec = [wave_nc,norm_flnc,unc]   # group it together
speck = rebin_spec(spec,w)  # rebin to the same wavelength as 1256.  Is this effectively smoothing????
spec2 = [wave_cl,norm_flcl,unc2]   # group it together
speck2 = rebin_spec(spec2,w)


wl_bin_nc = speck[0].value
fl_bin_nc = speck[1].value
wl_bin_cl = speck2[0].value
fl_bin_cl = speck2[1].value

# Normalize SED of 1256 to peak flux
norm_1256 = f/max(f)

plt.plot(wl_bin_nc,fl_bin_nc)
plt.plot(wl_bin_cl,fl_bin_cl)
plt.plot(w,norm_1256)

# TODO: Add legend, axis labels and all the good stuff!



