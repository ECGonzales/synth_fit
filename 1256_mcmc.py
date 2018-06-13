from astrodbkit import astrodb
import astropy.units as q
import pickle
import mcmc_fit.mcmc_fit as mc
import numpy as np
import time
start_time = time.time()
from synth_fit.Calc_chisquare_modifiedEG import *
#from synth_fit.calc_chisq import *


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
                      param_lims=[('teff', 1600, 3000, 50), ('logg', 3.5, 5.5, 0.5)], fill_holes=False, bands=[],
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


# ---------------------------------------------------------------------------------------------------------------------
#  ---------------------- Run the Chi-squared individually, b/c of code issues with normalizing -----------------------
# ---------------------------------------------------------------------------------------------------------------------
# Entire SED
test_all(w, f, e, mg2, ['teff', 'logg', 'f_sed', 'k_zz'], smooth=False, resolution=None, shortname='J1256-0224')
# Output: (array([  2.40000000e+03,   4.00000000e+00,   2.00000000e+00, 0.00000000e+00]), 0.17684231505573014) Constant
# Output for norm region: (array([  2.40000000e+03,   5.50000000e+00,   2.00000000e+00,
#          4.00000000e+00]), 286.8504874757507, 0.016133323255104089)
# params, chisqr, reduced chisqr

test_all(w,f,e,mg,['teff', 'logg'], smooth=False, resolution=None, shortname='J1256-0224')
# Output: (array([ 1650.,     4.]), 0.014132433189028338)
# Output for norm region: (array([ 2500.,     5.]), 85.216968270177901, 0.0047928553582777221)

# FIRE spectrum only
test_all(wavelength, flux, err, mg2, ['teff', 'logg', 'f_sed', 'k_zz'], smooth=False, resolution=None,
         shortname='J1256-0224_Fireonly')  # This won't run for some reason.

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

# ---------------------------------------------------------------------------------------------
# Pull out the models in order to plot them up. The is for the normalized by a constant method
# ---------------------------------------------------------------------------------------------
# ma_db = astrodb.Database(model_atmosphere_db)
# ms_nocloud = ma_db.query("select wavelength, flux from marley_saumon where logg=4.0 and teff=2400 and f_sed=2 and k_zz=0", fmt='dict')
# bt_cloud = ma_db.query("select wavelength, flux from bt_settl_2013 where logg=4.0 and teff=1650", fmt='dict')
#
# # pull out the wavelength and flux
# flux_nc = ms_nocloud[0]['flux']
# wave_nc = ms_nocloud[0]['wavelength']
# flux_cl = bt_cloud[0]['flux']
# wave_cl = bt_cloud[0]['wavelength']
#
# # smooth to same resolution
# unc = np.ones(len(flux_nc))
# spec = [wave_nc, flux_nc, unc]   # group it together
# speck = rebin_spec(spec, w)
#
# wl_bin_nc = speck[0].value
# fl_bin_nc = speck[1].value
# #
# # Normalize via that constant for the models (you don't apply this to the object)
# mult1 = f * fl_bin_nc / (e ** 2)
# bad = np.isnan(mult1)
# mult = np.sum(mult1[~bad])
# sq1 = fl_bin_nc * fl_bin_nc / (e ** 2)
# mult2 = float(sum(sq1[~bad]))
# ck = mult / mult2
# norm_mod_flux_nc = fl_bin_nc * ck
#
# # Normalizing to the max flux. (This isn't what happened in the code (did the constant), so this isn't vaild per say)
# norm_flnc=flux_nc/max(flux_nc)
# norm_flcl=flux_cl/max(flux_cl)
#
# # create an uncertainty array for the model
# unc = np.ones(len(norm_flnc))
# unc2 = np.ones(len(norm_flcl))
#
# spec = [wave_nc, norm_flnc, unc]   # group it together
# speck = rebin_spec(spec, w)  # rebin to the same wavelength as 1256.  Is this effectively smoothing????
# spec2 = [wave_cl, norm_flcl, unc2]   # group it together
# speck2 = rebin_spec(spec2, w)
#
#
# wl_bin_nc = speck[0].value
# fl_bin_nc = speck[1].value
# wl_bin_cl = speck2[0].value
# fl_bin_cl = speck2[1].value
#
# # Normalize SED of 1256 to peak flux
# norm_1256 = f/max(f)
#
# plt.plot(wl_bin_nc, fl_bin_nc)
# plt.plot(wl_bin_cl, fl_bin_cl)
# plt.plot(w, norm_1256)

# ---------------------------------------------------------------------------------------------------------------------
# This is for the normalized by across a region method. I chose to normalized everything near the Y band peak for J1256,
# which was over the 0.98-0.988 micron region. This needs to be updated in the Chi-Squared when running on other objects
# ----------------------------------------------------------------------------------------------------------------------
# Pull out models where normalized to max flux.
ma_db = astrodb.Database(model_atmosphere_db)
ms_norm = ma_db.query("select wavelength, flux from marley_saumon where logg=5.5 and teff=2400 and f_sed=2 and k_zz=4",
                      fmt='dict')
bt_norm = ma_db.query("select wavelength, flux from bt_settl_2013 where logg=5.0 and teff=2500", fmt='dict')

# Re-Read in SED of 1256-0224 so that there is not units and can manipulate to plot.
SED_path = '../Atmospheres_paper/Data/Gaia1256-0224 (L3.5sd) SED.txt'  # Move up to Modules and into folder
w, f, e = np.loadtxt(SED_path, delimiter=" ", unpack=True)

# pull out the wavelength and flux
flux_ms = ms_norm[0]['flux']
wave_ms = ms_norm[0]['wavelength']
flux_bt = bt_norm[0]['flux']
wave_bt = bt_norm[0]['wavelength']

# normalize all of the models and 1256 over the 0.98-0.988 region
norm_region_ms = np.where(np.logical_and(wave_ms >= 0.98, wave_ms <= 0.988))
flux_ms_norm = flux_ms / np.average(flux_ms[norm_region_ms])
norm_region_bt = np.where(np.logical_and(wave_bt >= 0.98, wave_bt <= 0.988))
flux_bt_norm = flux_bt / np.average(flux_bt[norm_region_bt])
norm_region_1256 = np.where(np.logical_and(w >= 0.98, w <= 0.988))
flux_1256_norm = f / np.average(f[norm_region_1256])

# create an uncertainty array for the model since they don't have uncertainties. This is needed for the rebinning.
unc_ms = np.ones(len(flux_ms_norm))
unc_bt = np.ones(len(flux_bt_norm))

spec_ms = [wave_ms, flux_ms_norm, unc_ms]   # group it together
speck_ms = rebin_spec(spec_ms, w)  # rebin to the same wavelength as 1256.  This effectively smoothing.
spec_bt = [wave_bt, flux_bt_norm, unc_bt]   # group it together
speck_bt = rebin_spec(spec_bt, w)

# define new names for wavelength and flux from the rebinned speck objects
wl_bin_ms = speck_ms[0].value
fl_bin_ms = speck_ms[1].value
wl_bin_bt = speck_bt[0].value
fl_bin_bt = speck_bt[1].value

# Create the final plots
fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.set_size_inches(10, 6.45)
for axis in ['top', 'bottom', 'left', 'right']:  # Thicken the frame
    ax1.spines[axis].set_linewidth(1.1)
ax1.tick_params(axis='both', labelsize=20, length=8, width=1.1)
plt.xlabel('Wavelength ($\mu$m)', fontsize=25)
plt.ylabel('Normalized Flux ($F_\lambda$)', fontsize=25)
plt.xlim([0.6, 2.42])
plt.ylim([-0.01, 1.5])

# Add the data and models
plt.plot(wl_bin_bt, fl_bin_bt, color="#3FB22B")
plt.plot(wl_bin_ms, fl_bin_ms, color="#B21987")
plt.plot(w, flux_1256_norm, color='k')
plt.tight_layout()

# Add Labels
ax1.annotate('BT-Settl 2013 \n$T_\mathrm{eff}$: 2500 K, log $g$: 5.0', xy=(1.75, 1.35), color='#3FB22B', fontsize=12)
ax1.annotate('Marley & Saumon 2008 \n$T_\mathrm{eff}$: 2400 K, log $g$: 5.5, f$_\mathrm{sed}$: 2, $K_\mathrm{zz}$: 4',
             xy=(1.75, 1.2), color='#B21987', fontsize=12)
ax1.annotate('J1256-0224 \n$T_\mathrm{eff}$: 2307$\pm$ 71 K, log $g$: 5.37$\pm$0.01', xy=(1.75, 1.05), color='k', fontsize=12)


plt.savefig('../Atmospheres_paper/Plots/Model_comparison.pdf')

