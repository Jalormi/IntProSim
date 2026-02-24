import numpy as np
from numpy import exp
from pulsar_sim.utils import convolve_funcs


# Creates a Gaussian shaped model
def make_gaussian(x, amp, phase, sigma):
	gauss = amp * exp( -(x-phase)**2./(2.*sigma**2.) )
	return ( gauss )

# Creates a Lorentzian shaped profile
def make_lorentzian( x, x0, a, gam ):
	lorentzian = a * gam**2 / ( gam**2 + ( x - x0 )**2)
	return(lorentzian)

# Thin scattering screen approximation
def pbf1(x, tscat):
	return ( exp(-x/tscat) )

# Thick scattering screen approximation (Eq. 13 in Williamson 1972)
def pbf2(x, tscat):
	x += 1e-20
	amp = np.sqrt( (np.pi * tscat) /(4 * x**3) )
	expo = (np.pi**2 * tscat) / (16 * x)
	return ( amp * exp(-expo) )

# Truncated screen approximation (Eq. 3 in Geyer et al. 2017)
def pbf3(x, tscat):
	x += 1e-20
	amp = 1./np.sqrt(np.pi * tscat * x)
	# suboptimal. But a workaround to avoid the zero bin peak
	amp[0] = amp[1]+amp[2]
	return ( amp * ( exp(-x/tscat) ) )

# Continuous media approximation (Eq 21 in Williamson 1972)
def pbf4(x, tscat):
	x += 1e-20
	amp = np.sqrt( (np.pi**5 * tscat**3) / (8*x**5) )
	expo = (np.pi**2 * tscat) / (4.*x)
	return( amp * exp(-expo) )

def apply_dm_smearing(profile, dm, freq_MHz, chanwidth_MHz, period_s):
    """
    Apply symmetric boxcar DM smearing without shifting peak.
    Works for frequency-scrunched single profile.

    profile        : 1D numpy array
    dm             : dispersion measure (pc cm^-3)
    freq_MHz       : observing frequency (MHz)
    chanwidth_MHz  : channel width (MHz)
    period_s       : pulsar period (seconds)
    """

    nbin = len(profile)

    # ---- DM smearing time in seconds ----
    # 8.3e6 gives microseconds
    delta_t_us = 8.3e6 * dm * chanwidth_MHz / (freq_MHz**3)

    delta_t_s = delta_t_us * 1e-6

    # ---- Convert to bins ----
    smear_bins = delta_t_s / period_s * nbin

    smear_bins = int(np.round(smear_bins))

    if smear_bins < 0.5:
        return profile

    # ---- Make symmetric boxcar ----
    kernel = np.ones(smear_bins)
    kernel /= kernel.sum()

    # ---- Convolve WITHOUT changing length ----
    smeared = np.convolve(profile, kernel, mode='same')

    return smeared

# Generates the pulse profile
def make_profile(nbin, snr, freq, dm, tscat, alpha, p0, spidx, fluxd, pbf, mjd, psr, simpfile, chanwidth, noscat, nodm, onlydm, onlyscat):

	"""
	Generates the profile based on the input params. One can define
	as many gaussian components as required in the file. The profile
	also gets scatter-broadened, if asked.
	"""

	x = np.arange(nbin, dtype=float)
	g_amp = []; g_phase = []; g_sigma = []
	g_wsc = []; g_hsc = []; g_psc = []
	l_amp = []; l_phase = []; l_gamma = []
	l_gsc = []; l_hsc = []; l_psc = []

	prof = np.zeros(nbin)
	dip_condition=False; param_dip = []
	# Reads the data of profile gaussian components from the file
	for line in open(simpfile):

		if line.startswith("gauss"):
			data = line.split()
			g_phase.append(float(data[1])*nbin)
			g_amp.append(float(data[2]))
			g_sigma.append((float(data[3]) / 2.35482)*nbin)
			g_psc.append(float(data[4]))
			g_hsc.append(float(data[5]))
			g_wsc.append(float(data[6]))

		if line.startswith("lorentz"):
			data = line.split()
			l_phase.append(float(data[1])*nbin)
			l_amp.append(float(data[2]))
			l_gamma.append((float(data[3]))*nbin)
			l_psc.append(float(data[4]))
			l_hsc.append(float(data[5]))
			l_gsc.append(float(data[6]))

		if line.startswith("edip"):
			dip_condition=True
			data = line.split()
			param_dip.append(str(data[4]))
			param_dip = param_dip[0]

	# Makes the combined profile of all the gaussian/lorentzian components
	for i in range (np.size(g_amp)):
		gstr = 'gauss'+str(i+1)
		if (dip_condition == True):
			if (param_dip.startswith('gauss')):
				a_f = 0;
				dd = np.loadtxt(psr+"_amplitudes.txt", dtype=float,usecols=0)
				ap = np.loadtxt(psr+"_amplitudes.txt", dtype=float,usecols=1)
				xp = np.where(np.isclose(dd,mjd))[0][0]
				p_f = g_phase[i]*(freq/1000.)**g_psc[i]
				w_f = g_sigma[i]*(freq/1000.)**g_wsc[i]
				if (gstr == param_dip):
					if (np.size(ap) >1):
						a_f = (ap[xp] + g_amp[i])*(1.0/freq)**g_hsc[i]
					if (np.size(ap) ==1):
						a_f = (g_amp[i]+ap)*(1.0/freq)**g_hsc[i]
				if (gstr != param_dip):
					a_f = g_amp[i]*(1.0/freq)**g_hsc[i]
				prof += make_gaussian(x, a_f, p_f, w_f)

		if (dip_condition == False) or (param_dip == 'dm'):
			p_f = g_phase[i]*(freq/1000.)**g_psc[i]
			w_f = g_sigma[i]*(freq/1000.)**g_wsc[i]
			a_f = g_amp[i]*(1.0/freq)**g_hsc[i]
			prof += make_gaussian(x, a_f, p_f, w_f)

	for i in range (np.size(l_amp)):
		p_f = l_phase[i]*(freq/1000.)**l_psc[i]
		w_f = l_gamma[i]*(freq/1000.)**l_gsc[i]
		a_f = l_amp[i]*(1.0/freq)**l_hsc[i]
		prof += make_gaussian(x, a_f, p_f, w_f)

	# Scales the flux density of the pulse profile
	prof /= np.sum(prof)
	fscale = fluxd*((freq/1400)**spidx)
	prof *= fscale

	# ----------------------------------
	# APPLY BROADENING EFFECTS
	# ----------------------------------

	# ---- CASE 1: Only DM ----
	if onlydm:

		prof = apply_dm_smearing(
			profile=prof,
			dm=dm,
			freq_MHz=freq,
			chanwidth_MHz=chanwidth,
			period_s=p0
		)

	# ---- CASE 2: Only Scattering ----
	elif onlyscat:

		if not noscat:
			tscat_f = tscat * ((freq / 1000.0) ** alpha)
			tscat_bins = nbin * (tscat_f / p0)

			scmodel = "pbf%d" % pbf
			scatfunc = globals()[scmodel]
			scatprof = scatfunc(x=x, tscat=tscat_bins)

			prof = convolve_funcs(prof, scatprof)

	# ---- CASE 3: Neither ----
	elif noscat and nodm:

		pass   # leave intrinsic profile unchanged

	# ---- CASE 4: Both (default behaviour) ----
	else:

		# Scattering
		if not noscat:
			tscat_f = tscat * ((freq / 1000.0) ** alpha)
			tscat_bins = nbin * (tscat_f / p0)

			scmodel = "pbf%d" % pbf
			scatfunc = globals()[scmodel]
			scatprof = scatfunc(x=x, tscat=tscat_bins)

			prof = convolve_funcs(prof, scatprof)

		# DM smearing
		if not nodm:
			prof = apply_dm_smearing(
				profile=prof,
				dm=dm,
				freq_MHz=freq,
				chanwidth_MHz=chanwidth,
				period_s=p0
			)

	# Adds noise based on the peak S/N from input
	noise = np.random.normal(0,np.max(prof)/snr,nbin)

	return (prof+noise)
