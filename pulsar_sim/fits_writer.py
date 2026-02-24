# pulsar_sim/fits_writer.py
import os
import numpy as np
import psrchive
from pulsar_sim.io import readparfile
from pulsar_sim.profiles import make_profile


def makefits(parfile, nsub, tsub, npol, nchan, nbin, freq_low, freq_high, snr, ar_mjd, dDM, tscat, alpha, spidx, fluxd, pbf, simpfile, noscat, nodm, onlydm, onlyscat, verbose=True):
	"""
	Produces the simulated archive in psrfits format based on the input arguments.
	"""

	# Reads name, position and DM of the pulsar
	ar_psr, ar_raj, ar_decj, ref_dm= readparfile(parfile)
	ar_dm = float(ref_dm) + dDM

	# Creates the frequency array for sub-bands
	freq_c = freq_low + (freq_high - freq_low)/2
	bw = freq_high - freq_low
	chanwidth = bw / nchan
	freqArr = np.linspace(freq_low + (chanwidth/2.0), freq_low + bw - (chanwidth/2.0), nchan)
	tscat_cfr = tscat * (1000./freq_c)**(-1.*alpha)
	# Writes out the DM tscat and alpha values of the epoch
	f = open(ar_psr+"_dmalptscat.dm",'a')
	f.write("%.6f %.10f %.10f %.6f\n" % (ar_mjd, ar_dm, tscat, alpha))
	f.close()

	# Creates an instance of psrfits file. Only possible with ASP backend.
	ar = psrchive.Archive_new_Archive('ASP')
	# reshaping the data matrix
	ar.resize(nsub,npol,nchan,nbin)
	# Setting the DM of the epoch
	ar.set_dispersion_measure(float(ar_dm))
	# Setting source name
	ar.set_source(ar_psr)
	# Setting the pulsar position
	ar.set_coordinates(psrchive.sky_coord(ar_raj + ar_decj))
	# Setting the centre frequency of observation
	ar.set_centre_frequency(freq_c)
	# Setting the bandwidth of the archive
	ar.set_bandwidth(bw)
	# Setting the telescope. Resets ASP to GMRT.
	ar.set_telescope('DE601')

	# Setting MJD of the archive
	start_MJD = psrchive.MJD(ar_mjd)
	epoch = start_MJD+tsub/2.0

	# setting time of subints
	for subint in ar:
		subint.set_epoch(epoch)
		subint.set_duration(tsub)
		epoch += tsub
		# Setting the frequencies of the sub-bands
		for ichan in range(nchan):
			subint.set_centre_frequency(ichan,freqArr[ichan])

	# Setting the par file. It also creates the predictor table
	ar.set_ephemeris(parfile)
	# setting the archive as dedispersed
	ar.set_dedispersed(True)
	ar.dedisperse()

	# Getting the folding period fo converting the tscat to nbins
	p0 = ar.get_Integration(0).get_folding_period() * 1000.0

	if verbose:
		print(f"Folding period (P0): {p0:.6f} ms")
	
    # Getting the bandshape. Currently disabled as a better way is sought.
	#chanamps = bandshape(nchan)

	# Main loop that writes the data for each channel into the profiles of archive
	for isub in range(nsub):
		sub = ar.get_Integration(isub)

		for ipol in range(npol):

			for ichan in range(nchan):
				prof = sub.get_Profile(ipol,ichan)
				prof.get_amps()[:] = make_profile(nbin, snr, freqArr[ichan], ar_dm, tscat, alpha, p0, spidx, fluxd, pbf, ar_mjd, ar_psr, simpfile, chanwidth, noscat, nodm, onlydm, onlyscat)
				#prof += chanamps[ichan] # Correcting for the bandshape

	# de-dedispersing the data and resetting the DM to reference DM.
	# This imitates the DM variation.
	ar.dededisperse()
	ar.set_dispersion_measure(float(ref_dm))

	# Writing the archive out and correcting for an important flag
	outfile = ar_psr+"."+str("%.3f"%ar_mjd)+"_"+str("%.0f" % freq_high)+"MHz."+ar.get_telescope()+".ar"
	ar.set_filename(outfile)
	ar.unload(outfile)
	os.popen("psredit -c aux:dmc=0 -m %s" % outfile).read()
	print ("Written out %s" % outfile)
	# Done.