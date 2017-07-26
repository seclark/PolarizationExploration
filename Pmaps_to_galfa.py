from __future__ import division
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import healpy as hp

#path = "/Volumes/DataDavy/Planck/"
#path = "/disks/jansky/a/users/goldston/susan/Planck/"
#infile1 = "HFI_SkyMap_353_2048_R2.02_full.fits"
path="/disks/jansky/a/users/goldston/susan/Planck/"
#infile1="mapPsI-ForSusanClark_fwhm80_ns128_AngSt1.fits"
#infile1="mapS-ForSusanClark_fwhm80_ns128_AngSt1.fits"
infile1="faraday.fits"
PfileQU = path + infile1 

print("loading", PfileQU)

Nside = 128
Npix = 12*Nside**2

#Pdata = fits.getdata(PfileQU) # this works for P
Pdata = hp.fitsfunc.read_map(PfileQU) # for faraday?


#Gfile = '/Volumes/DataDavy/GALFA/DR2/FullSkyWide/GALFA_HI_W_S0900_V-090.9kms.fits'
Gfile = '/disks/jansky/a/users/goldston/zheng/151019_NHImaps_SRcorr/data/GNHImaps_SRCORR_final/NHImaps/GALFA-HI_NHI_VLSR-90+90kms.fits'
ghdu = fits.open(Gfile)
#plt.figure(figsize=(20,10))
#plt.imshow(np.sum(ghdu[0].data[:, 100:900, 100:5000], 0))
gwcs = wcs.WCS(Gfile)
xax = np.linspace(1,ghdu[0].header['NAXIS1'], ghdu[0].header['NAXIS1'] ).reshape(ghdu[0].header['NAXIS1'], 1)
yax = np.linspace(1,ghdu[0].header['NAXIS2'], ghdu[0].header['NAXIS2'] ).reshape(1,ghdu[0].header['NAXIS2'])
test = gwcs.all_pix2world(xax, yax, 1)
RA = test[0]
Dec = test[1]
c = SkyCoord(ra=RA*u.degree, dec=Dec*u.degree, frame='icrs')
cg = c.galactic
tt = np.asarray(cg.l.rad)
pp = np.pi/2-np.asarray(cg.b.rad)
#tt = np.asarray(c.ra.rad)
#pp = np.pi/2-np.asarray(c.dec.rad)


Pproj = hp.pixelfunc.get_interp_val(Pdata.T ,pp, tt, nest=False)



ghdu[0].data = Pproj
#outname = path + "mapS-ForSusanClark_fwhm80_ns128_AngSt1" + "_projected_GALFAallsky_RING_T_Gal.fits"
outname = path + "faraday" + "_projected_GALFAallsky_RING_T_Gal.fits"
ghdu.writeto(outname)



