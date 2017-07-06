from __future__ import division
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import healpy as hp

#path = "/Volumes/DataDavy/Planck/"
path = "/disks/jansky/a/users/goldston/susan/Planck/"
infile1 = "HFI_SkyMap_353_2048_R2.02_full.fits"
PfileQU = path + infile1 

hdulistQU = fits.open(PfileQU)
hdulistQU.info()
tbdata = hdulistQU[1].data

print(hdulistQU[1].header)

Tdata = tbdata.field('I_STOKES').reshape(49152*1024)
Qdata = tbdata.field('Q_STOKES').reshape(49152*1024)
Udata = tbdata.field('U_STOKES').reshape(49152*1024)
#polfrac = np.sqrt(Qdata**2 + Udata**2)/Tdata

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
#taugalfa = hp.pixelfunc.get_interp_val(tau ,pp, tt, nest=True)

TQUcube = np.zeros([3, ghdu[0].header['NAXIS2'], ghdu[0].header['NAXIS1']])
TQUcube[0, :, :] = hp.pixelfunc.get_interp_val(Tdata ,pp, tt, nest=False)
TQUcube[1, :, :] = hp.pixelfunc.get_interp_val(Qdata ,pp, tt, nest=False)
TQUcube[2, :, :] = hp.pixelfunc.get_interp_val(Udata ,pp, tt, nest=False)

#planckTproj = hp.pixelfunc.get_interp_val(Tdata.T ,pp, tt, nest=False)
#planckQproj = hp.pixelfunc.get_interp_val(Qdata.T ,pp, tt, nest=False)
#planckUproj = hp.pixelfunc.get_interp_val(Udata.T ,pp, tt, nest=False)

ghdu[0].data = TQUcube
ghdu[0].header["NAXIS3"] = 3
outname = path + "HFI_SkyMap_353_2048_R2.02_full" + "_TQUprojected_GALFAallsky.fits"
ghdu.writeto(outname)

