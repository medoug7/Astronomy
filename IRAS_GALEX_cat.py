# Get IRAS WISE and Galex FUV data

from astroquery.ipac.irsa import Irsa
import astropy.units as u
from astropy.coordinates import SkyCoord as coord
from astroquery.mast import Catalogs
import numpy as np



# function to get W4 (22um) magnitude info from the ALLWISE source catalogue
# Input: position
# Returns: Name, W4 mag, W4 err, angular dist
def get_W4(Ra, Dec, verb=False):
    # query
    table = Irsa.query_region(coord(ra = Ra, dec = Dec,
                            unit=(u.deg, u.deg)), catalog='allwise_p3as_psd', radius=2 * u.arcmin)
    names = table.colnames
    
    if verb:
        # printing
        print(names[0],'\t\t', names[20],'\t',names[21],'\t', names[-2])
        print(table[0][0], '\t', table[0][20],'\t+\t',table[0][21],'\t\t', table[0][-2], 'arcsec')

    return table[0][0], table[0][20], table[0][21], table[0][-2]


# get 22um fluxes and luminosity from position and distance
# Returns W4 Flux, Flux err, Lum, Lum err
def get_W4flux(Ra, Dec, dist):
    # conversion to AB mag is done as explained in 
    # https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
    # with cosmological h = 0.7 correction applied
    
    # monochromatic color
    lam = 22 * u.um # 22 um
    
    if type(Ra) == list:
        W4 = [get_W4(r, Dec[n])[1] + 6.62 - 5*np.log10(0.73) for n,r in enumerate(Ra)]
        W4_err = [get_W4(r, Dec[n])[2] for n,r in enumerate(Ra)]
        
        for i in range(len(W4_err)):
            if W4_err[i] == np.ma.core.MaskedConstant:
                W4_err[i] = 0.5
        
        Flux = [(w*u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz, u.spectral_density(lam)) for w in W4]
        Flux_err = [f*(np.log(10)/2.5)*abs(W4_err[n]) for n,f in enumerate(Flux)]
        
        L = [(4*np.pi*F*lam.to(u.Hz, equivalencies=u.spectral())*(dist[n]*u.Mpc.to(u.cm)*u.cm)**2) for n,F in enumerate(Flux)]
        L_err = [(4*np.pi*F*lam.to(u.Hz, equivalencies=u.spectral())*(dist[n]*u.Mpc.to(u.cm)*u.cm)**2) for n,F in enumerate(Flux_err)]
        
    else:
        W4 = (get_W4(Ra, Dec)[1] + 6.62 - 5*np.log10(0.73))
        W4_err = get_W4(Ra, Dec)[2]
        if W4_err == np.ma.core.MaskedConstant:
            W4_err = 0.05

        Flux = (W4*u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz, u.spectral_density(lam))
        Flux_err = Flux*(np.log(10)/2.5)*abs(W4_err)

        L = 4*np.pi*Flux*lam.to(u.Hz, equivalencies=u.spectral())*(dist*u.Mpc.to(u.cm)*u.cm)**2
        L_err = 4*np.pi*Flux_err*lam.to(u.Hz, equivalencies=u.spectral())*(dist*u.Mpc.to(u.cm)*u.cm)**2

    # flux density in erg/s/cm²/Hz, Luminosity in erg/s
    return Flux, Flux_err, L, L_err
    
    
    
    
# get GALEX FUV magnitudes from position
def get_FUV(Ra, Dec, verb=False):    
    if type(Ra) == list:
        M = []
        M_err = []
        pos = [str(ra)+' '+str(Dec[n]) for n,ra in enumerate(Ra)]
        for n,p in enumerate(pos):
            try:
                l = len(M)
                cat = Catalogs.query_object(p, catalog="Galex")
                for i in range(len(cat)):
                    if type(cat[i]['fuv_mag']) != np.ma.core.MaskedConstant:
                        M.append(cat[i]['fuv_mag'])
                        M_err.append(cat[i]['fuv_magerr'])
                        break
                    else:
                        pass
                
            except:
                M.append('nan')
                M_err.append('nan')
                
            if l == len(M):
                M.append('nan')
                M_err.append('nan')
            
            if verb:
                print(n+1, '-', round(100*(n+1)/len(pos)), '% done')
            
    else:
        pos = str(Ra)+' '+str(Dec)
        try:
            cat = Catalogs.query_object(pos, catalog="Galex")[:3]
            for i in range(3):
                if type(cat[i]['fuv_mag']) != np.ma.core.MaskedConstant:
                    M = cat[i]['fuv_mag']
                    M_err = cat[i]['fuv_magerr']
                    break
                else:
                    pass
        except:
            M = 'nan'
            M_err = 'nan'
            
    return M, M_err


# get FUV flux and luminosities from magnitudes and distances
def get_FUVflux(mag, mag_err, dist):
    # conversion to AB mag is done as explained in 
    # https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
    
    # color band
    lam = (1350+1750)/2 * u.AA # 1350-1750 angstroms
        
    if type(mag) == list:
        for i in range(len(mag)):
            if mag[i] == 'nan':
                mag[i] = 999
                mag_err[i] = 0
        
        Flux = [((m - 5*np.log10(0.73))*u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz, u.spectral_density(lam)) for m in mag]
        Flux_err = [f*(np.log(10)/2.5)*abs(mag_err[n]) for n,f in enumerate(Flux)]
        
        L = [(4*np.pi*F*lam.to(u.Hz, equivalencies=u.spectral())*(dist[n]*u.Mpc.to(u.cm)*u.cm)**2) for n,F in enumerate(Flux)]
        L_err = [(4*np.pi*F*lam.to(u.Hz, equivalencies=u.spectral())*(dist[n]*u.Mpc.to(u.cm)*u.cm)**2) for n,F in enumerate(Flux_err)]
        
    else:
        
        if mag == 'nan':
            mag = 999
            mag_err = 0

        Flux = ((mag - 5*np.log10(0.73))*u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz, u.spectral_density(lam))
        Flux_err = Flux*(np.log(10)/2.5)*abs(mag_err)


        L = 4*np.pi*Flux*lam.to(u.Hz, equivalencies=u.spectral())*(dist*u.Mpc.to(u.cm)*u.cm)**2
        L_err = 4*np.pi*Flux_err*lam.to(u.Hz, equivalencies=u.spectral())*(dist*u.Mpc.to(u.cm)*u.cm)**2

    # flux density in erg/s/cm²/Hz, Luminosity in erg/s
    return Flux, Flux_err, L, L_err
    
    
