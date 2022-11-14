########################
# Tools to search the Nasa Sloan Atlas
# (both versions 0.1.2 and 1.0.1)
#
# v0.1.2 includes some emission line measurements
########################

import numpy as np
h, c, H0 = 0.73, 3*10**5, 73 # km/s

info = ['IAUNAME_1','NSAID_1', 'RA_1', 'DEC_1', 'Z_1', 'ABSMAG','SERSIC_N_1', 'PETROTH50', 'MASS', 'VDISP', 'D4000']

# search by NSAID (data table and NSAID)
def search_n(data, NSAID):
    d = 0
    for i in range(len(data['NSAID'])):
        if data['NSAID'][i] == int(NSAID):
            d = i
            #print(d)
            if d != 0:
                break
    if d == 0:
        return NSAID
    return d

# search by NSAID_0.1.2 for the crossed table  (data table and NSAID)
def search_cross_n(data, NSAID, v=0):
    d = 0
    
    if v==0:
        for i in range(len(data['NSAID_1'])):
            if data['NSAID_1'][i] == int(NSAID):
                d = i
                #print(d)
                if d != 0:
                    break
        if d == 0:
            return NSAID
        return d
    elif v==1:
        for i in range(len(data['NSAID_2'])):
            if data['NSAID_2'][i] == int(NSAID):
                d = i
                #print(d)
                if d != 0:
                    break
        if d == 0:
            return NSAID
        return d
		

# search by IAU name (data table and NAME)
def search_N(data, NAME):
    d = 0
    for i in range(len(data['IAUNAME'])):
        if data['IAUNAME'][i] == NAME:
            d = i
            #print(d)
            if d != 0:
                break
    if d == 0:
        return NAME
    return d


# return galaxy properties from NSA_0.1.2 (data table, row id)
def prin_gal_prop(data, d, verb=True):
    D, D_err = 0, 0
    info = ['IAUNAME','NSAID', 'RA', 'DEC', 'Z', 'ZDIST', 'ABSMAG','SERSIC_N', 'PETROTH50',  'MASS', 'VDISP', 'D4000']
    if verb:
        try:
            for i in info:
                if i == 'ZDIST':
                    D = data[d][i]*c/H0
                    print(i, ':', data[d][i], '-', D, 'Mpc')
                elif i == 'ABSMAG':
                    print(i, 'g-band:', data[d][i][1],'/ FUV-band:', data[d][i][-1])
                elif i == 'PETROTH50':
                    print(i, ':', data[d][i], 'arcsec -', np.pi/180*data[d][i]*D/3.6, 'kpc')
                elif i == 'MASS':
                    print(i, ':', np.log10(data[d][i]/h**2))
                elif i == 'VDISP':
                    if data[d][i] < 70:
                        print(i, ':', data[d][i], 'km/s, value below resolution!')
                    else:
                        print(i, ':', data[d][i], 'km/s')

                else:
                    print(i, ':', data[d][i])
        
        except:
            print('Error in id:', d)
        

    l = [data[d]['IAUNAME'], data[d]['NSAID'], 
         np.round(data[d]['RA'],6), np.round(data[d]['DEC'],6), 
         np.round(data[d]['Z'],4), np.round(data[d]['ZDIST']*c/H0,3), np.round(data[d]['ABSMAG'][1],4), 
         np.round(data[d]['SERSIC_N'],3), np.round(data[d]['PETROTH50'],3), np.round(np.pi/180*data[d]['PETROTH50']*data[d]['ZDIST']*c/H0/3.6,3),
         np.round(np.log10(data[d]['MASS']/h**2),3), np.round(data[d]['VDISP'],3), np.round(data[d]['D4000'],4)]
         
    info = ['IAUNAME','NSAID', 'RA', 'DEC', 'Z', 'ZDIST', 'ABSMAG','SERSIC_N', 'PETROTH50', 'r1/2',  'MASS', 'VDISP', 'D4000']
    line_dict = {i:info[i] for i in range(len(info))} # dict ST : NSA_ID 
    
    return l, line_dict

# Return line measurements from NSA_0.1.2 (data table, row id)
def prin_lines(data, d, verb=True):
    line_info = ['IAUNAME','NSAID',
                 'HAFLUX','HAFLUXERR','HAVMEAS', 'HAVMERR','HAEW', 'HAEWERR', 
                 'HBFLUX','HBFLUXERR', 'HBVMEAS', 'HBVMERR','HBEW', 'HBEWERR', 
                 'O3FLUX','O3FLUXERR', 'O3MEAS','O3VMERR',
                 'O2FLUX','O2FLUXERR', 'O2MEAS','O2VMERR',
                 'O1FLUX','O1FLUXERR', 'O1VMEAS','O1VMERR',
                 'N2FLUX', 'N2FLUXERR', 'N2VMEAS','N2VMERR',
                 'S2FLUX', 'S2FLUXERR', 'S2VMEAS','S2VMERR', 'S2RATIO']
    
    line_dict = {i:line_info[i] for i in range(len(line_info))} # dict ST : NSA_ID 
    
    if verb:
        for i in line_info:
            if i == 'IAUNAME' or i == 'NSAID':
                print(i, ':', data[d][i])

            if i == 'HAFLUX':
                print('-', i, ': (', data[d][i], '+', data[d]['HAFLUXERR'], ')*10^-17 erg/cm^2/s')
                print('FWHM: (', data[d]['HAVMEAS'], '+', data[d]['HAVMERR'], ') km/s')
                print('EW: (', data[d]['HAEW'], '+', data[d]['HAEWERR'], ') A')
            if i == 'HBFLUX':
                print('-', i, ': (', data[d][i], '+', data[d]['HBFLUXERR'], ')*10^-17 erg/cm^2/s')
                print('FWHM: (', data[d]['HBVMEAS'], '+', data[d]['HBVMERR'], ') km/s')
                print('EW: (', data[d]['HBEW'], '+', data[d]['HBEWERR'], ') A')
            if i == 'O3FLUX':
                print('-', i, ': (', data[d]['O3FLUX'], '+', data[d]['O3FLUXERR'], ')*10^-17 erg/cm^2/s')
                print('FWHM: (', data[d]['O3VMEAS'], '+', data[d]['O3VMERR'], ') km/s')
            if i == 'O2FLUX':
                print('-', i, ': (', data[d]['O2FLUX'], '+', data[d]['O2FLUXERR'], ')*10^-17 erg/cm^2/s')
                print('FWHM: (', data[d]['O2VMEAS'], '+', data[d]['O2VMERR'], ') km/s')
            if i == 'O1FLUX':
                print('-', i, ': (', data[d]['O1FLUX'], '+', data[d]['O3FLUXERR'], ')*10^-17 erg/cm^2/s')
                print('FWHM: (', data[d]['O1VMEAS'], '+', data[d]['O1VMERR'], ') km/s')
            if i == 'N2FLUX':
                print('-', i, ': (', data[d][i], '+', data[d]['N2FLUXERR'], ')*10^-17 erg/cm^2/s')
                print('FWHM: (', data[d]['N2VMEAS'], '+', data[d]['N2VMERR'], ') km/s')
            if i == 'S2FLUX':
                print('-', i, ': (', data[d][i], '+', data[d]['S2FLUXERR'], ')*10^-17 erg/cm^2/s')
                print('FWHM: (', data[d]['S2VMEAS'], '+', data[d]['S2VMERR'], ') km/s')
                print('ratio [SII]6731 / [SII]6716 :', data[d]['S2RATIO'])


    l = [data[d]['IAUNAME'], data[d]['NSAID'],
         np.round(data[d]['HAFLUX'],2), np.round(data[d]['HAFLUXERR'],1), np.round(data[d]['HAVMEAS'],1), np.round(data[d]['HAVMERR'],1), np.round(data[d]['HAEW'],1), np.round(data[d]['HAEWERR'],1),
         np.round(data[d]['HBFLUX'],2), np.round(data[d]['HBFLUXERR'],1), np.round(data[d]['HBVMEAS'],1), np.round(data[d]['HBVMERR'],1), np.round(data[d]['HBEW'],1), np.round(data[d]['HBEWERR'],1),
         np.round(data[d]['O3FLUX'],2), np.round(data[d]['O3FLUXERR'],1), np.round(data[d]['O3VMEAS'],1), np.round(data[d]['O3VMERR'],1),
         np.round(data[d]['O2FLUX'],2), np.round(data[d]['O2FLUXERR'],1), np.round(data[d]['O2VMEAS'],1), np.round(data[d]['O2VMERR'],1),
         np.round(data[d]['O1FLUX'],2), np.round(data[d]['O1FLUXERR'],1), np.round(data[d]['O1VMEAS'],1), np.round(data[d]['O1VMERR'],1),
         np.round(data[d]['N2FLUX'],2), np.round(data[d]['N2FLUXERR'],1), np.round(data[d]['N2VMEAS'],1), np.round(data[d]['N2VMERR'],1),
         np.round(data[d]['S2FLUX'],2), np.round(data[d]['S2FLUXERR'],1), np.round(data[d]['S2VMEAS'],1), np.round(data[d]['S2VMERR'],1), np.round(data[d]['S2RATIO'],3)]
    
    return l, line_dict


# return galaxy properties from crossed NSA table (data table, row id)
def prin_gal_prop_cross(data, d, verb=True):
    D, D_err = 0, 0
    
    if verb:
        for i in info:
            if i == 'Z_1':
                D = data[d][i]*c/H0
                print(i, ':', data[d][i], '-', D, 'Mpc')
            elif i == 'ABSMAG':
                print(i, 'g-band:', data[d][i][1],'/ FUV-band:', data[d][i][-1])
            elif i == 'PETROTH50':
                print(i, ':', data[d][i], 'arcsec -', np.pi/180*data[d][i]*D/3.6, 'kpc')
            elif i == 'MASS':
                print(i, ':', np.log10(data[d][i]/h**2))
            elif i == 'VDISP':
                if data[d][i] < 70:
                    print(i, ':', data[d][i], 'km/s, value below resolution!')
                else:
                    print(i, ':', data[d][i], 'km/s')

            else:
                print(i, ':', data[d][i])
        

    l = [data[d]['IAUNAME_1'], data[d]['NSAID_1'], 
         np.round(data[d]['RA_1'],6), np.round(data[d]['DEC_1'],6), 
         np.round(data[d]['Z_1'],4), 
         np.round(data[d]['SERSIC_N_1'],3), np.round(data[d]['PETROTH50'],3), 
         np.round(np.log10(data[d]['MASS']/h**2),3), np.round(data[d]['VDISP'],3), np.round(data[d]['D4000'],4)]
    
    return l
    

# return galaxy properties from NSA_1.0.1 (data table, row id)
def prin_gal2(data, d):
    info2 = ['IAUNAME','NSAID', 'RA', 'DEC', 'Z']
    if type(d) == str:
        print('No matches for', d)
    
    else:
        D, D_err = 0, 0
        for i in info2:
            if i == 'ZDIST':
                D = data[d][i]*c/H0
                print(i, ':', data[d][i], '-', D, 'Mpc')
            elif i == 'SERSIC_ABSMAG':
                print(i, 'g-band:', data[d][i][1],'/ FUV-band:', data[d][i][-1])
            elif i == 'ELPETRO_MASS':
                print(i, ':', np.log10(data[d][i]/h**2))
            elif i == 'ELPETRO_TH50':
                print(i, ':', data[d][i][1],'arcsec -', np.pi/180*data[d][i][1]*D/3.6, 'kpc')
            else:
                print(i, ':', data[d][i])
    
    l = [data[d]['IAUNAME'], data[d]['NSAID'], np.round(data[d]['RA'],6), np.round(data[d]['DEC'],6), np.round(data[d]['Z'],4), np.round(data[d]['ZDIST']*c/H0,3), np.round(data[d]['SERSIC_N'],3), np.round(data[d]['ELPETRO_TH50_R'],3), np.round(np.log10(data[d]['ELPETRO_MASS']/h**2),3)] 
    return l

