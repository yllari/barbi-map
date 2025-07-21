import numpy as np
from astropy.io import fits 

def get_redd_ext(l:np.ndarray|float, 
        b:np.ndarray|float, 
        r:np.ndarray|float, 
        fits_file,
    ):
    ''' Insert galactic coordinates and return reddening and extinction from Barbillon 2025+ (in prep.) differential reddening map
    The map value and error is obtained by calculating the mean and standard deviation reddening (extinction) for stars inside the
    bin after randomly sampling their reddening (extinction) and position with a gaussian as wide as their error. Since they are 
    phisically limited by 0, no negative values are introduced

    :param l: galactic latitude, from 0 to 360 with 0 pointing towards the galactic center
    :param b: galactic longitude, from -90 to 90 with 0 pointing towards the galactic center
    :param r: maximum radius
    :param E_map: Differential reddening that encodes E_{r+1} - E_r where r identifies increasing radius bins:
    :param Ag_map: Differential extinction map encoded as E_map

    :return: array-like with cumulative reddening and same shape as l, array-like with cumulative extinction and same shape as l
    '''
    def __read_fits_barbi(filename):
        ''' Read processed cumulative redd and extinction map from Barbillon's Gaia data and returns reddening and extinction maps 
        as well as discretization parameters engraved in the header. 
        
        The extinction, reddening (and their error) maps have BINSB X BINSL X BINSR discretization
        NaN values in red or ext represent no stars have been found. Use this information accordingly
        inf values represent zero standard deviation, i.e, no stars or same value in all random samplings. Use this information accordingly

        :param filename:
        :return: red map, ext map, red err map, ext err map, delta_l, delta_b, delta_r
        '''
        with fits.open(filename) as f:
            header = f[0].header
            E_map = f[1].data
            Ag_map = f[2].data

            E_map_err = f[3].data
            Ag_map_err = f[4].data

        # Extracting all information about the discretization
        l_lim0 = header['LIML0']
        l_lim1 = header['LIML1']
        l_bins = header['BINSL']

        delta_l = (l_lim1 - l_lim0)/l_bins

        b_lim0 = header['LIMB0']
        b_lim1 = header['LIMB1']
        b_bins = header['BINSB']

        delta_b = (b_lim1 - b_lim0)/b_bins


        r_lim0 = header['LIMR0']
        r_lim1 = header['LIMR1']
        r_bins = header['BINSR']

        delta_r = (r_lim1 - r_lim0)/r_bins

        discret_deltas = [delta_b, delta_l, delta_r]
        discret_lim0 = [b_lim0, l_lim0, r_lim0] # 0, 0, 0 in default map

        return E_map, E_map_err, Ag_map, Ag_map_err, discret_lim0, discret_deltas


    E_map, E_map_err, Ag_map, Ag_map_err, discret_lim0, discret_deltas = __read_fits_barbi(fits_file)
    b_lim0 , l_lim0, r_lim0 = discret_lim0
    delta_b, delta_l, delta_r = discret_deltas

    # Just to have consistent notation with Barbillon
    # this does not waste memory, its a reference (?)
    theta = b
    phi = l
    delta_theta = delta_b
    delta_phi = delta_l
    theta_lim0 = b_lim0 
    phi_lim0 = l_lim0 


    if type(r) is np.ndarray:
        r_indx = ( (r - r_lim0)/delta_r).astype(int)
    else:
        r_indx = int( (r - r_lim0)/delta_r)

    if type(theta) is np.ndarray:
        theta_indx = ( (theta - theta_lim0)/delta_theta).astype(int)
    else:
        theta_indx = int( (theta - theta_lim0)/delta_theta)

    if type(phi) is np.ndarray:
        phi_indx = ( (phi - phi_lim0)/delta_phi).astype(int)
    else:
        phi_indx = int( (phi - phi_lim0)/delta_phi) 

    E_val = E_map[theta_indx, phi_indx, r_indx]
    A_val = Ag_map[theta_indx, phi_indx, r_indx]
    E_val_err = E_map_err[theta_indx, phi_indx, r_indx]
    A_val_err = Ag_map_err[theta_indx, phi_indx, r_indx]

    return E_val, A_val, E_val_err, A_val_err