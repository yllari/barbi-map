import numpy as np
from astropy.io import fits 
from astropy.table import Table
import matplotlib.pyplot as plt
import vaex

class DiscretizeReddMapMC():
    ''' Class to generate the discretized map with mean and std values given
    by MC sampling from the error coming from individual stars after binning.
    This process includes poissonan effects and uncertaintites.

    Attributes:
        cumul_data: information about all stars with extinction data
        outfile: 
        rel_err_lim:  relative error limit
        n_bins_b: number of bins to use in b coordinate 
        n_bins_l: number of bins to use in l coordinate 
        n_bins_r: number of bins to use in r coordinate 
        r_lim_min: lower limit of distance to be imposed in calculate_map
        r_lim: upper limit of distance to be imposed in calculate_map
        verbose: whether to do plots and prints, defaults False

    
    Methods:
        calculate_map: calculates the discretized reddening and extinction map and 
        stores it in H_tot_(red,ag) and H_tot_(red,ag)_err
        save_map: saves the map to outfile with information in the header about the
        discretization


    
    '''

    def __init__(self, file_info, outfile, rel_err_lim=0.5, 
                 n_bins_b=180, n_bins_l=360, n_bins_r=15, 
                 r_lim_min=0, r_lim=2.5,
                 verbose=False):
        '''
        :param file_info: direction to the fits data with all information about the cumulative reddening 
            The format should include: dist -- distance in pc
                                        dist_err -- distance error in pc
                                        b -- longitude from -90 to 90 
                                        l -- latitude from 0 to 360
                                        E_bprp -- Reddening in Gaia bands 
                                        E_bprp_err -- Reddening error 
                                        Ag -- Extinction in Gaia G band
                                        Ag_err -- Extinction error in G band
        :param outfile: where to write the out fit file
        :param rel_err_lim: relative error limit to be imposed to extinction and reddening error/value < rel_err_lim 
        :param n_bins_b: number of bins to use in b coordinate 
        :param n_bins_l: number of bins to use in l coordinate 
        :param n_bins_r: number of bins to use in distance coordinate
        :param r_lim_min: lower limit of distance to be imposed
        :param r_lim: upper limit of distance to be imposed
        :param verbose: whether to do plots and prints, defaults False
        '''
        self.cumul_data = fits.open(file_info)[1].data
        self.outfile = outfile 
        self.cumul_data['b'] = self.cumul_data['b'] # -90 to 90
        self.cumul_data['dist'] /= 1000 # to kpc
        self.cumul_data['dist_err'] /= 1000 # to kpc

        self.rel_err_lim = rel_err_lim

        self.r_lim_min = r_lim_min
        self.r_lim = r_lim
        self.n_bins_b = n_bins_b
        self.n_bins_l = n_bins_l
        self.n_bins_r = n_bins_r

        self.verbose = verbose

        if self.verbose:
            columns = self.cumul_data.columns
            print('COLUMNS', columns)
            print('LEN', len(self.cumul_data))

    def calculate_map(self, niters:int, avoid_ind_neg_vals:bool = False):
        ''' Calculate the map by randomly sampling the values of reddening, extinction and 
        distance as a clipped gaussian (-2,2) of the errors. Then, the mean an std values are given 
        for both reddening and extinction. 

        WARNING: INDIVIDUAL VALUES OF THIS PROCEDURE ARE NOT CORRECT, SINCE EXTINCTION AND REDDENING ARE 
        NOT INDEP. VARIABLES, BUT WE ARE CONSIDERING AVERAGE (STD) STATISTICAL VALUES IN EACH BIN.

        :param niters: number of MC extractions for each star, around 100 is okay
        :param avoid_ind_neg_vals: whether clip negative values of reddening to zero during the extractions. 
                                    At the end, all average values inside the binning are fixed to >=0, but this
                                    fixex INDIVIDUAL stars behaviour (not very important, around ~ 10 stars in all sample suffer this)
        '''

        # Cut in relative error for reddening and extinction
        self.cumul_data = self.cumul_data[(self.cumul_data['E_bprp_err']/self.cumul_data['E_bprp'] < self.rel_err_lim) ]
        self.cumul_data = self.cumul_data[(self.cumul_data['Ag_err']/self.cumul_data['Ag'] < self.rel_err_lim) ]

        # filter by distance previous to histogram, just to clean a bit
        filter_bydist = self.cumul_data['dist'] < r_lim + 0.5

        self.cumul_data = self.cumul_data[filter_bydist]

        x = self.cumul_data['b']
        y = self.cumul_data['l']
        z = self.cumul_data['dist'] 
        z_err = self.cumul_data['dist_err']

        red = np.array(self.cumul_data['E_bprp'])
        red_err = np.array(self.cumul_data['E_bprp_err'])
        ag = np.array(self.cumul_data['Ag'])
        ag_err = np.array(self.cumul_data['Ag_err'])


        df = vaex.from_astropy_table(Table(self.cumul_data))
        df = df[['b','l','dist','E_bprp', 'Ag']]
        df.rename('E_bprp', 'red')
        df.rename('Ag', 'ag')

        df_list = [df]
        for i in range(niters):
            if self.verbose:
                print('MC extraction', i)
            z_= z + np.clip(np.random.normal(size=z.shape), -2,2)*z_err
            red_ = red + np.clip(np.random.normal(size=red.shape), -2,2)*red_err
            ag_ = ag + np.clip(np.random.normal(size=ag.shape), -2,2)*ag_err

            if avoid_ind_neg_vals:
                red_[red_ < 0] = 0
                ag_[ag_ < 0] = 0
                # z outside are not taken into account either way

            values = [x, y, z_, red_, ag_]
            columns = ['b', 'l', 'dist', 'red', 'ag']
            out_dict = dict(zip(columns, values))
            df_list.append( vaex.from_dict(out_dict))


        # Vaex
        df = vaex.concat(df_list)
        self.H_tot_red = df.mean(df.red, binby=[df.b, df.l, df.dist], limits=[[-90, 90], [0,360], [self.r_lim_min, self.r_lim]], shape=(self.n_bins_b, self.n_bins_l, self.n_bins_r))
        self.H_tot_ag = df.mean(df.ag, binby=[df.b, df.l, df.dist], limits=[[-90, 90], [0,360], [self.r_lim_min, self.r_lim]], shape=(self.n_bins_b, self.n_bins_l, self.n_bins_r))

        self.H_tot_red_err = df.std(df.red, binby=[df.b, df.l, df.dist], limits=[[-90, 90], [0,360], [self.r_lim_min, self.r_lim]], shape=(self.n_bins_b, self.n_bins_l, self.n_bins_r))
        self.H_tot_ag_err = df.std(df.ag, binby=[df.b, df.l, df.dist], limits=[[-90 , 90], [0,360], [self.r_lim_min, self.r_lim]], shape=(self.n_bins_b, self.n_bins_l, self.n_bins_r))


        # Ensure only positive values
        self.H_tot_red[self.H_tot_red < 0 ] = 0
        self.H_tot_red[self.H_tot_ag < 0 ] = 0

        # No variation means no variation in stars, cannot say anything
        # about the error
        self.H_tot_red_err[self.H_tot_red_err == 0] = np.inf
        self.H_tot_ag_err[self.H_tot_ag_err == 0] = np.inf

        if self.verbose:
            fig, ax = plt.subplots(1, 2)
            ax[0].imshow(self.H_tot_red[:,:,-1], extent=[0,360,-90,90] , interpolation=None)
            ax[0].set_title('Reddening at max dist')
            ax[1].imshow(self.H_tot_ag_err[:,:,-1], extent=[0,360,-90,90], interpolation=None)
            ax[1].set_title('Reddening error at max dist')
            ax[0].set_xlabel('l')
            ax[1].set_xlabel('l')

            ax[0].set_ylabel('b')
            ax[1].set_ylabel('b')
            plt.show()


    def save_map(self):
        '''Save the map in fits format

        The fits structure is:
        - primaryHDU, general information 
        - ImageHDU, Cumulative reddening map 
        - ImageHDU, Cumulative extinction map 
        - ImageHDU, Cumulative reddening error map 
        - ImageHDU, Cumulative extinction error map 

        The primary HDU stores the number of bins as well as the limits used 
        for the discretization. 
        
        '''
        if self.H_tot_red is None:
            print('You should run calculate_map first!!')
        else:
            # Save mask and engrave all info on header
            primary_HDU = fits.PrimaryHDU()
            hdr = primary_HDU.header
            hdr['BINSB'] = self.n_bins_b
            hdr.comments['BINSB'] = 'Number fo pixels in galactic latitude, b'
            hdr['BINSL'] = self.n_bins_l
            hdr.comments['BINSL'] = 'Number of pixels in galactic longitude, l'
            hdr['BINSR'] = self.n_bins_r
            hdr.comments['BINSR'] = 'Number fo pixels in distance to the sun, r'


            hdr['LIMB0'] = -90
            hdr.comments['LIMB0'] = 'Starting angle for galactic latitude, b'
            hdr['LIMB1'] = 90
            hdr.comments['LIMB1'] = 'Ending angle for galactic latitude, b'

            hdr['LIML0'] = 0
            hdr.comments['LIML0'] = 'Starting angle for galactic longitude, l'
            hdr['LIML1'] = 360
            hdr.comments['LIML1'] = 'Ending angle for galactic longitude, l'

            hdr['LIMR0'] = self.r_lim_min
            hdr.comments['LIMR0'] = 'Starting radius for distance to the sun, r'
            hdr['LIMR1'] = self.r_lim
            hdr.comments['LIMR1'] = 'Ending radius for distance to the sun, r'

            tab_1 = fits.ImageHDU(self.H_tot_red)
            tab_2 = fits.ImageHDU(self.H_tot_ag)
            tab_3 = fits.ImageHDU(self.H_tot_red_err)
            tab_4 = fits.ImageHDU(self.H_tot_ag_err)

            hdu_list = fits.HDUList([
                primary_HDU,
                tab_1,
                tab_2,
                tab_3,
                tab_4
            ])

            hdu_list.writeto(self.outfile, overwrite=True)



if __name__ == '__main__':
    file_dir = 'map/Yllari_selected_data.fits'
    outfile = 'cumul_red_ag_barbillon.fits'
    rel_err_lim = 0.5
    n_bins_r = 15
    n_bins_l = int(360)
    n_bins_b = int(180)
    r_lim = 2.5 # 2.5 kpc limit to minimize the effect of high extinction and low complet. in non-penetrated dust clouds
    r_lim_min = 0 

    niters = 10
    disc_map = DiscretizeReddMapMC(file_dir, outfile, 
                        rel_err_lim, 
                        n_bins_b, n_bins_l, n_bins_r, 
                        r_lim_min, r_lim, verbose=True)
    
    disc_map.calculate_map(niters)
    disc_map.save_map()