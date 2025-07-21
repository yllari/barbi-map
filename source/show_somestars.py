from calc_ext import get_redd_ext
import numpy as np
import vaex
import os
import matplotlib.pyplot as plt




if __name__ == '__main__':
    # ------------------------- Load all data from ED23 Halo sample (approx l-b-uniformly distributed) ---------------------
    df = vaex.open(os.path.join('samples', 'halo_l22.hdf5'))
    print(df.column_names)
    df_bayes = vaex.open(os.path.join('samples','halo_bayes.hdf5'))[['source_id', 'E_BP_RP', 'AG']]
    df_bayes.rename('E_BP_RP', 'E_BP_RP_bayes')
    df_bayes.rename('AG', 'AG_bayes')


    df_total = df.join(df_bayes, how='left', on='source_id')
    df_total.rename('E_BP_RP', 'E_BP_RP_l22')
    df_total.rename('AG', 'AG_l22')

    l = np.array(df_total['l'].values)
    b = np.array(df_total['b'].values) 
    dist = np.array(df_total['distance'].values)
    red, ag, red_err, ag_err = get_redd_ext(l,b,dist, 'cumul_red_ag_barbillon.fits')
    print(len(df_bayes['E_BP_RP_bayes'].values))


    # Values with reddening (filters nan or extremely high values)
    print('TOTAL', len(df_total))
    print('WITH RED in barbi', np.sum( (red < 5000) ))
    print('with red in bayes', np.sum((df_bayes['E_BP_RP_bayes'].values < 5000)))
    print('with red in l22', np.sum((df_total['E_BP_RP_l22'].values < 5000)))


    red_l22 = np.array(df_total['E_BP_RP_l22'].values)
    red_bayes = np.array(df_total['E_BP_RP_bayes'].values)
    ag_l22 = np.array(df_total['AG_l22'].values)
    ag_bayes = np.array(df_total['AG_bayes'].values)



    from matplotlib.colors import LogNorm
    plt.hist2d(red, ag, bins=100, range=[[0,3], [0,6]], norm=LogNorm())
    plt.xlabel('Reddening')
    plt.ylabel('AG')
    plt.annotate('This must be an almost-linear with almost no dispers. relation', (0.1, 5))
    plt.savefig('plots/red_ag.png', dpi=300)


    # ------------------------------------------ Comparison between bayestar, l22, and Barbi -----------------
    fig, axs = plt.subplots(2,3)

    upper_lim = 1.5 # Just for plotting
    fig.suptitle('Redd comparison')

    # Barbillon vs l22
    axs[0,0].set_title('Barbi - l22')
    axs[0,0].hist2d(red, red_l22, bins=100, range=[[0, upper_lim], [0,upper_lim]]       , norm=LogNorm() )
    axs[0,0].plot([0,upper_lim], [0, upper_lim], label='1-to-1', color='black')
    axs[0,0].legend()
    axs[0,0].set_xlabel('Barbi')
    axs[0,0].set_ylabel('l22')

    # Barbillon vs Bayestar
    axs[0,1].set_title('Barbi - bayes')
    axs[0,1].hist2d(red, red_bayes, bins=100, range=[[0, upper_lim], [0, upper_lim]]    , norm=LogNorm() )
    axs[0,1].plot([0,upper_lim], [0, upper_lim], color='black')
    inset = axs[0,1].inset_axes([0.7,0.1, 0.4,0.4])
    inset.hist2d(red, red_bayes, bins=30, range=[[0, 0.15], [0, 0.15]]    , norm=LogNorm() )
    inset.plot([0,0.15], [0, 0.15], color='black')
    axs[0,1].indicate_inset_zoom(inset, edgecolor="black")
    axs[0,1].set_xlabel('Barbi')
    axs[0,1].set_ylabel('bayes')

    # l22 vs Bayestar 
    axs[0,2].set_title('l22 - bayes')
    axs[0,2].hist2d(red_l22, red_bayes, bins=100, range=[[0, upper_lim], [0, upper_lim]], norm=LogNorm() )
    axs[0,2].plot([0,upper_lim], [0, upper_lim], color='black')
    axs[0,2].set_xlabel('l22')
    axs[0,2].set_ylabel('bayes')

    def plot_hist_med_qs(data1, data2, ax, lim):
        diff = data1 - data2
        ax.hist(diff, bins=100, range=[-lim, lim] )
        perc_25 = np.nanpercentile(diff, 25)
        median = np.nanmedian(diff)
        perc_75 = np.nanpercentile(diff, 75)
        ax.axvline(perc_25, label=f"25% {'{:.2f}'.format(perc_25)}", color='red', linestyle='--')
        ax.axvline(median, label=f"50% {'{:.2f}'.format(median)}", color='red')
        ax.axvline(perc_75, label=f"75% {'{:.2f}'.format(perc_75)}", color='red', linestyle='--')
        ax.legend()
        return ax

    plot_hist_med_qs(red, red_l22, axs[1,0], upper_lim/5)
    axs[1,0].set_xlabel('Barbi - l22')
    plot_hist_med_qs(red, red_bayes, axs[1,1], upper_lim/5)
    axs[1,1].set_xlabel('Barbi - bayes')
    plot_hist_med_qs(red_l22, red_bayes, axs[1,2], upper_lim/5)
    axs[1,2].set_xlabel('l22 - bayes')
    fig.set_size_inches((12,6))

    fig.savefig('plots/maps_comparison.png', dpi=300)
