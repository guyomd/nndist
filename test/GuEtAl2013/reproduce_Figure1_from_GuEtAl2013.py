import numpy as np
from matplotlib import pyplot as plt



def plotTRsmooth(T, R, colmap='jet', nbins=20, add_markers=False, 
                 TRconst=None, figname=None, bounds = None):
    """
    Produces a smoothed plot of rescaled distance R vs. rescaled time T for a set of space-time-magnitude
    distances obtained using Zaliapin (2007)'s approach. Smoothing is done using a 2-D Gaussian kernel density estimate.

    :param colmap: str, Name of colormap (any  Matplotlib colormap, e.g. 'Spectral_r', 'binary', 'jet', 'Blues', 'viridis').
    :param nbins: int, Number of (hexagonal) bins along X and Y axes for density plot.
    :param add_markers: bool, Specify whether or not to overlay markers for each nearest-neighbor distance
    :param TRconst: float, If specified, add a "T.R=const" line on plot. Value corresponds to its slope in a log-log plot.
    :param figname: str, optional figure name

    :returns matplotlib.pyplot.Figure instance
    """
    from scipy.stats import gaussian_kde
    xlabel = r'Rescaled time, $\log_{10} T$'
    ylabel = r'Rescaled distance, $\log_{10} R$'
    title = 'Distribution of time and space components'
    in0 = np.where((T != 0) & (R != 0))[0]
    print(f'>> {len(R) - len(in0)} events removed due to R or T equal to 0')

    print(min(T))
    print(max(T))
    print(min(R))
    print(max(R))
    x = np.log10(T[in0])
    y = np.log10(R[in0])
    iev = np.isfinite(x) & np.isfinite(y)
    x = x[iev]
    y = y[iev]
    data = np.vstack([x, y])
    kde = gaussian_kde(data)
    if bounds is None:
        bounds = [x.min(), x.max(), y.min(), y.max()]
    print(bounds)

    # evaluate on a regular grid
    xgrid = np.linspace(x.min(), x.max(), nbins)
    ygrid = np.linspace(y.min(), y.max(), nbins)
    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    Z = kde.evaluate(np.vstack([Xgrid.ravel(), Ygrid.ravel()]))

    # Plot the result as an image
    h = plt.figure()
    cMap = plt.cm.get_cmap(colmap)
    """
    plt.imshow(Z.reshape(Xgrid.shape),
        origin='lower', aspect='auto',
        extent=bounds,
        cmap=cMap)
    """
    plt.contourf(Z.reshape(Xgrid.shape),
        origin='lower', #aspect='auto',
        extent=[x.min(), x.max(), y.min(), y.max()],
        cmap=cMap)
    plt.xlim((bounds[0], bounds[1]))
    plt.ylim((bounds[2], bounds[3]))
    hax = plt.gca()
    cbar = plt.colorbar()
    cbar.set_label('Density')

    # If required, add markers:
    if add_markers:
        plt.plot(x, y, '.k', markersize=1)

    # If required, add power-law trend:
    if TRconst is not None:
        Tval = hax.get_xlim()
        Rval = hax.get_ylim()
        hax.plot(Tval, np.subtract(np.log10(TRconst),Tval), 'k--', linewidth=1)
        hax.set_xlim(Tval)
        hax.set_ylim(Rval)
    hax.set_xlabel(xlabel)
    hax.set_ylabel(ylabel)
    hax.set_title(title)

    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
    else:
        h.show()
    return h


def plot_histogram(eta, figname=None, density=False):
    """
    Produce histogram of the nearest-neighbor space-time-magnitude distance SELF.ETA

    :param figname: str, optional figure name
    :returns matplotlib.pyplot.Figure instance
    """
    in0 = np.where((eta != 0) & np.isfinite(eta))[0]
    print(f'>> {len(eta) - len(in0)} events removed due to Nearest-Neighbor distance equal to 0')
    etaLog = np.log10(eta[in0])
    emin = etaLog.min()
    emax = etaLog.max()
    etaLog_xbins = np.linspace(emin, emax, num=50)
    etaLog_ycnts, _ = np.histogram(etaLog, etaLog_xbins, density=density)

    h = plt.figure()
    plt.bar(etaLog_xbins[:-1] + 0.5 * (etaLog_xbins[1:] - etaLog_xbins[:-1]),
            etaLog_ycnts)
    plt.xlabel(r'$\log_{10} \;{\eta}^{\star}$')
    if density:
        plt.ylabel('PDF')
    else:
        plt.ylabel('Counts')
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
    else:
        h.show()
    return h





if __name__ == "__main__":

    configfile = 'config_GuEtAl2013.txt'
    print(f'>> Loading parameters from file "{configfile}"')
    with open(configfile, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line.startswith('#'):
                pass
            else:
                key, value = line.split(':')
                if key == 'parameter_w':
                    w_prm = float(value)
                elif key == 'parameter_p':
                    p_prm = float(value)
                elif key == 'parameter_q':
                    q_prm = float(value)
                elif key == 'parameter_eta0':
                    eta0_prm = float(value)
                elif key == 'parameter_alpha0':
                    alpha0_prm = float(value)
                elif key == 'fractal_dimension':
                    df_prm = float(value)
    print(f'>> Parameters used for nearest-neighbor distance computation:')
    print(f'-- w = {w_prm}')
    print(f'-- p = {p_prm}')
    print(f'-- q = {q_prm}')
    print(f'-- df = {df_prm}')


    results = np.loadtxt('output_GuEtAl2013.txt', delimiter=";", skiprows=1)
    dates = results[1:, 0]
    lats = results[1:, 1]
    lons = results[1:, 2]
    mags = results[1:, 3]
    deps = results[1:, 4]
    nnd = results[1:, 5]
    R = results[1:, 6]
    T = results[1:, 7]
    index_anc = results[1:, 8]
    mag_anc = results[1:, 9]

    # Convert dates and distances into the original units (days and meters):
    b_Gu2013 = 1.09
    df_Gu2013 = 2.0 # see legend of Figure %8 in Gu et al (2013)
    dates *= 365.25 * 24 * 60* 60  

    # Convert T and R values to appropriate units (seconds, meters):
    Tc = T * 365.25 * 24 * 60* 60 * np.power(10, q_prm * w_prm* mag_anc) * np.power(10, -b_Gu2013 * 0.5 * mag_anc) 
    # For R --> convert value using the (df, b) parameters of Gu et al. (2013):
    Rc = np.power(np.power(R * np.power(10, p_prm * w_prm * mag_anc), 1 / df_prm) * 1000, 
                df_Gu2013) * np.power(10, -b_Gu2013 * 0.5 * mag_anc)
    eta = Tc * Rc


    # Make T-R density plot:
    plotTRsmooth(Tc, 
                 Rc, 
                 colmap='gray_r', 
                 nbins=100, 
                 add_markers=False, 
                 TRconst=1E8, 
                 figname='TRdensity_GuEtAl2013.png',
                 bounds=(-1, 7, 0, 10))
    plot_histogram(Tc * Rc, 
                   figname='ETAhisto_GuEtAl2013.png', 
                   density=False)