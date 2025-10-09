import numpy as np
import pygmt




if __name__ == "__main__":

    #configfile = 'config_ZBZ2020_alpha0.txt'  # Configuration with parameter alpha0 = 0.0
    configfile = 'config_ZBZ2020_alpha0.1.txt'  # Configuration with parameter alpha0 = 0.1

    #resultsfile = 'output_declust_ZBZ2020_alpha0.txt'  # Results obtained with parameter alpha0 = 0.0
    resultsfile = 'output_declust_ZBZ2020_alpha0.1.txt'  # Results obtained with parameter alpha0 = 0.1


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


    results = np.loadtxt(resultsfile, delimiter=";", skiprows=1)
    dates = results[:, 0]
    lats = results[:, 1]
    lons = results[:, 2]
    mags = results[:, 3]
    deps = results[:, 4]
    nnd = results[:, 5]
    R = results[:, 6]
    T = results[:, 7]
    index_anc = results[:, 8]
    mag_anc = results[:, 9]
    p_bgnd = results[:, 10]; 
    is_bgnd = results[:, 11]; 
    norm_prox = results[:, 12]; 
    avg_nn_distance = results[:, 13]
    ib = is_bgnd == 1


    # Compute proportion and number of background events:
    nall = results.shape[0]
    nbg = is_bgnd.sum()
    nbg_w = p_bgnd.sum()
    prop_bgnd_w = nbg_w / nall
    prop_bgnd_n = nbg / nall 
    print(f'>> Number of background events:\n\t- counts in realization: {nbg}\n\t- in weights: {nbg_w}')
    print(f'>> Proportion of background events:\n\t- counts  in realization: {prop_bgnd_n * 100:.2f}%\n\t- in weights: {prop_bgnd_w * 100:.2f}%')

    # Map declustered catalogue realization
    fig = pygmt.Figure()
    fig.basemap(projection='X10c',
                region=[0, 600, 0, 600],
                frame=['WSne', 'xaf100+lx-coordinate, km', 'yaf100+ly-coordinate, km'])
    fig.plot(x=lons[ib],
             y=lats[ib],
             style='p0.05c',
             fill='black')
    fig.savefig('MapZBZ2020.png', dpi=300)

    # Plots X vs. Time (declustered catalogue realization)
    fig = pygmt.Figure()
    fig.basemap(projection='X15c/8c',
                region=[1, 22, 0, 600],
                frame=['WSne', 'xaf1+lTime, yrs', 'yaf100+ly-coordinate, km'])
    """
    fig.plot(x=dates,
             y=lons,
             style='p0.03c',
             fill='grey')
    """
    fig.plot(x=dates[ib],
             y=lats[ib],
             style='p0.05c',
             fill='black')
    fig.savefig('TX_ZBZ2020.png', dpi=300)
