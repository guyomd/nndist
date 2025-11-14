import numpy as np
import pygmt

# Load values:
A = np.loadtxt('p_vs_alpha0.txt', delimiter=';', skiprows=1)
ia = np.argsort(A[:,0])
ALPHA0 = A[ia, 0]  # alpha0 values
BZ10_med = A[ia, 1] # p-values
BZ10_Q2p5 = A[ia, 2]
BZ10_Q97p5 = A[ia, 3]
BZ100_med = A[ia, 4]
BZ100_Q2p5 = A[ia, 5]
BZ100_Q97p5 = A[ia, 6]
KS_med = A[ia, 7]
KS_Q2p5 = A[ia, 8]
KS_Q97p5 = A[ia, 9]

# Make plots:
pygmt.config(MAP_GRID_PEN_PRIMARY="0.5p,black,dotted",
             MAP_GRID_PEN_SECONDARY="0.5p,black,dotted")
fig = pygmt.Figure()
with fig.subplot(3, 1, 
                 sharex=True,
                 subsize=('9c','5c'),
                 clearance=('s1c', 'n1c'),
                 autolabel=True):
    
    # Brown-Zhao test for 10 intervals:
    with fig.set_panel(0):
        fig.basemap(projection='X?/X?l',
                    region=[-1.0, 1.0, 0.005, 1.0],
                    frame=['xafg0.2+l@~\\141@~@-0@-', 'ya2f2g2+lp-value', 'WSne+t@:11:Brown-Zhao test (10 intervals)@::'])
        uncert_polyg = [(x, y) for x, y in zip(ALPHA0, BZ10_Q97p5)]
        uncert_polyg += [(x, y) for x, y in zip(np.flipud(ALPHA0), np.flipud(BZ10_Q2p5))]
        uncert_polyg = np.array(uncert_polyg)

        fig.plot(x=uncert_polyg[:, 0],
                y=uncert_polyg[:, 1],
                close=True,
                fill='gray',
                transparency=40)
        fig.plot(x=ALPHA0,
                 y=BZ10_Q2p5,
                 style='c0.1c',
                 pen='0.5p,black,solid',
                 fill='black')
        fig.plot(x=ALPHA0,
                 y=BZ10_Q2p5,
                 pen='0.5p,black,solid',
                 label='2.5%')
        fig.plot(x=ALPHA0,
                 y=BZ10_Q97p5,
                 style='c0.1c',
                 pen='0.5p,black,solid',
                 fill='black')
        fig.plot(x=ALPHA0,
                 y=BZ10_Q97p5,
                 pen='0.5p,black,solid',
                 label='97.5%')
        fig.plot(x=ALPHA0, 
                y=BZ10_med,
                style='c0.1c',
                pen='0.5p,red,solid',
                fill='red')
        fig.plot(x=ALPHA0, 
                y=BZ10_med,
                pen='1p,red,solid',
                label='median')
    
    # Brown-Zhao test for 100 intervals:
    with fig.set_panel(1):
        fig.basemap(projection='X?/X?l',
                    region=[-1.0, 1.0, 0.005, 1.0],
                    frame=['xafg0.2+l@~\\141@~@-0@-', 'ya2f2g2+lp-value', 'WSne+t@:11:Brown-Zhao test (100 intervals)@::'])
        uncert_polyg = [(x, y) for x, y in zip(ALPHA0, BZ100_Q97p5)]
        uncert_polyg += [(x, y) for x, y in zip(np.flipud(ALPHA0), np.flipud(BZ100_Q2p5))]
        uncert_polyg = np.array(uncert_polyg)
        fig.plot(x=uncert_polyg[:, 0],
                y=uncert_polyg[:, 1],
                close=True,
                fill='gray',
                transparency=40)
        fig.plot(x=ALPHA0,
                 y=BZ100_Q2p5,
                 style='c0.1c',
                 pen='0.5p,black,solid',
                 fill='black')
        fig.plot(x=ALPHA0,
                 y=BZ100_Q2p5,
                 pen='0.5p,black,solid',
                 label='2.5%')
        fig.plot(x=ALPHA0,
                 y=BZ100_Q97p5,
                 style='c0.1c',
                 pen='0.5p,black,solid',
                 fill='black')
        fig.plot(x=ALPHA0,
                 y=BZ100_Q97p5,
                 pen='0.5p,black,solid',
                 label='97.5%')
        fig.plot(x=ALPHA0, 
                y=BZ100_med,
                style='c0.1c',
                pen='0.5p,red,solid',
                fill='red')
        fig.plot(x=ALPHA0, 
                y=BZ100_med,
                pen='1p,red,solid',
                label='median')

    # Kolmogorov-Smirnow test:
    with fig.set_panel(2):
        fig.basemap(projection='X?/X?l',
                    region=[-1.0, 1.0, 0.005, 1.0],
                    frame=['xafg0.2+l@~\\141@~@-0@-', 'ya2f2g2+lp-value', 'WSne+t@:11:Kolmogorov-Smirnov test@::'])
        uncert_polyg = [(x, y) for x, y in zip(ALPHA0, KS_Q97p5)]
        uncert_polyg += [(x, y) for x, y in zip(np.flipud(ALPHA0), np.flipud(KS_Q2p5))]
        uncert_polyg = np.array(uncert_polyg)
        fig.plot(x=uncert_polyg[:, 0],
                y=uncert_polyg[:, 1],
                close=True,
                fill='gray',
                transparency=40)
        fig.plot(x=ALPHA0,
                 y=KS_Q2p5,
                 style='c0.1c',
                 pen='0.5p,black,solid',
                 fill='black')
        fig.plot(x=ALPHA0,
                 y=KS_Q2p5,
                 pen='0.5p,black,solid',
                 label='2.5%')
        fig.plot(x=ALPHA0,
                 y=KS_Q97p5,
                 style='c0.1c',
                 pen='0.5p,black,solid',
                 fill='black')
        fig.plot(x=ALPHA0,
                 y=KS_Q97p5,
                 pen='0.5p,black,solid',
                 label='97.5%')
        fig.plot(x=ALPHA0, 
                y=KS_med,
                style='c0.1c',
                pen='0.5p,red,solid',
                fill='red')
        fig.plot(x=ALPHA0, 
                y=KS_med,
                pen='1p,red,solid',
                label='median')
        

fig.savefig('tests.png', dpi=300)

