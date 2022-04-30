def calculate_Tg(result_fname: str, make_plot: bool = True) -> int:
    '''Method to calculate glass transition temperature based on the
    result file obtained from TgMeasurement Procedure

    Parameters:
        result_fname (str): Name of the result file from TgMeasurement 
                            Procedure
        make_plot (bool): Whether to make a plot to visualize the fitting

    Returns:
        Tg (int): Glass transition temperature of the system
    '''

    from scipy import optimize
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    def piecewise_linear(x, x0, y0, k1, k2):
        return np.piecewise(
            x, [x < x0],
            [lambda x: k1 * x + y0 - k1 * x0, lambda x: k2 * x + y0 - k2 * x0])

    # Read result file (should have first column as temperature,
    # and second column as the corresponding density)
    df = pd.read_csv(result_fname,
                     header=1,
                     usecols=[1, 2],
                     names=['Temp', 'Rho'],
                     delim_whitespace=True)
    x = np.array(df['Temp'])
    y = np.array(df['Rho'])

    p0 = (250, 0.85, -0.0001, 0.0001)
    p, e = optimize.curve_fit(piecewise_linear, x, y, p0)
    xd = np.linspace(x.min(), x.max(), 100)

    if make_plot:
        plt.rcParams.update({'font.size': 12})
        plt.plot(x, y, "ko", fillstyle='none', markersize=8)
        plt.plot(xd, piecewise_linear(xd, *p), 'r')
        plt.ylabel('Density $(g/cm^3)$')
        plt.xlabel('Temperature ($K$)')
        plt.savefig('temp_vs_density.png', dpi=300)

    Tg = p[0]
    print('Glass transition temperature of the system is', Tg)

    return Tg
