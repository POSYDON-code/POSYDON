__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
]


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def func_1(x,y,z):
    # Spheroid centered at [x,y]=[-1,1]
    radius = 1; a = 1.2; b = 2; c = 0.8
    return  ( (x+1)/a )**2 + ( (y-1)/b )**2 + ( z/c )**2 - radius**2

def func_2(x,y,z):
    # Modified Rosenbrock function
    a = 1; b = 3; c = abs(z) ;
    rosenbrock_func = (a-x)**2 + b*(y-x**2)**2
    return rosenbrock_func + c*(x-z)**2

def func_3(x,y,z):
    return np.sin( (x-0.25) ) + (y)*z

def func_4(x,y,z):
    return np.where( z < 0,  abs((x)**3-(y)**3)*np.exp(-z), 0 )

def func_5(x,y,z):
    # Spheroid centered at [x,y,z]=[0.5,0.5,0.85]
    radius = 0.2; a = 0.6; b = 0.2; c = 0.3;
    return ((x-0.5)/a )**2+( (y-0.5)/b )**2+( (z-0.85)/c )**2 - radius**2

def analytic_classification_3D( X, Y, Z, f_splits=[1, 2.5,-0.5, 1.8, 0.3]):
    """Unique classes 1,2,3,4,6,8 for a total of 6 classes. This nicely matches
    the colormap Dark2 for plotting."""
    ret_f1 = func_1(X,Y,Z)
    ret_f2 = func_2(X,Y,Z)
    ret_f3 = func_3(X,Y,Z)
    ret_f4 = func_4(X,Y,Z)
    ret_f5 = func_5(X,Y,Z)

    where_above_1 = ret_f1 > f_splits[0]
    where_above_2 = ret_f2 > f_splits[1]
    where_above_3 = ret_f3 > f_splits[2]
    where_above_4 = ret_f4 > f_splits[3]
    where_above_5 = ret_f5 > f_splits[4]

    # First layer - purple, yellow
    base_class = np.where( where_above_3, 3, 6)

    # Second layer - gray
    base_class[ where_above_2 ] = 8

    # Third layer - green
    combined_val = np.logical_and( where_above_1, where_above_2 ) # intersect with gray
    base_class[ combined_val ] = 1

    # Fourth layer - pink
    base_class[ where_above_4 ] = 4

    # Fifth layer - orange
    base_class[ np.invert(where_above_5) ] = 2
    return base_class

def analytic_regression_3D(x,y,z):
    # Something simple to start with
    return np.sin( np.sqrt(x**2 + y**2 + z**2) * 2*np.pi/0.75 )


def get_raw_output_3D(X,Y,Z):
    """Get the raw output for classification and regression for
    the 3D synthetic data set. Unique classes include the
    following [1,2,3,4,6,8].

    Returns
    -------
    classification_data : ndarray
        Class data (unique integers) in 3d classification space.
    regression_data : ndarray
        Regression output frmo 3d function.
    """
    if isinstance( X, (float,int) ):
        X = np.array([X]); Y=np.array([Y]); Z=np.array([Z])
    elif isinstance( X, list ):
        X= np.array(X); Y=np.array(Y); Z=np.array(Z)

    classification_data = analytic_classification_3D(X,Y,Z)
    regression_data = analytic_regression_3D(X,Y,Z)
    return classification_data, regression_data


def get_output_3D(X,Y,Z):
    """Takes input points (X,Y,Z) in 3D and returns a DataFrame with the
    corresponding 3D classification and regression functions.

    Returns
    -------
    data_frame : DataFrame
        Pandas dataframe including inputs, class, and outputs.
    """
    if isinstance( X, (float,int) ):
        X = np.array([X]); Y=np.array([Y]); Z=np.array([Z])
    elif isinstance( X, list ):
        X= np.array(X); Y=np.array(Y); Z=np.array(Z)

    classification_data, regression_data = get_raw_output_3D(X,Y,Z)
    if X.ndim > 1 or Y.ndim > 1 or Z.ndim > 1:
        classification_data = classification_data.flatten()
        regression_data = regression_data.flatten()
        X=X.flatten(); Y=Y.flatten(); Z=Z.flatten()
    data_frame = pd.DataFrame()
    data_frame["input_1"] = X
    data_frame["input_2"] = Y
    data_frame["input_3"] = Z
    data_frame["class"] = classification_data.flatten()
    data_frame["output_1"] = regression_data.flatten()
    return data_frame

# EXAMPLE
# x_vals = np.linspace(-1,1,15)
# X, Y, Z = np.meshgrid( x_vals, x_vals+0.01, x_vals-0.02, indexing='ij')
# pd_data = get_output_3D( X,Y,Z )
# print( pd_data )



# -----------------------------------------------------------------------------
###########################
##        PLOTTING       ##
###########################

if (__name__ == "main"):

    # SLICE PLOTS
    all_Z_vals = [-1, -0.5, 0, 0.5, 1]
    fig, (subs_xy, subs_xz) = plt.subplots(nrows = 2, ncols = len(all_Z_vals),
                             figsize=( 4*len(all_Z_vals), 7.5), dpi=50 )

    for ct, zed_val in enumerate(all_Z_vals):
        N = 90
        # X-Y slice
        X, Y = np.meshgrid( np.linspace(-1,1,N), np.linspace(-1,1,N) )
        Z = np.ones( X.shape ) * zed_val
        f_out = analytic_classification_3D( X, Y, Z )

        image = subs_xy[ct].pcolor( X, Y, f_out,
                                     cmap="Dark2", vmin=0.5, vmax=8.5 )
        subs_xy[ct].set_title("z = {:.2f}".format(zed_val), fontsize=14 )

        # X-Z slice
        X, Z = np.meshgrid( np.linspace(-1,1,N), np.linspace(-1,1,N) )
        Y = np.ones( X.shape ) * zed_val
        f_out_2 = analytic_classification_3D( X, Y, Z )
        image = subs_xz[ct].pcolor( X, Z, f_out_2,
                                   cmap="Dark2", vmin=0.5, vmax=8.5 )
        subs_xz[ct].set_title("y = {:.2f}".format(zed_val), fontsize=14 )

        # horizontal lines
        subs_xy[ct].hlines( all_Z_vals, xmin=-1, xmax=1, color="w", linestyle='--', alpha =0.35 )

    # split subplot into two exes and change one to cbar
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    for sub_plot_ax in [subs_xy[ct], subs_xz[ct]]:
        divider = make_axes_locatable( sub_plot_ax )
        cax2 = divider.append_axes("right", size="5%", pad=0.1)
        fig.colorbar(image, cax=cax2)

    subs_xy[0].set_ylabel("Y", fontsize=15, rotation=0)
    subs_xz[0].set_ylabel("Z", fontsize=15, rotation=0)
    for i in range(len(all_Z_vals)):
        subs_xy[i].set_xlabel("X", fontsize=15);
        subs_xz[i].set_xlabel("X", fontsize=15);

    #fig.colorbar(image, ax=subs_xy[ct])
    plt.subplots_adjust(hspace=0.35)
    plt.show()



    # 3D Plots
    from matplotlib.colormaps import get_cmap

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    n = 45
    x_vals = np.linspace(-1,1,n)
    y_vals = np.linspace(-1.01,1,n)
    z_vals = np.linspace(-1,1.01,n)
    X, Y, Z = np.meshgrid( x_vals, y_vals, z_vals, indexing='ij')
    test_f_out = analytic_classification_3D( X, Y, Z )

    # Plotting all points for a given class
    # 4, 1, 8, 3, 6, 2
    for ct, number in enumerate( [4,6,2] ):
        cmap = get_cmap('Dark2')
        rgba = cmap( number/8.5 )
        where_good_pts = np.where( test_f_out == number )
        ax.scatter( X[where_good_pts], Y[where_good_pts], Z[where_good_pts],
                   color=rgba, alpha=0.55, s=10 )

    # Plotting all points in a z slice
    slice_zed_val = 0
    slice_val = np.argmin( abs(slice_zed_val-z_vals) )  # which zed val
    # for ct, number in enumerate( np.unique(test_f_out) ):
    #     cmap = get_cmap('Dark2')
    #     rgba = cmap( number/8.5 )
    #     where_good_pts = np.where( test_f_out[:,:,slice_val] == number )
    #     ax.scatter( X[:,:,slice_val][where_good_pts],
    #                 Y[:,:,slice_val][where_good_pts],
    #                 Z[:,:,slice_val][where_good_pts],
    #                 color=rgba, alpha=0.85, s=10 )

    ax.set_xlim3d(1,-1); ax.set_ylim3d(1,-1); ax.set_zlim3d(-1,1)
    ax.set_xlabel('X', fontsize=15); ax.set_ylabel('Y', fontsize=15); ax.set_zlabel('Z', fontsize=15)
    ax.view_init(elev=89.5, azim=90)
    plt.show()





    # Regression

    x_vals = np.linspace(-1,1,15)
    X, Y, Z = np.meshgrid( x_vals, x_vals+0.01, x_vals-0.02, indexing='ij')
    pd_data = get_output_3D( X,Y,Z )

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter_plot = ax.scatter( X, Y, Z, c=pd_data["output_1"].values, alpha=1, s=10 )
    ax.set_xlim3d(1,-1); ax.set_ylim3d(1,-1); ax.set_zlim3d(-1,1)
    ax.set_xlabel('X', fontsize=15); ax.set_ylabel('Y', fontsize=15); ax.set_zlabel('Z', fontsize=15)
    ax.view_init(elev=89.5, azim=90)

    fig.colorbar(scatter_plot)
    plt.show()
