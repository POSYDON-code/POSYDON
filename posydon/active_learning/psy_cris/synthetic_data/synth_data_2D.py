__authors__ = [
    "Kyle Akira Rocha <kylerocha2024@u.northwestern.edu>",
]


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if (__name__ == "main"):
    save_figs = True   # plots will be saved
    save_data = True   # saving sampled data to csv
else:
    save_figs = False
    save_data = False

############# Sample Params #############
Nx = 10; Ny = 10
file_name = "synth_data.dat"


############# Cls / Regr Geometry #############
def cls_curve(x,y):
    vals = regr_func(x,y)
    curve_bool = ( y - (-0.1)*(x**3) - 1 ) > -1  # x^3 curve
    line = np.where( curve_bool, 2, vals )

    bound_1 = 0.4 ; bound_2 = -0.28
    cls1 = np.where( vals > bound_1, 1, line ) # banana
    cls2 = np.where( (cls1 < bound_1) == (cls1 > bound_2), 0, cls1 ) # background
    cls3 = np.where( vals < bound_2, -1, cls2 ) # triangle
    return cls3

def regr_func(x,y):
    eps = .9    # ~ 1/r padding
    a = 1.  ;  b = 1.   # ~ ellipse
    return  ( np.tanh((x)**3+(y)**2) )/( np.sqrt( (x/a)**2+(y/b)**2) + eps)
################################################

x = np.linspace(-3,3, 140)
y = np.linspace(-3,3, 140)
X,Y = np.meshgrid(x,y)

# %matplotlib notebook
# %matplotlib inline
if (__name__ == "main"):
    print("Making analytic plot...")
    fig, subs = plt.subplots( nrows=1, ncols=2, dpi=120, figsize=(10,4),
                                gridspec_kw={'width_ratios': [0.85, 1]})
    subs[0].set_title("Classification")
    subs[0].pcolormesh( X, Y, cls_curve(X,Y) )

    subs[1].set_title("Regression")
    pcm = subs[1].pcolormesh(X,Y, regr_func(X,Y), cmap="PiYG") #PiYG
    fig.colorbar( pcm )
    my_cont = subs[1].contour(X,Y, regr_func(X,Y), colors='k',
                          levels=[-0.4, -0.3, -0.15, 0, 0.15, 0.3, 0.4], )
    my_cont.clabel( fmt='%1.2f', )

    for i in range(2):
        subs[i].set_xlabel(r"$x$", fontsize=12); subs[i].set_ylabel(r"$y$", fontsize=12)
    if save_figs: plt.savefig("cls_regr_original.png")
    plt.show()


############# Sampling #############
# x,y vals to query with classifcationa & regression functions
x_vals = np.linspace(-3,3, int(Nx))
y_vals = np.linspace(-3,3, int(Ny))
if (__name__=="main") and save_data:
    print("Sampling {0} points...".format(len(x_vals)*len(y_vals)) )

points = []
for i in x_vals:
    for j in y_vals:
        points.append( np.array([i,j]) )
points = np.array(points)

results = cls_curve( points.T[0], points.T[1] )
unique_classes = np.unique( results )
mapping = {-1:"A", 0:"B", 1:"C", 2:"D"}
str_classification = [ mapping[val] for val in results.astype(int)]

df = pd.DataFrame()
df["input_1"] = points.T[0]
df["input_2"] = points.T[1]
df["class"] = str_classification
df["output_1"] = regr_func( points.T[0], points.T[1] )

if save_data:
    df.to_csv( file_name, index = False );
    print("Saved to '{0}'.".format(file_name))


if (__name__ == "main"):
    print("Making sampled points plot...")
    fig, subs = plt.subplots( nrows=1, ncols=2, figsize=(8,4), dpi=120 )
    subs[0].set_title("Sampled Points")
    subs[0].scatter( df["input_1"], df["input_2"], c = results, s=40, )

    subs[1].set_title("Overlay")
    subs[1].contourf( X, Y, cls_curve(X,Y), levels = 20, vmin=-1, vmax=2, alpha=1. )
    #subs[1].pcolormesh( X, Y, cls_curve(X,Y), levels = 10 )
    subs[1].scatter( df["input_1"], df["input_2"], c = results, s=40, edgecolors='w'  )

    for i in range(2):
        subs[i].set_xlim(-3.1,3.1); subs[i].set_ylim(-3.1,3.1)
        subs[i].set_xlabel(r"$x$", fontsize=12); subs[i].set_ylabel(r"$y$", fontsize=12)
    if save_figs: plt.savefig("cls_regr_sampled.png")
    plt.show()


def get_raw_output_2D(x,y):
    """Get the raw output for classification and regression functions
    for the 2D synthetic data set. Original class data to strings given by
    the following relation {-1:"A", 0:"B", 1:"C", 2:"D"}. """
    if isinstance( x, (float,int) ):
        x = np.array([x]); y=np.array([y])
    elif isinstance(x, list):
        x= np.array(x); y=np.array(y)

    classification_output = cls_curve(x,y)
    regression_output = regr_func(x,y)
    return classification_output, regression_output

# For queries
def get_output_2D(x,y):
    """For a set of query points (x,y) in the range (-3,3)
    return a DataFrame with the inputs and outputs for
    classificaiton and regression.
    """
    if isinstance( x, (float,int) ):
        x = np.array([x]); y=np.array([y])
    elif isinstance(x, list):
        x= np.array(x); y=np.array(y)

    cls_results, regr_out = get_raw_output_2D(x,y)
    if x.ndim > 1 or y.ndim > 1:
        cls_results = cls_results.flatten()
        regr_out = regr_out.flatten()
        x = x.flatten(); y = y.flatten()
    str_results = np.array( [mapping[val] for val in cls_results.astype(int)] )
    data_frame = pd.DataFrame()
    data_frame["input_1"] = x
    data_frame["input_2"] = y
    data_frame["class"] = str_results
    data_frame["output_1"] = regr_out
    return data_frame
