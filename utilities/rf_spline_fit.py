"""
Contains a function to fit on-axis rf data using a spline. The on-axis data 
should be in two columns, where the first column is the z position and the 
second column is the field data. Writes the fitted data as well as the 
1st, 2nd, and 3rd derivatives to file. Can be called within another python 
script or using the CLI.

For details on using the CLI, call 'python rf_spline_fit.py -h'.

It is expected that the output file will be used with btype 105 in Impact-T.

If the rf_spline_fit function is called from another python script, the input 
dictionary should have the following format:

    fit_dict = {
        "input":"rfdata.in",
        "output":"rfdata1",
        "type":"rf_cavity",
        "s":1e-9,
        "k":5
        }
        
    rf_spline_fit(fit_dict)

The spline_derivative_array function and the default s and k values were 
taken from

https://github.com/ChristopherMayes/openPMD-beamphysics

"""

import numpy as np
import pandas as pd
import argparse
from scipy.interpolate import UnivariateSpline

def spline_derivative_array(z, fz, s=1e-9, k=5):
    # K : degree of smoothing spline
    # s : smoothing factor, s means interpolate through all points
    # Make spline and derivatives
    S = UnivariateSpline(z, fz, k=k, s=s)
    fz1 = S.derivative()(z)
    fz2 = S.derivative(2)(z)
    fz3 = S.derivative(3)(z)
    
    a = np.array([fz, fz1, fz2, fz3]).T
    
    return a

def rf_spline_fit(fit_dict):
    df_input = pd.read_csv(fit_dict['input'],
                           header=None,
                           names=['z','Bz'],
                           sep='\s+')
    
    # Normalize Fields
    df_input.Bz = df_input.Bz/df_input.max().Bz
    
    # fit to spline
    dfield_spline = spline_derivative_array(df_input.z,
                                            df_input.Bz,
                                            s=fit_dict['s'],
                                            k=fit_dict['k'])
    
    # make add header/footer for use with Impact-T btype 105
    if fit_dict['type'].lower()=='rf_cavity':
        header = (
        f'{dfield_spline.shape[0]}   {df_input.min().z}'
        f'{df_input.max().z}   {df_input.max().z}')
        footer = (
            f'1 0.0d0 0.0d0 0.0d0\n'
            f'0.0d0 0.0d0 0.0d0 0.0d0\n')     
    elif(fit_dict['type'].lower()=='solenoid'):
        header = (
            f'1 0.0d0 0.0d0 0.0d0\n'
            f'0.0d0 0.0d0 0.0d0 0.0d0\n'
            f'{dfield_spline.shape[0]}   {df_input.min().z}   '
            f'{df_input.max().z}   {df_input.max().z}')
        footer = ''
    elif(fit_dict['type']==""):
        header = (
        f'{dfield_spline.shape[0]}   {df_input.min().z}   '   
        f'{df_input.max().z}   {df_input.max().z}')   
        footer=''
    else:
        raise Exception("Type of field not specified")
        
    np.savetxt(fit_dict['output'], 
               dfield_spline, 
               fmt='%0.16e',
               header=header,
               footer=footer,
               comments='')
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Applies Spline Fit to on-axis RF data")
    parser.add_argument("-o", 
                        "--output", 
                        type=str, 
                        default="rfdata.out",
                        help=("File name of the output file of the spline fit"
                        "and its derivatives. Default is rfdata.out."))
    parser.add_argument("-i", 
                        "--input", 
                        type=str, 
                        default="rfdata.in",
                        help="File name of the on-axis field data. "
                        "Default is rfdata.in.")
    parser.add_argument("-s", 
                        type=float, 
                        default=1e-9,
                        help=("Positive smoothing factor used to choose the "
                        "number of knots. See scipy.interpolate.UnivariateSpline. "
                        "Default is 1e-9."))
    parser.add_argument("-k", 
                        type=float, 
                        default=5,
                        help=("Degree of the smoothing spline. "
                        "See scipy.interpolate.UnivariateSpline. "
                        "Default is 5"))
    parser.add_argument("-t",
                        "--type",
                        default='',
                        choices=['rf_cavity', 'solenoid',''], 
                        help=("Type of fields being fit. "
                        "Either rf_cavity, solenoid, or "
                        "an empty string '' for neither. Default is ''."))
    args = parser.parse_args()
    
    rf_spline_fit(vars(args)) 