import os
import warnings
import pprint
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import numpy as np
import pickle
import pandexo.engine.justplotit as jpi

from pandexo.engine.load_modes import SetDefaultModes

import matplotlib.pyplot as plt
'''
def run_pandexo(output,spectra):
    
    """
    manual changes are done to catalog.py and config.py in pandexo to make this code work
    """

    exo_dict = jdi.load_exo_dict()
    exo_dict['observation']['sat_unit'] = "%"
    exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well 
    exo_dict['observation']['noccultations'] = 10 #number of transits 
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later 
                                                 #without having to redo the calcualtion
    exo_dict['observation']['fraction'] = 1.0    #fraction of time in transit versus out = in/out
    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath 
    
    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 8.0               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 4000             #in K 
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0  

    #exo_dict['planet']['type'] = 'constant'
    #exo_dict['planet']['depth'] = 0.01                      #other options include "um","nm" ,"Angs", "secs" (for phase curves)
    exo_dict['planet']['exopath'] = spectra
    exo_dict['planet']['w_unit'] = "um"
    exo_dict['planet']['f_unit'] = "fp/f*"#'rp^2/r*^2'
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0 
    
    print('Starting TEST run')
    jdi.run_pandexo(exo_dict, ['NIRSpec G140M'], output_file=output)
    print('SUCCESS') 

def read_output():
    
    
    out = pickle.load(open("ref_trans_sim2.p", 'rb'))
    
    print out["FinalSpectrum"]["wave"]
    print out["FinalSpectrum"]["spectrum"]
    
    return 
    m1,n1,x1,y1,e1 = jpi.jwst_1d_spec(out, R=100, num_tran=100, model=True, x_range=[1,5])
    
    plt.plot(m1,n1)
    plt.errorbar(x1,y1,e1)
    plt.show()
    
'''


def try_pandexo_tutorial(spectra):
    
    exo_dict = jdi.load_exo_dict()
    exo_dict['observation']['sat_unit'] = "%"
    exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
    exo_dict['observation']['noccultations'] = 2 #number of transits
    exo_dict['observation']['R'] = None          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                                 #without having to redo the calcualtion
    #exo_dict['observation']['fraction'] = 1.0    #fraction of time in transit versus out = in/out
    
    exo_dict['observation']['baseline'] = 4.0*60.0*60.0    #time spent observing out of transit, make sure to speciy units
    exo_dict['observation']['baseline_unit'] = 'total'    #fraction of time in transit versus out = in/out 
    
    exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                                 #to a wavelength dependent noise floor solution (units are ppm)
    exo_dict['star']['type'] = 'phoenix'        #phoenix or user (if you have your own)
    exo_dict['star']['mag'] = 8.0               #magnitude of the system
    exo_dict['star']['ref_wave'] = 1.25         #For J mag = 1.25, H = 1.6, K =2.22.. etc (all in micron)
    exo_dict['star']['temp'] = 5500             #in K
    exo_dict['star']['metal'] = 0.0             # as log Fe/H
    exo_dict['star']['logg'] = 4.0              #log surface gravity cgs        
        
    exo_dict['planet']['exopath'] = spectra
    exo_dict['planet']['w_unit'] = 'um'                      #other options include "um", "Angs", "secs" (for phase curves)
    exo_dict['planet']['f_unit'] = 'rp^2/r*^2'               #other options are 'fp/f*'
    exo_dict['planet']['transit_duration'] = 2.0*60.0*60.0   #transit duration in seconds


    inst_dict = jdi.load_mode_dict('NIRSpec Prism')
    inst_dict["configuration"]["detector"]["subarray"] = 'sub512'
    inst_dict["configuration"]["detector"]["readmode"] = 'nrs'
    
    out = jdi.run_pandexo(exo_dict, ['MIRI LRS'], output_file="ref_trans_sim_c.p")
    
    
    m1,n1,x1,y1,e1 = jpi.jwst_1d_spec(out, R=100, num_tran=100, model=True, x_range=[1,5])
    
    plt.plot(m1,n1)
    plt.errorbar(x1,y1,e1)
    plt.show()


    """
    Works:
        NIRCam F444W
        NIRSpec G395M
        NIRSpec G395H (works but.. big noise bar?)
        NIRSpec G235H(works but.. big noise bar?)
        NIRSpec G140H
        MIRI LRS
        NIRISS SOSS (but takes a while...)
        
    
    Not Work
        NIRSpec Prism
        NIRSpec G235M
        NIRCam F322W
        NIRSpec G140M
    
    
    jwst ~ 350
    etc3D ~1170
    
    """
    
def read_output():
    
    
    out = pickle.load(open("ref_trans_sim2.p", 'rb'))
    
    print out["FinalSpectrum"]["wave"]
    print out["FinalSpectrum"]["spectrum"]
    
   
    m1 = out["OriginalInput"]["model_wave"]
    n1 = out["OriginalInput"]["model_spec"]
    
 
    #m1,n1,x1,y1,e1 = jpi.jwst_1d_spec(out, R=100, num_tran=100, model=True, x_range=[1,5])
    
    plt.plot(m1,n1)
    #plt.errorbar(x1,y1,e1)
    plt.show()


def test_instrument(inst):
    
    
    
    info = SetDefaultModes(inst)

    print info.pick()


if __name__ == "__main__":
    
    spectra = 'bio_trans1.txt'
    try_pandexo_tutorial(spectra)
    
    #test_instrument("NIRSpec G395M")
    #test_instrument("NIRSpec G235M")
    #read_output()
    
    
    """
    output = "ref_trans_sim2.p"
    spectra = 'ref_trans.txt'
    #run_pandexo(output,spectra)
    read_output()
    """
    """
    output = "bio_trans_sim.p"
    spectra = 'bio_trans.txt'
    run_pandexo(output,spectra)    
    """