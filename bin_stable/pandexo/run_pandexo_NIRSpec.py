
import warnings
import pprint
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi # THIS IS THE HOLY GRAIL OF PANDEXO
import numpy as np

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
    jdi.run_pandexo(exo_dict, ['MIRI LRS'], output_file=output)
    print('SUCCESS') 


 
if __name__ == "__main__":
    output = "ref_trans_sim.p"
    spectra = 'ref_trans.txt'
    run_pandexo(output,spectra)
    
    output = "bio_trans_sim.p"
    spectra = 'bio_trans.txt'
    run_pandexo(output,spectra)    
    