"""

data processing and playing with pendexo results 


Horizontal Bin: 100ppm

data should be a combination between NIRSpec and MIRI
detection could indicate if the data is from NIRSpec of MIRI.


3 sigma detection: signal1-signal2 > 3*(signal1_error+signal2_error)


"""

"""
['FinalSpectrum', 'OriginalInput', 'timing_div', 'warnings_div', 
'PandeiaOutTrans', 'input_div', 'warning', 'RawData', 'timing', 'input']

FinalSpectrum 
['spectrum_w_rand', 'spectrum', 'error_w_floor', 'wave']

OriginalInput 
['model_spec', 'model_wave']

timing_div warnings_div PandeiaOutTrans 
['sub_reports', 'information', 'warnings', 'transform', '2d', 'scalar', '1d', 'input']

input_div warning 
['Num Groups Reset?', 'Group Number Too Low?', 'Group Number Too High?', 'Saturated?', 'Non linear?', '% full well high?']

RawData 
['var_in', 'rn[out,in]', 'electrons_out', 'e_rate_in', 'wave', 'electrons_in', 'error_no_floor', 'bkg[out,in]', 'e_rate_out', 'var_out']

timing 
['Seconds per Frame', 'Number of Transits', 'Time/Integration incl reset (sec)', 'Observing Efficiency (%)', 'Num Integrations Out of Transit', 'Num Integrations In Transit', 'APT: Num Integrations per Occultation', 'APT: Num Groups per Integration', 'Transit+Baseline, no overhead (hrs)', 'Transit Duration']

input 
['Target Mag', 'Readmode', 'Instrument', 'Disperser', 'Filter', 'Calculation Type', 'Mode', 'Saturation Level (electons)', 'Aperture', 'Subarray', 'Primary/Secondary']
"""



import pickle
import pprint
import numpy as np
import matplotlib.pyplot as plt

import pandexo.engine.justdoit as jdi
import pandexo.engine.justplotit as jpi

def load_pandexo(instrument,filename):


    if type(instrument) == str:

        pkl_file = open(filename, 'rb')
        
        data = pickle.load(pkl_file)
        
        wave_ori = data["OriginalInput"]["model_wave"]
        spec_ori = data["OriginalInput"]["model_spec"]#*10**6
        
        wave = data["FinalSpectrum"]["wave"]
        spec = data["FinalSpectrum"]["spectrum"]
        error = data["FinalSpectrum"]["error_w_floor"]
    
        return wave_ori,spec_ori,wave,spec,error

    if len(instrument) == 1:
    
        pkl_file = open(filename, 'rb')
        
        data = pickle.load(pkl_file)
        
        print data[0]['NIRSpec G140M']["FinalSpectrum"].keys()
        
        return
        
        wave_ori = data[0][instrument[0]]["OriginalInput"]["model_wave"]
        spec_ori = data[0][instrument[0]]["OriginalInput"]["model_spec"]#*10**6
        
        wave = data[0][instrument[0]]["FinalSpectrum"]["wave"]
        spec = data[0][instrument[0]]["FinalSpectrum"]["spectrum"]
        error = data[0][instrument[0]]["FinalSpectrum"]["error_w_floor"]
    
        return wave_ori,spec_ori,wave,spec,error

    else:
        
        pkl_file1 = open(filename, 'rb')
        
        data = pickle.load(pkl_file1)
        
        
        wave_ori = []
        spec_ori = []
        wave = []
        spec = []
        error = []
        
        for i,inst in enumerate(instrument):
        
            wave_ori = np.concatenate([wave_ori,data[i][inst]["OriginalInput"]["model_wave"]])
            spec_ori = np.concatenate([wave_ori,data[i][inst]["OriginalInput"]["model_spec"]])
            
            wave = np.concatenate([wave_ori,data[i][inst]["FinalSpectrum"]["wave"]])
            spec = np.concatenate([wave_ori,data[i][inst]["FinalSpectrum"]["spectrum"]])
            error = np.concatenate([wave_ori,data[i][inst]["FinalSpectrum"]["error_w_floor"]])
        
        return wave_ori,spec_ori,wave,spec,error

def play_with_pandexo_output():
    
    file1 = "ETC.p"
    file2 = "ETC-2.p"
        
    wave_ori1,spec_ori1,wave1,spec1,error1 = load_pandexo(file1)
    wave_ori2,spec_ori2,wave2,spec2,error2 = load_pandexo(file2)
    
    
    
    
    plt.plot(wave_ori1,spec_ori1)
    plt.plot(wave_ori2,spec_ori2)
    #plt.errorbar(wave1,spec1,yerr=error1)
    
    
    plt.show()

def display_pandexo_output():

    file1 = "Test_multi_instrument.p"
    file1 = "Test_NIRSpec_instrument2.p"
    
    instrument = ""#['NIRSpec Prism']#,'MIRI LRS']
    #instrument = ['NIRSpec G140M']#, 'NIRSpec G235M','NIRSpec G395M']
    #run.run_pandexo(file1)
    
    pkl_file = open(file1, 'rb')
    out = pickle.load(pkl_file)
    

    """
    num_tran = 100
    R = 100

    ntran_old = out['timing']['Number of Transits']
    to = out['timing']["Num Integrations Out of Transit"]
    ti = out['timing']["Num Integrations In Transit"]
    #remove any nans 
    y = out['FinalSpectrum']['spectrum_w_rand']
    x = out['FinalSpectrum']['wave'][~np.isnan(y)]
    err = out['FinalSpectrum']['error_w_floor'][~np.isnan(y)]
    y = y[~np.isnan(y)]

    new_wave = jpi.bin_wave_to_R(x, R)


    print out['RawData']['electrons_in']
    
    out = jpi.uniform_tophat_sum(new_wave,x, out['RawData']['electrons_out']*num_tran/ntran_old)
    inn = jpi.uniform_tophat_sum(new_wave,x, out['RawData']['electrons_in']*num_tran/ntran_old)
    return
    
    vout = jpi.uniform_tophat_sum(new_wave,x, out['RawData']['var_out']*num_tran/ntran_old)
    vin = jpi.uniform_tophat_sum(new_wave,x, out['RawData']['var_in']*num_tran/ntran_old)
    var_tot = (to/ti/out)**2.0 * vin + (inn*to/ti/out**2.0)**2.0 * vout
    if dict['input']['Primary/Secondary']=='fp/f*':
        fac = -1.0
    else:
        fac = 1.0
    rand_noise = np.sqrt((var_tot))*(np.random.randn(len(new_wave)))
    raw_spec = (out/to-inn/ti)/(out/to)       
    sim_spec = fac*(raw_spec + rand_noise )
    x = new_wave
    y = sim_spec
    err = np.sqrt(var_tot)
    """
    x,y, e = jpi.jwst_1d_spec(out, R=50, num_tran=100, model=True, x_range=[.8,11.28])
    
    #wave_ori1,spec_ori1,wave1,spec1,error1 = load_pandexo(instrument, file1)
    #plt.plot(wave_ori1,spec_ori1)
    #plt.errorbar(wave1,spec1,yerr=error1)
    
    
    plt.show()

def display_compare():

    file1 = "ref_trans_sim.p"
    file2 = "bio_trans_sim_c.p"
    
    Bins = 100
    Tran = 100
    

    out1 = pickle.load(open(file1, 'rb'))
    out2 = pickle.load(open(file2, 'rb'))
    
    m1,n1,x1,y1,e1 = jpi.jwst_1d_spec(out1, R=Bins, num_tran=Tran, model=True, x_range=[1,13])
    m2,n2,x2,y2,e2 = jpi.jwst_1d_spec(out2, R=Bins, num_tran=Tran, model=True, x_range=[1,13])
    
    
    plt.plot(m1,n1, label="base_ref")
    plt.errorbar(x1,y1,e1,fmt="o",markersize='5', label="base_obs")
    
    plt.plot(m2,n2, label="bio_ref")
    plt.errorbar(x2,y2,e2,fmt="*",markersize='5', label="bio_obs")
    
    plt.legend()
    plt.title("Simulated Exoplanet Atmosphere Observed Spectra using Pandexo ETC \nwith %s bins and %s transits"%(Bins,Tran))
    plt.xlabel('Wavelength [microns]')
    plt.ylabel('fp/f*')
    
    plt.show()
    """
    """
    #detection = y2-y1# +3*(e1+e2)
    x1 = np.array(x1[0])
    
    y1 = np.array(y1[0])
    y2 = np.array(y2[0])
    e1 = np.array(e1[0])
    e2 = np.array(e2[0])
    
    detection = y2-y1# +3*(e1+e2)
    error = 3*(e1+e2)
    
    """
    plt.errorbar(x1,detection,error)

    plt.title("Simulated Exoplanet Atmosphere Detection with %s bins and %s transits"%(Bins,Tran))
    plt.xlabel('Wavelength [microns]')
    plt.ylabel('fp/f*')    
    plt.show()
    """
    
    detected_x = []
    detected_y1 = []
    detected_y2 = []
    for i,bin in enumerate(detection-error):
        if bin > 0:
            print x1[i],"detected"
            
            detected_x.append(x1[i])
            detected_y1.append(y1[i])
            detected_y2.append(y2[i])


    plt.title("Simulated Exoplanet Atmosphere Detection with %s bins and %s transits at 3 Sigma"%(Bins,Tran))
    plt.xlabel('Wavelength [microns]')
    plt.ylabel('fp/f*') 

    plt.plot(m1,n1)
    plt.plot(m2,n2)   
    plt.plot(detected_x,detected_y1,".")
    plt.plot(detected_x,detected_y2,".")

    
    
    plt.show()
    # next, find local maximum...
    
    
if __name__ == "__main__":
    
    display_compare()
    
    
    
    