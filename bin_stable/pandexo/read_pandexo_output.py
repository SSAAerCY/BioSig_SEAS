import cPickle as pickle
import matplotlib.pyplot as plt

filename = "ETC.p"


data = pickle.load( open( filename, "rb" ) )

keys = data.keys()

for i in keys:
    print i
    try:
        print data[i].keys()
    except:
        pass


w_rand = data["FinalSpectrum"]["spectrum_w_rand"]
spectrum = data["FinalSpectrum"]["spectrum"]
error = data["FinalSpectrum"]["error_w_floor"]
wave = data["FinalSpectrum"]["wave"]



wave1     = data["OriginalInput"]["model_wave"]
                
spectrum1 = data["OriginalInput"]["model_spec"]

plt.plot(wave1,spectrum1)
plt.errorbar(wave,spectrum, yerr=error)
#plt.plot(wave,spectrum)

plt.show()

