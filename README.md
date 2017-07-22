# BioSig_SEAS


This is the Simulator for Exoplanet Atmosphere Spectra for biosignature studies

Code Download:

git clone https://github.com/azariven/BioSig_SEAS.githttps://github.com/azariven/BioSig_SEAS.git


All tested exectution files are under /bin_stable
To run the simulation, please check /bin_stable/Main/

For code development, please put testing execution files under /bin
All code and modules should follow be two directory deep by default.

Recent Work adding CIA. To Run, 


# Calculate the CIA Spectra for H2-H2
Data download, you only need to do this once

cd BioSig_SEAS/bin_stable/data_download
python download_HITRAN_CIA.py

cd BioSig_SEAS/bin_stable/spectra
python display_CIA_spectra.py
