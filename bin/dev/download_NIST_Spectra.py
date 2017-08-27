"""
Download NIST Spectra from NIST Chemistry Workbook

Column in ID 

('SMILES', 'text'), 
('Compatible_SMILES', 'text'), 
('Database_Number', 'text'), 
('Molecular_Weight', 'text'), 
('Formula', 'text'), 
('IUPAC_chemical_name', 'text'), 
('Number_of_Atoms', 'text'), 
('Non_H_atoms', 'text'), 
('InChI', 'text'), 
('InChIKey', 'text'), 
('Boiling_Point', 'text'), 
('Produced_by_Life', 'text')
"""


import os
import sys
import urllib2
import httplib
import mechanize
import numpy as np
from bs4 import BeautifulSoup
from openpyxl import Workbook, load_workbook

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils as utils
import SEAS_Utils.common_utils.db_management2 as dbm
from SEAS_Utils.common_utils.DIRs import molecule_info, NIST_Spectra
from SEAS_Utils.common_utils.data_downloader import link_loader
from SEAS_Utils.common_utils.timer import simple_timer
from SEAS_Utils.common_utils.data_saver import save_txt

def get_spectra(link):
    
    soup = BeautifulSoup(link_loader(link).read(), "html.parser")
    
    
    if "Spectrum not found" in str(soup):
        return "N",""
    else:
        return "Y",str(soup)
    

def download_NIST():

    kwargs = {  "dir"        :molecule_info,
                "db_name"    :"Molecule_DB2.db",
                "user"       :"Azariven",
                "DEBUG"      :False,
                "REMOVE"     :True,
                "BACKUP"     :False,
                "OVERWRITE"  :True}

    cross_db = dbm.database(**kwargs)
    cross_db.access_db()
    table_name = "ID"

    #cmd = "SELECT InChIKey from %s WHERE Formula='C2H2'"%table_name
    cmd = "SELECT InChIKey,Formula from %s"%table_name
    
    cross_db.c.execute(cmd)
    data = np.array(cross_db.c.fetchall()).T
    
    inchikey_list = data[0]
    formula_list = data[1]
    
    print inchikey_list
    print formula_list
    

    
    Timer = simple_timer(4)
    counter = 0
    total = len(inchikey_list)
    print total
    sys.exit()
    for inchikey in inchikey_list:
    
        if counter%100 == 0:
            print Timer.progress(counter,total)
    
        link = "http://webbook.nist.gov/cgi/cbook.cgi?InChI=%s&Units=SI"%inchikey


    #data_link = "http://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C74862&amp;Index=0&amp;Type=IR"

    
        #print inchikey, get_cas(link)

        counter+=1
    
    
def download_NIST2():


    inputfilename = os.path.join(molecule_info,"ASM4.3_Lite.xlsx")
    inputsheetname = "Database"

    WBI = load_workbook(inputfilename)
    data_sheet_ranges = WBI[inputsheetname]
    
    Data = []
    Timer = simple_timer(4)
    num = 2
    total = 16420
    while True:

        Smile            = utils.to_str(data_sheet_ranges['A%d'%num].value)
        if Smile == None:
            break
        
        cas      = utils.to_str(data_sheet_ranges['M%d'%num].value)
        if cas == None:
            num+=1
            continue
        
        num+=1
        print num,
        
        
        for i in range(4):
            Index = i
            
            data_link = "http://webbook.nist.gov/cgi/cbook.cgi?JCAMP=%s&amp;Index=%s&amp;Type=IR"%(cas,Index)
            
            result, spectra_data = get_spectra(data_link)
            
            if Index == 0:
                data_sheet_ranges['N%d'%num].value = result
            
            if result == "N":
                print ""
                break
            else:
                name = "%s_%s.jdx"%(cas,Index)
                path = os.path.join(NIST_Spectra,cas)
                print "Saving %s"%name
                save_txt(path, name, spectra_data, extension=".jdx", check=True)

        
    WBI.save(inputfilename) 
    
    


if __name__ == "__main__":
    
    
    download_NIST2()
    
    
    
    
    