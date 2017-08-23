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
from bs4 import BeautifulSoup
from openpyxl import Workbook, load_workbook

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))



import SEAS_Utils as utils
import SEAS_Utils.common_utils.db_management2 as dbm
from SEAS_Utils.common_utils.DIRs import molecule_info
from SEAS_Utils.common_utils.data_downloader import link_loader



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

    cmd = "SELECT InChIKey from %s WHERE Formula='C2H2'"%table_name

    cross_db.c.execute(cmd)
    inchikey = cross_db.c.fetchall()[0][0]
    
    link = "http://webbook.nist.gov/cgi/cbook.cgi?InChI=%s&Units=SI"%inchikey

    print link


if __name__ == "__main__":
    
    
    download_NIST()
    
    
    
    
    