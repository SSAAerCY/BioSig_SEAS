"""

Production code to get cas registry number for each molecule in the asm list

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
from SEAS_Utils.common_utils.DIRs import molecule_info
from SEAS_Utils.common_utils.data_downloader import link_loader
from SEAS_Utils.common_utils.timer import simple_timer


    
def get_cas(link):
    
    
    soup = BeautifulSoup(link_loader(link).read(), "html.parser")
    
    try:
        for i in soup.findAll("a"):
            if "Download the identifier in a file" in i.text:
                cas = i["href"].split("=")[-1]

                return cas
        return None
    except:
        return None   
    
def main():    


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
        
        if num%100 == 0:
            print Timer.progress(num,total)
            WBI.save(inputfilename)
        
        inchikey = utils.to_str(data_sheet_ranges['J%d'%num].value).split("=")[-1]

        link = "http://webbook.nist.gov/cgi/cbook.cgi?InChI=%s&Units=SI"%inchikey
        
        result = get_cas(link)
        
        data_sheet_ranges['M%d'%num].value = str(result)
        
        num+=1
        
    WBI.save(inputfilename)

if __name__ == "__main__":
    

    main()
    