"""

Testing ground for visualizing the excel data from result

"""

import os
import sys
import numpy as np


from openpyxl import load_workbook, Workbook

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils as utils
from SEAS_Utils.common_utils.data_loader import NIST_Smile_List
from SEAS_Utils.common_utils.DIRs import *




def visualizer():
    
    
    name = "Second_Stage_NIST_Detectable_Abundance_Earth_Grey_Cloud.xlsx"
    path = os.path.join(Result_DIR,"Second_Stage_Simulation")
    
    filename = os.path.join(path,name)
    
    sheet = "General Result"

    WBI = load_workbook(filename) 

    inputsheetdata = WBI["EA_10000_0.2"]
    
    
    for i in range(534):
        num = i+2
        
        
        large = 0
        medium = 0
        small = 0
        
        for j in range(8):
            col = chr(ord("G")+j)
            value = inputsheetdata["%s%s"%(col,num)].value
            
            try:
                value = float(value)
                if value <= 1:
                    large +=1
                elif value <= 100:
                    medium +=1
                elif value <= 10000:
                    small +=1
            except:
                pass
        
            
        preference1 = large*5+medium*2+small   
        inputsheetdata["%s%s"%("O",num)].value = preference1
        
        preference2 = large*2+medium*1.5+small  
        inputsheetdata["%s%s"%("P",num)].value = preference2

    WBI.save(filename)
    

if __name__ == "__main__":
    
    visualizer()
    
    
    
    
    