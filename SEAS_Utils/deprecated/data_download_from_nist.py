# -*- coding: utf-8 -*-

"""

This module is used for acquiring data from NIST through web scraping


The previous version uses only the inchikey to dig through nist. 
Although inchikey are unique, it's not well documented and spectra might be missing
We motivate a multi-layer approach to scrap nist using all naming methods 
by piping the not found molecules through another round of scraping using another alias
until a molecule is found or the list of alias is exhausted. 

There are two approach. One is to scrape all the molecule using one list, another is to 
scrape each molecule using all its alias. For simplicity we choose the former method

webpage = urllib2.urlopen('http://en.wikipedia.org/wiki/Main_Page')
soup = BeautifulSoup(webpage)
for anchor in soup.find_all('a'):
    print(anchor.get('href', '/'))
    
V 2.0
"""
from bs4 import BeautifulSoup
import requests
import numpy as np
import time
import csv
import sys


def data_process(filename,spliter=""):
    """
    process data files into numpy arrays that can be plotted.
    input files are arbitrary columns.
    """
    start = time.time()
    ndata = []
    reader = csv.reader(open(filename, 'rU'))
    for line in reader:
        data = np.array(line,ndmin=2)  # convert data into numpy array 
        if ndata == []:     
            ndata = data
        else:
            ndata = np.concatenate((ndata,data))      # data stitching
    ndata = np.transpose(ndata)    

    print "Data Reading Time Elapse: ",time.time()-start
    return ndata

def list_of_names(filename):
    
    data = data_process(filename,",")
    names = data[0][1:]
    
    # getting the molecules
    new_name = []
    formula = []
    for i in names:
        name = i.split("\xca")
        if name == ['']:
            break
        new_name.append(name[0])
        formula.append(name[1])
    
    return new_name


def inchi_names(filename):
    data = []
    for line in csv.reader(open(filename,"r")):
        data.append(line[0])
    return data

def checks(data,keyword):
    links = []
    for link in data:
        if keyword in link.get_text():
            datafile = link.a.get('href')
            links.append("http://webbook.nist.gov%s"%datafile)
    return links


def function1(filename):
    respond1 = requests.get(filename)
    soup1 = BeautifulSoup(respond1.text)
    data1 = soup1.find_all('li')
    datalink1 = []
    for link1 in data1:
        if "Other data available" in link1.get_text():
            deeperdata = link1.ul.find_all("li")
            datalink1 = checks(deeperdata,"IR Spectrum")
    return datalink1

def function2(datalink1):
    datalink2 = []
    for i in datalink1:
        respond2 = requests.get(i)
        soup2 = BeautifulSoup(respond2.text)
        data2 = soup2.find_all('ul')
        for link2 in data2:
            try:
                datafile2 = link2.li.a.get('href')
                if "Type=IR-SPEC" in datafile2:
                    datalink2.append("http://webbook.nist.gov%s"%datafile2)
            except AttributeError:
                pass
    return datalink2
    
def function3_1(datalink1,tracker,notgas,T):
    """
    one spectra
    """
    respond3 = requests.get(datalink1[0])
    soup3 = BeautifulSoup(respond3.text)
    data3 = soup3.find_all('p')
    datalink3 = checks(data3,"Download")
    
    dataphase = soup3.find_all("table")
    temp = []
    
    trackersave,notgassave = [],[]
    for t in tracker:
        trackersave.append(t)
    for n in notgas:
        notgassave.append(n)
    
    for i in BeautifulSoup(str(dataphase)).find_all("tr"):
        if "State" in str(i) and "States" not in str(i):
            if ("gas" in str(i)) or ("GAS" in str(i)):
                print i
                temp.append(True)
            else:
                print i
                tracker.append(str(i))
                notgas.append(T)
                temp.append(False)
    
    if True in temp:
        return True,datalink3,trackersave,notgassave
    else:
        return False,datalink3,tracker,notgas

def function3_2(datalink2,tracker,notgas,T):
    """
    multiple spectra
    """
    # getting datafile from final data page
    datalink3 = []
    temp = []
    trackersave,notgassave = [],[]
    for t in tracker:
        trackersave.append(t)
    for n in notgas:
        notgassave.append(n)
        
    for i in datalink2:
        respond3 = requests.get(i)
        soup3 = BeautifulSoup(respond3.text)
        data3 = soup3.find_all('p')
        datalink3.append(checks(data3,"Download"))
        
        #print data3
        dataphase = soup3.find_all("table")
        for j in BeautifulSoup(str(dataphase)).find_all("tr"):
            if "State" in str(j) and "States" not in str(j):
                if ("gas" in str(j)) or ("GAS" in str(j)):
                    print j
                    temp.append(True)
                else:
                    print j
                    tracker.append(str(i))
                    notgas.append(T)
                    temp.append(False)
    
    if True in temp:
        return True,datalink3,trackersave,notgassave
    else:
        return False,datalink3,tracker,notgas

if __name__ =="__main__":
    
    #filename = "../Spectra_Data/Light_Gas.csv"
    #Target = list_of_names(filename)
    
    # all the non-gas
    #a = ['QMMFVYPAHWMCMS-UHFFFAOYSA-N', 'TUDWMIUPYRKEFN-UHFFFAOYSA-N', 'VOVUARRWDCVURC-UHFFFAOYSA-N', 'CYNYIHKIEHGYOZ-UHFFFAOYSA-N', 'FXXACINHVKSMDR-UHFFFAOYSA-N', 'GBBZLMLLFVFKJM-UHFFFAOYSA-N', 'HFEHLDPGIKPNKL-UHFFFAOYSA-N', 'ILAHWRKJUDSMFH-UHFFFAOYSA-N', 'JTWWWQGSFTWWDL-UHFFFAOYSA-N', 'KGBXLFKZBHKPEV-UHFFFAOYSA-N', 'MBNVSWHUJDDZRH-UHFFFAOYSA-N', 'RENMDAKOXSCIGH-UHFFFAOYSA-N', 'UMGDCJDMYOKAJW-UHFFFAOYSA-N', 'WEYOKDYZYYMRSQ-UHFFFAOYSA-N', 'WXEHBUMAEPOYKP-UHFFFAOYSA-N', 'WXRGABKACDFXMG-UHFFFAOYSA-N', 'XSQUKJJJFZCRTK-UHFFFAOYSA-N', 'ZDHXKXAHOVTTAH-UHFFFAOYSA-N', 'AEMRFAOFKBGASW-UHFFFAOYSA-N', 'ANGGPYSFTXVERY-UHFFFAOYSA-N', 'AUQDITHEDVOTCU-UHFFFAOYSA-N', 'BEBCJVAWIBVWNZ-UHFFFAOYSA-N', 'CUONGYYJJVDODC-UHFFFAOYSA-N', 'DOHQOGRRQASQAR-UHFFFAOYSA-L', 'FDNAPBUWERUEDA-UHFFFAOYSA-N', 'FOCAUTSVDIKZOP-UHFFFAOYSA-N', 'GFZMFCVDDFHSJK-UHFFFAOYSA-N', 'GYXHHICIFZSKKZ-UHFFFAOYSA-N', 'HBNYJWAFDZLWRS-UHFFFAOYSA-N', 'HCSDJECSMANTCX-UHFFFAOYSA-N', 'HEDKQVNHJZBFQR-UHFFFAOYSA-N', 'IHZAEIHJPNTART-UHFFFAOYSA-N', 'IJOOHPMOJXWVHK-UHFFFAOYSA-N', 'IWYRWIUNAVNFPE-UHFFFAOYSA-N', 'JLUFWMXJHAVVNN-UHFFFAOYSA-N', 'JMYVMOUINOAAPA-UHFFFAOYSA-N', 'JVQIKJMSUIMUDI-UHFFFAOYSA-N', 'OIFAHDAXIUURLN-UHFFFAOYSA-N', 'OYUNTGBISCIYPW-UHFFFAOYSA-N', 'PAVZHTXVORCEHP-UHFFFAOYSA-N', 'PXAJQJMDEXJWFB-UHFFFAOYSA-N', 'PZFBULOUMNPBFA-UHFFFAOYSA-N', 'QLNJFJADRCOGBJ-UHFFFAOYSA-N', 'RAXXELZNTBOGNW-UHFFFAOYSA-N', 'STZZWJCGRKXEFF-UHFFFAOYSA-N', 'SYZRZLUNWVNNNV-UHFFFAOYSA-N', 'TZGPACAKMCUCKX-UHFFFAOYSA-N', 'UCXUKTLCVSGCNR-UHFFFAOYSA-N', 'VAWRKITUPUFMHV-UHFFFAOYSA-N', 'VEZXCJBBBCKRPI-UHFFFAOYSA-N', 'VXIVSQZSERGHQP-UHFFFAOYSA-N', 'WEGOLYBUWCMMMY-UHFFFAOYSA-N', 'WFCLYEAZTHWNEH-UHFFFAOYSA-N', 'WQYSXVGEZYESBR-UHFFFAOYSA-N', 'WTKZEGDFNFYCGP-UHFFFAOYSA-N', 'XGEGHDBEHXKFPX-UHFFFAOYSA-N', 'YGKHJWTVMIMEPQ-UHFFFAOYSA-N', 'YOWQWFMSQCOSBA-UHFFFAOYSA-N', 'BTVWZWFKMIUSGS-UHFFFAOYSA-N', 'BUZZUHJODKQYTF-UHFFFAOYSA-N', 'CEBKHWWANWSNTI-UHFFFAOYSA-N', 'CQFQAARMEJVWAL-UHFFFAOYSA-N', 'DGJMPUGMZIKDRO-UHFFFAOYSA-N', 'DLDJFQGPPSQZKI-UHFFFAOYSA-N', 'DNSISZSEWVHGLH-UHFFFAOYSA-N', 'FIYBYNHDEOSJPL-UHFFFAOYSA-N', 'FQPSGWSUVKBHSU-UHFFFAOYSA-N', 'FTAHXMZRJCZXDL-UHFFFAOYSA-N', 'FYOWZTWVYZOZSI-UHFFFAOYSA-N', 'GCSJLQSCSDMKTP-UHFFFAOYSA-N', 'GSCLMSFRWBPUSK-UHFFFAOYSA-N', 'HHENWUYAOCBSAE-UHFFFAOYSA-N', 'HIZVCIIORGCREW-UHFFFAOYSA-N', 'HSJKGGMUJITCBW-UHFFFAOYSA-N', 'HTWIZMNMTWYQRN-UHFFFAOYSA-N', 'IAHFWCOBPZCAEA-UHFFFAOYSA-N', 'JFEVIPGMXQNRRF-UHFFFAOYSA-N', 'JLSJEUQOXVVCPN-UHFFFAOYSA-N', 'JOHLTMWXHJLNDE-UHFFFAOYSA-N', 'JOYRKODLDBILNP-UHFFFAOYSA-N', 'KMTRUDSVKNLOMY-UHFFFAOYSA-N', 'KPSSIOMAKSHJJG-UHFFFAOYSA-N', 'KYQCOXFCLRTKLS-UHFFFAOYSA-N', 'LCTONWCANYUPML-UHFFFAOYSA-N', 'LHYKTQVFLKHQSR-UHFFFAOYSA-N', 'MCTWTZJPVLRJOU-UHFFFAOYSA-N', 'NXMXETCTWNXSFG-UHFFFAOYSA-N', 'OIBMEBLCOQCFIT-UHFFFAOYSA-N', 'OPCJOXGBLDJWRM-UHFFFAOYSA-N', 'ORMSTDJYMPIZAO-UHFFFAOYSA-N', 'PBMFSQRYOILNGV-UHFFFAOYSA-N', 'PDQAZBWRQCGBEV-UHFFFAOYSA-N', 'PNECSTWRDNQOLT-UHFFFAOYSA-N', 'QUKGLNCXGVWCJX-UHFFFAOYSA-N', 'RFFFKMOABOFIDF-UHFFFAOYSA-N', 'RXKJFZQQPQGTFL-UHFFFAOYSA-N', 'SHEINYPABNPRPM-UHFFFAOYSA-N', 'SKDFWEPBABSFMG-UHFFFAOYSA-N', 'LCTONWCANYUPML-UHFFFAOYSA-N', 'VAYTZRYEBVHVLE-UHFFFAOYSA-N', 'VDBCTDQYMZSQFQ-UHFFFAOYSA-N', 'VJHTZTZXOKVQRN-UHFFFAOYSA-N', 'VLERMNIUDRUNQO-UHFFFAOYSA-N', 'VLIDBBNDBSNADN-UHFFFAOYSA-N', 'VUIWJRYTWUGOOF-UHFFFAOYSA-N', 'WASQWSOJHCZDFK-UHFFFAOYSA-N', 'WCXXISMIJBRDQK-UHFFFAOYSA-N', 'WDXYVJKNSMILOQ-UHFFFAOYSA-N', 'WFKAJVHLWXSISD-UHFFFAOYSA-N', 'WGJCBBASTRWVJL-UHFFFAOYSA-N', 'XSDCTSITJJJDPY-UHFFFAOYSA-N', 'XXZVBSPSFVERCB-UHFFFAOYSA-N', 'YLJJAVFOBDSYAN-UHFFFAOYSA-N', 'YVWPNDBYAAEZBF-UHFFFAOYSA-N', 'ZQDPJFUHLCOCRG-AATRIKPKSA-N', 'ZXCYIJGIGSDJQQ-UHFFFAOYSA-N']
    #b = ['<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C557686&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1330, 10% CS2 FOR 1330-450 CM<sup>-1</sup>)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C624737&Units=SI&Type=IR-SPEC&Index=2#IR-SPEC', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C556569&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1330, 10% CS2 FOR 1330-450 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (NUJOL MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 2.6-7.5, 10% CS2 FOR 7.5-26)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (MINERAL OIL MULL)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C18283937&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">VAPOR AT 30 mmHg PRESSURE; $$ RESEARCH PURITY</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL, FLUOROLUBE FOR 3800-1333 AND NUJOL FOR 1333-450 CM <sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (5% CCl4 FOR 3800-1300, 4% CS2 FOR 1300-650, AND 5% CCl4 FOR 650-250)<br/>\nMORE THAN 99.5% PURE</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C558178&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (0.4 mg / 650 mg KBr DISC)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C109773&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL), FLUOROLUBE FOR 3800-1330 CM<sup>-1</sup>, NUJOL FOR 1330-460 CM<sup>-1</sup><br/>\nSPECTRAL FEATURE AT 1336 CM<sup>-1</sup> IS AN ARTIFACT CAUSED BY A CHANGE IN SOLVENT<br/>\nSPECTRAL FEATURE AT 2008 CM<sup>-1</sup> IS AN ARTIFACT CAUSED BY GRATING CHANGE</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (2% IN CCl4 FOR 700-500 CM<sup>-1</sup>)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C79118&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C109820&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (MINERAL OIL MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1340, 10% CS2 FOR 1340-400 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 5000-1330, 10% CS2 FOR 1330-625 CM<sup>-1</sup>)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C353548&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C75774&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C765344&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (9% CCl4 FOR 3800-1300, 10% CS2 FOR 1300-650, AND 9% CCl4 FOR 650-250)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% IN CCl4 FOR 3800-1333, 10% IN CS2 FOR 1333-450, AND 10% IN C6H12 FOR 450-200 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1330, 10% CS2 FOR 1330-440 CM<sup>-1</sup>)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C920376&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID; $$ RESEARCH PURITY</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1340, 10% CS2 FOR 1340-400) VS SOLVENT</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C21020246&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)<br/>\nSPECTRAL FEATURE AT 1993 CM<sup>-1</sup> IS PROBABLY AN ARTIFACT CAUSED BY GRATING CHANGE</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C3018120&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1340, 10% CS2 FOR 1340-400 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (11% CCl4 FOR 3800-1300, 2% CS2 FOR 1300-650, AND 11% CCl4 FOR 650-250)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">VAPOR AT 24.5 mmHg PRESSURE; $$ PURITY 95%</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 2-7.6, SATURATED-LESS THAN 10% -NaCl ADDED CS2 FOR 7.4-16)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 2.6-7.5, 10% CS2 FOR 7.5-25 MICRON)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1330, 10% CS2 FOR 1330-450 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (5% CCl4 FOR 3800-1335, 5% CS2 FOR 1335-400)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (OIL MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1300, 10% CS2 FOR 1300-400 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)<br/>\nSPECTRAL FEATURE AT 1992 CM<sup>-1</sup> IS AN ARTIFACT CAUSED BY GRATING CHANGE</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% IN CCl4 FOR 3800-1333, AND 10% IN CS2 FOR 1333-440 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL, FLUOROLUBE FOR 3800-1300 AND NUJOL FOR 1300-250 CM <sup>-1</sup>)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C694053&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL), FLUOROLUBE FOR 3800-1330 CM<sup>-1</sup>, NUJOL FOR 1330-400 CM<sup>-1</sup><br/>\nSPECTRAL CONTAMINATION DUE TO AN UNKNOWN AROUND 400 CM<sup>-1</sup></td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C754052&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 2.6-7.5, 10% CS2 FOR 7.5-22 MICRON)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (CCl4 FOR 2.0-7.5, CS2 FOR 7.5-16.0 MICRON) VS SOLVENT</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C2567148&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID(MINERAL OIL MULL)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (0.7 mg / 500 mg KBr DISC)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C51796&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C96491&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1340, 10% CS2 FOR 1340-450 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 2.6-7.5, 10% CS2 FOR 7.5-24 MICRON)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C127173&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (0.7 mg / 500 mg KBr DISC)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (2% CCl4 FOR 5000-1330, 2% CS2 FOR 1330-625 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 2.5-7.5, 10% CS2 FOR 7.5-16 MICRONS)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL, FLUOROLUBE FOR 3800-1333 AND NUJOL FOR 1333-400 CM <sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% IN CCl4 FOR 3800-1300, 2% IN CS2 FOR 1300-600, AND 10% IN CCl4 FOR 600-250 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (KBr PELLET)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (1% KBr DISC)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% IN CCl4 FOR 3800-1340 AND 10% IN CS2 FOR 1340-450 CM<sup>-1</sup>)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C127173&Units=SI&Type=IR-SPEC&Index=0#IR-SPEC', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID (NEAT)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (1 mg / 650 mg KBr DISC)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (KBr PELLET)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (1 mg /650 mg KBr DISC)<br/>\nPURITY - ANALYTICAL<br/>\nSPECTRAL CONTAMINATION DUE TO H2O AROUND 3500 CM<sup>-1</sup></td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 2.5-7.5, 10% CS2 FOR 7.5-16 MICRONS)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION UNDER NITROGEN (10% CCl4 FOR 5000-1330, 10% CS2 FOR 1330-625 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (CS2 FOR 3800-400 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">LIQUID</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (SPLIT MULL)<br/>\nSPECTRAL FEATURE AT 1342 CM<sup>-1</sup> IS AFFECTED BY A CHANGE IN GRATING</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (KBr PELLET)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLID (BETWEEN SALTS)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1330, 10% CS2 FOR 1330-400 CM<sup>-1</sup>)</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1300, 2% CS2 FOR 1300-620, AND 10% CCl4 FOR 620-240) VS SOLVENT</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-900, 2% CS2 FOR 900-640, AND 10% CCl4 FOR 650-250)<br/>\nMORE THAN 99% PURE</td></tr>', '<tr><th align="left" valign="top">State</th><td align="left" valign="top">SOLUTION (10% CCl4 FOR 3800-1330, 10% CS2 FOR 1330-460 CM<sup>-1</sup>)</td></tr>', 'http://webbook.nist.gov/cgi/cbook.cgi?ID=C616239&Units=SI&Type=IR-SPEC&Index=1#IR-SPEC']


    
    
    filename = "Inchikey.txt"
    Target = inchi_names(filename)
    
    good = 0
    start = time.time()
    one, three = 0,0
    count = 0
    gas = 0
    notgas = []
    tracker = []
    has_spectra = []
    for T in Target:
        #filename = "http://webbook.nist.gov/cgi/cbook.cgi?Name=%s&Units=SI"%T
        count+=1
        print count,gas,good,one,three
        
        filename = "http://webbook.nist.gov/cgi/cbook.cgi?InChI=%s&Units=SI"%T
        
        datalink1 = function1(filename)
        if datalink1 == []:
            print T,":\nNo Spectra Available 1"
            one +=1
            continue
        
        datalink2 = function2(datalink1)
        
        if len(datalink2) == 0:
            gases,datalink3,tracker,notgas = function3_1(datalink1,tracker,notgas,T)
        else:
            gases,datalink3,tracker,notgas = function3_2(datalink2,tracker,notgas,T)
        if gases == True:
            gas+=1
            gases = False

        if datalink3 == []:
            print T,"\nNo Spectra Available 3"
            three +=1
            continue
        
        print T
        for i in datalink3:
            if i == []:
                continue
            else:
                print str(i)      
        good+=1
        if gases:
            has_spectra.append(T)
        
        
    
    total = len(Target)
    print "%s/%s"%(str(good),str(total))
    print "Gas: ",gas,"/",total
    print "Total Time: ",time.time()-start
    print notgas
    print tracker
    
    print len(has_spectra)
    sys.stdout = open("nist_has_spectra_gasonly.txt","w")
    for i in has_spectra:
        print i