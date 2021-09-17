import sys, os, re
from math import sqrt
from statistics import mean, stdev
from array import array
from ROOT import *
import ctypes

if len(sys.argv) != 3 :
    print("execute as: \npython time_characterization.py path/to/folder/containing/CalTemp_measN_amulet_files path/to/oscilloscope/file.txt")
    raise Exception("you have to pass the correct arguments")

def GetDeltas(t):
    return [j-i for i, j in zip(t[:-1], t[1:])]

def GetFiles(path):
    files = []
    folder = os.fsencode(path)
    for file in os.listdir(folder):
        filename = os.fsdecode(file)
        if "CalTemp" in filename:
            files.append(path+filename)
    return files

def GetDigitizerPeriods(files):
    periods_up, periods_dwn, periods_upErr, periods_dwnErr = {}, {}, {}, {}
    for file in files:
        meas_num = int(re.findall(r'\d+', file)[0])
        f = TFile.Open(file)
        tree = f.Get("amulet")
        period_dwn, period_up = [], []
        for entry in range(0,tree.GetEntries()):
            tree.GetEntry(entry)
            time_ups = getattr(tree,"decoded.ch0timeups")
            time_dwns = getattr(tree,"decoded.ch0timedwns")
            period_dwn += GetDeltas(time_ups)
            period_up += GetDeltas(time_dwns)
        periods_up[meas_num] = mean(period_up)
        periods_dwn[meas_num] = mean(period_dwn)
        periods_upErr[meas_num] = stdev(period_up)/sqrt(len(period_up))
        periods_dwnErr[meas_num] = stdev(period_dwn)/sqrt(len(period_dwn))
        f.Close()
    return {"up":periods_up,"dwn":periods_dwn,"upErr":periods_upErr, "dwnErr":periods_dwnErr}

def GetOscilloscopePeriod(filename):
    periods = {}
    with open(filename) as f:
        for line in f.readlines():
            data = line.split()
            periods[int(data[0])] = (float(data[1])*1e-9,float(data[2])*1e-9)
    return periods

if __name__ == "__main__":
    
    digi_periods = GetDigitizerPeriods(GetFiles(sys.argv[1]))
    osc_periods = GetOscilloscopePeriod(sys.argv[2])

    digiups = digi_periods["up"]
    digidwns = digi_periods["dwn"]
    digiupserr = digi_periods["upErr"]
    digidwnserr = digi_periods["dwnErr"]

    osc, oscerr, digi, digierr, digiup,  digidwn, digiuperr, digidwnerr, diff, differr = array('d'), array('d'), array('d'), array('d'),array('d'), array('d'),array('d'), array('d'),array('d'), array('d')

    for i in range(1,len(digiups)+1):
        osc.append(osc_periods[i][0])
        oscerr.append(osc_periods[i][1])
        digiup.append(digiups[i])
        digidwn.append(digidwns[i])
        digiuperr.append(digiupserr[i])
        digidwnerr.append(digidwnserr[i])
        digi.append(mean([digiups[i],digidwns[i]]))
        digierr.append(sqrt(digiupserr[i]**2+digidwnserr[i]**2))
        diff.append(osc[i-1]-digi[i-1])
        differr.append(sqrt(oscerr[i-1]**2+digierr[i-1]**2))
    gStyle.SetOptFit(1111)

    c = TCanvas()
    g = TGraphErrors(len(osc), osc, digiup, oscerr, digiuperr)
    axisTitles=";oscilloscope T[s];digitizer T[s]"
    gtitle = "digitizer_characterization_up"
    g.SetNameTitle(gtitle,gtitle+axisTitles)
    g.SetMarkerSize(1)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(kBlue)
    g.Fit("pol1")
    g.Draw()

    c1 = TCanvas()
    g1 = TGraphErrors(len(osc), osc, digidwn, oscerr, digidwnerr)
    axisTitles=";oscilloscope T[s];digitizer T[s]"
    gtitle = "digitizer_characterization_dwn"
    g1.SetNameTitle(gtitle,gtitle+axisTitles)
    g1.SetMarkerSize(1)
    g1.SetMarkerStyle(20)
    g1.SetMarkerColor(kBlue)
    g1.Fit("pol1")
    g1.Draw()
    
    digi.pop(12)
    digierr.pop(12)
    osc.pop(12)
    oscerr.pop(12)
    diff.pop(12)
    differr.pop(12)
    
    c2 = TCanvas()
    g2 = TGraphErrors(len(osc), osc, digi, oscerr, digierr)
    axisTitles=";oscilloscope T[s];digitizer T[s]"
    gtitle = "digitizer_characterization"
    g2.SetNameTitle(gtitle,gtitle+axisTitles)
    g2.SetMarkerSize(1.3)
    g2.SetMarkerStyle(20)
    g2.SetMarkerColor(kBlue)
    g2.Fit("pol1")
    g2.Draw("ap")
    c2.SetTicks()
    c2.SetGrid()
    c2.Update()
    
    
    c3 = TCanvas()
    g3 = TGraphErrors(len(osc), osc, diff, oscerr, differr)
    axisTitles=";oscilloscope T[s];residuals [s]"
    gtitle = "digitizer_characterization_residuals"
    g3.SetNameTitle(gtitle,gtitle+axisTitles)
    g3.SetMarkerSize(1.3)
    g3.SetMarkerStyle(20)
    g3.SetMarkerColor(kBlue)
    g3.Fit("pol0")
    g3.Draw()
    c3.SetTicks()
    c3.SetGrid()
    c3.Update()

    import pdb; pdb.set_trace()

