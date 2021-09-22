from ROOT import *
import sys, os, re
from array import array
from statistics import mean, stdev

def GetFilesList(path):
    files = []
    folder = os.fsencode(path)
    for file in os.listdir(folder):
        filename = os.fsdecode(file)
        if "min" in filename and "max" in filename:
            min_ = int(re.findall(r'\d+', filename)[0])
            max_ = int(re.findall(r'\d+', filename)[1])
            files.append((min_,max_,GetFilesList(filename+"/")))
        if "run" in filename and "ALL" in filename:
            files.append(path+filename)
    return files

def GetVar(tree, varName):
    res = []
    for entry in tree:
        exec('globals()["var"]=entry.'+varName)
        res.append(var)
    return res



if __name__ == "__main__":

    gInterpreter.Declare("typedef struct {double sqFall=0, sqFallErr=0; double sqRise=0, sqRiseErr=0; double sqWdt = sqFall - sqRise; double sqWdtErr = TMath::Sqrt(sqFall*sqFall + sqRise*sqRise);} Square_Signal;")
    gInterpreter.Declare("typedef struct {Square_Signal start; Square_Signal stop; double dtFall = stop.sqFall - start.sqFall; double dtFallErr = TMath::Sqrt(stop.sqFallErr*stop.sqFallErr + start.sqFallErr*start.sqFallErr); double dtRise = stop.sqRise - start.sqRise; double dtRiseErr = TMath::Sqrt(stop.sqRiseErr*stop.sqRiseErr + start.sqRiseErr*start.sqRiseErr); } Decay_Event;") 
    
    runName = sys.argv[1]

    with open("ThOpt_run"+runName+".txt", "w+") as txt:
        acc, up, dwn, min_, max_ = array('d'), array('d'), array('d'), array('d'), array('d')
        for j,el in enumerate(GetFilesList(".")):
            minTh = float(el[0])
            maxTh = float(el[1])
            min_.append(minTh)
            max_.append(maxTh)
            fileList = el[2]
            for file in fileList:
                if runName in file:
                    f = TFile.Open(file)
                    a = f.Get("accepted").GetTitle()
                    u = f.Get("up decay").GetTitle()
                    d = f.Get("dwn decay").GetTitle()
                    f.Close()
                    aud = ((acc,a),(up,u),(dwn,d))
                    for res,cut in aud:
                        cutseffs = cut.split()
                        for i,name in enumerate(cutseffs):
                            if name == "cumulative":
                                numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
                                rx = re.compile(numeric_const_pattern, re.VERBOSE)
                                cumulative_eff = float(rx.findall(cutseffs[i+1])[0])
                                res.append(cumulative_eff)
                    txt.write(str(minTh)+"\t"+str(maxTh)+"\t"+str(acc[j])+"\t"+str(up[j])+"\t"+str(dwn[j])+"\n")

    c = TCanvas()
    g = TGraph(len(min_), min_, acc)
    axisTitles=";lower threshold limit [%pulse width]; #varepsilon [%]"
    gtitle = "cumulative efficiency accepted"
    g.SetNameTitle(gtitle,gtitle+axisTitles)
    g.SetMarkerSize(1)
    g.SetMarkerStyle(20)
    g.SetMarkerColor(kBlue)
    g.Draw("ap")
    c.SetTicks()
    c.SetGrid()
    c.Update()

    c1 = TCanvas()
    g1 = TGraph(len(min_), min_, up)
    gtitle = "cumulative efficiency up decay"
    g1.SetNameTitle(gtitle,gtitle+axisTitles)
    g1.SetMarkerSize(1)
    g1.SetMarkerStyle(20)
    g1.SetMarkerColor(kBlue)
    g1.Draw("ap")
    c1.SetTicks()
    c1.SetGrid()
    c1.Update()
    
    c2 = TCanvas()
    g2 = TGraph(len(min_), min_, dwn)
    gtitle = "cumulative efficiency dwn decay"
    g2.SetNameTitle(gtitle,gtitle+axisTitles)
    g2.SetMarkerSize(1)
    g2.SetMarkerStyle(20)
    g2.SetMarkerColor(kBlue)
    g2.Draw("ap")
    c2.SetTicks()
    c2.SetGrid()
    c2.Update() 

    fall, rise, startWdt, stopWdt, min_, max_ = array('d'), array('d'), array('d'), array('d'), array('d'), array('d')
    efall, erise, estartWdt, estopWdt = array('d'), array('d'), array('d'), array('d')
    for j,el in enumerate(GetFilesList(".")):
        minTh = float(el[0])
        maxTh = float(el[1])
        min_.append(minTh)
        max_.append(maxTh)
        fileList = el[2]
        for file in fileList:
            if runName in file:
                fi = TFile.Open(file)
                tr = fi.Get("DwnDecays")
                if not tr:
                    continue
                f = GetVar(tr,"DwnDecay.dtFall")
                r = GetVar(tr,"DwnDecay.dtRise")
                ws = GetVar(tr,"DwnDecay.stop.sqWdt")
                wr = GetVar(tr,"DwnDecay.start.sqWdt")
                fi.Close()
                aud = ((fall,efall,f),(rise,erise,r),(stopWdt,estopWdt,ws),(startWdt,estartWdt,wr))
                for res,eres,cut in aud:
                    res.append(mean(cut))
                    eres.append(stdev(cut))

    c3 = TCanvas()
    g3 = TGraph(len(min_), min_, fall)
    axisTitles=";lower threshold limit [%pulse width]; time [s]"
    gtitle = "fall"
    g3.SetNameTitle(gtitle,gtitle+axisTitles)
    g3.SetMarkerSize(1)
    g3.SetMarkerStyle(20)
    g3.SetMarkerColor(kBlue)
    g3.Draw("ap")
    c3.SetTicks()
    c3.SetGrid()
    c3.Update()

    c4 = TCanvas()
    g4 = TGraph(len(min_), min_, efall)
    gtitle = "stdev fall"
    g4.SetNameTitle(gtitle,gtitle+axisTitles)
    g4.SetMarkerSize(1)
    g4.SetMarkerStyle(20)
    g4.SetMarkerColor(kBlue)
    g4.Draw("ap")
    c4.SetTicks()
    c4.SetGrid()
    c4.Update()
    
    c5 = TCanvas()
    g5 = TGraph(len(min_), min_, rise)
    gtitle = "rise"
    g5.SetNameTitle(gtitle,gtitle+axisTitles)
    g5.SetMarkerSize(1)
    g5.SetMarkerStyle(20)
    g5.SetMarkerColor(kBlue)
    g5.Draw("ap")
    c5.SetTicks()
    c5.SetGrid()
    c5.Update()
    
    c6 = TCanvas()
    g6 = TGraph(len(min_), min_, erise)
    gtitle = "stdev rise"
    g6.SetNameTitle(gtitle,gtitle+axisTitles)
    g6.SetMarkerSize(1)
    g6.SetMarkerStyle(20)
    g6.SetMarkerColor(kBlue)
    g6.Draw("ap")
    c6.SetTicks()
    c6.SetGrid()
    c6.Update()
    
    
    c7 = TCanvas()
    g7 = TGraph(len(min_), min_, startWdt)
    gtitle = "startWdt"
    g7.SetNameTitle(gtitle,gtitle+axisTitles)
    g7.SetMarkerSize(1)
    g7.SetMarkerStyle(20)
    g7.SetMarkerColor(kBlue)
    g7.Draw("ap")
    c7.SetTicks()
    c7.SetGrid()
    c7.Update()


    c8 = TCanvas()
    g8 = TGraph(len(min_), min_, stopWdt)
    gtitle = "stopWdt"
    g8.SetNameTitle(gtitle,gtitle+axisTitles)
    g8.SetMarkerSize(1)
    g8.SetMarkerStyle(20)
    g8.SetMarkerColor(kBlue)
    g8.Draw("ap")
    c8.SetTicks()
    c8.SetGrid()
    c8.Update()

    c9 = TCanvas()
    g9 = TGraph(len(min_), min_, estartWdt)
    gtitle = "stdev startWdt"
    g9.SetNameTitle(gtitle,gtitle+axisTitles)
    g9.SetMarkerSize(1)
    g9.SetMarkerStyle(20)
    g9.SetMarkerColor(kBlue)
    g9.Draw("ap")
    c9.SetTicks()
    c9.SetGrid()
    c9.Update()


    c10 = TCanvas()
    g10 = TGraph(len(min_), min_, estopWdt)
    gtitle = "stdev stopWdt"
    g10.SetNameTitle(gtitle,gtitle+axisTitles)
    g10.SetMarkerSize(1)
    g10.SetMarkerStyle(20)
    g10.SetMarkerColor(kBlue)
    g10.Draw("ap")
    c10.SetTicks()
    c10.SetGrid()
    c10.Update()
    
    import pdb; pdb.set_trace()




                            


