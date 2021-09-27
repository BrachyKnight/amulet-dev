from ROOT import *
import sys, os, re
from array import array

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


if __name__ == "__main__":

    if len(sys.argv) != 3:
        raise Exception("input error")

    runName = sys.argv[1]
    path = sys.argv[2]

    with open("ThOpt_run"+runName+".txt", "w+") as txt:
        acc, up, dwn, min_, max_ = array('d'), array('d'), array('d'), array('d'), array('d')
        for j,el in enumerate(GetFilesList(path)):
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
    axisTitles=";lower threshold limit [%pulse height]; #varepsilon [%]"
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

        

    import pdb; pdb.set_trace()




                            


