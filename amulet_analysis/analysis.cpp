/*
Developed using ROOT 6.22/07 (guaranteed to work with version 6.22/07)
COMPILE WITH:
ON UBUNTU:
c++ -Wall -o lifetimev2 lifetimev2.cpp amulet.cc amulet.h `root-config --cflags --glibs`
*/
#include <TApplication.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <chrono>


#include <ROOT/RDataFrame.hxx>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TBrowser.h>
#include <TMath.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TString.h>
#include <TObject.h>
#include <TVirtualPad.h>
#include <TSystem.h>

#include "amulet.h"

using ROOT::RDF::RNode;
using ROOT::RDataFrame;
using ROOT::VecOps::RVec;

#define ROOT_ANALYSIS 1

using namespace std;
typedef RVec<double> rvecdouble_t;

void CompleteAnalysis( amulet* amu, string rootOut, string outDir, bool stabilityAnalysis  ){
	
	//amu->SnapshotAssociatedDataFramesInRootFile();
	
	//names
	string rootOutName = rootOut;
	rootOutName = rootOut.substr(rootOutName.find_last_of("/")+1);
	rootOutName.erase(rootOutName.end()-14, rootOutName.end());//assuming it is called somethingLifetimes.root it translates it into something
	string UpDecayName = rootOutName + "Up";
	string DwnDecayName = rootOutName + "Dwn";
	string outDecayCanvasName = rootOutName + "DecayFit";
	string combinedDecayName = rootOutName + "CombinedFit";
	vector<string> names = {"", outDecayCanvasName, outDir, UpDecayName, DwnDecayName, combinedDecayName};
	
	
	//parameters of the experiment
	double StartSignalWidth = 10e-6;
	double CFDSignalsWidth = 300e-9; //all data under this bin are not reliable
	double fitLimits[3];
	if(((TString)rootOutName).Contains("run1")){
		cout<<endl<<"RUN1: nessun materiale"<<endl;
		int nBinsHistoDecay = 300;//216; //60; //213;
		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
		std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run2")){
		cout<<endl<<"RUN2: alluminio + legno (da rigettare)"<<endl;
		int nBinsHistoDecay = 300;//216; //60; //213;
		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
		std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run3")){
		int nBinsHistoDecay = 300;//216; //60; //213;
		cout<<endl<<"RUN3: solo alluminio"<<endl;
		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
		std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run4")){
		cout<<endl<<"RUN4: sale"<<endl;
		int nBinsHistoDecay = 286;//244 consigliato per non normalizzata con due tau;
		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
		std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run5")){
		cout<<endl<<"RUN5: carbonio, again"<<endl;
		int nBinsHistoDecay = 100;
		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
		std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run6")){
                cout<<endl<<"RUN6: studi di fondo"<<endl;
                int nBinsHistoDecay = 200;
                auto init = std::initializer_list<double>({25e-6, 41e-6, (double)nBinsHistoDecay});
                std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run7")){
                cout<<endl<<"RUN7: carbonio, again!"<<endl;
                int nBinsHistoDecay = 100;
       		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
                std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run8")){
                cout<<endl<<"RUN8 campo magnetico spento"<<endl;
		int nBinsHistoDecay = 100;
		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
		std::copy(init.begin(), init.end(), fitLimits);
	}else if(((TString)rootOutName).Contains("run9")){
                cout<<endl<<"RUN9 campo magnetico spento"<<endl;
		int nBinsHistoDecay = 40;
		auto init = std::initializer_list<double>({CFDSignalsWidth*(1+0.7), (StartSignalWidth+CFDSignalsWidth)*(1-0.05), (double)nBinsHistoDecay});
		std::copy(init.begin(), init.end(), fitLimits);
	}else{
		cout<<"RUN NOT IMPLEMENTED"<<endl;
		throw std::runtime_error("you must implement binning also for this run");
	}
	
	amu->Perform_Lifetime_Fit(fitLimits, names, true, "RSE");
/*	
	//--------------------------------SIGNAL WIDTH and START and STOP signal analysis-----------------------------------------
	bool STARTandSTOP_SignalsAnalysis = false;
	if(STARTandSTOP_SignalsAnalysis){
		//--- SIGNAL WIDTH analysis ---
		int nBinsSTARTWidth = 2000;
		double lowerRangeSTARTWidth = 0.10e-6;
		double upperRangeSTARTWidth = 0.30e-6;
		int nBinsSTOPWidth = 2000;
		double lowerRangeSTOPWidth = 0.10e-6;
		double upperRangeSTOPWidth = 0.30e-6;
		amu->SignalWidth_Analysis({(double)nBinsSTARTWidth, lowerRangeSTARTWidth, upperRangeSTARTWidth, (double)nBinsSTOPWidth, lowerRangeSTOPWidth, upperRangeSTOPWidth}, "STARTandSTOP_SignalsAnalysis/WIDTH");

		//--- START and STOP signal analysis ---
		int nBinsSTART = fitLimits[2];
		double lowerRangeSTART = 0.;
		double upperRangeSTART = StartSignalWidth + 0.1*StartSignalWidth;
		//STOP signal should be everytime in the same position
		int nBinsSTOP = 2000;
		double lowerRangeSTOP = 9.8e-6;
		double upperRangeSTOP = 10.9e-6;
		amu->StartAndStop_analysis({(double)nBinsSTART, lowerRangeSTART, upperRangeSTART, (double)nBinsSTOP, lowerRangeSTOP, upperRangeSTOP}, "STARTandSTOP_SignalsAnalysis/POSITION");
	}
	
	if(stabilityAnalysis){
		bool stability1D = false;
		if(stability1D){
			//amu->SetDoPulls(false);
			names = {"stabilityAnalysis/1D", outDecayCanvasName, outDir, UpDecayName, DwnDecayName, combinedDecayName};
			int N_iterations = 800-30;
			int nMin = 30; 
			int nMax = 800;
			double lowerLimit = fitLimits[0];
			double upperLimit = fitLimits[1];
			vector<double> iter_settings = {(double)N_iterations, (double)nMin, (double)nMax, lowerLimit, upperLimit};
			vector<vector<double>> listToDraw = {};//{{lowerLimit, upperLimit, 10}, {lowerLimit, upperLimit, 208}}; //a list of lists of the ones you want to draw: {{lowLim, upLim, nBin}, ...}. if {{1}} draw everything.
			
			auto t1 = std::chrono::high_resolution_clock::now();
			amu->Stability_Analysis1DnBins(names, iter_settings, listToDraw, "EP");
			auto t2 = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::minutes>( t2 - t1 ).count();
			auto mean_duration_per_iteration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count()/N_iterations;
			std::cout<<"1D N_iters = "<<N_iterations<<" completed in in "<<duration<<" minutes"<<endl<<"needed mean of "<<mean_duration_per_iteration<<" ms/iteration"<<endl<<endl;
		}
		
		bool stability2D = false;
		if(stability2D){
			amu->SetDoPulls(false); //FIXME: here segmentation violation if set to true
			names = {"stabilityAnalysis/2D", outDecayCanvasName, outDir, UpDecayName, DwnDecayName, combinedDecayName};
			int N_iterationsBin = (350-100)/2; //125; //if you put a value larger than nMax-nMin it will do it increasing one bin at a time
			int nMin = 100; 
			int nMax = 350;
			int N_iterationsRange = 70; //70;
			double lowerLimitMax = fitLimits[0]*(1-0.3);
			double upperLimitMax = fitLimits[1]*(1+0.03);
			double rangeMin = 8e-6;
			vector<double> iter_settings = {(double)N_iterationsBin, (double)nMin, (double)nMax, (double)N_iterationsRange, lowerLimitMax, upperLimitMax, rangeMin};
			vector<vector<double>> listToDraw = {}; //a list of lists of the ones you want to draw: {{lowLim, upLim, nBin}, ...}. if {{1}} draw everything.
			
			auto t1 = std::chrono::high_resolution_clock::now();
			amu->Stability_Analysis2DnBinsSymRange(names, iter_settings, listToDraw, "SURF7 Z");
			auto t2 = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::minutes>( t2 - t1 ).count();
			auto mean_duration_per_iteration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count()/(N_iterationsBin*N_iterationsRange);
			std::cout<<"2D N_iters = "<<N_iterationsBin*N_iterationsRange<<" completed in in "<<duration<<" minutes"<<endl<<"needed mean of "<<mean_duration_per_iteration<<" ms/iteration"<<endl<<endl;
		}
		
		bool stability3D = false;
		if(stability3D){
			amu->SetDoPulls(false);
			names = {"stabilityAnalysis/3D", outDecayCanvasName, outDir, UpDecayName, DwnDecayName, combinedDecayName};
			int N_iterationsBin = 10; //if you put a value larger than nMax-nMin it will do it increasing one bin at a time
			int nMin = 100; 
			int nMax = 350;
			int N_iterationsLowerLimit = 10; //if you put -1 or a value larger than nMax-nMin it will do it increasing one bin at a time
			double lowerLimitMin = fitLimits[0]*(1-0.6);
			double lowerLimitMax = fitLimits[0]*(1+0.6);
			int N_iterationsUpperLimit = 10;
			double upperLimitMin = fitLimits[1]*(1-0.35);
			double upperLimitMax = fitLimits[1]*(1+0.1);
			vector<double> iter_settings = {(double)N_iterationsBin, (double)nMin, (double)nMax, (double)N_iterationsLowerLimit, lowerLimitMin, lowerLimitMax, (double)N_iterationsUpperLimit, upperLimitMin, upperLimitMax};
			vector<vector<double>> listToDraw = {}; //a list of lists of the ones you want to draw: {{lowLim, upLim, nBin}, ...}. if {{1}} draw everything.
			gStyle->SetCanvasPreferGL(true);

			auto t1 = std::chrono::high_resolution_clock::now();
			amu->Stability_Analysis3DnBinsLowUpLims(names, iter_settings, listToDraw, "BOX2 Z");
			auto t2 = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::minutes>( t2 - t1 ).count();
			auto mean_duration_per_iteration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count()/(N_iterationsBin*N_iterationsLowerLimit*N_iterationsUpperLimit);
			std::cout<<"3D N_iters = "<<N_iterationsBin*N_iterationsLowerLimit*N_iterationsUpperLimit<<" completed in in "<<duration<<" minutes"<<endl<<"needed mean of "<<mean_duration_per_iteration<<" ms/iteration"<<endl<<endl;
		}
	}*/
}

int main(int argc, char** argv)
{
	//set names from input line
	if ( argc != 3 && argc != 4 ){
		std::cout<<"\nINPUT ERROR"<<endl
				 <<"needs:\noutputDir \nrootInputFile\n"<<endl
				 <<"as an example \n(ubuntu):\n./lifetime results/ results/run1meas1Timings.root"<<endl;
		return 1;
	}
	bool setRootBatchMode = (argc == 4) ? true : false;
	string  outDir      = argv[1],
	       	rootIn      = argv[2],
		rootOut     = rootIn ,
		rootOutName = rootOut.substr(rootOut.find_last_of("/")+1);
	rootOut = outDir + rootOutName;
	rootOut.replace(rootOut.end()-12, rootOut.end(), "Lifetimes.root"); //assuming it is called somethingTimings.root it translates it into somethingLifetimes.root
	rootOutName.replace(rootOutName.end()-12, rootOutName.end(), "Lifetimes");

	//open file and create RDataFrame
	//ROOT::EnableImplicitMT();
	ROOT::RDataFrame df_histUp("UpDecays", rootIn.c_str());
	ROOT::RDataFrame df_histDwn("DwnDecays", rootIn.c_str());

	#ifdef ROOT_ANALYSIS
	//SET ROOT ENVIRONMENT
	gStyle->SetOptTitle(true);
	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(110);
	(setRootBatchMode) ? gROOT->SetBatch(true) : gROOT->SetBatch(false);
	TApplication *myApp = (setRootBatchMode) ? NULL : new TApplication("ROOT example", &argc, argv);

	// !!! FIT FUNCTION !!!
	bool DoPulls = false;
	//here you choose the function and you can also fix parameters
	// !!! MEGA WARNING !!!
	//devi scrivere la funzione utilizzando esattamente gli stessi caratteri perche il confronto
	//viene fatto comparando delle stringhe e sono stati previsti solo alcuni prototipi
	//si veda la funzione partialFit nella classe amulet

	//CASO BASIC
	string myFormula = "[N0]*exp(-x/[#tau])+[B]"; //caso base
	//CASO BASIC NORMALIZZATO
	//string myFormula = "(exp(-x/[#tau])+[b])/([#tau]*(exp(-[lowerLimit]/[#tau])-exp(-[upperLimit]/[#tau]))-[lowerLimit]*[b]+[upperLimit]*[b])";//caso base normalizzato
	//CASO BASIC SOLO FONDO
	//string myFormula = "[B]";

	//CASO DUE MATERIALI DIVERSI CON RAPPORTO
	bool DoNormalize = false;
	//string myFormula = "([N0]/([f]+1))*([f]*exp(-x/[#tau])+exp(-x/[#tau1])) + [B]";
	//CASO DUE MAT NORMALIZZATO
	//string myFormula = "([f]*exp(-x/[#tau])+exp(-x/[#tau1])+[b])/(-[lowerLimit]*[b] + [upperLimit]*[b] + (exp(-[lowerLimit]/[#tau]) - exp(-[upperLimit]/[#tau]))*[f]*[#tau] + (e^(-[lowerLimit]/[#tau1]) - exp(-[upperLimit]/[#tau1]))*[#tau1])";
	//fissando #tau al valore misurato senza alluminio e #tau1 al valore noto dall'alluminio
	//string myFormula = "([f]*exp(-x/[2.11927e-06])+exp(-x/[0.880e-6])+[b])/(-[lowerLimit]*[b] + [upperLimit]*[b] + (exp(-[lowerLimit]/[2.11927e-06]) - exp(-[upperLimit]/[2.11927e-06]))*[f]*[2.11927e-06] + (e^(-[lowerLimit]/[0.880e-6]) - exp(-[upperLimit]/[0.880e-6]))*[0.880e-6])";
	//fissando #tau1 al valore dell'alluminio e lasciando #tau libero
	//string myFormula = "([f]*exp(-x/[#tau])+exp(-x/[0.880e-6])+[b])/(-[lowerLimit]*[b] + [upperLimit]*[b] + (exp(-[lowerLimit]/[#tau]) - exp(-[upperLimit]/[#tau]))*[f]*[#tau] + (e^(-[lowerLimit]/[0.880e-6]) - exp(-[upperLimit]/[0.880e-6]))*[0.880e-6])";
	//fissando #tau al valore del vuoto e #tau1 al valore del carbonio
	//string myFormula = "([f]*exp(-x/[2.197e-06])+exp(-x/[2.026e-6])+[b])/(-[lowerLimit]*[b] + [upperLimit]*[b] + (exp(-[lowerLimit]/[2.197e-06]) - exp(-[upperLimit]/[2.197e-06]))*[f]*[2.197e-06] + (e^(-[lowerLimit]/[2.026e-6]) - exp(-[upperLimit]/[2.026e-6]))*[2.026e-6])";
	//fissando #tau al valore del vuoto e #tau1 libero
	//string myFormula = "([f]*exp(-x/[2.197e-6])+exp(-x/[#tau1])+[b])/(-[lowerLimit]*[b] + [upperLimit]*[b] + (exp(-[lowerLimit]/[2.197e-6]) - exp(-[upperLimit]/[2.197e-6]))*[f]*[2.197e-6] + (e^(-[lowerLimit]/[#tau1]) - exp(-[upperLimit]/[#tau1]))*[#tau1])";


	//OLD
	//string myFormula = "([N0]/([f]+1))*([f]*exp(-x/[#tau])+exp(-x/[#tau1])) + [B]"; //con rapporto //tutto libero
	//string myFormula = "[N0]*exp(-x/[2.197e-6])+[N1]*exp(-x/[#tau1])+[B]"; //senza rapporto //per stimare [#tau1] dato [#tau]
	//string myFormula = "[N0]*exp(-x/[2.137e-6])+[N1]*exp(-x/[#tau1])+[B]"; //senza rapporto //per stimare [#tau1] dato [#tau] misurato in precedenza
	//string myFormula = "([N0]/([f]+1))*([f]*exp(-x/[2.197e-6])+exp(-x/[#tau1])) + [B]"; //con rapporto //per stimare [#tau1] dato [#tau]
	//string myFormula = "([N0]/([f]+1))*([f]*exp(-x/[2.197e-6])+exp(-x/[2.026e-6])) + [B]"; //per stimare [f] assumendo i tau libero e del carbonio
	//string myFormula = "([N0]/([f]+1))*([f]*exp(-x/[2.197e-6])+exp(-x/[1.26e-6])) + [B]"; //per stimare [f] assumendo i tau libero e del'alluminio media pesata calcolata rob
	//string myFormula = "([N0]/([f]+1))*([f]*exp(-x/[2.197e-6])+exp(-x/[0.88e-6])) + [B]"; //per stimare [f] assumendo i tau libero e del'alluminio
	//string myFormula = "[N0]*exp(-x/[#tau])+[N1]*exp(-x/[#tau1])+[B]"; //senza rapporto //tutto libero


	//create amulet object (contains all analysis methods for muon experment with digitizer)
	string rootOutResults = outDir + rootOutName + "Analyzed.root";
	auto amu = new amulet(df_histUp, df_histDwn, rootOutResults.c_str(),"UPDATE");
	amu->SetDoNormalize(DoNormalize);
	amu->SetDoPulls(DoPulls);
	amu->SetFunction(myFormula);
	if(((TString)rootOutName).Contains("run1"))
		amu->SetTotalRate(0.0428968); //0.042Hz corrisponde a circa 2.57/minuto
	else
		cout<<"total rate for this run not yet evaluated"<<endl;
	CompleteAnalysis(amu, rootOut, outDir, true);


/*
	//-------------------New Signal-Width-Based Filters (F2)-----------------------
	//define cuts
	bool signalWidthBasedFilters = false;
	if(signalWidthBasedFilters){
		strMap_t UpFilters;
		strMap_t DwnFilters;
		if(((TString)rootOutName).Contains("run1")){
			cout<<endl<<"aplying F2 filters for run 1:"<<endl;
			string UpSTOPposCut = "( (UpSTOPSquareMidPos>=10.02e-6 && UpSTOPSquareMidPos<=10.05e-6) ||  (UpSTOPSquareMidPos>=10.745e-6 && UpSTOPSquareMidPos <= 10.768e-6))";
			string DwnSTOPposCut = "( (DwnSTOPSquareMidPos<=10.042e-6 && DwnSTOPSquareMidPos>=10.02e-6) ||  (DwnSTOPSquareMidPos>=10.74e-6 && DwnSTOPSquareMidPos <= 10.765e-6))";
			UpFilters = { {"Up STOP-position-based filter", UpSTOPposCut} };
			DwnFilters = { {"Dwn STOP-position-based filter", DwnSTOPposCut} };
		}else if(((TString)rootOutName).Contains("run3")){
			cout<<"filters for run3 not yet implemented"<<endl;
		}else{
			cout<<"filters for this run not yet implemented"<<endl;
		}

		if ( !(UpFilters.size() == 0 && DwnFilters.size() == 0) ){
			string rootOutF2Results = outDir + rootOutName + "Analyzed_F2.root";
			auto amuFiltered = new amulet(amu->ApplyFiltersAndPrintReport(UpFilters, DwnFilters, rootOutF2Results.c_str(), "UPDATE"));
			amuFiltered->SetDoNormalize(DoNormalize);
			amuFiltered->SetDoPulls(DoPulls);
			amuFiltered->SetFunction(myFormula);
			//CompleteAnalysis(amuFiltered, outDir + rootOutName + "_F2Lifetimes.root", outDir, false);
			delete amuFiltered;
		}
	}
*/
	delete amu;
	if (!setRootBatchMode) myApp->Run();
	#endif
	return 0;
}
