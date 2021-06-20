/* 
 * Author: Massimo Girola
 * Date: June 2021
 * Developed using ROOT 6.24/00
 * c++ -Wall -o lifetime lifetime.cpp `root-config --cflags --glibs`
 */
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

#include <TMath.h>
#include <TString.h>
#include <TSystem.h> //gInterpreter->Declare()
#include <TChain.h>

using ROOT::RDF::RNode;
using ROOT::RDataFrame;
using ROOT::VecOps::RVec;
using std::vector;
using std::runtime_error;
using std::cout;
using std::endl;
using std::string;

//struct representing a square signal: falling edge and rising edge 
typedef struct {
	double sqFall, sqFallErr;
	double sqRise, sqRiseErr;
	double sqWdt = sqFall - sqRise;
	double sqWdtErr = TMath::Sqrt(sqFall*sqFall + sqRise*sqRise);
} Square_Signal;

//struct representing a decay event: two square signals start and stop
typedef struct {
	Square_Signal start;
	Square_Signal stop;
	double dtFall = stop.sqFall - start.sqFall;
	double dtFallErr = TMath::Sqrt(stop.sqFallErr*stop.sqFallErr + start.sqFallErr*start.sqFallErr);
	double dtRise = stop.sqRise - start.sqRise;
	double dtRiseErr = TMath::Sqrt(stop.sqRiseErr*stop.sqRiseErr + start.sqRiseErr*start.sqRiseErr);
} Decay_Event;

//evaluate the decay events for UP decays
Decay_Event GetSameChSquare( RVec<double> vFall, RVec<double> vFallErr, RVec<double> vRise, RVec<double> vRiseErr  ){
	if (!(vFall.size() == 2 && vRise.size() == 2 && vFallErr.size() == 2 && vRiseErr.size() == 2) ){ //exactly two signals in the same pulse	
		cout<<"size of vector is "<<vRise.size()<<endl;  
		throw std::runtime_error( "\n check filters, which signals should I choose?" ); 
	}
	Square_Signal start{ vFall[0], vFallErr[0], vRise[0], vRiseErr[0] };
	Square_Signal stop { vFall[1], vFallErr[1], vRise[1], vRiseErr[1] };
	Decay_Event decay{ start, stop };
	if (decay.dtFall < 0 || decay.dtRise < 0)
		throw std::runtime_error( "\n DeltaT<0 does not make any sense" ); //consecutive signals in the same pulse
	return decay;
}

//evaluate the decay events for DOWN decays
Decay_Event GetDiffChSquare(	RVec<double> v0Fall, RVec<double> v0FallErr,  RVec<double> v1Fall, RVec<double> v1FallErr,
	       			RVec<double> v0Rise, RVec<double> v0RiseErr,  RVec<double> v1Rise, RVec<double> v1RiseErr ){
	if (	!(  (v1Fall.size()   ==v0Fall.size()   ) && (v1Rise.size()   ==v0Rise.size()   )
		&&  (v1FallErr.size()==v0FallErr.size()) && (v1RiseErr.size()==v0RiseErr.size()) 
		&&  (v1Rise.size()   == 1 )            )
	   ){   //exactly one signal in each pulse (one for each channel)
		cout<<"size of vectors are "<<v0Fall.size()<<" and "<<v1Fall.size()<<endl;  
		throw std::runtime_error( "\ncheck filters, which signals should I choose?" ); 
	}
	Square_Signal stop { v1Fall[0], v1FallErr[0], v1Rise[0], v1RiseErr[0] };
	Square_Signal start{ v0Fall[0], v0FallErr[0], v0Rise[0], v0RiseErr[0] };
	Decay_Event decay{ start, stop };
	if (decay.dtFall < 2*abs(decay.dtRiseErr) || decay.dtRise < 2*abs(decay.dtFallErr))
		cout<<"\n !!!!!!!WARNING!! DeltaT<0 IN DOWN DECAYS, MEANING???????"<<endl; //in principle it could since we are in two different chs, however I check it is compatible with the error
	return decay;
}

int main(int argc, char** argv)
{
	//set names from input line
	if ( !(argc >= 4) ){ 
		std::cout<<"\nINPUT ERROR"<<endl
		<<"needs:\n./path/to/executable/ path/to/output/directory/ \npath/to/root/input/files/directory/"<<endl
		<<"RUNnumber MEASnumbers (dont specify for all)"<<endl
		<<"examples: (for saving output into ../DAQprocessed/ and taking data from ../DAQpreproccesed/"
		<<"./executable ../DAQprocessed/ ../DAQpreproccesed/ 1 2 //for run1meas2"<<endl
		<<"./executable ../DAQprocessed/ ../DAQpreprocessed/ 4 6a //for run4meas6a"<<endl
		<<"./executable ../DAQprocessed/ ../DAQpreprocessed/ 9 1 2 3 //for run9meas1 run9meas2 run9meas3"<<endl
		<<"./executable ../DAQprocessed/ ../DAQpreprocessed/ 1 //for all run1"<<endl;
		return 1;
	}

	//declare the struct to the CLING compiler
	gInterpreter->Declare("typedef struct {double sqFall, sqFallErr; double sqRise, sqRiseErr; double sqWdt = sqFall - sqRise; double sqWdtErr = TMath::Sqrt(sqFall*sqFall + sqRise*sqRise);} Square_Signal;");
	gInterpreter->Declare("typedef struct {Square_Signal start; Square_Signal stop; double dtFall = stop.sqFall - start.sqFall; double dtFallErr = TMath::Sqrt(stop.sqFallErr*stop.sqFallErr + start.sqFallErr*start.sqFallErr); double dtRise = stop.sqRise - start.sqRise; double dtRiseErr = TMath::Sqrt(stop.sqRiseErr*stop.sqRiseErr + start.sqRiseErr*start.sqRiseErr); } Decay_Event;");
	
	//setting I\O file names
	string outDir = argv[1];
	string inDir  = argv[2];
	string runN = argv[3];
	string OutFileName = outDir+"run"+runN+"meas";
	std::vector<string> measNs;	
	std::cout<<"I will process the following files:"<<std::endl;
	if( argc == 4 ){
		measNs.push_back("*");
		printf("%s/run%smeas*_amulet\n",inDir.c_str(),runN.c_str());
		OutFileName += "_ALL";
	}else{
		for(int i = 4; i<argc; i++)
			measNs.push_back(argv[i]);
		for(auto measN : measNs){
			printf("%s/run%smeas%s_amulet\n",inDir.c_str(),runN.c_str(),measN.c_str());
			OutFileName += "_" + measN;
		}
	}
	OutFileName += ".root";
	std::cout<<"and save the results in "<<OutFileName<<std::endl<<std::endl;

	//open file(s)
	TChain myChain("amulet");
	for(auto measN : measNs)
		myChain.Add( (inDir+"run"+runN+"meas"+measN+"_amulet.root").c_str() );
	
	//use all cores of your machine MultiThreading to speed up analysis
	//(read RDataFrame documentation for MT warnings)
	ROOT::EnableImplicitMT();
	
	//create df
	ROOT::RDataFrame d(myChain);

	//define cuts
	string ch0NupNdwn = "(ch0Nup==ch0Ndwn)";
	string ch1NupNdwn = "(ch1Nup==ch1Ndwn)";
	string UpDecayCut = "(ch0Nup==2)";
	string DwnDecayCut ="((ch0Nup==ch1Nup) && (ch0Nup==1))";
		
	//apply cuts (filters) and print report
	auto upDecays = d.Filter(ch0NupNdwn, "Ch0 Nup == Ndwn").Filter(UpDecayCut, "UpDecay (ch0Nup==2)");
	auto dwnDecays =d.Filter(ch1NupNdwn, "Ch1 Nup == Ndwn").Filter(ch0NupNdwn, "Ch0 Nup == Ndwn").Filter(DwnDecayCut,"DwnDecay (ch0Nup==ch1Nup==1)");
	upDecays.Report()->Print();
	dwnDecays.Report()->Print();

	//evaluate lifetimes from digital signals defining new columns in the dataframe
	auto dtUpDec  = upDecays.Define("decayUP"   , GetSameChSquare, {"ch0timedwns","ch0timedwnsErr","ch0timeups","ch0timeupsErr" } );
	auto dtDwnDec = dwnDecays.Define("decayDOWN", GetDiffChSquare, {"ch0timeups", "ch0timeupsErr", "ch1timeups", "ch1timeupsErr","ch0timedwns", "ch0timedwnsErr", "ch1timedwns", "ch1timedwnsErr"} );
	
	//save dataframes in external file
	RNode df_histUp = *dtUpDec.Snapshot("UpDecays", OutFileName.c_str(), {"measN","idx","decayUP"} );
	ROOT::RDF::RSnapshotOptions optsDwnSnapshot;
	optsDwnSnapshot.fMode = "UPDATE"; //to write the tree in the same file
	RNode df_histDwn = *dtDwnDec.Snapshot("DwnDecays", OutFileName.c_str(), {"measN","idx","decayDOWN"},  optsDwnSnapshot );



	cout<<"\nDATAFRAME CREATED IN FILE "<<OutFileName<<endl<<endl;
	
	return 0;
}
