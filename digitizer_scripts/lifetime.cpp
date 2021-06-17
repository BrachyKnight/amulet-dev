/*
Developed using ROOT 6.22/07 (guaranteed to work with version 6.22/07)
COMPILE WITH:
ON UBUNTU:
c++ -Wall -o lifetime lifetime.cpp amulet.cc amulet.h `root-config --cflags --glibs`
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


using ROOT::RDF::RNode;
using ROOT::RDataFrame;
using ROOT::VecOps::RVec;

#define ROOT_ANALYSIS 1

using namespace std;
typedef RVec<double> rvecdouble_t;

std::pair<double,double> DeltaTandErrSameCh( rvecdouble_t v, rvecdouble_t verr ){
	if (v.size() != 2){	
		cout<<"size of vector is "<<v.size()<<endl;  
		throw std::runtime_error( "\n check filters, which signals should I choose?" ); 
	}
	double DeltaT = v[1]-v[0];
	if (DeltaT < 0)
		throw std::runtime_error( "\n DeltaT<0 does not make any sense" ); //here there must be an error because these are consecutive square waves of the same waveform
	double DeltaT_err = sqrt(verr[1]*verr[1] + verr[0]*verr[0]);
	return make_pair( DeltaT, DeltaT_err );
}
std::pair<double,double> DeltaTandErrDifferentChs( rvecdouble_t v0, rvecdouble_t v0err,  rvecdouble_t v1, rvecdouble_t v1err ){
	if (!( (v1.size()==v0.size()) == 1 )){	
		cout<<"size of vectors are "<<v0.size()<<" and "<<v1.size()<<endl;  
		throw std::runtime_error( "\ncheck filters, which signals should I choose?" ); 
	}
	double DeltaT = v1[0]-v0[0];
	if (DeltaT < 0)
		cout<<"\nWARNING: DeltaT="<<DeltaT<<"<0, are columns in the right order?\n If yes you should treat these events\n"; //here it could happen? what are these signals?
	double DeltaT_err = sqrt(v0err[0]*v0err[0] + v1err[0]*v1err[0]);
	return make_pair( DeltaT, DeltaT_err );
}

std::pair<double,double> PositionMidpointStartStopSquareWaves(rvecdouble_t vUp, rvecdouble_t vDwn, rvecdouble_t v1Up, rvecdouble_t v1Dwn){
//returns an ordered pair with the mean position btw the times of the rising and the falling edges of the (START, STOP) square waves
	if ( vUp.size() == 2 && vDwn.size() == 2 ){
		return make_pair((vUp[0] + vDwn[0])/2., (vUp[1] + vDwn[1])/2.);
	}
	else if ( vUp.size() == 1 && vDwn.size() == 1 && v1Up.size() == 1 && v1Dwn.size() == 1 ){
		return make_pair((vUp[0] + vDwn[0])/2., (v1Up[0] + v1Dwn[0])/2.);
	}
	else{
		cout<<"size of vectors are:"<<endl
		    <<"first signal: "<<vUp.size()<<" and "<<vDwn.size()<<endl
		    <<"second signal: "<<v1Up.size()<<" and "<<v1Dwn.size()<<endl;
		throw std::runtime_error( "\n check filters, which signals should I choose?" );
	}
}

std::pair<double,double> WidthStartStopSquareWaves(rvecdouble_t vUp, rvecdouble_t vDwn, rvecdouble_t v1Up, rvecdouble_t v1Dwn){ 
//returns an ordered pair with the mean position btw the times of the rising and the falling edges of the (START, STOP) square waves
	if ( vUp.size() == 2 && vDwn.size() == 2 ){
		return make_pair((vUp[0] - vDwn[0]), (vUp[1] - vDwn[1]));
	}
	else if ( vUp.size() == 1 && vDwn.size() == 1 && v1Up.size() == 1 && v1Dwn.size() == 1 ){
		return make_pair((vUp[0] - vDwn[0]), (v1Up[0] - v1Dwn[0]));
	}
	else
	{
		cout<<"size of vectors are:"<<endl
		    <<"first signal: "<<vUp.size()<<" and "<<vDwn.size()<<endl
		    <<"second signal: "<<v1Up.size()<<" and "<<v1Dwn.size()<<endl;
		throw std::runtime_error( "\n check filters, which signals should I choose?" );
	}
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
	string outDir      = argv[1],
	       rootIn      = argv[2],
		   rootOut     = rootIn ,
		   rootOutName = rootOut.substr(rootOut.find_last_of("/")+1);
	rootOut = outDir + rootOutName;
	rootOut.replace(rootOut.end()-12, rootOut.end(), "Lifetimes.root"); //assuming it is called somethingTimings.root it translates it into somethingLifetimes.root
	rootOutName.replace(rootOutName.end()-12, rootOutName.end(), "Lifetimes");
	
	//open file and create RDataFrame
	//ROOT::EnableImplicitMT();       
	ROOT::RDataFrame d("Timings", rootIn.c_str());

	//define cuts
	string ch0NupNdwn = "(ch0Nup==ch0Ndwn)";
	string ch1NupNdwn = "(ch1Nup==ch1Ndwn)";
	string UpDecayCut = "(ch0Nup==2)";
	string DwnDecayCut ="((ch0Nup==ch1Nup)==1)";
		
	//define filters and print report
	auto upDecays = d.Filter(ch0NupNdwn, "Ch0 Nup == Ndwn").Filter(UpDecayCut, "UpDecay (ch0Nup==2)");
	auto dwnDecays =d.Filter(ch1NupNdwn, "Ch1 Nup == Ndwn").Filter(ch0NupNdwn, "Ch0 Nup == Ndwn").Filter(DwnDecayCut,"DwnDecay (ch0Nup==ch1Nup==1)");
	upDecays.Report()->Print();
	dwnDecays.Report()->Print();

	//evaluate lifetimes from digital signals defining new columns in the dataframe
	auto dtUpDec  = upDecays.Define("DeltaTandErrSameCh"  , DeltaTandErrSameCh                                  , {"ch0timeups", "ch0timeupsErr" } )
	                        .Define("Updt"                , [](std::pair<double,double> dt ){return dt.first;  }, {"DeltaTandErrSameCh"          } )
	                        .Define("UpdtErr"             , [](std::pair<double,double> dt ){return dt.second; }, {"DeltaTandErrSameCh"          } )
	                        .Define("STARTandSTOPpos"     , PositionMidpointStartStopSquareWaves                , {"ch0timeups", "ch0timedwns", "ch1timeups", "ch1timedwns"} )
	                        .Define("UpSTARTSquareMidPos" , [](std::pair<double,double> pos){return pos.first; }, {"STARTandSTOPpos"             } )
	                        .Define("UpSTOPSquareMidPos"  , [](std::pair<double,double> pos){return pos.second;}, {"STARTandSTOPpos"             } )
	                        .Define("wdtStartStop"        , WidthStartStopSquareWaves                           , {"ch0timeups", "ch0timedwns", "ch1timeups", "ch1timedwns"} )
	                        .Define("UpwdtSTART"          , [](std::pair<double,double> wdt){return wdt.first; }, {"wdtStartStop"                } )
	                        .Define("UpwdtSTOP"           , [](std::pair<double,double> wdt){return wdt.second;}, {"wdtStartStop"                } );

	auto dtDwnDec =dwnDecays.Define("DeltaTandErrDifferentChs", DeltaTandErrDifferentChs                        , {"ch0timeups", "ch0timeupsErr", "ch1timeups", "ch1timeupsErr"} )
	                        .Define("Dwndt"               , [](std::pair<double,double> dt){return dt.first;  } , {"DeltaTandErrDifferentChs"    } )
	                        .Define("DwndtErr"            , [](std::pair<double,double> dt){return dt.second; } , {"DeltaTandErrDifferentChs"    } )
	                        .Define("STARTandSTOPpos"     , PositionMidpointStartStopSquareWaves                , {"ch0timeups", "ch0timedwns", "ch1timeups", "ch1timedwns"} )
	                        .Define("DwnSTARTSquareMidPos", [](std::pair<double,double> pos){return pos.first; }, {"STARTandSTOPpos"             } )
	                        .Define("DwnSTOPSquareMidPos" , [](std::pair<double,double> pos){return pos.second;}, {"STARTandSTOPpos"             } )
	                        .Define("wdtStartStop"        , WidthStartStopSquareWaves                           , {"ch0timeups", "ch0timedwns", "ch1timeups", "ch1timedwns"} )
	                        .Define("DwnwdtSTART"         , [](std::pair<double,double> wdt){return wdt.first; }, {"wdtStartStop"                } )
	                        .Define("DwnwdtSTOP"          , [](std::pair<double,double> wdt){return wdt.second;}, {"wdtStartStop"                } );
	
	//save dataframes in external file
	RNode df_histUp = *dtUpDec.Snapshot("DeltaT UpDecays", rootOut.c_str(), {"Updt", "UpdtErr", "UpSTARTSquareMidPos", "UpSTOPSquareMidPos", "UpwdtSTART", "UpwdtSTOP"} );
	ROOT::RDF::RSnapshotOptions optsDwnSnapshot;
	optsDwnSnapshot.fMode = "UPDATE"; //to write the tree in the same file
	RNode df_histDwn = *dtDwnDec.Snapshot("DeltaT DwnDecays", rootOut.c_str(), {"Dwndt", "DwndtErr", "DwnSTARTSquareMidPos", "DwnSTOPSquareMidPos", "DwnwdtSTART", "DwnwdtSTOP"}, optsDwnSnapshot );
	cout<<"\nDATAFRAME CREATED IN FILE "<<rootOut.c_str()<<endl<<endl;
	
	return 0;
}
