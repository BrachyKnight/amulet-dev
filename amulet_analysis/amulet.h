//muon analysis class for laboratory particle physics UNIMIB
//amulet = Analysis MUon LifETime 
//Developed using ROOT 6.22/07 (guaranteed to work with version 6.22/07)

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <stdexcept>

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
#include <TList.h>

#ifndef amulet_h
#define amulet_h

using ROOT::RDataFrame;
using ROOT::RDF::RNode;
using ROOT::RDF::RSnapshotOptions;

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::abs;
using std::map;
using std::pair;

typedef pair<string,double> pairstrdoub_t;
typedef vector<pairstrdoub_t> vecWithName_t;
typedef vector<pair<string, vecWithName_t >> nestedVecWithName_t;
typedef vector<TH1D*> vecHisto_t;
typedef vector<vector<TH1D*>> nestedVecHisto_t;
typedef map<string, string> strMap_t;


class amulet
{
public:
    //constructors
    amulet();
    amulet(RNode up, RNode dwn, const char* OutFilePath, const char* OutFileOption );
    amulet(RDataFrame up, RDataFrame dwn, const char* OutFilePath, const char* OutFileOption );
    amulet(const amulet &old_amu); //copy constr
    ~amulet();
    
    //output file methods
    inline const bool IsFileOK() { return ( !(_OutFile==NULL) && !(_OutFile->IsZombie()) ); };
    inline const short int ChangeOption( const char* OutFileOption  ) { if(_OutFile!=NULL) return _OutFile->ReOpen(OutFileOption); else return -1; }; //returns 0 in case the mode was successfully modified, 1 in case the mode did not change and -1 in case of failure.
    const bool SetOutputFile(const char* rootOutName, const char* rootFileOpt);
    inline TFile* GetOutputFile() { return _OutFile; };
    const bool MakeAndChangeCurrentDir( const char* subdirName );
    
	//utility methods
	inline const char* GetOption() {if(_OutFile!=NULL) return _OutFile->GetOption(); else return "Failed to retrieve option: _OutFile==NULL";}
	inline const string GetFunction() { return _function; };
	inline void SetFunction(string function ) { _function = function; };
	inline const short int GetVerboseLevel() { return _verbosity; };
	inline void SetVerboseLevel(short int verbosity ) { _verbosity = verbosity; };
	inline double GetTotalRate() { return _totalRate; };
	inline void SetTotalRate(double rate ) { _totalRate = rate; };
	inline bool IsSetTotalRate() { if( _totalRate == 0 ) return false; else return true; }
	inline const bool IsDoNormalize() { return _DoNormalize; };
	inline void SetDoNormalize(bool DoNormalize) { _DoNormalize = DoNormalize; };
	inline const bool IsDoPulls() { return _DoPulls; };
	inline void SetDoPulls(bool DoPulls) { _DoPulls = DoPulls; };
	
    //filter and snapshot methods
    amulet ApplyFiltersAndPrintReport( strMap_t UpFilterMap, strMap_t DwnFilterMap, const char* NewOutFileName, const char* option );
    void SnapshotAssociatedDataFramesInRootFile();
    
    //analysis methods
    nestedVecWithName_t Perform_Lifetime_Fit(double* fitLims,  vector<string> names, bool drawPlots, const char* fitOpts);
    void Stability_Analysis1DnBins(vector<string> names, vector<double> iter_settings,  vector<vector<double>> drawList, TString DrawOpts = "EP", bool drawErrorsOnPlot = true, bool plotOnlyValids = true);
    void Stability_Analysis2DnBinsSymRange(vector<string> names, vector<double> iter_settings,  vector<vector<double>> drawList, TString DrawOpts = "SURF2", bool drawErrorsOnPlot = false, bool plotOnlyValids = true);
    void Stability_Analysis3DnBinsLowUpLims(vector<string> names, vector<double> iter_settings,  vector<vector<double>> drawList, TString DrawOpts = "BOX2 Z", bool drawErrorsOnPlot = false, bool plotOnlyValids = true);
    void SignalWidth_Analysis(std::vector<double> binnings, const char* outDirNameInFile);
	void StartAndStop_analysis(std::vector<double> binnings, const char* outDirNameInFile);

private:
    RNode _dfUp;
    RNode _dfDwn;
    TFile* _OutFile;
    TString _OutFileName;
    string _function = "[N0]*exp(-x/[#tau])+[B]";
    double _totalRate = 0;
    bool _DoNormalize = false;
    short int _verbosity = 1;
    bool _DoPulls = true;

    template<class T = ROOT::RDF::RResultPtr<TH1D>> vecWithName_t partialFit(T histo, const double fitRange[2],  string name, const char* fitOpts, bool drawPlots, vector<TList*> &histoDrawLists);
    template<class T = ROOT::RDF::RResultPtr<TH1D>> vecWithName_t partialPulls(T histo, const TF1* FitFunction, string name, bool drawPlots, TList *&histoPullToDrawList);
    void drawHistoList(TList* DrawList, TString CanvasName, string outDecayCanvasDir="");
    bool drawThisOne( double* fitLims, vector<vector<double>> drawList);
    template<class T = nestedVecHisto_t> void drawNestedVecHistos ( T histoMatrix, vector<TCanvas*> canvases, const char* outDirNameInFile, TString DrawOpts = "P");
    Int_t countpads(TVirtualPad *pad);
};

#endif
