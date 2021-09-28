/* 
 * Author: Massimo Girola
 * Date: June 2021
 * Developed using ROOT 6.24/00
 * c++ -Wall -o amulet amulet.cpp `root-config --cflags --glibs`
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
#include <map>
#include <utility>

#include <ROOT/RDataFrame.hxx>

#include <TMath.h>
#include <TString.h>
#include <TSystem.h> //gInterpreter->Declare()
#include <TChain.h>
#include <TNamed.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>

using ROOT::RDF::RNode, ROOT::RDataFrame, ROOT::VecOps::RVec;
using std::vector, std::pair, std::runtime_error, std::cout, std::endl, std::string, std::map, std::make_pair, std::to_string;

//struct representing a square signal: falling edge and rising edge 
typedef struct {
	double sqFall=0, sqFallErr=0;
	double sqRise=0, sqRiseErr=0;
	double sqWdt = sqRise - sqFall;
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

void WriteWidthHistos( RNode dfUp, RNode dfDwn, TFile* f ){
	double wdt = 0.5e-9;
	double min=0, max=0;
	
	min = dfDwn.Min("Decay.start.sqWdt").GetValue();
	max = dfDwn.Max("Decay.start.sqWdt").GetValue();
	TH1D startWdtDwn = *dfDwn.Histo1D<double>({"Start_Width_Dwn",
						"Start Width;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min, max},
						"Decay.start.sqWdt");
	
	min = dfDwn.Min("Decay.stop.sqWdt").GetValue();
	max = dfDwn.Max("Decay.stop.sqWdt").GetValue();
	TH1D stopWdtDwn = *dfDwn.Histo1D<double>({"Stop_Width_Dwn",
						"Stop Width;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min, max},
						"Decay.stop.sqWdt");
	
	min = dfUp.Min("Decay.start.sqWdt").GetValue();
	max = dfUp.Max("Decay.start.sqWdt").GetValue();
	TH1D startWdtUp = *dfUp.Histo1D<double>({"Start_Width_Up",
						"Start Width;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min, max},
						"Decay.start.sqWdt");
	
	min = dfUp.Min("Decay.stop.sqWdt").GetValue();
	max = dfUp.Max("Decay.stop.sqWdt").GetValue();
	TH1D stopWdtUp = *dfUp.Histo1D<double>({"Stop_Width_Up",
						"Stop Width;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min, max},
						"Decay.stop.sqWdt");
	
	vector<TH1D> wdth = {startWdtDwn,stopWdtDwn,startWdtUp,stopWdtUp};
	for( const auto & hist : wdth )
		hist.Write("",TObject::kOverwrite);
}


void WriteStopSignalFallTimeHistos( RNode dfUp, RNode dfDwn, TFile* f ){
	double wdt = 1.5e-9;
	double min=0, max=0;
	
	min = dfDwn.Min("Decay.stop.sqFall").GetValue();
	max = dfDwn.Max("Decay.stop.sqFall").GetValue();
	TH1D stopDwn = *dfDwn.Histo1D<double>({"Dwn_Stop_Time",
						"Stop Time Dwn Decays;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min-0.1*(max-min), max+0.1*(max-min)},
						"Decay.stop.sqFall");
	
	min = dfUp.Min("Decay.stop.sqFall").GetValue();
	max = dfUp.Max("Decay.stop.sqFall").GetValue();
	TH1D stopUp = *dfUp.Histo1D<double>({"Up_Stop_Time",
						"Stop Width Up Decays;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min-0.1*(max-min), max+0.1*(max-min)},
						"Decay.stop.sqFall");
	
	vector<TH1D> wdth = {stopDwn,stopUp};
	for( const auto & hist : wdth )
		hist.Write("",TObject::kOverwrite);
}

vector<TNamed> GetEffString(ROOT::RDF::RCutFlowReport CutReport){
	vector<ROOT::RDF::TCutInfo> CutInfos(CutReport.begin(),CutReport.end());
	vector<TNamed> results;
	const auto allEntries = CutInfos.empty() ? 0ULL : CutInfos.begin()->GetAll();
	for (auto &&ci : CutInfos) {
		const auto &name = ci.GetName();
		const auto pass = ci.GetPass();
		const auto all = ci.GetAll();
		const auto eff = ci.GetEff();
		const auto cumulativeEff = 100.f * float(pass) / float(allEntries);
		string effStr;
		effStr += name + " pass = "+ pass +"all = " + all +  "-- eff = " + eff + "cumulative eff = " + cumulativeEff+" ";
		results.push_back(TNamed(name,effStr));
	}
	return results;
}

map<const char*,bool> GetRunBasedStopCuts( short int run, double sqFall ){
	bool up=false, dwn=false;
	switch (run){
		case 1:
		up =  (sqFall<9.9e-6)    || (sqFall>10.4e-6 && sqFall<10.62e-6);
		dwn = up;
		break;
		case 3:
		up =  (sqFall<9.9e-6)    || (sqFall>10.4e-6 && sqFall<10.62e-6);
		dwn = up;
		break;
		case 4:
		up =  (sqFall<8.75e-6) || (sqFall>8.9e-6 && sqFall<9.185e-6) || (sqFall>11.1e-6 && sqFall<11.47e-6);
		dwn = up;
		break;
		case 5:
		up = sqFall<12.052e-6;
		dwn = up;
		break;
		case 6:
		up = sqFall<85.89e-6;
		dwn = up;
		break;
		case 7:
		up = sqFall<21.3765e-6;
		dwn = up;
		break;
		case 8:
		up = (sqFall<17.36e-6) || (sqFall>19.4e-6 && sqFall<19.65e-6) || (sqFall>21.1e-6 && sqFall<21.37e-6);
		dwn = up;
		break;
		case 9:
		up =  (sqFall<12.05e-6) || (sqFall>16.8e-6 && sqFall<17.07e-6) || (sqFall>21e-6 && sqFall<21.374e-6);
		dwn = up;
		break;
		case 10:
		up =  (sqFall<14.21e-6) || (sqFall>17e-6 && sqFall<17.364e-6) || (sqFall>19e-6 && sqFall<19.65e-6);
		dwn = up;
		break;
		default:
		cout<<"Warning: sqFall not implemented for run "<<run<<endl;
		dwn = true;
		up = true;
	}
	return {{"dwn",!dwn}, {"up",!up}, {"antidwn",dwn}, {"antiup",up}};
}

map<const char*, RNode> ApplyRunBasedStopCuts(TFile* f, RNode dfUp, RNode dfDwn){
	auto dwnDecays=dfDwn.Filter([](short int run, double sqFall)->bool{
					return GetRunBasedStopCuts(run,sqFall)["dwn"];}, 
					{"runN","Decay.stop.sqFall"} , "StopTimeFilterDwnDecay");
	auto upDecays = dfUp.Filter([](short int run, double sqFall)->bool{
					return GetRunBasedStopCuts(run,sqFall)["up"];}, 
					{"runN","Decay.stop.sqFall"}  , "StopTimeFilterUpDecay ");
	auto dwnDecaysBad=dfDwn.Filter([](short int run, double sqFall)->bool{
					return GetRunBasedStopCuts(run,sqFall)["antidwn"];}, 
					{"runN","Decay.stop.sqFall"} , "StopTimeFilterDwnDecay_BadEvents");
	auto upDecaysBad = dfUp.Filter([](short int run, double sqFall)->bool{
					return GetRunBasedStopCuts(run,sqFall)["antiup"];}, 
					{"runN","Decay.stop.sqFall"}  , "StopTimeFilterUpDecay_BadEvents");
	auto upReport = upDecays.Report();
	auto dwnReport = dwnDecays.Report();
	dwnReport->Print();
	upReport->Print();
	dwnDecaysBad.Report()->Print();
	upDecaysBad.Report()->Print();
	auto upEffStr = GetEffString( (*upReport) );
	auto dwnEffStr = GetEffString( (*dwnReport) );
	for( const auto & name : upEffStr )
		name.Write("",TObject::kOverwrite);
	for( const auto & name : dwnEffStr )
		name.Write("",TObject::kOverwrite);
	return {make_pair<const char*, RNode>("dwn",dwnDecays),
		make_pair<const char*, RNode>("up",upDecays),
		make_pair<const char*, RNode>("upbad",upDecaysBad),
		make_pair<const char*, RNode>("dwnbad",dwnDecaysBad)};
}

void WriteDecayHistos(RNode dfUp, RNode dfDwn, TFile* RootOut, double binWdt, TString name_prefix){
	const char* var = "Decay.dtFall";
	
	double max = dfUp.Max(var).GetValue();
	TH1D hup = *dfUp.Histo1D<double>(	{name_prefix+"_UpDecay",name_prefix+"UpDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,5)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	hup.Write("",TObject::kOverwrite);
	
	max = dfDwn.Max(var).GetValue();
	TH1D hdwn = *dfDwn.Histo1D<double>(	{name_prefix+"_DwnDecay",name_prefix+"DwnDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,5)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	hdwn.Write("",TObject::kOverwrite);
}

map<const char*, TH1D> GetDecayHistos(RNode df, double binWdt, TString name_prefix){
	const char* var = "dt";
	
	double max = df.Max(var).GetValue();
	TH1D htot = *df.Histo1D<double>(	{name_prefix+"Decay",name_prefix+"Decay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	
	max = df.Filter("topology==1").Max(var).GetValue();
	TH1D hup = *df.Filter("topology==1")
		      .Histo1D<double>(	{name_prefix+"UpDecay",name_prefix+"UpDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	
	max = df.Filter("topology==0").Max(var).GetValue();
	TH1D hdwn = *df.Filter("topology==0")
		       .Histo1D<double>( {name_prefix+"DwnDecay",name_prefix+"DwnDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);

	return {{"up",hup},{"dwn",hdwn},{"htot",htot}};
}

void WriteDecayTree(RNode dfUp, RNode dfDwn, TFile* f){
	std::string_view var = "Decay.dtFall";
	auto ups = dfUp.Take<double>(var);
	auto dwns= dfDwn.Take<double>(var);
	double dt;
	bool topology; 
	TTree tree("decays","decays");
	tree.Branch("dt",&dt,"dt/D");
	tree.Branch("topology",&topology,"topology/O");
	for( const auto & updt : ups ){
		dt = updt;
		topology = true;
		tree.Fill();
	}
	for( const auto & dwndt : dwns ){
		dt = dwndt;
		topology = false;
		tree.Fill();
	}
	tree.Write();
}

class AmuletFitCore {
	public:
		AmuletFitCore(TH1D h, double rmin, double rmax, TString fitopt):
		_rmin(rmin), _rmax(rmax)
		{
			if( !fitopt.Contains("S") )
				fitopt+"S";
			_fitopt = fitopt;
			_h = h;
			_status = false;
		};
		~AmuletFitCore(){ delete _func; };
		TFitResult AmuletFit(string myformula){
			auto func = new TF1("decay_law",myformula.c_str(),_rmin,_rmax,"NL");
			auto tau_idx = func->GetParNumber("#tau");
			func->SetParameter(tau_idx, 2.197e-6);
			func->SetParameter("N",_h.GetBinContent(_h.FindBin(_rmin))*2.197e-6);
			func->SetParameter("B",_h.GetBinContent(_h.FindBin(_rmax)));
			TFitResultPtr fitres = _h.Fit(func,_fitopt);	
			_status = fitres->IsValid() && fitres->Status()==0 && fitres->HasMinosError(tau_idx);
			_func = func;
			return *fitres;
		};
		inline const bool GetStatus(){return _status;};
		inline const short int GetParIdx(const char* parname){ return (_func) ?  _func->GetParNumber(parname) : -1 ;}
	private:
		TH1D _h;
		Option_t *_fitopt;
		double _rmin, _rmax;
		bool _status;
		TF1* _func;
};

void WriteBinNumberStabilityFixedRange(RNode decaydf, double rmin, double rmax, int N_iter, TString RootOut){
	double time_res_from_run11 = 2.179e-9;
	double binWdtMin = time_res_from_run11/2;
	double binWdtMax = 1.5e-6;
	double step = (binWdtMax-binWdtMin)/static_cast<long double>(N_iter);
	auto gtau  = TGraphAsymmErrors(N_iter);
	auto gchi2 = TGraph(N_iter);
	TString name="BinWidthStability";
	TString title="Bin width stability; Bin Width [s]; #tau_{#mu} [s]";
	gtau.SetNameTitle(name+"_tau",title);
	gchi2.SetNameTitle(name+"_chi2","#chi2 bin width stability; Bin Width [s]; #frac{#chi2}{NDF}");
	for( int i = 0; i<N_iter; i++ ){
		double binWdt = binWdtMin + static_cast<long double>(i)*step;
		int NbinsInFitRange = (rmax-rmin)/binWdt;
		cout<<"iter: "<<i<<" of "<<N_iter<<"\twdt: "<<binWdt<<"\tN bins in fit range: "<<NbinsInFitRange<<endl;
		auto hs = GetDecayHistos(decaydf, binWdt, to_string(binWdt));
		//auto hup = hs["up"], hdwn = hs["dwn"];
		auto htot = hs["htot"];
		//string myformula ="[N]*exp(-x/[#tau])+[b]";
		//string myformula ="([N]/[#tau])*exp(-x/[#tau])+[b]";
		string myformula ="[N]*ROOT::Math::exponential_pdf(x,1./[#tau])+[B]";
		//fitres->Print("V");
		auto fitcore = AmuletFitCore(htot,rmin,rmax,"LEQRS");
		TFitResult fitres = fitcore.AmuletFit(myformula);
		auto tau_idx = fitcore.GetParIdx("#tau");
		if( fitcore.GetStatus() ){
			gtau.SetPoint(i, binWdt, fitres.Parameter(tau_idx));
			gtau.SetPointError(i, 0, 0, abs(fitres.LowerError(tau_idx)), fitres.UpperError(tau_idx));
			gchi2.SetPoint(i, binWdt, fitres.Chi2()/fitres.Ndf());
		}else{
			cout<<"Fit for bin wdt = "<<binWdt<<" has failed, "<<
			"status: "<<fitres.Status()<<" valid: "<<fitres.IsValid()<<endl;
		}
	}
	TFile f(RootOut,"UPDATE");
	if(f.IsZombie()){
		cout<<"Failed to open file "<<RootOut<<endl;
	}else{
		gtau.Write("",TObject::kOverwrite);
		gchi2.Write("",TObject::kOverwrite);
		f.Close();
	}
}

int main(int argc, char** argv)
{
	//set names from input line
	if ( !(argc >= 3) ){ 
		std::cout<<"\nINPUT ERROR"<<endl
		<<"needs:\n./path/to/executable/ path/to/output/root/file.root path/to/input/root/file.root [...add as many input files as you wish]"<<endl
		<<"./executable ../DAQresults/muLifetimeC.root ../DAQprocessed/run1meas_ALL.root ../DAQprocessed/run5meas_ALL.root ../DAQprocessed/run7meas_ALL.root"<<endl
		<<"./executable ../DAQprocessed/muLifetimeAl.root ../DAQprocessed/run3meas_ALL.root"<<endl;
		return 1;
	}

	//declare the struct to the CLING compiler
	gInterpreter->Declare("typedef struct {double sqFall=0, sqFallErr=0; double sqRise=0, sqRiseErr=0; double sqWdt = sqRise - sqFall; double sqWdtErr = TMath::Sqrt(sqFall*sqFall + sqRise*sqRise);} Square_Signal;");
	gInterpreter->Declare("typedef struct {Square_Signal start; Square_Signal stop; double dtFall = stop.sqFall - start.sqFall; double dtFallErr = TMath::Sqrt(stop.sqFallErr*stop.sqFallErr + start.sqFallErr*start.sqFallErr); double dtRise = stop.sqRise - start.sqRise; double dtRiseErr = TMath::Sqrt(stop.sqRiseErr*stop.sqRiseErr + start.sqRiseErr*start.sqRiseErr); } Decay_Event;");
	gSystem->Load("libMathMore");
	
	//use all cores of your machine MultiThreading to speed up analysis
	//(read RDataFrame documentation for MT warnings)
	ROOT::EnableImplicitMT();
		
	//setting I\O file names
	cout<<"I will process the following files:"<<endl;
	TString RootOut = argv[1];
	vector<TString> files;
	for(int i = 2; i<argc; i++)
		files.push_back(argv[i]);
	TChain chainup("UpDecays");
	TChain chaindwn("DwnDecays");
	for(const auto & file : files){
		cout<<file<<endl;
		chainup.Add(file);
		chaindwn.Add(file);
	}
	
	//create df
	ROOT::RDataFrame dfUp(chainup);
	ROOT::RDataFrame dfDwn(chaindwn);

	TFile* outFile = TFile::Open(RootOut,"RECREATE");
	WriteWidthHistos(dfUp, dfDwn, outFile);
	WriteStopSignalFallTimeHistos(dfUp, dfDwn, outFile);
	double binWdt = 2.179e-9*2;
	WriteDecayHistos(dfUp, dfDwn, outFile, binWdt, "");
	auto res = ApplyRunBasedStopCuts(outFile, dfUp, dfDwn);
	WriteDecayHistos(res.at("upbad"), res.at("dwnbad"), outFile, binWdt, "bad");
	WriteDecayHistos(res.at("up"), res.at("dwn"), outFile, binWdt, "filtered");
	WriteDecayTree(dfUp, dfDwn, outFile);
	outFile->Close();
	RDataFrame decaydf("decays",RootOut);
	double rmin = 0.6e-6, rmax = 9.1e-6;
	int N_iters = 10;	
	WriteBinNumberStabilityFixedRange(decaydf, rmin, rmax, N_iters, RootOut);


	cout<<"File "<<RootOut<<" RECREATED"<<endl<<endl;
	
	return 0;
}
/* idea da provare per funzione: \TODO capire bene se la normalizzazione e corretta
		std::ostringstream rminstr, rmaxstr, nbinsstr;
		rminstr<<rmin;
		rmaxstr<<rmax;
		nbinsstr<<NbinsInFitRange;
		string myformula = "(2./("+nbinsstr.str()+"*[#tau]*(exp(-"+rminstr.str()+"/[#tau])-exp(-"+rmaxstr.str()+"/[#tau]))))*exp(-x/[#tau])+[B]";
*/
