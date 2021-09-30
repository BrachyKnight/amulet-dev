/* 
 * AMULET: Analysis MUon LifETime
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
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Math/PdfFunc.h> //exponential_pdf() //uniform_pdf()

using ROOT::RDF::RNode, ROOT::RDataFrame, ROOT::VecOps::RVec, ROOT::Math::exponential_pdf, ROOT::Math::uniform_pdf;
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

enum class RunConfiguration : short int {
	kRunNotPresent = 0,
	kCarbon = 1,
	kAl = 2,
	kNaCl = 3,
	kBkg = 4,
	kBoff = 5,
	kBon = 6,
};

long double MEASTIME;

double NormExpRange( double t, double lambda, double rmin, double rmax){
	return exponential_pdf(t,lambda)/(exp(-lambda*rmin)-exp(-lambda*rmax));
}

//core class for AMULET: Analysis MUon LifETime.
//use AmuletFitCore in order to fit histogram h with one of the available functions
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
			_binWdt = h.GetXaxis()->GetBinWidth(h.FindBin(rmin+(rmax-rmin)/2.));
		};
		~AmuletFitCore(){ delete _func; };
		enum class FuncType : short int{
			kNormExpB = 1,
			kRoberto1 = 2,
			kRoberto2 = 3,
			kConst = 4,
			kUniform = 5,
			kNormExpNormB = 6,
			kMassimo1 = 7,
			kDefault = 1,
		};
		TFitResult AmuletFit(FuncType fType){
			ChooseFormula(fType);
			auto tau_idx = _func->GetParNumber("#tau");
			_func->SetParameter(tau_idx, 2.197e-6);
			TFitResultPtr fitres = _h.Fit(_func,_fitopt);	
			_status = fitres->IsValid() && fitres->Status()==0 && fitres->HasMinosError(tau_idx);
			return *fitres;
		};
		inline const bool GetStatus(){return _status;};
		inline const short int GetParIdx(const char* parname){ return (_func) ?  _func->GetParNumber(parname) : -1 ;}
		inline TH1D* GetHistoClone(const char* newname=""){return (TH1D*)(_h.Clone(newname));};
	private:
		TH1D _h;
		Option_t *_fitopt;
		double _rmin, _rmax, _binWdt;
		bool _status;
		TF1* _func;

		void ChooseFormula( FuncType fType ){
			const long double rmi = _rmin, rma = _rmax, binw = _binWdt;
			const long double meastime = MEASTIME;
			const long double nbins = (long double)static_cast<int>(((rma-rmi)/binw));
			const long double nentr = (long double)static_cast<long long int>((_h.Integral(_h.FindBin(_rmin),_h.FindBin(_rmax))));
			//list of the available functions
			switch (fType) {
				case FuncType::kNormExpB: //exponential pdf + const bkg
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau = p[0];
								double B = p[1];
								double N = p[2];
								double t = x[0];
								double lambda = 1./tau;
								double signal = N*exponential_pdf(t,lambda);
								double bkg = B;
								return signal + bkg; 
							}, rmi, rma, 3 );
					_func->SetParNames("#tau","b","N");
					_func->SetParameter("b",_h.GetBinContent(_h.FindBin(_rmax)));
					_func->SetParameter("N",_h.GetBinContent(_h.FindBin(_rmin))*(2.197e-6));
				}
				break;
				case FuncType::kNormExpNormB: //exponential pdf + uniform pdf (const bkg)
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau = p[0];
								double B = p[1];
								double N = p[2];
								double t = x[0];
								double lambda = 1./tau;
								double signal = N*exponential_pdf(t,lambda);
								double bkg = B;
								return signal + bkg; 
							}, rmi, rma, 3, "NL" );
					_func->SetParNames("#tau","B","N");
					_func->SetParameter("B",_h.GetBinContent(_h.FindBin(_rmax))*(_rmax-_rmin));
					_func->SetParameter("N",_h.GetBinContent(_h.FindBin(_rmin))*(2.197e-6));
				}
				break;
				case FuncType::kRoberto1: //normalized pdf (two parameters) [B represents N_bkg/nBins]
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau = p[0];
								double B = p[1];
								double t = x[0];
								double lambda = 1./tau;
								double N = ((nentr-B*nbins)*binw)/(exp(-rmi*lambda)-exp(-rma*lambda));
								double signal = N*exponential_pdf(t,lambda); 
								double bkg = B;
								return (t<=rma && t>=rmi) ? signal+bkg : 0.;
							}, _rmin, _rmax, 2, "NL" );
					_func->SetParNames("#tau","b");
					_func->SetParameter("b",_h.GetBinContent(_h.FindBin(_rmax)));
				}
				break;
				case FuncType::kRoberto2: //normalized pdf (two parameters) [B represents N_bkg in fit range]
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau = p[0];
								double B = p[1]/nbins;
								double t = x[0];
								double lambda = 1./tau;
								double N = ((nentr-B*nbins)*binw)/(exp(-rmi*lambda)-exp(-rma*lambda));
								double signal = N*exponential_pdf(t,lambda); 
								double bkg = B;
								return (t<=rma && t>=rmi) ? signal+bkg : 0.;
							}, _rmin, _rmax, 2, "NL" );
					_func->SetParNames("#tau","B");
					_func->SetParameter("B",_h.GetBinContent(_h.FindBin(_rmax))*nbins);
				}
				break;
				case FuncType::kMassimo1:
				{
					_func = new TF1("decay_law",
							[meastime,rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau = p[0];
								double rate_cc = p[1];
								double t = x[0];
								double lambda = 1./tau;
								double A = ((rma-rmi)*meastime)/nbins;
								double rate_tot = nentr/meastime;
								double signal = NormExpRange(t, lambda, rmi, rma); 
								double bkg = uniform_pdf(t,rmi,rma);
								double pdf = A*( (rate_tot-rate_cc)*signal + rate_cc*bkg  );
								return (t<=rma && t>=rmi) ? pdf : 0.;
							}, _rmin, _rmax, 2, "NL" );
					_func->SetParNames("#tau","rcc");
					_func->SetParameter("rcc",0.004);
				}
				break;
				case FuncType::kConst: //normalized pdf (two parameters)
				{	
					_func = new TF1("decay_law",
							[rmi,rma](double*x, double *p)->double{	
								double B = p[1];
								double t = x[0];
								return (t<rma && t>rmi) ? B : 0.;
							}, _rmin, _rmax, 1, "NL" );
					_func->SetParameter("b",_h.GetBinContent(_h.FindBin(_rmax)));
				}
				break;
				case FuncType::kUniform: //normalized pdf (two parameters)
				{
					_func = new TF1("decay_law",
							[rmi,rma](double*x, double *p)->double{	
								double B = p[1];
								double t = x[0];
								return B*uniform_pdf(t,rmi,rma);
							}, _rmin, _rmax, 1, "NL" );
					_func->SetParameter("B",_h.GetBinContent(_h.FindBin(_rmax)));
				}
				break;
			}
			//if( !((TString)_fitopt).Contains("Q") ) cout<<"using: "<<formula<<endl;
		};
};

void WriteWidthHistos( RNode dfUp, RNode dfDwn, TFile* f ){
	double wdt = 0.5e-9;
	double min=0, max=0;
	
	min = dfDwn.Min<double>("Decay.start.sqWdt").GetValue();
	max = dfDwn.Max<double>("Decay.start.sqWdt").GetValue();
	TH1D startWdtDwn = *dfDwn.Histo1D<double>({"Start_Width_Dwn",
						"Start Width;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min, max},
						"Decay.start.sqWdt");
	
	min = dfDwn.Min<double>("Decay.stop.sqWdt").GetValue();
	max = dfDwn.Max<double>("Decay.stop.sqWdt").GetValue();
	TH1D stopWdtDwn = *dfDwn.Histo1D<double>({"Stop_Width_Dwn",
						"Stop Width;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min, max},
						"Decay.stop.sqWdt");
	
	min = dfUp.Min<double>("Decay.start.sqWdt").GetValue();
	max = dfUp.Max<double>("Decay.start.sqWdt").GetValue();
	TH1D startWdtUp = *dfUp.Histo1D<double>({"Start_Width_Up",
						"Start Width;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min, max},
						"Decay.start.sqWdt");
	
	min = dfUp.Min<double>("Decay.stop.sqWdt").GetValue();
	max = dfUp.Max<double>("Decay.stop.sqWdt").GetValue();
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
	
	min = dfDwn.Min<double>("Decay.stop.sqFall").GetValue();
	max = dfDwn.Max<double>("Decay.stop.sqFall").GetValue();
	TH1D stopDwn = *dfDwn.Histo1D<double>({"Dwn_Stop_Time",
						"Stop Time Dwn Decays;time [s];counts per "+TString(to_string(wdt*1e9).substr(0,5))+" ns",
						static_cast<int>((max-min)/wdt), min-0.1*(max-min), max+0.1*(max-min)},
						"Decay.stop.sqFall");
	
	min = dfUp.Min<double>("Decay.stop.sqFall").GetValue();
	max = dfUp.Max<double>("Decay.stop.sqFall").GetValue();
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
	
	double max = dfUp.Max<double>(var).GetValue();
	TH1D hup = *dfUp.Histo1D<double>(	{name_prefix+"_UpDecay",name_prefix+"UpDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,5)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	hup.Write("",TObject::kOverwrite);
	
	max = dfDwn.Max<double>(var).GetValue();
	TH1D hdwn = *dfDwn.Histo1D<double>(	{name_prefix+"_DwnDecay",name_prefix+"DwnDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,5)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	hdwn.Write("",TObject::kOverwrite);
}

map<const char*, TH1D> GetDecayHistos(RNode df, double binWdt, TString name_prefix){
	const char* var = "dt";
	
	double max = df.Max<double>(var).GetValue();
	TH1D htot = *df.Histo1D<double>(	{name_prefix+"Decay",name_prefix+"Decay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	
	max = df.Filter("topology==1").Max<double>(var).GetValue();
	TH1D hup = *df.Filter("topology==1")
		      .Histo1D<double>(	{name_prefix+"UpDecay",name_prefix+"UpDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	
	max = df.Filter("topology==0").Max<double>(var).GetValue();
	TH1D hdwn = *df.Filter("topology==0")
		       .Histo1D<double>( {name_prefix+"DwnDecay",name_prefix+"DwnDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);

	return {{"up",hup},{"dwn",hdwn},{"htot",htot}};
}

void WriteDecayTree(RNode dfUp, RNode dfDwn, TFile* f){
	std::string_view var = "Decay.dtFall";
	auto ups = dfUp.Take<double>(var);
	auto dwns= dfDwn.Take<double>(var);
	std::string_view startWdtVar = "Decay.start.sqWdt";
	std::string_view stopWdtVar = "Decay.stop.sqWdt";
	vector<double> up_start_wdt = *dfUp.Take<double>(startWdtVar);
	vector<double> dw_start_wdt = *dfDwn.Take<double>(startWdtVar);
	vector<double> up_stop_wdt = *dfUp.Take<double>(stopWdtVar);
	vector<double> dw_stop_wdt = *dfDwn.Take<double>(stopWdtVar);
	vector<short int> upMeasN = *dfUp.Take<short int>("measN");
	vector<short int> dwMeasN = *dfDwn.Take<short int>("measN");
	vector<unsigned long long> upidx = *dfUp.Take<unsigned long long>("idx");
	vector<unsigned long long> dwidx = *dfDwn.Take<unsigned long long>("idx");
	vector<short int> upRunN = *dfUp.Take<short int>("runN");
	vector<short int> dwRunN = *dfDwn.Take<short int>("runN");
	double dt, startWdt, stopWdt;
	short int runN, measN;
       	unsigned long long idx;
	bool topology; 
	TTree tree("decays","decays");
	tree.Branch("dt",&dt,"dt/D");
	tree.Branch("topology",&topology,"topology/O");
	tree.Branch("startWdt",&startWdt,"startWdt/D");
	tree.Branch("stopWdt", &stopWdt, "stopWdt/D" );
	tree.Branch("measN",&measN,"measN/B");
	tree.Branch("idx", &idx, "idx/l" );
	tree.Branch("runN",&runN, "runN/B");
	{ 	long int i = 0;
		for( const auto & updt : ups ){
			dt = updt;
			startWdt = up_start_wdt[i];
			stopWdt = up_stop_wdt[i];
			measN = upMeasN[i];
			idx = upidx[i];
			runN = upRunN[i];
			topology = true;
			i++;
			tree.Fill();
		}
	}
	{	long int i = 0;
		for( const auto & dwndt : dwns ){
			dt = dwndt;
			startWdt = dw_start_wdt[i];
			stopWdt = dw_stop_wdt[i];
			measN = dwMeasN[i];
			idx = dwidx[i];
			runN = dwRunN[i];
			topology = false;
			i++;
			tree.Fill();
		}
	}
	tree.Write("",TObject::kOverwrite);
}

//generate a list of n exponentially spaced numbers between first and last
vector<long double> ExpList(long double first, long double last, unsigned long int n){
    vector<long double> vector(n); 
    long double m = (long double) 1 / (n - 1);
    long double quotient = pow(last / first, m);

    vector[0] = first;

    for (unsigned long int i = 1; i < n; i++)
        vector[i] = vector[i - 1] * quotient;

    return vector;
}

void WriteBinNumberStabilityFixedRange(RNode decaydf, double rmin, double rmax, int N_iter, TString RootOut, AmuletFitCore::FuncType func_type){
	long double time_res_from_run11 = 2.179e-9;
	long double binWdtMin = time_res_from_run11/10;
	long double binWdtMax = 1e-6*2; //non ha senso mettere piu di cosi
	//long double step = (binWdtMax-binWdtMin)/static_cast<long double>(N_iter);
	auto exp_spaced = ExpList(binWdtMin, binWdtMax, N_iter);
	auto gtau  = TGraphAsymmErrors(N_iter);
	auto gchi2 = TGraph(N_iter);
	TString name="BinWidthStability";
	TString title="Bin width stability;Bin Width [ns];#tau_{#mu} [#mus]";
	gtau.SetNameTitle(name+"_tau",title);
	gchi2.SetNameTitle(name+"_chi2","#chi2 bin width stability;Bin Width [ns];#frac{#chi2}{NDF}");
	for( int i = 0; i<N_iter; i++ ){
		//long double binWdt = binWdtMin + static_cast<long double>(i)*step;
		long double binWdt = exp_spaced[i];
		int NbinsInFitRange = (rmax-rmin)/binWdt;
		cout<<"iter: "<<i<<" of "<<N_iter<<"\twdt: "<<binWdt<<"\tN bins in fit range: "<<NbinsInFitRange<<endl;
		auto hs = GetDecayHistos(decaydf, binWdt, to_string(binWdt));
		//auto hup = hs["up"], hdwn = hs["dwn"];
		auto htot = hs["htot"];
		auto fitcore = AmuletFitCore(htot,rmin,rmax,"LERSN0Q");
		TFitResult fitres = fitcore.AmuletFit(func_type);
		//fitres.Print("V");
		auto tau_idx = fitcore.GetParIdx("#tau");
		if( fitcore.GetStatus() ){
			gtau.SetPoint(i, binWdt*1e9, fitres.Parameter(tau_idx)*1e6);
			gtau.SetPointError(i, 0, 0, abs(fitres.LowerError(tau_idx))*1e6, fitres.UpperError(tau_idx)*1e6);
			gchi2.SetPoint(i, binWdt*1e9, fitres.Chi2()/fitres.Ndf());
		}else{
			cout<<"Fit for bin wdt = "<<binWdt<<" has FAILED, "<<
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

void WriteRangeStabilityFixedBins(RNode decaydf, double binWdt, int NitMin, int NitMax, TString RootOut, AmuletFitCore::FuncType func_type){
	long double stopWdtUp = decaydf.Filter("topology==1").Mean<double>("stopWdt").GetValue();
	long double stopWdtDw = decaydf.Filter("topology==0").Mean<double>("stopWdt").GetValue();
	long double startWdt = decaydf.Mean<double>("startWdt").GetValue();
	long double stopWdt = std::min(stopWdtUp,stopWdtDw);
	long double rmin_min = stopWdt + startWdt - 500e-9;
	long double rmin_max =  11e-6;
	long double rmax_min = rmin_min + 500e-9;
	long double rmax_max =  11.5e-6;
	cout<<"Range stability: "<<endl;
	cout<<"rmin min = "<<rmin_min<<"\t rmin max = "<<rmin_max<<endl;
	cout<<"rmax min = "<<rmax_min<<"\t rmax max = "<<rmax_max<<endl<<endl;
	long double stepMin = (rmin_max-rmin_min)/static_cast<long double>(NitMin);
	long double stepMax = (rmax_max-rmax_min)/static_cast<long double>(NitMax);
	auto gtau  = TGraph2DErrors();
	auto gchi2 = TGraph2D();
	TString name="RangeStability";
	TString title="Range stability;rmin [#mus];rmax [#mus]; #tau_{#mu} [#mus]";
	gtau.SetNameTitle(name+"_tau",title);
	gchi2.SetNameTitle(name+"_chi2","#chi2 range stability;rmin [#mus];rmax [#mus]; #frac{#chi2}{NDF}");
	int Npoint = 0;
	for( int i = 0; i<NitMin; i++ ){
		long double rmin = rmin_min + static_cast<long double>(i)*stepMin;
		for( int j = 0; j<NitMax; j++ ){
			long double rmax = rmax_min + static_cast<long double>(j)*stepMax;
			int NbinsInRange = abs(rmax-rmin)/binWdt;
			if( rmax < rmin || NbinsInRange<6 )
				continue;
			cout<<"iter: ("<<i<<","<<j<<") of ("<<NitMin<<","<<NitMax<<")"<<"\trmin: "<<rmin<<"\trmax: "<<rmax
			    <<"\tN bins in fit range: "<<NbinsInRange<<endl;
			auto hs = GetDecayHistos(decaydf, binWdt, "rmin"+to_string(rmin*1e-6)+"rmax"+to_string(rmax*1e-6)+"#mus");
			//auto hup = hs["up"], hdwn = hs["dwn"];
			auto htot = hs["htot"];
			auto fitcore = AmuletFitCore(htot,rmin,rmax,"LERSN0Q");
			TFitResult fitres = fitcore.AmuletFit(func_type);
			//fitres.Print("V");
			auto tau_idx = fitcore.GetParIdx("#tau");
			if( fitcore.GetStatus() ){
				gtau.AddPoint(rmin*1e6, rmax*1e6, fitres.Parameter(tau_idx)*1e6);
				Npoint++;
				gtau.SetPointError(Npoint, 0, 0, fitres.Error(tau_idx)*1e6);
				gchi2.AddPoint(rmin*1e6, rmax*1e6, fitres.Chi2()/fitres.Ndf());
			}else{
				cout<<"Fit for rmin = "<<rmin<<" and rmax = "<<rmax<<" has FAILED, "<<
				"status: "<<fitres.Status()<<" valid: "<<fitres.IsValid()<<endl;
			}
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

//stability analysis: opt "B" for bin width fit stab, opt: "R" for range stab.
void StabilityAnalysis(RNode decaydf, Option_t *opt, TString RootOut, double rmin, double rmax, double binWdt, AmuletFitCore::FuncType func_type){
	if( TString(opt).Contains("B")){
		int N_iters = 100;
		WriteBinNumberStabilityFixedRange(decaydf, rmin, rmax, N_iters, RootOut, func_type);
	}	
	if( TString(opt).Contains("R") ){
		int N_iters_min = 100;
		int N_iters_max = 100;
		WriteRangeStabilityFixedBins( decaydf, binWdt, N_iters_min, N_iters_max, RootOut, func_type);

	}
}

void WriteLifetimeFit(RNode decaydf, TString RootOut, double binWdt, double rmin, double rmax, AmuletFitCore::FuncType func_type, const char* name_prefix = ""){
		gStyle->SetOptFit(1111);
		auto hs = GetDecayHistos(decaydf, binWdt, name_prefix);
		//auto hup = hs["up"], hdwn = hs["dwn"];
		auto htot = hs["htot"];
		auto fitcore = AmuletFitCore(htot,rmin,rmax,"LERS");
		TFitResult fitres = fitcore.AmuletFit(func_type);
		fitres.Print("V");
		auto f = TFile(RootOut,"UPDATE");
		auto c = TCanvas("decay_fit","decay_fit",500,500);
		c.cd();
		fitcore.GetHistoClone()->Draw();
		c.Write("",TObject::kOverwrite);
		f.Close();
}

RunConfiguration ChooseConfig(TString RootOut){
	if      ( RootOut.Contains("muLifetimeC.root")    )
	       return RunConfiguration::kCarbon;
	else if ( RootOut.Contains("muLifetimeAl.root")   )	
		return RunConfiguration::kAl;
	else if ( RootOut.Contains("muLifetimeNaCl.root") )
		return RunConfiguration::kNaCl;
	else if ( RootOut.Contains("muLifetimeBKG.root")  )
		return RunConfiguration::kBkg;
	else if ( RootOut.Contains("muLifetimeBon.root")  )
		return RunConfiguration::kBon;
	else if ( RootOut.Contains("muLifetimeBoff")      )
		return RunConfiguration::kBoff;
	else
		return RunConfiguration::kRunNotPresent;
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
	
	double binWdt = 0, rmin = 0, rmax = 0;
	AmuletFitCore::FuncType fType = AmuletFitCore::FuncType::kDefault;
	auto measconf = ChooseConfig(RootOut);
	switch( measconf ){
		case RunConfiguration::kCarbon:
			binWdt = 20e-9;  //opt
			rmin = 500e-9;   //opt?
			rmax = 10200e-9; //opt?
			MEASTIME = 3.0464e06; //sigma?
			fType = AmuletFitCore::FuncType::kNormExpB;
		break;	
		case RunConfiguration::kAl:
			binWdt = 40e-9; //opt?
			rmin = 0.5e-6;
			rmax = 9.1e-6;
			fType = AmuletFitCore::FuncType::kDefault;
		break;	
		case RunConfiguration::kNaCl:
			binWdt = 30e-9; //opt??
			rmin = 0.5e-6;
			rmax = 9.1e-6;
			fType = AmuletFitCore::FuncType::kDefault;
		break;	
		case RunConfiguration::kBkg:
			binWdt = 1.2e-7; //quella che sceglie in automatico il ttree
			rmin = 30e-6; //scelta guardando plot
			rmax = 38e-6; //scelta guardando plot
			MEASTIME = 1194420.;
			fType = AmuletFitCore::FuncType::kConst;
		break;	
		case RunConfiguration::kBon:
			binWdt = 60e-9; //opt?
			rmin = 0.5e-6;
			rmax = 9.1e-6;
			fType = AmuletFitCore::FuncType::kDefault;
		break;	
		case RunConfiguration::kBoff:
			binWdt = 60e-9; //opt?
			rmin = 0.5e-6;
			rmax = 9.1e-6;
			fType = AmuletFitCore::FuncType::kDefault;
		break;
		case RunConfiguration::kRunNotPresent:
			binWdt = 0;
			rmin = 0;
			rmax = 0;
			throw std::runtime_error("Configuration for "+RootOut+" not present");
		break;
	}
	
	//create df
	ROOT::RDataFrame dfUp(chainup);
	ROOT::RDataFrame dfDwn(chaindwn);

	TFile* outFile = TFile::Open(RootOut,"UPDATE");
	WriteWidthHistos(dfUp, dfDwn, outFile);
	WriteStopSignalFallTimeHistos(dfUp, dfDwn, outFile);
	WriteDecayHistos(dfUp, dfDwn, outFile, binWdt, "");
	auto res = ApplyRunBasedStopCuts(outFile, dfUp, dfDwn);
	WriteDecayHistos(res.at("upbad"), res.at("dwnbad"), outFile, binWdt, "bad");
	WriteDecayHistos(res.at("up"), res.at("dwn"), outFile, binWdt, "filtered");
	WriteDecayTree(res.at("up"), res.at("dwn"), outFile);
	outFile->Close();
	RDataFrame decaydf("decays",RootOut);

	
	if(measconf != RunConfiguration::kBkg)
		StabilityAnalysis(decaydf, "B", RootOut, rmin, rmax, binWdt, fType );

	WriteLifetimeFit(decaydf, RootOut, binWdt, rmin, rmax, fType);
	
	cout<<"File "<<RootOut<<" UPDATED"<<endl<<endl;
	gROOT->EndOfProcessCleanups();
	cout<<"term"<<endl;	
	
	return 0;
}
