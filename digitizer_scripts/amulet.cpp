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
#include <TAxis.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TGraph2DErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <Math/PdfFunc.h> //exponential_pdf() //uniform_pdf()

using ROOT::RDF::RNode, ROOT::RDataFrame, ROOT::VecOps::RVec, ROOT::Math::exponential_pdf;
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
long double TAUMATERIAL;
long double TAUP = 2.19703421e-6; 
long double TAU_carbon = 2025e-9;
long double TAU_Al = 864.6e-9; //880e-9??????
long double TAU_NaCl = 696e-9;
double CHARGERATIORANGECORRECTION = 1.;

double NormExpRange( double x, double lambda, double rmin, double rmax){
	return exponential_pdf(x,lambda)/(exp(-lambda*rmin)-exp(-lambda*rmax));
}
double uniform_pdf( double x, double rmin, double rmax ){
	return ( ( rmax!=rmin && x<=rmax && x>=rmin ) ? 1.0/(rmax-rmin) : 0.0 );
}

//core class for AMULET: Analysis MUon LifETime.
//use AmuletFitCore in order to fit histogram h with one of the available functions
class AmuletFitCore {
	public:
		AmuletFitCore(TH1D h, double rmin, double rmax, TString fitopt):
		_rmin(rmin), _rmax(rmax)
		{
			if( !fitopt.Contains("S") )
				fitopt+="S";
			_fitopt = fitopt;
			_h = h;
			_status = false;
			_binWdt = h.GetXaxis()->GetBinWidth(h.FindBin(rmin+(rmax-rmin)/2.));
		};
		~AmuletFitCore(){ delete _func; };
		enum class FuncType : short int{
			kNormExpB = 1,
			kConst = 2,
			kUniform = 3,
			kNormExpNormB = 4,
			kMassimo = 5,
			kMassimoFreeN = 6,
			kMassimoFraction = 7,
			kMaterial = 8,
		};
		TFitResult AmuletFit(FuncType fType){
			ChooseFormula(fType);
			auto tau_idx = _func->GetParNumber("#tau");
			if( tau_idx != -1 )
				_func->SetParameter(tau_idx, 2.197e-6);
			TFitResultPtr fitres = _h.Fit(_func, _fitopt);
			_status = (tau_idx != -1) ? fitres->IsValid() && fitres->Status()==0 && fitres->HasMinosError(tau_idx) : fitres->IsValid() && fitres->Status()==0;
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
			const double rmi = _rmin, rma = _rmax, binw = _binWdt;
			const long double meastime = MEASTIME;
			const double nbins = static_cast<int>((rma-rmi)/binw);
			const double nentr = static_cast<double>((_h.Integral(_h.FindBin(_rmin),_h.FindBin(_rmax))));
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
				case FuncType::kMassimo:
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau = p[0];
								double nbkg = p[1];
								double t = x[0];
								double lambda = 1./tau;
								double signal = NormExpRange(t, lambda, rmi, rma); 
								double bkg = uniform_pdf(t,rmi,rma);
								double pdf = binw*( (nentr-nbkg)*signal + nbkg*bkg  );
								return (t<=rma && t>=rmi) ? pdf : 0.;
							}, _rmin, _rmax, 2, "NL" );
					_func->SetParNames("#tau","N_{bkg}");
					_func->SetParameter("N_{bkg}",0.004*meastime);
				}
				break;
				case FuncType::kMassimoFreeN:
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau = p[0];
								double nbkg = p[1];
								double ndecay = p[2];
								double t = x[0];
								double lambda = 1./tau;
								double signal = NormExpRange(t, lambda, rmi, rma); 
								double bkg = uniform_pdf(t,rmi,rma);
								double pdf = binw*( ndecay*signal + nbkg*bkg  );
								return (t<=rma && t>=rmi) ? pdf : 0.;
							}, _rmin, _rmax, 3, "NL" );
					_func->SetParNames("#tau","N_{bkg}","N_{decay}");
					_func->SetParameter("N_{bkg}",0.004*meastime);
					_func->SetParameter("N_{decay}",0.004*meastime);
				}
				break;
				case FuncType::kMassimoFraction:
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{
								double tau_p = p[0], lambda_p = 1./tau_p;
								double tau_m = p[1], lambda_m = 1./tau_m;
								double nbkg = p[2], frac_p = p[3];
								double t = x[0];
								double signal = (1-frac_p)*NormExpRange(t, lambda_m, rmi, rma) + frac_p*NormExpRange(t, lambda_p, rmi, rma);
								double bkg = uniform_pdf(t,rmi,rma);
								double pdf = binw*( (nentr-nbkg)*signal + nbkg*bkg  );
								return (t<=rma && t>=rmi) ? pdf : 0.;
							}, _rmin, _rmax, 4, "NL" );
					_func->SetParNames("#tau_{#mu^{+}}","#tau_{#mu^{-}}","N_{bkg}","f_{#mu^{+}}");
					_func->SetParameter("N_{bkg}",0.004*meastime);
					_func->FixParameter(0, TAUP);
					_func->FixParameter(1, TAUMATERIAL);
					_func->SetParameter(3, 0.5);
				}
				break;
				case FuncType::kMaterial:
				{
					_func = new TF1("decay_law",
							[rmi,rma,binw,nbins,nentr](double*x, double *p)->double{	
								double tau_p  = p[0], lambda_p = 1./tau_p;
								double tau_C  = p[1], lambda_C = 1./tau_C;
								double tau_M  = p[2], lambda_M = 1./tau_M;
								double nbkg   = p[3];
								double f_p    = p[4];
								double f_Mat  = p[5];
								double f_C    = 1.-(f_p+f_Mat);
								double t      = x[0];
								double mu_p   = NormExpRange(t, lambda_p, rmi, rma);
								double mu_C   = NormExpRange(t, lambda_C, rmi, rma);
								double mu_Mat = NormExpRange(t, lambda_M, rmi, rma);
								double signal = f_p*mu_p + f_Mat*mu_Mat + f_C*mu_C;
								double bkg    = uniform_pdf(t,rmi,rma);
								double pdf    = binw*( (nentr-nbkg)*signal + nbkg*bkg  );
								return (t<=rma && t>=rmi) ? pdf : 0.;
							}, _rmin, _rmax, 6, "NL" );
					_func->SetParNames("#tau_{#mu^{+}}","#tau_{#mu^{-}}_{C}","#tau_{#mu^{-}}_{mat}","N_{bkg}","f_{#mu^{+}}","f_{#mu^{-}}_{mat}");
					_func->FixParameter(0, TAUP);
					_func->FixParameter(1, TAU_carbon);
					_func->SetParameter(2, TAUMATERIAL);
					_func->SetParameter(3, 0.004*meastime);
					_func->FixParameter(4, 0.53488 + CHARGERATIORANGECORRECTION);
					_func->SetParameter(5, 0.1);
					//_func->FixParameter(6,0.3);

				}
				break;
				case FuncType::kConst: //normalized pdf (two parameters)
				{	
					_func = new TF1("bkg_law",
							[rmi,rma](double*x, double *p)->double{	
								double B = p[0];
								double t = x[0];
								return (t<=rma && t>=rmi) ? B : 0.;
							}, _rmin, _rmax, 1, "NL" );
					_func->SetParNames("b");
					_func->SetParameter("b",_h.GetBinContent(_h.FindBin(_rmin+(_rmax-_rmin)/2.)));
				}
				break;
				case FuncType::kUniform: //normalized pdf (two parameters)
				{
					_func = new TF1("bkg_law",
							[rmi,rma](double*x, double *p)->double{	
								double B = p[0];
								double t = x[0];
								return (t<=rma && t>=rmi) ? B*uniform_pdf(t,rmi,rma) : 0.;
							}, _rmin, _rmax, 1, "NL" );
					_func->SetParNames("B");
					_func->SetParameter("B",_h.GetBinContent(_h.FindBin(_rmax))*(rma-rmi));
				}
				break;
			}
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
	TH1D htot = *df.Histo1D<double>(	{name_prefix+"Decay",name_prefix+"Decay;Time [s];Events per "+to_string(binWdt*1e06).substr(0,5)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	
	max = df.Filter("topology==1").Max<double>(var).GetValue();
	TH1D hup = *df.Filter("topology==1")
		      .Histo1D<double>(	{name_prefix+"UpDecay",name_prefix+"UpDecay;Time [s];Events per "+to_string(binWdt*1e06).substr(0,5)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	
	max = df.Filter("topology==0").Max<double>(var).GetValue();
	TH1D hdwn = *df.Filter("topology==0")
		       .Histo1D<double>( {name_prefix+"DwnDecay",name_prefix+"DwnDecay;Time [s];Events per "+to_string(binWdt*1e06).substr(0,5)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);

	return {{"up",hup},{"dwn",hdwn},{"htot",htot}};
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

void ExportTxt(RNode decaydf, const char* var, TString RootOut){
	auto v = decaydf.Take<double>(var);
	auto TxtOut = RootOut.ReplaceAll(".root",".txt");
	std::ofstream outxt ((const char*)TxtOut);
	if (outxt.is_open()){
		for ( const auto & d : v )
			outxt<<d<<endl;
	}
       	else cout << "Unable to open file" <<endl;
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
	int runN, measN;
       	unsigned long long idx;
	bool topology; 
	TTree tree("decays","decays");
	tree.Branch("dt",&dt,"dt/D");
	tree.Branch("topology",&topology,"topology/O");
	tree.Branch("startWdt",&startWdt,"startWdt/D");
	tree.Branch("stopWdt", &stopWdt, "stopWdt/D" );
	tree.Branch("measN",&measN,"measN/I");
	tree.Branch("idx", &idx, "idx/l" );
	tree.Branch("runN",&runN, "runN/I");
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
	auto exp_spaced = ExpList(binWdtMin, binWdtMax, N_iter);
	auto gtau  = TGraphAsymmErrors();
	auto gchi2 = TGraph();
	auto gprob = TGraph();
	TString name="fType"+to_string(static_cast<short int>(func_type))+"BinWdtStab";
	TString title="Bin width stability;Bin Width [ns];#tau_{#mu} [#mus]";
	gtau.SetNameTitle(name+"_tau",title);
	gchi2.SetNameTitle(name+"_chi2","#chi2 bin width stability;Bin Width [ns];#frac{#chi2}{NDF}");
	gprob.SetNameTitle(name+"_prob","prob bin width stability;Bin Width [ns];prob [%]");
	vector<int> nBinsVec(exp_spaced.size());
	std::transform(exp_spaced.begin(),exp_spaced.end(),nBinsVec.begin(),[rmax,rmin](double bwdt)->int{return (rmax-rmin)/bwdt;});
	auto last = std::unique(nBinsVec.begin(), nBinsVec.end());
    	nBinsVec.erase(last, nBinsVec.end());
	int i = 0;
	for( const auto & NbinsInFitRange : nBinsVec ){
		double binWdt = (rmax-rmin)/NbinsInFitRange;
		cout<<"bin width: "<<binWdt<<" n bins in range: "<<NbinsInFitRange;
		auto hs = GetDecayHistos(decaydf, binWdt, to_string(binWdt));
		//auto hup = hs["up"], hdwn = hs["dwn"];
		auto htot = hs["htot"];
		auto fitcore = AmuletFitCore(htot,rmin,rmax,"LERSN0Q");
		TFitResult fitres = fitcore.AmuletFit(func_type);
		//fitres.Print("V");
		const char* tau_name = "#tau";
		int tau_idx = fitcore.GetParIdx(tau_name);
		if ( tau_idx < 0 ){
			tau_name = "#tau_{#mu^{-}}_{mat}";
			tau_idx = fitcore.GetParIdx(tau_name);
			gtau.GetYaxis()->SetTitle(TString(tau_name)+" [#mus]");
		}
		if( fitcore.GetStatus() ){
			gtau.AddPoint(binWdt*1e9, fitres.Parameter(tau_idx)*1e6);
			gtau.SetPointError(i, 0, 0, abs(fitres.LowerError(tau_idx))*1e6, fitres.UpperError(tau_idx)*1e6);
			gchi2.AddPoint(binWdt*1e9, fitres.Chi2()/fitres.Ndf());
			gprob.AddPoint(binWdt*1e9, fitres.Prob()*100);
			i++;
			cout<<" fit success! prob: "<<fitres.Prob()*100;
			if (fitres.Prob()*100 > 5) cout<<" >5% !!!!!";
			cout<<endl;
		}else{
			cout<<" Fit FAILED, status: "<<fitres.Status()<<" valid: "<<fitres.IsValid()<<endl;
		}
	}
	TFile f(RootOut,"UPDATE");
	if(f.IsZombie()){
		cout<<"Failed to open file "<<RootOut<<endl;
	}else{
		gtau.Write("",TObject::kOverwrite);
		gchi2.Write("",TObject::kOverwrite);
		gprob.Write("",TObject::kOverwrite);
		f.Close();
	}
}

void WriteRangeStabilityFixedBins(RNode decaydf, double binWdt, int NitMin, int NitMax, TString RootOut, AmuletFitCore::FuncType func_type, vector<double> limits = {}){
	long double stopWdtUp = decaydf.Filter("topology==1").Mean<double>("stopWdt").GetValue();
	long double stopWdtDw = decaydf.Filter("topology==0").Mean<double>("stopWdt").GetValue();
	long double startWdt = decaydf.Mean<double>("startWdt").GetValue();
	long double stopWdt = std::min(stopWdtUp,stopWdtDw);
	long double rmin_min = stopWdt + startWdt - 300e-9;
	long double rmin_max =  11e-6;
	long double rmax_min = rmin_min + 500e-9;
	long double rmax_max =  11.5e-6;
	TString name="fType"+to_string(static_cast<short int>(func_type))+"RangeStab";
	if (limits.size() == 4){
		rmin_min = limits[0];
		rmin_max = limits[1];
		rmax_min = limits[2];
		rmax_max = limits[3];
		for ( const auto & limit : limits)
			name += "_"+to_string(limit*1e06);
	}
	cout<<"Range stability: "<<endl;
	cout<<"rmin min = "<<rmin_min<<"\t rmin max = "<<rmin_max<<endl;
	cout<<"rmax min = "<<rmax_min<<"\t rmax max = "<<rmax_max<<endl<<endl;
	long double stepMin = (rmin_max-rmin_min)/static_cast<long double>(NitMin);
	long double stepMax = (rmax_max-rmax_min)/static_cast<long double>(NitMax);
	auto gtau  = TGraph2DErrors();
	auto gchi2 = TGraph2D();
	TString title="Range stability;t_{min} [#mus];t_{max} [#mus]; #tau_{#mu} [#mus]";
	gtau.SetNameTitle(name+"Tau",title);
	gchi2.SetNameTitle(name+"Chi2","#chi2 range stability;rmin [#mus];rmax [#mus]; #frac{#chi2}{NDF}");
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
			const char* tau_name = "#tau";
			int tau_idx = fitcore.GetParIdx("#tau");
			if ( tau_idx < 0 ){
				tau_name = "#tau_{#mu^{-}}_{mat}";
				tau_idx = fitcore.GetParIdx(tau_name);
				gtau.GetZaxis()->SetTitle(((TString)tau_name)+" [#mus]");
			}
			if( fitcore.GetStatus() ){
				gtau.AddPoint(rmin*1e6, rmax*1e6, fitres.Parameter(tau_idx)*1e6);
				gtau.SetPointError(Npoint, 0, 0, fitres.Error(tau_idx)*1e6);
				gchi2.AddPoint(rmin*1e6, rmax*1e6, fitres.Chi2()/fitres.Ndf());
				Npoint++;
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


void WriteLifetimeFit(RNode decaydf, TString RootOut, double binWdt, double rmin, double rmax, AmuletFitCore::FuncType func_type, const char* topology = "htot", const char* name_prefix = ""){
	
	cout<<endl<<endl<<endl;	
	cout<<"LIFETIME FIT for "<<RootOut<<endl;
	cout<<"\tfType: "<<static_cast<short int>(func_type)<<endl;
	cout<<"\tbWdt: "<<binWdt<<endl;
	cout<<"\trmin: "<<rmin<<"\n\trmax: "<<rmax<<endl;
	cout<<"\tn bins: "<<(rmax-rmin)/binWdt<<endl;

	gStyle->SetOptFit(1111);
	gStyle->SetOptStat(10);
	gStyle->SetOptTitle(0);
	gROOT->ForceStyle();
	
	auto hs = GetDecayHistos(decaydf, binWdt, name_prefix);
	auto hfit = hs[topology]; //topology = "htot", "up", "dwn"
	auto fitcore = AmuletFitCore(hfit,rmin,rmax,"LERS");
	TFitResult fitres = fitcore.AmuletFit(func_type);
	fitres.Print("V");
	cout<<"Chi2 prob (GoF): "<<fitres.Prob()*100<<endl;
	auto f = TFile(RootOut,"UPDATE");
	TString name = "DecayFit_fType" + to_string(static_cast<short int>(func_type));
	if ((string)topology != "htot"){
		TString top(topology);
		top.ToUpper();
		name += "_topology"+top;
	}
	auto c = TCanvas(name,"decay_fit",1173,590);
	c.cd();
	auto h = fitcore.GetHistoClone();
	auto func = (TF1*)h->GetListOfFunctions()->At(0);
	func->SetLineWidth(2);
	double tmin, tmax;
	func->GetRange(tmin,tmax);
	h->GetYaxis()->SetTitleOffset(0.91);
	//h->GetYaxis()->SetNdivisions(505);
	h->GetYaxis()->SetTitleSize(0.05);
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetTitleOffset(0.95);
	h->GetXaxis()->SetNdivisions(520);
	h->GetXaxis()->SetTitleSize(0.05);
	h->GetXaxis()->SetLabelSize(0.05);
	h->SetFillStyle(3003);
	h->SetFillColor(4);
	h->SetMarkerColor(kBlack);
	h->SetMarkerStyle(7);
	h->Sumw2();
	auto hc = (TH1D*)h->Clone();
	hc->Draw("hist LF2");
	h->Draw("e1 same");
	hc->SetLineWidth(0);
	//h->DrawClone("same hist");
	/*double bkg_idx = func->GetParNumber("N_{bkg}");
	if( bkg_idx != -1 ){ //to show background component explicitly
		double Nbkg = func->GetParameter(bkg_idx);
		double eNbkg = func->GetParError(bkg_idx);
		double biw  = h->GetBinWidth(h->FindBin(tmin+(tmax-tmin)/2.));
		auto bkg = TF1("Bkg Fit (casual coinc)",[tmin,tmax,biw](double *x, double *p){return p[0]*biw*uniform_pdf(x[0],tmin,tmax);},rmin,rmax,1,"NL");
		bkg.SetParameter(0,Nbkg);
		bkg.SetParError(0,eNbkg);
		bkg.SetLineColor(kViolet);
		bkg.SetLineStyle(kDashed);
		bkg.DrawClone("same");
		func->SetTitle("S+B Fit");
		func->DrawClone("same");
		h->SetTitle("Data");
		c.BuildLegend(0.55,0.45,0.76,0.62,"","e1 x0 l p");
	}*/
	c.SetGrid();
	c.SetTicks();
	c.Modified();
	c.Update();
	c.SaveAs((RootOut.ReplaceAll(".root",name+".pdf")).ReplaceAll("DAQresults","/DAQresults/plots/"),"pdf");
	c.Write("",TObject::kOverwrite);
	
	f.Close();
}

void WriteAsymmetry(RNode df, TString RootOut, double binWdt, vector<short int> runs = {}, TString name_prefix = "Asymm"){
	auto f = TFile(RootOut,"UPDATE");
	const char* var = "dt";
	double maxup = df.Filter("topology==1").Max<double>(var).GetValue();
	double maxdw = df.Filter("topology==0").Max<double>(var).GetValue();
	double max = std::max(maxup,maxdw);
	TH1D hup = *df.Filter("topology==1").Histo1D<double>(	{name_prefix+"UpDecay",name_prefix+"UpDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
					static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	TH1D hdwn = *df.Filter("topology==0").Histo1D<double>( {name_prefix+"DwnDecay",name_prefix+"DwnDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
					static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);

	auto hasymmetry = hup.GetAsymmetry(&hdwn);
	hasymmetry->Write("",TObject::kOverwrite);
	delete hasymmetry;
	
	TString name_prefix_old = name_prefix;
	for( const auto & run : runs ){
		name_prefix += "run"+to_string(run)+"_";
		hup = *df.Filter("topology==1")
			 .Filter([run](int runN){return runN==run;},{"runN"})
			 .Histo1D<double>({name_prefix+"UpDecay",name_prefix+"UpDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
		hdwn = *df.Filter("topology==0")
			  .Filter([run](int runN){return runN==run;},{"runN"})
			  .Histo1D<double>( {name_prefix+"DwnDecay",name_prefix+"DwnDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
						static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);

		hasymmetry = hup.GetAsymmetry(&hdwn);
		hasymmetry->Write("",TObject::kOverwrite);
		delete hasymmetry;
		name_prefix = name_prefix_old;	
	}
	
	f.Close();

}

//stability analysis: opt "B" for bin width fit stab, opt: "R" for range stab.
void StabilityAnalysis(RNode decaydf, Option_t *opt, TString RootOut, double rmin, double rmax, double binWdt, AmuletFitCore::FuncType func_type){
	if( TString(opt).Contains("B")){
		int N_iters = 1000;
		WriteBinNumberStabilityFixedRange(decaydf, rmin, rmax, N_iters, RootOut, func_type);
	}	
	if( TString(opt).Contains("R") ){
		int N_iters_min = 100;
		int N_iters_max = 100;
		WriteRangeStabilityFixedBins( decaydf, binWdt, N_iters_min, N_iters_max, RootOut, func_type);

	}
}

int main(int argc, char** argv){
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
	
	gStyle->SetPalette(kViridis);
	gStyle->SetNumberContours(999);
	
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
	vector<AmuletFitCore::FuncType> fTypes;
	AmuletFitCore::FuncType fStabType = AmuletFitCore::FuncType::kMassimo;
	vector<short int> runs;
	vector<const char*> topologies;
	auto measconf = ChooseConfig(RootOut);
	switch( measconf ){
		case RunConfiguration::kCarbon:
			binWdt = 29.321e-9;  //OPT
			rmin = 500e-9;       //OPT
			rmax = 10e-6;        //OPT
			MEASTIME = 3.0464e06;
			TAUMATERIAL = TAU_carbon;
			runs = {1,5,7};
			fTypes.push_back(AmuletFitCore::FuncType::kMassimo);
			fTypes.push_back(AmuletFitCore::FuncType::kMassimoFraction);
			fStabType = AmuletFitCore::FuncType::kMassimo;
		break;	
		case RunConfiguration::kAl:
			binWdt = 38.15e-9; //opt?
			rmin = 500e-9;  //TODO
			rmax = 10e-6;   //TODO
			MEASTIME = 1775700;
			TAUMATERIAL = TAU_Al;
			topologies = {"htot"/*,"up","dwn"*/};
			CHARGERATIORANGECORRECTION = 0.054516; 
			fTypes.push_back(AmuletFitCore::FuncType::kMaterial);
			fStabType = AmuletFitCore::FuncType::kMaterial;
		break;	
		case RunConfiguration::kNaCl:
			binWdt = 18.18e-9; //opt?
			rmin = 0.5e-6;  //TODO
			rmax = 8.7e-6;  //TODO
			MEASTIME = 2338100;
			TAUMATERIAL = TAU_NaCl;
			topologies = {"htot"/*,"up","dwn"*/};
			CHARGERATIORANGECORRECTION = 0.058460; 
			fTypes.push_back(AmuletFitCore::FuncType::kMaterial);
			fStabType = AmuletFitCore::FuncType::kMaterial;
		break;	
		case RunConfiguration::kBkg:
			binWdt = 1.2e-7; //quella che sceglie in automatico il ttree
			rmin = 30e-6;    //scelta guardando plot
			rmax = 38e-6;    //scelta guardando plot
			MEASTIME = 1194420.;
			TAUMATERIAL = 0;
			fTypes.push_back(AmuletFitCore::FuncType::kUniform);
			fTypes.push_back(AmuletFitCore::FuncType::kConst);
		break;	
		case RunConfiguration::kBon:
			binWdt = 60e-9; //TODO
			rmin = 0.5e-6;  //TODO
			rmax = 9.1e-6;  //TODO
			MEASTIME = 3115900;
			TAUMATERIAL = TAU_carbon;
			fTypes.push_back(AmuletFitCore::FuncType::kMassimo);
			fTypes.push_back(AmuletFitCore::FuncType::kMassimoFraction);
		break;	
		case RunConfiguration::kBoff:
			binWdt = 60e-9; //TODO
			rmin = 0.5e-6;  //TODO
			rmax = 9.1e-6;	//TODO
			MEASTIME = 2935200;
			TAUMATERIAL = TAU_carbon;
			runs = {10,8};
			fTypes.push_back(AmuletFitCore::FuncType::kMassimo);
			fTypes.push_back(AmuletFitCore::FuncType::kMassimoFraction);
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
	
	//ExportTxt(decaydf, "dt", RootOut);
	if(measconf != RunConfiguration::kBkg)
		StabilityAnalysis(decaydf, "BR", RootOut, rmin, rmax, binWdt, fStabType );
	for( const auto & fType : fTypes )
		if( topologies.size() == 0 )
			WriteLifetimeFit(decaydf, RootOut, binWdt, rmin, rmax, fType);
		else 
			for( const auto & top : topologies )
				WriteLifetimeFit(decaydf, RootOut, binWdt, rmin, rmax, fType, top);
	
	WriteAsymmetry(decaydf, RootOut, 6.5e-7/3, runs);
		
	
	cout<<"File "<<RootOut<<" UPDATED"<<endl<<endl;
	
	return 0;
}
