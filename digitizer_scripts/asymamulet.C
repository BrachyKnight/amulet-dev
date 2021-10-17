
using ROOT::RDF::RNode, ROOT::RDataFrame, ROOT::VecOps::RVec, ROOT::Math::exponential_pdf;
using std::vector, std::pair, std::runtime_error, std::cout, std::endl, std::string, std::map, std::make_pair, std::to_string;

TH1* GetDFAsymmetry(RNode df, double binWdt, double max = 0, TString global_cut = "", TString name_prefix = "Asymm"){
	const char* var = "dt";
	if (max==0){
		double maxup = df.Filter("topology==1").Max<double>(var).GetValue();
		double maxdw = df.Filter("topology==0").Max<double>(var).GetValue();
		max = std::max(maxup,maxdw);
	}
	RNode cutdf = df;
	if( global_cut != ""  )
		cutdf = df.Filter((const char*)global_cut);
	TH1D hup = *cutdf.Filter("topology==1").Histo1D<double>(	{name_prefix+"UpDecay",name_prefix+"UpDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
					static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);
	TH1D hdwn = *cutdf.Filter("topology==0").Histo1D<double>( {name_prefix+"DwnDecay",name_prefix+"DwnDecay;time [s];events per "+to_string(binWdt*1e06).substr(0,4)+" #mus",
					static_cast<int>((max+0.1*max)/binWdt),0,max+0.1*max},	var);

	auto hasymmetry = hup.GetAsymmetry(&hdwn);
	
	return hasymmetry;	
}

RNode GetDecayTreeDF(TString RootOut){
	TFile* outFile = TFile::Open(RootOut,"UPDATE");
	RDataFrame decaydf("decays",RootOut);
	outFile->Close();
	delete outFile;
	return decaydf;
}

vector<RNode> GetMagneticDF(){
	RNode BonDF = GetDecayTreeDF("../DAQresults/muLifetimeBon.root");
	RNode BoffDF = GetDecayTreeDF("../DAQresults/muLifetimeBoff.root");
	return {BonDF, BoffDF};
}

TH1* GetFourierTransform(TH1* asymmetry){
  	
	TH1* hfft = new TH1D(*dynamic_cast<TH1D*>(asymmetry));

	asymmetry->FFT(hfft, "MAG R2C EX");

	double	fftSize = asymmetry->GetNbinsX();
	double timeStep = asymmetry->GetBinWidth(2);
	
	// here is where I "rescale" the x-axis
	// the fftSize-1 is for an overflow/underflow bin offset that biases the frequency
	TH1D* hFFT = new TH1D("HistFFT", "HistFFT"
			, fftSize / 2
			, 0
			, hfft->GetNbinsX() / (fftSize*timeStep));

	// copy FFT histogram to new, scaled x-axis histogram 
	for(int j = 0; j < hFFT->GetNbinsX(); j++)
		hFFT->SetBinContent(j, hfft->GetBinContent(j));
	return hfft;
}

TH1* GetBKGSubtractedAsymmetry(){
	auto res = GetMagneticDF();
	auto sig = GetDFAsymmetry(res[0], (10e-6-500e-9)/25, 10e-6);
	auto bkg = GetDFAsymmetry(res[1], (10e-6-500e-9)/25, 10e-6);

	auto asym = (TH1*)sig->Clone("asym");
	asym->Add(bkg, -1.);
	delete sig;
	delete bkg;
	return asym;
}

TH1* GetSignalAsymmetry(){
	auto res = GetMagneticDF();
	auto sig = GetDFAsymmetry(res[0], (10e-6-500e-9)/25, 10e-6);

	return sig;
}

TH1* GetBkgAsymmetry(){
	auto res = GetMagneticDF();
	auto bkg = GetDFAsymmetry(res[1], (10e-6-500e-9)/25, 10e-6);

	return bkg;
}

void asymamulet(){

	gStyle->SetOptFit(1111);
	gStyle->SetOptTitle(0);
	TCanvas* c0 = new TCanvas("Bkg","Bkg",500,500);
	auto bkg = GetBkgAsymmetry();
	bkg->Draw("e1 p");
	TCanvas* c1 = new TCanvas("BkgSub","BkgSub",500,500);
	auto sigsub =  GetBKGSubtractedAsymmetry();
	sigsub->Draw("e1 p");
	TCanvas* c2 = new TCanvas();
	auto sig = GetSignalAsymmetry();
	sig->Draw("e1 p");
	
	
	auto f = new TF1("func",[](double *x, double *p){
			double Nd = 2518;
			double Nu = 1580;
			double A = Nu + Nd;
			double B = Nu - Nd;
			double bu = 68.06;
			double bd = 66.16;
			double C = bu - bd;
			double D = bu + bd;
			double xi = p[0];
			double Xi = xi/6;
			double omega = p[1];
			double t = x[0];
			double tau = 2.033e-6;//(2.192e-6*(1/pow(6.97e-8,2))+1.953e-6*(1/pow(4.52e-8,2)))/(1/pow(6.97e-8,2)+1/pow(4.52e-8,2));
			double NUM = A*Xi*cos(omega*t) + B + C*exp(t/tau);
			double DEN = B*Xi*cos(omega*t) + A + D*exp(t/tau);
			return NUM/DEN;;
			}, 2e-6, 10e-6, 2,"NL");
	f->SetParNames("#xi","#omega");
	f->SetParameters(0.35,2e06);

	sig->Fit(f, "R", "", 2e-6, 10e-6);
	
	c0->cd();
	
	auto f1 = new TF1("f1",[](double *x, double *p){
			double Nd = 4267;
			double Nu = 2562;
			double A = Nu + Nd;
			double B = Nu - Nd;
			double bu = 119;
			double bd = 132;
			double C = bu - bd;
			double D = bu + bd;
			double t = x[0];
			double tau = 2.01e-6; 
			double NUM = B + C*exp(t/tau);
			double DEN = A + D*exp(t/tau);
			return NUM/DEN + p[0];;
			}, 2e-6, 10e-6, 1,"NL");
	f1->Draw("same");
	f1->SetParNames("K");
	bkg->Fit(f1,"R","",2e-6,10e-6);
	double chi2 = bkg->Chisquare(f1, "R");
	cout<<"chi2 = "<<chi2<<endl;
	cout<<"prob = "<<TMath::Prob(chi2, 21)<<endl;

}

int main(){
	asymamulet();
}
