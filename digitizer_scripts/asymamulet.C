
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

void asymamulet(){

	auto res = GetMagneticDF();
	auto sig = GetDFAsymmetry(res[0], (10e-6-500e-9)/25, 10e-6);
	auto bkg = GetDFAsymmetry(res[1], (10e-6-500e-9)/25, 10e-6);

	auto asym.Clone("asym");

}

int main(){
	asymamulet();
}
