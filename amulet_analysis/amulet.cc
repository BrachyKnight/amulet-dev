//muon analysis class for laboratory particle physics UNIMIB
//amulet = Analysis MUon LifETime 
//Developed using ROOT 6.22/07 (guaranteed to work with version 6.22/07)

#include "amulet.h"

//constructors
amulet::amulet()
   	: _dfUp(RDataFrame(0)), _dfDwn(RDataFrame(0)), _OutFile(NULL)
   	{ }

amulet::amulet(RNode up, RNode dwn, const char* OutFilePath, const char* OutFileOption )
   	: _dfUp(up), _dfDwn(dwn){
   		_OutFile = new TFile( OutFilePath, OutFileOption );
   		const bool fileDirSuccess = (_OutFile==NULL) ? gROOT->cd() : _OutFile->cd();
   		if(_OutFile!=NULL){ //extract file name without extension and remove "Lifetimes" if it is present
   			TString name = (TString)(((TString*)_OutFile->GetName())->Tokenize('/')->Last())->GetName();
   			if( name.Contains("Lifetimes") )
   				name.ReplaceAll("Lifetimes", "");
   			if( name.Contains("Analyzed") )
   				name.ReplaceAll("analyzed", "");
   			if ( name.EndsWith(".root") )
   				name.ReplaceAll(".root", "");
   			_OutFileName = name;
   		}
   		if(_verbosity>0){
   			(!fileDirSuccess) ? cout<<"Warning: directory has not changed, current directory:"<<endl : cout<<"Current directory:"<<endl;
   			gDirectory->pwd();
   		}
   	}

amulet::amulet(RDataFrame up, RDataFrame dwn, const char* OutFilePath, const char* OutFileOption )
   	: _dfUp(up), _dfDwn(dwn){
   		_OutFile = new TFile( OutFilePath, OutFileOption );
   		bool fileDirSuccess = (_OutFile==NULL) ? gROOT->cd() : _OutFile->cd();
   		if(_OutFile!=NULL){ //extract file name without extension and remove "Lifetimes" if it is present
   			TString name = (TString)(((TString*)_OutFile->GetName())->Tokenize('/')->Last())->GetName();
   			if( name.Contains("Lifetimes") )
   				name.ReplaceAll("Lifetimes", "");
   			if( name.Contains("Analyzed") )
   				name.ReplaceAll("analyzed", "");
   			if ( name.EndsWith(".root") )
   				name.ReplaceAll(".root", "");
   			_OutFileName = name;
   		}
   		if(_verbosity>0){
   			(!fileDirSuccess) ? cout<<"Warning: directory has not changed, current directory:"<<endl : cout<<"Current directory:"<<endl;
   			gDirectory->pwd();
   		}
   	}
   
amulet::amulet(const amulet &old_amu) //cpy constr forwards to master constructor
	: amulet(old_amu._dfUp, old_amu._dfDwn, old_amu._OutFile->GetName(), old_amu._OutFile->GetOption())
	{ }

amulet::~amulet(){
		_OutFile->Close();
		delete _OutFile;  
	}

const bool amulet::SetOutputFile(const char* rootOutName, const char* rootFileOpt){
	if (_OutFile == NULL ){
		_OutFile = new TFile( rootOutName, rootFileOpt );
		TString name = (TString)(((TString*)_OutFile->GetName())->Tokenize('/')->Last())->GetName();
		if( name.Contains("Lifetimes") )
			name.ReplaceAll("Lifetimes", "");
		if ( name.EndsWith(".root") )
			name.ReplaceAll(".root", "");
		_OutFileName = name;
		return true;
	}else{
		 if(_verbosity>0) cout<<"File associated to amulet object already present"<<endl;
		 return false;
	}
}

const bool amulet::MakeAndChangeCurrentDir( const char* subdirName ){
	bool OK = true;
	TDirectory* dir = NULL;
	if( IsFileOK() ){
		if(!(subdirName && !subdirName[0])){ //this should mean if (subdirName != "") ....
			dir = ( _OutFile->GetDirectory(subdirName)==nullptr ) ? _OutFile->mkdir(subdirName) : _OutFile->GetDirectory(subdirName);
			if (dir==NULL)
				OK=false;
			else
				_OutFile->GetDirectory(subdirName)->cd();
		}else _OutFile->cd();
	}else OK=false;
	if( _verbosity>0 || !OK ){
		if (!OK) throw std::runtime_error(" Warning: impossible to create directory ");//cout<<" Warning: impossible to create directory "<<subdirName<<endl;
		cout<<"Current directory:"<<endl;
		gDirectory->pwd();
	}
	return OK;
}

amulet amulet::ApplyFiltersAndPrintReport(  strMap_t UpFilter, strMap_t DwnFilter, const char* NewOutFileName, const char* option  ){
	//what follows is a trick to make Filter operation recursevly for all the filters in map
	//notice that the same can be done also for the define operation
   	auto latestDFup = std::make_unique<RNode>(_dfUp);
   	auto latestDFdwn = std::make_unique<RNode>(_dfDwn);
   	for (auto const& filter : UpFilter)
		latestDFup = std::make_unique<RNode>(latestDFup->Filter(filter.second, filter.first));                 
   	for (auto const& filter : DwnFilter)
		latestDFdwn = std::make_unique<RNode>(latestDFdwn->Filter(filter.second, filter.first));                                                                             
                                           
	cout<<"\nFILTERS:\nfilters report dfUp"<<endl;
   	auto UpFilterReport = latestDFup->Report();
	UpFilterReport->Print();
	cout<<"filters report dfDwn"<<endl;
	auto DwnFilterReport = latestDFdwn->Report();
	DwnFilterReport->Print();
	cout<<endl;
	return amulet( *latestDFup,  *latestDFdwn, NewOutFileName, option);
}

void amulet::SnapshotAssociatedDataFramesInRootFile(){
	if(IsFileOK()){
		string name = _OutFile->GetName();
		_OutFile->Close();
		delete _OutFile;
		RSnapshotOptions optsSnapshot;
		optsSnapshot.fMode = "UPDATE";
		optsSnapshot.fLazy = false;
		optsSnapshot.fOverwriteIfExists = true;
		_dfUp.Snapshot<double, double, double, double, double, double>("DeltaTUpDecays", name.c_str(),  _dfUp.GetColumnNames(), optsSnapshot);
		_dfDwn.Snapshot<double, double, double, double, double, double>("DeltaTDwnDecays", name.c_str(), _dfDwn.GetColumnNames(), optsSnapshot);
		if(_verbosity>1) cout<<"\nDATAFRAME CREATED IN FILE "<<name<<endl<<endl;
		_OutFile = new TFile( name.c_str() , "UPDATE" );
	} else cout<<"ERROR in SnapshotAssociatedDataFramesInRootFile()"<<endl;
}

//analysis methods
static TF1* fixParameters( TF1* myfunc ) //if parameter name is a number then fix that parameter and set that number
{
	for ( int i = 0; i<myfunc->GetNpar(); i++ )
	{
		TString parName = myfunc->GetParName(i);
		if( parName.IsFloat() )
			myfunc->FixParameter( i, parName.Atof() );
	}
	return myfunc;
}

static TF1* fixParametersForNormilized( TF1* myfunc, double lowerLimit, double upperLimit ) //if parameter name is a number then fix that parameter and set that number
{
	for ( int i = 0; i<myfunc->GetNpar(); i++ )
	{
		TString parName = myfunc->GetParName(i);
		if( parName.IsFloat() )
			myfunc->FixParameter( i, parName.Atof() );
		else if( parName.Contains("upperLimit") )
			myfunc->FixParameter(i, upperLimit);
		else if( parName.Contains("lowerLimit") )
			myfunc->FixParameter(i, lowerLimit);
	}
	return myfunc;
}

template<class T>
vecWithName_t amulet::partialFit(T histo, const double fitRange[2],  string name, const char* fitOpts, bool drawPlots, vector<TList*> &histoDrawLists){
	bool WasBatch = gROOT->IsBatch();
	gROOT->SetBatch(kTRUE);
	//extract fit limits
	double lowerLimit = fitRange[0];
	double upperLimit = fitRange[1];

	//create fit function
	const char* functionName = (name+"_func").c_str();
	auto fitFunc = new TF1(functionName, _function.c_str(), lowerLimit, upperLimit );
	//automize parameter initialization based on the _function
	//notice that only some functions are provided so be careful!!!
	if(name != "asymmetry") //TODO fissare la normalizzazione perche terranova ha detto che e sbagliato lasciarla libera
	{
		if(fitFunc->GetNpar() == 3 && _function == "[N0]*exp(-x/[#tau])+[B]") //caso base
		{
			fitFunc->SetParameters(2e-6, histo->GetBinContent(histo->GetBin(lowerLimit)) , histo->GetBinContent(histo->GetBin(upperLimit)));
			//fitFunc->SetParLimits(0,1e-8,5e-5);
		}
		//"([N0]/([f]+1))*([f]*exp(-x/[#tau])+exp(-x/[#tau1])) + [B]" //caso con rapporto
		else if (fitFunc->GetNpar() == 5 && ((TString)_function).BeginsWith("([") && ((TString)_function).Contains("])+exp(-x/[") && ((TString)_function).EndsWith("])) + [B]"))
		{
			fitFunc->SetParameters(histo->GetMean()/2., histo->GetMean(), histo->GetBinContent(histo->FindBin(upperLimit)), histo->GetBinContent(histo->FindBin(lowerLimit)), 1.);
			fitFunc->SetParLimits(0, 0.8e-7,0.8e-5); //per Al 0.8e-6
			fitFunc->SetParLimits(1, 1e-8,1e-5);  //2.197e-6
			//fitFunc->SetParLimits(3,0,histo->GetBinContent(histo->FindBin(lowerLimit))*(1.3));
			//fitFunc->SetParLimits(4,0.4,1.4);
			fitFunc = fixParameters(fitFunc);
		}
		//"[N0]*exp(-x/[#tau])+[N1]*exp(-x/[#tau1])+[B]" //caso senza rapporto
		else if (fitFunc->GetNpar() == 5 && ((TString)_function).BeginsWith("[") && ((TString)_function).Contains("]*exp(-x/[") && ((TString)_function).EndsWith("])+[B]"))
		{
			fitFunc->SetParameters(histo->GetMean()/2., histo->GetMean(), histo->GetBinContent(histo->FindBin(upperLimit)), histo->GetBinContent(histo->FindBin(lowerLimit))/2., histo->GetBinContent(histo->FindBin(lowerLimit))/2.);
			fitFunc->SetParLimits(1, 1e-8,1e-5);  //2.197e-6
			fitFunc->SetParLimits(0, 0.8e-7,0.8e-5); //per Al 0.8e-6
			//fitFunc->SetParLimits(2,0,histo->GetBinContent(histo->GetBin(lowerLimit))*(1.3));
			fitFunc = fixParameters(fitFunc);
		}//"(exp(-x/[#tau])+[b])/([#tau]*(exp(-[lowerLimit]/[#tau])-exp(-[upperLimit]/[#tau]))-[lowerLimit]*[b]+[upperLimit]*[b])"
		else if(_DoNormalize && fitFunc->GetNpar() == 4 && ((TString)_function).BeginsWith("(exp(-x/") && ((TString)_function).Contains("])/([") && ((TString)_function).Contains("]))-[lowerLimit]*") && ((TString)_function).EndsWith("])") ){
			fitFunc->SetParameters(2.197e-6, histo->GetBinContent(histo->FindBin(upperLimit))/(histo->Integral(histo->FindBin(lowerLimit), histo->FindBin(upperLimit))));
			fitFunc = fixParametersForNormilized(fitFunc, lowerLimit, upperLimit);
		}//"([f]*exp(-x/[#tau])+exp(-x/[#tau1])+[b])/(-[lowerLimit]*[b] + [upperLimit]*[b] + (exp(-[lowerLimit]/[#tau]) - exp(-[upperLimit]/[#tau]))*[f]*[#tau] + (e^(-[lowerLimit]/[#tau1]) - exp(-[upperLimit]/[#tau1]))*[#tau1])"
		else if(_DoNormalize && fitFunc->GetNpar() == 6 && ((TString)_function).BeginsWith("([f]*exp(-x/") && ((TString)_function).Contains(")/(-[lowerLimit]*") && ((TString)_function).Contains("]))*[f]*[") && ((TString)_function).EndsWith("])")){
			fitFunc->SetParameters(2.45944e-06, 1.12796e-06, histo->GetBinContent(histo->FindBin(upperLimit))/(histo->Integral(histo->FindBin(lowerLimit), histo->FindBin(upperLimit))),1.);
			fitFunc = fixParametersForNormilized(fitFunc, lowerLimit, upperLimit);
		}else if(_function =="[B]"){
			fitFunc->SetParameter(0,histo->GetBinContent(histo->FindBin(35e-6)));
		}else
			throw std::runtime_error(((string)("Have you enabled/disabled _DoNormalize??????\nCASE NOT IMPLEMENTED, CHECK YOUR FUNCTION\nnPars="+std::to_string(fitFunc->GetNpar())+"\nformula="+fitFunc->GetExpFormula()+"\n it must be a string exactly equal to what expected!!!!")).c_str());
	}
	//fit and extract results
   	TFitResultPtr r = histo->Fit(fitFunc, fitOpts);

	const bool isvalid = r->IsValid();
	const int status = r->Status();
	const double chi2_reduced = r->Chi2()/r->Ndf();
	const double prob = r->Prob();
   	vecWithName_t fitResults = { {"fit validity", (double)isvalid}, {"fit status", (double)status}, /*{"#chi^{2}", r->Chi2()},*/ {"#chi^{2}/NDF", chi2_reduced}, {"fit prob", prob} };
	vector<double> params(r->GetParams(), r->GetParams()+r->NPar());
	vector<double> paramsErr = r->Errors();
	vector<string> parNames;
	vector<string> parNamesErr;
	for(unsigned int i = 0; i<r->NPar(); i++){
		parNames.push_back("par"+r->ParName(i));
		parNamesErr.push_back("err"+r->ParName(i));
	}

	vecWithName_t pars;
	vecWithName_t parsErr;
	assert(params.size() == parNames.size());
	for (size_t i = 0; i < parNames.size(); ++i){
    	pars.push_back( make_pair( parNames[i], params[i] ) );
		parsErr.push_back( make_pair( parNamesErr[i], paramsErr[i] ) );
	}
	//stima coincidenze casuali partendo dal parametro b del fit
	/*int idxB = 1;
   	if( _DoNormalize && name != "asymmetry" && drawPlots){
   		double b = r->Parameter(idxB), bErr = r->Error(idxB);
   		cout<<"-----------------------------------------"<<endl;
   		if(r->ParName(idxB)=="b"){
   			double bInRange = b*(upperLimit-lowerLimit), bInRangeErr = (upperLimit-lowerLimit)*bErr;
   			cout<<"from fit b_fit = "<<b<<"+-"<<bErr<<endl;
   			cout<<"[b_fit*(b-a)] = " <<bInRange<<"+-"<<bInRangeErr<<endl;
   			if( IsSetTotalRate() ){
   				cout<<"calcolo nel quaderno di lab dice che [b_fit*(b-a)] dovrebbe essere uguale a Rcc/Rtot"<<endl;
   				double Rtot = GetTotalRate();
   				double w = histo->GetBinWidth(1);
   				cout<<"Rtot = "<<Rtot<<endl;
   				cout<<"(b-a) = "<<(upperLimit-lowerLimit)<<endl;
   				cout<<"Bin Width w = "<<w<<endl;
   				cout<<"moltiplicando [b_fit*(b-a)]*Rtot quindi dovrei ottenere Rcc"<<endl;
   				cout<<"[b_fit*(b-a)]*Rtot = "<<bInRange*Rtot<<"+-"<<bInRangeErr*Rtot<<" (assuemendo Rtot senza errore....)"<<endl;
   				//cout<<"[b_fit*(b-a)]*Rtot/w = "<<bInRange*Rtot/w<<"+-"<<bInRangeErr*Rtot/w<<" (assuemendo Rtot senza errore....)"<<endl;
   				//cout<<"[b_fit*(b-a)]*Rtot/(b-a) = "<<bInRange*Rtot/(upperLimit-lowerLimit)<<"+-"<<bInRangeErr*Rtot/(upperLimit-lowerLimit)<<" (assuemendo Rtot senza errore....)"<<endl;
   				pars.push_back( { "Rcc stimato", bInRange*Rtot } ); 
   				parsErr.push_back( { "Rcc_error", bInRangeErr } );
   			}else cout<<"total rate not set"<<endl;
   		}else{
   		   	cout<<r->ParName(idxB)<<" = "<<b<<" +- "<<bErr<<endl;
   			cout<<"parameter of index "<<idxB<<" is not equal to \'b\' but instead is equal to "<<r->ParName(idxB)<<endl;
   		}
   		cout<<"-----------------------------------------"<<endl;
   	}*/
   	//print some other infos:
   	if(!_DoNormalize && name != "asymmetry" && drawPlots){
   		double I, Ierr;
   		I = histo->IntegralAndError(histo->FindBin(lowerLimit), histo->FindBin(upperLimit), Ierr);
   		double wdth = histo->GetBinWidth(histo->FindBin(lowerLimit));
   		cout<<"-----------------------------------------"<<endl;
   		cout<<"range: lower = "<<lowerLimit<<" upper = "<<upperLimit<<endl;
   		cout<<"Integral in range = "<<I<<"+-"<<Ierr<<"           "<<endl;
   		cout<<"bin width = "                       <<wdth<<endl;
   		cout<<"n bins in range = range/binWidth = "<<(upperLimit-lowerLimit)/wdth<<endl;
   		cout<<"n bins in range = (bin upper - bin lower) = "<<(histo->FindBin(upperLimit ) -histo->FindBin(lowerLimit))<<endl;
   		cout<<"-----------------------------------------"<<endl;
   	}
	fitResults.insert(fitResults.end(), pars.begin(), pars.end());
	fitResults.insert(fitResults.end(), parsErr.begin(), parsErr.end());
	
	if(drawPlots){
      	histo->SetTitle((name+"; #Deltat [s]; N events").c_str());
	  	histo->SetName(name.c_str());
      	if ( histoDrawLists.size() == 0 )
      		histoDrawLists.push_back( new TList() );
      	else if ( histoDrawLists.size() >= 3 )
      		cout<<"WARNING IN TLIST VECTOR IN FUNCTION partialFit FIT"<<endl;
      histoDrawLists[0]->Add(&(*histo)); 
	}
	else if (!_DoPulls)
		delete fitFunc;
	
	if(_DoPulls){
		if (drawPlots){
      		if ( histoDrawLists.size() == 1 )
      			histoDrawLists.push_back( new TList() );
      		else if ( histoDrawLists.size() == 0 || histoDrawLists.size() >= 3 )
      			cout<<"WARNING IN TLIST VECTOR IN FUNCTION partialFit PULLS"<<endl;
      	}else{
      		histoDrawLists.push_back( NULL );
		}
		vecWithName_t pull_results = partialPulls(histo, fitFunc, name, drawPlots, histoDrawLists[1]);
		if( pull_results.size() != 0 )
			fitResults.insert(fitResults.end(), pull_results.begin(), pull_results.end());
		else
			if(_verbosity>1) cout<<"Impossible to produce PULLS for "<<name<<endl;
	}
	gROOT->SetBatch(WasBatch);
	return fitResults;
}

template<class T = ROOT::RDF::RResultPtr<TH1D>> vecWithName_t amulet::partialPulls(T histo, const TF1* FitFunction, string name, bool drawPlots, TList *&histoPullToDrawList){
	TH1D* histoPull = new TH1D((name+"_pulls").c_str(), (name+"_pulls ;pull ;N").c_str(), histo->GetNbinsX()/10, 0, 0);
	double xmin=0., xmax=0.;
	FitFunction->GetRange(xmin,xmax);
	int binMin = std::max(1, histo->FindBin(xmin));
	int binMax = std::min(histo->FindBin(xmax),histo->GetNbinsX());
	if ( (binMax - binMin ) <= 1 )
		return {}; 	//returns empty vector if pulls are impossible to evaluate
	//is this correct? is correct bin width as X-error? Should I take into account error on the evaluated _function value that comes from error on parameters? remember that it is related (correlation) with the original histogram.
	//is in general the pull method defined like this correct?
	for( int i = binMin; i<=binMax; i++ )
		histoPull->Fill((histo->GetBinContent(i) - FitFunction->Eval( histo->GetXaxis()->GetBinCenter(i) ) )/sqrt( pow(histo->GetBinError(i),2) + pow(histo->GetXaxis()->GetBinWidth(i) ,2) ));
	TFitResultPtr r = histoPull->Fit( "gausn", "SLQ" );
	if(drawPlots && histoPullToDrawList != NULL){
		histoPullToDrawList->Add(histoPull);
		histoPullToDrawList->SetName("PULLS_list");
	}
	if( r->IsValid() && r->Status()==0)
		return { {"PULL_mean", r->Parameter(1)}, {"PULL_sigma", r->Parameter(2)} };
	else
		return {};
}

nestedVecWithName_t amulet::Perform_Lifetime_Fit(double* fitLims,  vector<string> names, bool drawPlots, const char* fitOpts){
   	nestedVecWithName_t fitResults;
	vecWithName_t partialFitResult;

	//extract fit limits
	const double lowerLimit = fitLims[0];
	const double upperLimit = fitLims[1];
	const double fitRange[2] = {lowerLimit, upperLimit};
	const int nBins = (int)fitLims[2];
	
	//extract names from input vector
	const char* outDirNameInFile = names[0].c_str();
	string outDecayCanvasDir  = names[2] + names[1]; //outDir + outDecayCanvasName
	const char* outDecayCanvasName = names[1].c_str();
	const char* UpDecayName = names[3].c_str();
	const char* DwnDecayName = names[4].c_str();
	const bool comb = (names.size() == 6);
	const char* combinedDecayName = (comb) ? names[5].c_str() : "NoComb";

        auto histUpDecay  = _dfUp.Histo1D<double>({UpDecayName, UpDecayName, nBins, lowerLimit*(0.5), upperLimit*(1+0.25)}, "decayUP.dtFall"); //TODO here upper limit changes with function range, probably wrong when performin stability analysis
        auto histDwnDecay = _dfDwn.Histo1D<double>({DwnDecayName, DwnDecayName, nBins, lowerLimit*(0.5), upperLimit*(1+0.25)}, "decayDOWN.dtFall");

	vector<TList*> histoDrawLists;
	
	if( !_DoNormalize ){
		//------------------------------------------UP DECAY FIT----------------------------------------;
	 	partialFitResult = partialFit(histUpDecay, fitRange,  UpDecayName, fitOpts, drawPlots, histoDrawLists);
	   	fitResults.push_back(make_pair(UpDecayName, partialFitResult));

		//--------------------------------------DOWN DECAY FIT------------------------------------------
	 	partialFitResult = partialFit(histDwnDecay, fitRange,  DwnDecayName, fitOpts, drawPlots, histoDrawLists);
	   	fitResults.push_back(make_pair(DwnDecayName, partialFitResult));

		//---------------------------------------COMBINED FIT-------------------------------------------
		if (comb){ //perform combined fit and asymmetry only if a name for combined is specified
			if (histDwnDecay->GetSumw2N() == 0) histDwnDecay->Sumw2(kTRUE);
			if (histUpDecay->GetSumw2N() == 0) histUpDecay->Sumw2(kTRUE);
			TH1D* histCombined = new TH1D(*histUpDecay); //initialize new histogram to the old one (same number of bins...)
			if(!histCombined->Add(&(*histUpDecay), &(*histDwnDecay))) cout<<"WARNING IN HISTOGRAM ADDITION"<<endl; //add the two
	 		partialFitResult = partialFit<TH1D*>(histCombined, fitRange,  combinedDecayName, fitOpts, drawPlots, histoDrawLists);
	   		fitResults.push_back(make_pair(combinedDecayName, partialFitResult));

			TH1D* asymmetry = dynamic_cast<TH1D*>(histUpDecay->GetAsymmetry(&(*histDwnDecay))); //TH1D asymmetry = (*histUpDecay - *histDwnDecay)/histCombined;
			string OldFunc = this->GetFunction();
			this->SetFunction("pol1");
			partialFitResult = partialFit<TH1D*>(asymmetry, fitRange, "asymmetry", fitOpts, drawPlots, histoDrawLists);
	   		this->SetFunction(OldFunc);
	   		fitResults.push_back(make_pair("asymmetry", partialFitResult));
		}
	}else{

		//------------------------------------------UP DECAY FIT----------------------------------------
		TH1D* histUpDecayTemp = (TH1D*)histUpDecay->Clone();
		histUpDecay->Scale(1./histUpDecay->Integral(histUpDecay->FindBin(lowerLimit), histUpDecay->FindBin(upperLimit)), "width");
	 	partialFitResult = partialFit(histUpDecay, fitRange,  UpDecayName, fitOpts, drawPlots, histoDrawLists);
	   	fitResults.push_back(make_pair(UpDecayName, partialFitResult));

		//--------------------------------------DOWN DECAY FIT------------------------------------------
	 	TH1D* histDwnDecayTemp = (TH1D*)histDwnDecay->Clone();
	 	histDwnDecay->Scale(1./histDwnDecay->Integral(histDwnDecay->FindBin(lowerLimit), histDwnDecay->FindBin(upperLimit)), "width");
	 	partialFitResult = partialFit(histDwnDecay, fitRange,  DwnDecayName, fitOpts, drawPlots, histoDrawLists);
	   	fitResults.push_back(make_pair(DwnDecayName, partialFitResult));
	   	
		if (comb){ //perform combined fit and asymmetry only if a name for combined is specified
			if (histDwnDecayTemp->GetSumw2N() == 0) histDwnDecayTemp->Sumw2(kTRUE);
			if (histUpDecayTemp->GetSumw2N() == 0) histUpDecayTemp->Sumw2(kTRUE);
			
			TH1D* histCombined = new TH1D(*histUpDecayTemp); //initialize new histogram to the old one (same number of bins...)
			if(!histCombined->Add(&(*histUpDecayTemp), &(*histDwnDecayTemp))) cout<<"WARNING IN HISTOGRAM ADDITION"<<endl; //add the two
			histCombined->Scale(1./histCombined->Integral(histCombined->FindBin(lowerLimit), histCombined->FindBin(upperLimit)), "width");
	 		partialFitResult = partialFit<TH1D*>(histCombined, fitRange,  combinedDecayName, fitOpts, drawPlots, histoDrawLists);
	   		fitResults.push_back(make_pair(combinedDecayName, partialFitResult));

			//WARNING: ASYMETTRY IS ALWAYS BETWEEN THE NOT-NORMALIZED ONES
			TH1D* asymmetry = dynamic_cast<TH1D*>(histUpDecayTemp->GetAsymmetry(&(*histDwnDecayTemp))); //TH1D asymmetry = (*histUpDecayTemp - *histDwnDecayTemp)/histCombined;
			string OldFunc = this->GetFunction();
			this->SetFunction("pol1");
			partialFitResult = partialFit<TH1D*>(asymmetry, fitRange, "asymmetry", fitOpts, drawPlots, histoDrawLists);
	   		this->SetFunction(OldFunc);
	   		fitResults.push_back(make_pair("asymmetry_btw_not_normalized", partialFitResult));
		}
		

	}

   	if(drawPlots){
   		if ( MakeAndChangeCurrentDir(outDirNameInFile) ){
   		   	for( auto histoList : histoDrawLists ){
   		   		drawHistoList(histoList, (TString)outDecayCanvasName, outDecayCanvasDir);
   		   		delete histoList;
   		   	}
			TDirectory* dir = gDirectory;
			dir->Save();
		}
   	}
   	return fitResults;
}

void amulet::drawHistoList(TList* DrawList, TString CanvasName, string outDir){
	if ( DrawList == NULL ) return;
	bool isPULL = ((TString)DrawList->GetName()).Contains("PULL");
	int nToDraw = DrawList->GetEntries();
	if(isPULL) CanvasName = CanvasName + "_PULL";
	TCanvas* c = new TCanvas(CanvasName,CanvasName,0,0,floor(sqrt(nToDraw))*600,floor(sqrt(nToDraw))*400);
	c->DivideSquare(nToDraw);
	unsigned int ifor = 0;
	for(const auto&& obj: *DrawList){
		ifor++;
		c->cd(ifor);
		if (isPULL)
			((TH1D*)obj)->DrawCopy();
		else {
			((TH1D*)obj)->DrawCopy("E");
			TString name = (TString)obj->GetName();
			if(((TH1D*)obj)->GetListOfFunctions()->GetSize() != 1)
				cout<<"WARNING: there should be only one TF1 associated to the histogram: it will be drawn the first one"<<endl;
			TF1* fitFunc = (TF1*)(((TH1D*)obj)->GetListOfFunctions()->At(0));
      		auto legend = (!(name.Contains("asymmetry"))) ? new TLegend(0.42,0.335,0.8452,0.506) : new TLegend(0.116267,0.877079,0.506706,0.996506);
      		legend->SetLineColor(kWhite);
      		if(fitFunc!=NULL) legend->AddEntry(fitFunc,fitFunc->GetExpFormula(),"l");
      		legend->AddEntry(name,name,"le");
      		legend->DrawClone();
		}
	}
	c->Modified();
	c->Update();
	c->Write("", TObject::kOverwrite);
	//if(!isPULL)
	//	c->SaveAs((outDir+".pdf").c_str(), ".pdf");
}

bool amulet::drawThisOne( double* fitLims, vector<vector<double>> drawList ){
	bool to_draw = false;
	if (drawList.size()>0){	
		if(drawList.size()==1 && drawList[0].size()==1 && drawList[0][0] == 1)
			return true;
		double binWidth = (fitLims[1] - fitLims[0])/fitLims[2];
		for(unsigned int ii = 0; ii<drawList.size(); ii++){
			if(drawList[ii].size() != 3){
				cout<<"ERROR in drawing list: i wont draw anything"<<endl;
				return false;
			}
			if(_verbosity>2) cout<<"looking for lower = "<<drawList[ii][0]<<" upper = "<<drawList[ii][1]<<" nBins = "<<drawList[ii][2]<<endl;
			to_draw += abs(drawList[ii][0]- fitLims[0])<binWidth && abs(drawList[ii][1]- fitLims[1])<binWidth && (int)drawList[ii][2] == (int)fitLims[2];
		}
		if(to_draw && _verbosity>2) cout<<"DRAWING lower = "<<fitLims[0]<<" upper = "<<fitLims[1]<<" nBins = "<<fitLims[2]<<endl;
	}
	return to_draw;
}

void amulet::Stability_Analysis1DnBins(vector<string> names, vector<double> iter_settings,  vector<vector<double>> drawList, TString DrawOpts, bool drawErrorsOnPlot, bool plotOnlyValids){
	//prepare variables to cycle on
   	if(_verbosity > 0) std::cout<<"\nPerforming stability analysis 1D nBins..."<<endl;
   	int N_iterationsBin = (int)iter_settings[0];
	int nMinBin = (int)iter_settings[1];
	int nMaxBin = (int)iter_settings[2];
	double lowerLimit = iter_settings[3];
	double upperLimit = iter_settings[4];
	int dN = ( (nMaxBin-nMinBin)>=N_iterationsBin && N_iterationsBin>0 ) ? (nMaxBin-nMinBin)/N_iterationsBin : 1;
	N_iterationsBin = (nMaxBin-nMinBin)/dN;
	int nBins = nMinBin;
	nMinBin = nMinBin - dN;
	
	TCanvas* c0 = new TCanvas(); //need this to not overwrite previous canvases created in other methods
	c0->SetBatch();
	
	//first fit to get structure of fit results
	const char* fitOpts = "QRSE";//(drawErrorsOnPlot) ? "QRSE" : "QRS";
	double fitLimsGuess[3] = {lowerLimit, upperLimit, (double)nBins+dN*N_iterationsBin/2.};
	nestedVecWithName_t fitResults = Perform_Lifetime_Fit(fitLimsGuess,  names, false, fitOpts);
	delete c0; //this is not needed anymore

	//create histogram structure
	vector<TCanvas*> canvasesSameType;
	nestedVecHisto_t stabilityHistos;
    for (auto decayPlot : fitResults){
		unsigned int nQuant = decayPlot.second.size();
		TString canvName = decayPlot.first + " Stability Plot 1DnBins";
		TCanvas* myC0 = new TCanvas(canvName,canvName,200,10,ceil(sqrt(nQuant))*700,ceil(sqrt(nQuant))*500);
		myC0->DivideSquare(nQuant);
		vecHisto_t tempHistList;
		unsigned int jfor = 0; //unfortunately range-based for loops with initializer list are supported from c++20
		for (auto quantity : decayPlot.second){ //from c++20 you could do for(int jfor = 0; auto quantity : decayPlot.second)
			myC0->cd(jfor+1);
			TString histName = decayPlot.first + "_" + quantity.first+"_stability1DnBins";
			TString histTitle = decayPlot.first +" "+quantity.first+"; nBins;" + quantity.first;
			TH1D* tempHist = new TH1D(histName, histTitle, N_iterationsBin, nMinBin+dN, nMaxBin);
			//tempHist->SetDirectory(0);
			tempHistList.push_back(tempHist);
			jfor++;
		}
		canvasesSameType.push_back(myC0);
		stabilityHistos.push_back(tempHistList);
	}

	//to print progress percentage in main loop
	int stepPrint = 10;
	int nextPrint = stepPrint;
	int nItMax = (nMaxBin-nMinBin)/dN;
	vector<vector<int>> vjErr, vjPar;
	//MAIN LOOP
	for( int nb = 1; nb<=N_iterationsBin; nb++ ){
		float percent = 100.*(nb)/(float)nItMax;
    	if ((_verbosity > 0) && percent >= nextPrint){
			cout << "\r" << std::setprecision(3) <<percent << "% in 1D stability"<<endl;
			nextPrint += stepPrint;
    	}
		//populate fitLims for this iteration
		double fitLims[3] = {lowerLimit, upperLimit, (double)nMinBin+nb*dN};
		if(_verbosity>2) cout<<"processing nBins = "<<nMinBin+nb*dN<<" for fit range ["<<lowerLimit<<","<<upperLimit<<"]"<<endl;
		
		//set canvas and histogram names for each try
		vector<string> namesN = names;
		namesN[0] = namesN[0]+"/partialPlots";
		namesN[1] = namesN[1]+std::to_string(nMinBin+nb*dN);
	  	for(unsigned int i = 3; i<names.size(); i++) 
			namesN[i] = namesN[i]+std::to_string(nMinBin+nb*dN);
		
		//decide wether to draw based on what specified in drawList
		bool to_draw = drawThisOne(fitLims, drawList);
		
		//perform fit for this fitLims
		fitResults = Perform_Lifetime_Fit(fitLims,  namesN, to_draw, fitOpts);
		//fill histograms
		if (drawErrorsOnPlot) vjErr.clear(), vjPar.clear();
		unsigned int ifor = 0; //needs this before c++20
		for(auto decayPlot : fitResults){
			vector<int> jErr, jPar;
			bool fitValid = (decayPlot.second[0].second == true && decayPlot.second[1].second == 0 );
			unsigned int jfor = 0; //needs this before c++20
			for (auto quantity : decayPlot.second){
				if(quantity.first == "fit validity" || quantity.first =="fit status" || fitValid){
					stabilityHistos[ifor][jfor]->SetBinContent(nb, quantity.second);
					stabilityHistos[ifor][jfor]->SetBinError(nb, 0);
				}
				if (drawErrorsOnPlot){
					if( ((TString)quantity.first).BeginsWith("err"))
						jErr.push_back(jfor);
					else if ( ((TString)quantity.first).BeginsWith("par"))
						jPar.push_back(jfor);
				}
				jfor++;
			}
			if (drawErrorsOnPlot){
				vjErr.push_back(jErr);
				vjPar.push_back(jPar);
			}
			ifor++;
		}
		//fill with errors
		if (drawErrorsOnPlot){
			unsigned int ifor = 0;
			for(auto decayPlot : fitResults){
				vector<int> jErr=vjErr[ifor], jPar=vjPar[ifor];
				if (jErr.size() == jPar.size()){
					bool fitValid = (decayPlot.second[0].second == true && decayPlot.second[1].second == 0 );
					unsigned int jfor = 0; 
					for (int j : jPar){
						if(fitValid)
							stabilityHistos[ifor][j]->SetBinError(nb, decayPlot.second[jErr[jfor]].second);
						jfor++;
					}
					ifor++;
				}
			}
		}
	}
	if(_verbosity > 0) cout<<"stability analysis complete"<<endl;
	const char* outDirNameInFile = names[0].c_str();
	drawNestedVecHistos(stabilityHistos, canvasesSameType, outDirNameInFile, DrawOpts);
}

void amulet::Stability_Analysis2DnBinsSymRange(vector<string> names, vector<double> iter_settings,  vector<vector<double>> drawList, TString DrawOpts, bool drawErrorsOnPlot, bool plotOnlyValids){
	//prepare variables to cycle on
   	if(_verbosity > 0) std::cout<<"\nPerforming stability analysis 2D nBins and symmetric range..."<<endl;
   	int N_iterationsBin = (int)iter_settings[0];
	int nMinBin = (int)iter_settings[1];
	int nMaxBin = (int)iter_settings[2];
	int dN = ( (nMaxBin-nMinBin)>=N_iterationsBin && N_iterationsBin>0 ) ? (nMaxBin-nMinBin)/N_iterationsBin : 1;
	N_iterationsBin = (nMaxBin-nMinBin)/dN;
	int nBins = nMinBin;
	nMinBin = nMinBin - dN;

	int N_iterationsRange = (int)iter_settings[3];
	double lowLim = iter_settings[4];
	double upLim = iter_settings[5];
	double centerLim = lowLim + (upLim - lowLim)/2.;
	double Rmin = iter_settings[6]; //start from this and go up
	double Rmax = upLim - lowLim;
	double dR = (Rmax-Rmin)/N_iterationsRange;
	Rmin = Rmin - dR;
	//double range = Rmin;
	
	TCanvas* c0 = new TCanvas(); //need this to not overwrite previous canvases created in other methods
	c0->SetBatch();
	
	//first fit to get structure of fit results
	double fitLimsGuess[3] = {lowLim, upLim, (double)nBins+dN*N_iterationsBin/2.};
	const char* fitOpts = "QRSE";//(drawErrorsOnPlot) ? "QRSE" : "QRS";
	nestedVecWithName_t fitResults = Perform_Lifetime_Fit(fitLimsGuess,  names, false, fitOpts);
	delete c0; //this is not needed anymore

	//create histogram structure
	vector<TCanvas*> canvasesSameType;
	vector<vector<TH2D*>> stabilityHistos;
    for (auto decayPlot : fitResults){
		unsigned int nQuant = decayPlot.second.size();
		TString canvName = decayPlot.first + " Stability Plot 2DnBinsSymRange";
		TCanvas* myC0 = new TCanvas(canvName,canvName,200,10,ceil(sqrt(nQuant))*700,ceil(sqrt(nQuant))*500);
		myC0->DivideSquare(nQuant);
		vector<TH2D*> tempHistList;
		unsigned jfor = 0; 
		for (auto quantity : decayPlot.second){
			myC0->cd(jfor+1);
			TString histName = decayPlot.first + "_" + quantity.first+"_stability2DnBinsSymRange";
			TString histTitle = decayPlot.first +" "+quantity.first+"; nBins; range;" + quantity.first;
			TH2D* tempHist = new TH2D(histName, histTitle, N_iterationsBin, nMinBin+dN, nMaxBin, N_iterationsRange, Rmin+dR, Rmax);
			//tempHist->SetDirectory(0);
			tempHistList.push_back(tempHist);
			jfor++;
		}
		canvasesSameType.push_back(myC0);
		stabilityHistos.push_back(tempHistList);
	}
	//to cout progress percentage in main loop
	float nIters = 0;
	const float stepPrint = 1.;
	float nextPrint = stepPrint;
	const float nItMax = N_iterationsBin*N_iterationsRange;
	vector<vector<int>> vjErr, vjPar;
	//main loop
	for( int nb = 1; nb<=N_iterationsBin; nb++ ){
		for( int nr = 1; nr<=N_iterationsRange; nr++ ){
			float percent = 100.*nIters/nItMax;
    		if ((_verbosity > 0) && percent >= nextPrint){
				cout << "\r" << std::setprecision(3) << percent << "% in 2D stability"<<endl;
				nextPrint += stepPrint;
    		}
			//populate fitLims for this iteration
			double lowLim = centerLim - (Rmin+nr*dR)/2.;
			double upLim = centerLim + (Rmin+nr*dR)/2.;
			double fitLims[3] = {lowLim, upLim, (double)nMinBin+nb*dN};
			if(_verbosity>2) cout<<"processing nBins = "<<nMinBin+nb*dN<<" for fit range ["<<lowLim<<","<<upLim<<"]"<<endl;
			
			//set canvas and histogram names for each try
			vector<string> namesN = names;
			namesN[0] = namesN[0]+"/partialPlots";
			namesN[1] = namesN[1]+std::to_string(nMinBin+nb*dN)+std::to_string(round(1e6*lowLim))+std::to_string(round(1e6*upLim));
			for(unsigned int i = 3; i<names.size(); i++) 
				namesN[i] = namesN[i]+std::to_string(nMinBin+nb*dN)+std::to_string(round(1e6*lowLim))+std::to_string(round(1e6*upLim));
			
			//decide wether to draw based on what specified in drawList
			bool to_draw = drawThisOne(fitLims, drawList);
			
			//perform fit for this fitLims
			fitResults = Perform_Lifetime_Fit(fitLims,  namesN, to_draw, fitOpts);
			
			//fill histograms
			if (drawErrorsOnPlot) vjErr.clear(), vjPar.clear();
			unsigned int ifor = 0; 
			for(auto decayPlot : fitResults){
				vector<int> jErr, jPar;
				bool fitValid = (plotOnlyValids) ? (decayPlot.second[0].second == true && decayPlot.second[1].second == 0 ) : true;
				unsigned int jfor = 0; 
				for (auto quantity : decayPlot.second){
					if(quantity.first == "fit validity" || quantity.first =="fit status" || fitValid){
						stabilityHistos[ifor][jfor]->SetBinContent(nb, nr, quantity.second);
						stabilityHistos[ifor][jfor]->SetBinError(nb, nr, 0);
					}
					if(drawErrorsOnPlot){
						if( ((TString)quantity.first).BeginsWith("err"))
							jErr.push_back(jfor);
						else if ( ((TString)quantity.first).BeginsWith("par"))
							jPar.push_back(jfor);
					}
					jfor++;
				}
				if(drawErrorsOnPlot){
					vjErr.push_back(jErr);
					vjPar.push_back(jPar);
				}
				ifor++;
			}
			//fill with errors
			if(drawErrorsOnPlot){
				unsigned int ifor = 0;
				for(auto decayPlot : fitResults){
					vector<int> jErr=vjErr[ifor], jPar=vjPar[ifor];
					if (jErr.size() == jPar.size()){
						bool fitValid = (decayPlot.second[0].second == true && decayPlot.second[1].second == 0 );
						unsigned int jfor = 0; 
						for (int j : jPar){
							if(fitValid)
								stabilityHistos[ifor][j]->SetBinError(nb, nr, decayPlot.second[jErr[jfor]].second);
							jfor++;
						}
						ifor++;
					}
				}
			}
			nIters++;
		}
	}
	if(_verbosity > 0) cout<<"stability analysis complete"<<endl;
	
	TString drawOpts = DrawOpts;
	const char* outDirNameInFile = names[0].c_str();
	drawNestedVecHistos(stabilityHistos, canvasesSameType, outDirNameInFile, DrawOpts);
}

void amulet::Stability_Analysis3DnBinsLowUpLims(vector<string> names, vector<double> iter_settings,  vector<vector<double>> drawList, TString DrawOpts, bool drawErrorsOnPlot, bool plotOnlyValids){
	//prepare variables to cycle on
   	if(_verbosity > 0) std::cout<<"\nPerforming stability analysis 3D nBins and lowLim and upLim..."<<endl;
   	int N_iterationsBin = (int)iter_settings[0];
	int nMinBin = (int)iter_settings[1];
	int nMaxBin = (int)iter_settings[2];
	int dN = ( (nMaxBin-nMinBin)>=N_iterationsBin && N_iterationsBin>0 ) ? (nMaxBin-nMinBin)/N_iterationsBin : 1;
	N_iterationsBin = ((nMaxBin-nMinBin)/dN);
	int nBins = nMinBin;
	nMinBin = nMinBin - dN;

	int N_iterationsLowLim = (int)iter_settings[3];
	double lowLimMin = iter_settings[4];
	double lowLimMax = iter_settings[5];
	double dLL = (lowLimMax-lowLimMin)/(double)N_iterationsLowLim;
	double lowLim = lowLimMin;
	lowLim = lowLim - dLL;

	int N_iterationsUpLim = (int)iter_settings[6];
	double upLimMin = iter_settings[7];
	double upLimMax = iter_settings[8];
	double dUL = (upLimMax-upLimMin)/(double)N_iterationsUpLim;
	double upLim = upLimMin;
	upLim = upLim - dUL;
	
	TCanvas* c0 = new TCanvas(); //need this to not overwrite previous canvases created in other methods
	c0->SetBatch();
	
	//first fit to get structure of fit results
	double fitLimsGuess[3] = {lowLim+dLL*N_iterationsLowLim/2., upLim+dUL*N_iterationsUpLim/2., (double)nBins+dN*N_iterationsBin/2.};
	const char* fitOpts = "QRSE";//(drawErrorsOnPlot) ? "QRSE" : "QRS";
	nestedVecWithName_t fitResults = Perform_Lifetime_Fit(fitLimsGuess,  names, false, fitOpts);
	delete c0; //this is not needed anymore

	//create histogram structure
	vector<TCanvas*> canvasesSameType;
	vector<vector<TH3D*>> stabilityHistos;
    for (auto decayPlot : fitResults){
		unsigned int nQuant = decayPlot.second.size();
		TString canvName = decayPlot.first + " Stability Plot 3D";
		TCanvas* myC0 = new TCanvas(canvName,canvName,200,10,ceil(sqrt(nQuant))*700,ceil(sqrt(nQuant))*500);
		myC0->DivideSquare(nQuant);
		vector<TH3D*> tempHistList;
		unsigned int jfor = 0; 
		for (auto quantity : decayPlot.second){
			myC0->cd(jfor+1);
			TString histName = decayPlot.first + "_" + quantity.first+"_stability3D";
			TString histTitle = decayPlot.first +" "+quantity.first+";nBins ;lowLim ;upLim";
			TH3D* tempHist = new TH3D(histName, histTitle, N_iterationsBin, nMinBin+dN, nMaxBin, N_iterationsLowLim, lowLimMin+dLL, lowLimMax, N_iterationsUpLim, upLimMin+dUL, upLimMax);
			//tempHist->SetDirectory(0);
			tempHistList.push_back(tempHist);
			jfor++;
		}
		canvasesSameType.push_back(myC0);
		stabilityHistos.push_back(tempHistList);
	}
	//to cout progress percentage in main loop
	int nIters = 0;
	float stepPrint = 1; //0.1;
	float nextPrint = stepPrint;
	float nItMax = N_iterationsBin*N_iterationsLowLim*N_iterationsUpLim;
	vector<vector<int>> vjErr, vjPar;
	//main loop
	for( int nb = 1; nb<=N_iterationsBin; nb++ ){
		//lowLim = lowLimMin + dLL;
		for( int nll = 1; nll<=N_iterationsLowLim; nll++ ){
			for( int nul = 1; nul<=N_iterationsUpLim; nul++ ){
				float percent = 100.*nIters/nItMax;
				if ((_verbosity > 0) && percent >= nextPrint){
					cout << "\r" << std::setprecision(3) << percent << "% in 3D stability"<<endl;
					nextPrint += stepPrint;
    			}
				double fitLims[3] = {lowLimMin+nll*dLL, upLimMin+nul*dUL, (double)nMinBin+nb*dN};
				if(_verbosity>2) cout<<"processing nBins = "<<nMinBin+nb*dN<<" for fit range ["<<(lowLimMin+nll*dLL)<<","<<(upLimMin+nul*dUL)<<"]"<<endl;
				
				//set canvas and histogram names for each try
				vector<string> namesN = names;
				namesN[0] = namesN[0]+"/partialPlots";
				namesN[1] = namesN[1]+std::to_string(nMinBin+nb*dN)+std::to_string(round(1e6*(lowLimMin+nll*dLL)))+std::to_string(round(1e6*(upLimMin+nul*dUL)));
				for(unsigned int i = 3; i<names.size(); i++)
					namesN[i] = namesN[i]+std::to_string(nMinBin+nb*dN)+std::to_string(round(1e6*(lowLimMin+nll*dLL)))+std::to_string(round(1e6*(upLimMin+nul*dUL)));
				
				//decide wether to draw based on what specified in drawList
				bool to_draw = drawThisOne(fitLims, drawList);
				
				//perform fit for this fitLims
				fitResults = Perform_Lifetime_Fit(fitLims,  namesN, to_draw, fitOpts);
				
				//fill histograms
				if (drawErrorsOnPlot) vjErr.clear(), vjPar.clear();
				unsigned int ifor = 0; 
				for(auto decayPlot : fitResults){
					vector<int> jErr, jPar;
					bool fitValid = (plotOnlyValids) ? (decayPlot.second[0].second == true && decayPlot.second[1].second == 0 ) : true;
					unsigned int jfor = 0; 
					for (auto quantity : decayPlot.second){
						if(quantity.first == "fit validity" || quantity.first =="fit status" || fitValid){
							stabilityHistos[ifor][jfor]->SetBinContent(nb, nll, nul, quantity.second);
							stabilityHistos[ifor][jfor]->SetBinError(nb, nll, nul, 0);
							//cout<<stabilityHistos[ifor][jfor]->GetXaxis()->FindBin(nBins)<<" "<<nb<<endl;
							//cout<<stabilityHistos[ifor][jfor]->GetYaxis()->FindBin(lowLim)<<" "<<nll<<endl; //FIXME something strange happens here because these two numbers are not equal
							//cout<<stabilityHistos[ifor][jfor]->GetYaxis()->GetBinCenter(stabilityHistos[ifor][jfor]->GetYaxis()->FindBin(lowLim))<<" "<<stabilityHistos[ifor][jfor]->GetYaxis()->GetBinCenter(nll)<<endl;
						}
						if (drawErrorsOnPlot){
							if( ((TString)quantity.first).BeginsWith("err"))
								jErr.push_back(jfor);
							else if ( ((TString)quantity.first).BeginsWith("par"))
								jPar.push_back(jfor);
						}
						jfor++;
					}
					if (drawErrorsOnPlot){
						vjErr.push_back(jErr);
						vjPar.push_back(jPar);
					}
					ifor++;
				}
				//fill with errors
				if (drawErrorsOnPlot){
					unsigned int ifor = 0; 
					for(auto decayPlot : fitResults){
						vector<int> jErr=vjErr[ifor], jPar=vjPar[ifor];
						if (jErr.size() == jPar.size()){
							bool fitValid = (decayPlot.second[0].second == true && decayPlot.second[1].second == 0 );
							unsigned int jfor = 0; 
							for (int j : jPar){
								if(fitValid)
									stabilityHistos[ifor][j]->SetBinError(nb, nll, nul, decayPlot.second[jErr[jfor]].second);
								jfor++;
							}
							ifor++;
						}
					}
				}
				nIters++;
				//upLim += dUL;
			}
			//cout<<lowLim<<" "<<lowLimMin+nll*dLL<<endl;
			//lowLim += dLL;
		}
		//nBins += dN;
	}
	if(_verbosity > 0) cout<<"stability analysis complete"<<endl;
	
	TString drawOpts = DrawOpts;
	const char* outDirNameInFile = names[0].c_str();
	drawNestedVecHistos(stabilityHistos, canvasesSameType, outDirNameInFile, DrawOpts);
}

template<class T = nestedVecHisto_t>
void amulet::drawNestedVecHistos ( T histoMatrix, vector<TCanvas*> canvases, const char* outDirNameInFile, TString DrawOpts ){
	bool dirOK = MakeAndChangeCurrentDir(outDirNameInFile);
	unsigned int ifor = 0; 
	for(auto ROW : histoMatrix ){
		canvases[ifor]->cd();
		unsigned int jfor = 0; 
		for (auto COL : ROW ){
			canvases[ifor]->cd(jfor+1);
			if(((TString)COL->ClassName()).Contains("1")){
				COL->SetMarkerSize(1);
				COL->SetMarkerStyle(33);
			}
			COL->DrawCopy(DrawOpts);
			jfor++;
		}
		canvases[ifor]->Modified();
		canvases[ifor]->Update();
		if(dirOK){
			canvases[ifor]->Write("", TObject::kOverwrite);
			TDirectory* dir = gDirectory;
			dir->Save();
		}
		ifor++;
	}
}

void amulet::SignalWidth_Analysis(std::vector<double> binnings, const char* outDirNameInFile){
	TString UpWidthSTART  = _OutFileName + "UpSTARTWidth" ,
		    UpWidthSTOP   = _OutFileName + "UpSTOPWidth"  ,
		    DwnWidthSTART = _OutFileName + "DwnSTARTWidth",
		    DwnWidthSTOP  = _OutFileName + "DwnSTOPWidth" ;

	int     nBinsSTARTWidth     = (int)binnings[0];
	double lowerRangeSTARTWidth =      binnings[1];
	double upperRangeSTARTWidth =      binnings[2];
	int    nBinsSTOPWidth       = (int)binnings[3];
	double lowerRangeSTOPWidth  =      binnings[4];
	double upperRangeSTOPWidth  =      binnings[5];
	
	if ( MakeAndChangeCurrentDir(outDirNameInFile) ){	
		(_dfUp.Histo1D<double> ({UpWidthSTART , UpWidthSTART , nBinsSTARTWidth, lowerRangeSTARTWidth, upperRangeSTARTWidth}, "UpwdtSTART") ) -> Write("", TObject::kOverwrite);
		(_dfUp.Histo1D<double> ({UpWidthSTOP  , UpWidthSTOP  , nBinsSTOPWidth , lowerRangeSTOPWidth , upperRangeSTOPWidth }, "UpwdtSTOP")  ) -> Write("", TObject::kOverwrite);
		(_dfDwn.Histo1D<double>({DwnWidthSTART, DwnWidthSTART, nBinsSTARTWidth, lowerRangeSTARTWidth, upperRangeSTARTWidth}, "DwnwdtSTART")) -> Write("", TObject::kOverwrite);
		(_dfDwn.Histo1D<double>({DwnWidthSTOP , DwnWidthSTOP , nBinsSTOPWidth , lowerRangeSTOPWidth , upperRangeSTOPWidth }, "DwnwdtSTOP") ) -> Write("", TObject::kOverwrite);
		if(_verbosity > 0) cout<<"SignalWidth_Analysis plots have been saved"<<endl;
		TDirectory* dir = gDirectory;
		dir->Save();
	}
}

void amulet::StartAndStop_analysis(std::vector<double> binnings, const char* outDirNameInFile){
	TString UpDecaySTART  = _OutFileName + "UpDecayMidpointSTART" ,
		    UpDecaySTOP   = _OutFileName + "UpDecayMidpointSTOP"  ,
		    DwnDecaySTART = _OutFileName + "DwnDecayMidpointSTART",
		    DwnDecaySTOP  = _OutFileName + "DwnDecayMidpointSTOP" ;
		   
	int    nBinsSTART      = (int)binnings[0];
	double lowerRangeSTART =      binnings[1];
	double upperRangeSTART =      binnings[2];
	int    nBinsSTOP       = (int)binnings[3];
	double lowerRangeSTOP  =      binnings[4];
	double upperRangeSTOP  =      binnings[5];
	
	if ( MakeAndChangeCurrentDir(outDirNameInFile) ){
		(_dfUp.Histo1D<double> ({UpDecaySTART , UpDecaySTART ,  nBinsSTART, lowerRangeSTART, upperRangeSTART}, "UpSTARTSquareMidPos" )) -> Write("", TObject::kOverwrite);
		(_dfUp.Histo1D<double> ({UpDecaySTOP  , UpDecaySTOP  ,  nBinsSTOP , lowerRangeSTOP , upperRangeSTOP }, "UpSTOPSquareMidPos"  )) -> Write("", TObject::kOverwrite);
		(_dfDwn.Histo1D<double>({DwnDecaySTART, DwnDecaySTART,  nBinsSTART, lowerRangeSTART, upperRangeSTART}, "DwnSTARTSquareMidPos")) -> Write("", TObject::kOverwrite);
		(_dfDwn.Histo1D<double>({DwnDecaySTOP , DwnDecaySTOP ,  nBinsSTOP , lowerRangeSTOP , upperRangeSTOP }, "DwnSTOPSquareMidPos" )) -> Write("", TObject::kOverwrite);
		if(_verbosity > 0) cout<<"StartAndStop_analysis plots have been saved"<<endl;
		TDirectory* dir = gDirectory;
		dir->Save();
	}
}

Int_t amulet::countpads(TVirtualPad *pad) {
   //count the number of pads in pad
   if (!pad) return 0;
   Int_t npads = 0;
   TObject *obj;
   TIter next(pad->GetListOfPrimitives());
   while ((obj = next())) {
      if (obj->InheritsFrom(TVirtualPad::Class())) npads++;
   }
   return npads;
}
