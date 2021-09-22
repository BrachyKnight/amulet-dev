/* 
 * Author: Massimo Girola
 * Date: June 2021
 * Developed using ROOT 6.24/00
 * c++ -Wall -o logicSignalAnalysis logicSignalAnalysis.cpp `root-config --cflags --glibs`
 */
 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <filesystem>

#include <ROOT/RDataFrame.hxx>
#include <TString.h>
#include <TPRegexp.h>
//STL usage
using std::cout, std::endl, std::make_pair, std::string, std::make_tuple, std::get, std::runtime_error;

//define custom type
typedef std::tuple<ROOT::VecOps::RVec<double>,ROOT::VecOps::RVec<double>,ROOT::VecOps::RVec<double>,ROOT::VecOps::RVec<double>> tuple4RVec_t;

std::pair<double,double> Center_and_Err( ROOT::VecOps::RVec<double> vec)
{
	return make_pair(Min(vec) + (Max(vec)-Min(vec))/2., (Max(vec)-Min(vec))/2.);
}

tuple4RVec_t analyze ( ROOT::VecOps::RVec<int> v, double sampfreq, double minTh, double maxTh ) 
{
	unsigned int nPts = v.size();
	if( nPts == 0 ){ //if the channel is empty simply return vectors with 0
		ROOT::VecOps::RVec<double> nullvec = {0};
		return make_tuple(nullvec,nullvec,nullvec,nullvec);
	}
	double deltaV = abs(Max(v)-Min(v));
	double vminTh = minTh*deltaV+Min(v), vmaxTh = maxTh*deltaV+Min(v);  //is considered only if between minTh and maxTh percent of pulse amplitude
	ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>> ups, dwns, timeups, timedwns;
	for( unsigned int i = 0; i<nPts; i++ )
	{
		ROOT::VecOps::RVec<double> auxv, auxt;
		bool is_btw = (v[i]<=vmaxTh)&&(v[i]>=vminTh);
		int j = i;
		while( is_btw )
		{
			auxv.push_back(v[j]);
			auxt.push_back(static_cast<double>(j)/sampfreq);
			is_btw = (v[j]<=vmaxTh)&&(v[j]>=vminTh);
			j++;
			i = j;
		}
		if( (auxv.size() == auxt.size()) && auxt.size()!=0 ) //cases where at least one point falls between the min threshold and the max threshold
		{
			if(auxv[0]>auxv[auxv.size()-1]) //check if first point is greater then last point (means wvf goes down)
				timedwns.push_back(auxt);
			else if (auxv[0]<auxv[auxv.size()-1])  //check if first point is less then last point (means wvf goes up)
				timeups.push_back(auxt);			
			else cout<<"error analyze"<<endl;
		}else if( (auxv.size() == auxt.size()) && auxt.size()==0 ){ //cases where signals jump straight from under min threshold to over max threshold
			if( i == 0 || i == v.size()-1) continue;
			if( v[i-1] < vminTh && v[i] < v[i+1] && v[i]>vmaxTh){ //jump from over max threshold value to under min threshold value
				auxt.push_back(static_cast<double>(i)/sampfreq);
				timeups.push_back(auxt);
			}else if( v[i-1] > vmaxTh && v[i] > v[i+1] && v[i]<vminTh ){ //jump from under min threshold value to over max threshold value
				auxt.push_back(static_cast<double>(i)/sampfreq);
				timedwns.push_back(auxt);
			}
		}else
			throw runtime_error("error in analyze() function");		
	}
	ROOT::VecOps::RVec<double> tup_center, tdwn_center, tup_cerr, tdwn_cerr; 
	for( unsigned int i = 0; i<timedwns.size(); i++ )
	{
		auto dwnstCE = Center_and_Err(timedwns[i]);
		tdwn_center.push_back( dwnstCE.first ), tdwn_cerr.push_back( dwnstCE.second );
	}
	for( unsigned int i = 0; i<timeups.size(); i++ )
	{
		auto upstCE = Center_and_Err(timeups[i]);
		tup_center.push_back( upstCE.first ), tup_cerr.push_back( upstCE.second );
	}	
	return make_tuple(tup_center, tdwn_center, tup_cerr, tdwn_cerr);
}

ROOT::VecOps::RVec<float> temporalize (ROOT::VecOps::RVec<int> ADC_ch, double sampfreq){
	ROOT::VecOps::RVec<float> aux;
	for(unsigned int i = 0; i < ADC_ch.size(); i++)
		aux.push_back(static_cast<float>(i)/(sampfreq));
	return aux;
}

int main(int argc, char** argv)
{
	double sampfreq = 1e9; //SET THE SAMPLE FREQUENCY OF YOUR DIGITIZER TO GET CORRECT RESULTS
	cout<<endl<<"ASSUMING DIGITIZER SAMPLE FREQUENCY EQUAL TO "<<sampfreq<<" Hz"<<endl;
	
	auto t0 = std::chrono::high_resolution_clock::now(); //to evaluate execution time
	if ( argc != 2 && argc != 4 ){
		std::cout<<"\nINPUT ERROR: argc = "<<argc<<endl
		<<"needs:\noutputDir \nrootInputFile\n"<<endl
		<<"optional arguments threshold values"<<endl
		<<"as an example \n./executable path/to/file/fileName.root minTh(optional) maxTh(optional)"<<endl;
		return 1;
	}
	double minTh=25, maxTh=75;
	bool ThOptimization = false;
	if ( argc == 4 ){
		minTh = std::stod(argv[2]);
		maxTh = std::stod(argv[3]);
		ThOptimization = true;
	}
	if(!(minTh<=100 && maxTh<=100 && minTh>=0 && maxTh>=0 && minTh<=maxTh))
		throw std::runtime_error("THRESHOLD VALUES ARE \% OF PULSE AMPLITUDE SO THEY MUST BE IN RANGE [0,100]");
	minTh/=100; maxTh/=100;
	cout<<"Using threshold values: \n minTh = "<<minTh<<" maxTh = "<<maxTh<<endl;
  
	//set I/O file names
	string	rootIn  = argv[1],
		rootOut = rootIn, //update existing TTree with new branches (columns)
		outDir = std::filesystem::path(rootIn).parent_path().u8string(); //c++17 way to extract path from filename

	//use all cores of your machine MultiThreading to speed up analysis (read RDataFrame documentation for MT warnings)
	ROOT::EnableImplicitMT();
	
	//tell snapshot to update root file
	ROOT::RDF::RSnapshotOptions opt;
	opt.fMode="UPDATE";
	opt.fOverwriteIfExists=true;
	
	//open TTree amulet using df
	ROOT::RDataFrame df_wvf("amulet", rootIn);
	//try except syntax used to detect if idx and temporalized column have already been evaluated	
	bool alreadyExecuted = false;
	try{
		//build an index to manage reshuffling on entries caused by MT and attribute meas number with a trick using regex
		df_wvf  .Define("idx","rdfentry_")
			.Define("measN", [rootIn](){ auto myRegex=TPRegexp("[1-9][0-9]*|0"); return static_cast<short int>(std::stoi(((TString)rootIn)(myRegex,((TString)rootIn).Index("meas")))); })
			.Define("ch0_wvf_time"	, [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return temporalize(ADC_ch, sampfreq); }, {"ch0_wvf_amp"} ) //add time point to wvf ch0
			.Define("ch1_wvf_time"	, [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return temporalize(ADC_ch, sampfreq); }, {"ch1_wvf_amp"} ) //add time point to wvf ch1
			.Snapshot("amulet_wvf",rootOut.c_str(), {"idx","measN","ch0_wvf_amp","ch0_wvf_time","ch1_wvf_amp","ch1_wvf_time"}, opt);
		//open file and delete old tree then rename the existing one (needed to be able to execute this script multiple times)	
		TFile myF(rootOut.c_str(),"UPDATE");
		gDirectory->Delete("amulet;*");
		auto *myT = myF.Get<TTree>("amulet_wvf");
		myT->SetNameTitle("amulet","amulet");
		gDirectory->Delete("amulet_wvf;*");
		myT->Write();
		myF.Close();
	}catch (std::exception& e){
		//if already executed we need to trick root into thinking that there is no friend tree
		//so that later it can redefine columns already present in the friend tree
		if( ((TString)e.what()).Contains("already present in TTree") ){
			alreadyExecuted = true;
			if(ThOptimization)
				df_wvf.Snapshot("amulet_aux",rootOut.c_str(), {"idx","measN","ch0_wvf_amp","ch1_wvf_amp"}, opt);
			else
				df_wvf.Snapshot("amulet_aux",rootOut.c_str(), {"idx","measN","ch0_wvf_amp","ch0_wvf_time","ch1_wvf_amp","ch1_wvf_time"}, opt);
		}else
			throw std::runtime_error(e.what());
	}

	//reopen new tree
	const char* name = (alreadyExecuted) ? "amulet_aux" : "amulet";
	ROOT::RDataFrame df_dec(name, rootIn);
	//define new columns with the decoded informations and create new tree
	df_dec	.Define("ch0Decoded"    , [sampfreq,minTh,maxTh]( ROOT::VecOps::RVec<int> ADC_ch ){ return analyze(ADC_ch, sampfreq, minTh, maxTh);     }, {"ch0_wvf_amp"} )
		.Define("ch1Decoded"    , [sampfreq,minTh,maxTh]( ROOT::VecOps::RVec<int> ADC_ch ){ return analyze(ADC_ch, sampfreq, minTh, maxTh);     }, {"ch1_wvf_amp"} )
		.Define("ch0Nup"        , []( tuple4RVec_t WVF_dec ){ return (int)get<0>(WVF_dec).size();  }, {"ch0Decoded" } )
		.Define("ch0Ndwn"       , []( tuple4RVec_t WVF_dec ){ return (int)get<1>(WVF_dec).size();  }, {"ch0Decoded" } ) 
		.Define("ch0timeups"    , []( tuple4RVec_t WVF_dec ){ return get<0>(WVF_dec);              }, {"ch0Decoded" } )
		.Define("ch0timedwns"   , []( tuple4RVec_t WVF_dec ){ return get<1>(WVF_dec);              }, {"ch0Decoded" } )
		.Define("ch0timeupsErr" , []( tuple4RVec_t WVF_dec ){ return get<2>(WVF_dec);              }, {"ch0Decoded" } )
		.Define("ch0timedwnsErr", []( tuple4RVec_t WVF_dec ){ return get<3>(WVF_dec);              }, {"ch0Decoded" } )
		.Define("ch1Nup"        , []( tuple4RVec_t WVF_dec ){ return (int)get<0>(WVF_dec).size();  }, {"ch1Decoded" } )
		.Define("ch1Ndwn"       , []( tuple4RVec_t WVF_dec ){ return (int)get<1>(WVF_dec).size();  }, {"ch1Decoded" } ) 
		.Define("ch1timeups"    , []( tuple4RVec_t WVF_dec ){ return get<0>(WVF_dec);              }, {"ch1Decoded" } )
		.Define("ch1timedwns"   , []( tuple4RVec_t WVF_dec ){ return get<1>(WVF_dec);              }, {"ch1Decoded" } )
		.Define("ch1timeupsErr" , []( tuple4RVec_t WVF_dec ){ return get<2>(WVF_dec);              }, {"ch1Decoded" } )
		.Define("ch1timedwnsErr", []( tuple4RVec_t WVF_dec ){ return get<3>(WVF_dec);              }, {"ch1Decoded" } )
		.Snapshot("decoded", rootOut.c_str(), { "idx",
							"ch0Nup","ch0Ndwn","ch0timeups","ch0timedwns","ch0timeupsErr","ch0timedwnsErr",
							"ch1Nup","ch1Ndwn","ch1timeups","ch1timedwns","ch1timeupsErr","ch1timedwnsErr" }, opt ); 
	

	// use BuildIndex to read the original tree and its (shuffled but indexed) friend
	TFile f(rootOut.c_str(),"UPDATE");
	auto *tmain = f.Get<TTree>("amulet");
	auto *tfriend = f.Get<TTree>("decoded");
	tfriend->BuildIndex("idx");
	tfriend->Write("",TObject::kOverwrite);
	tmain->AddFriend(tfriend);
	tmain->Write("",TObject::kOverwrite);
	if(alreadyExecuted) gDirectory->Delete("amulet_aux;*");
	f.Close();
	
	cout<<endl<<"\nDATAFRAME UPDATED IN FILE "<<rootOut.c_str()<<endl<<endl;
	
	auto t3 = std::chrono::high_resolution_clock::now();
	auto total_duration = std::chrono::duration_cast<std::chrono::seconds>( t3 - t0 ).count();
	std::cout<<"Total execution time "<<total_duration<<" seconds"<<endl<<endl;
	return 0;
}
