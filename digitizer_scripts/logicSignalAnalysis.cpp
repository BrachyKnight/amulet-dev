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

#include <ROOT/RDataFrame.hxx>
#include <TString.h>
#include <TPRegexp.h>

typedef std::tuple<ROOT::VecOps::RVec<double>,ROOT::VecOps::RVec<double>,ROOT::VecOps::RVec<double>,ROOT::VecOps::RVec<double>> tuple4RVec_t;

using namespace std;

std::pair<double,double> Center_and_Err( ROOT::VecOps::RVec<double> vec)
{
	return make_pair(Min(vec) + (Max(vec)-Min(vec))/2., (Max(vec)-Min(vec))/2.);
}

tuple4RVec_t analyze ( ROOT::VecOps::RVec<int> v, double sampfreq ) 
{
	unsigned int nPts = v.size();
	double vminTh = 0.2*Max(v), vmaxTh = 0.7*Max(v);  //is considered only if between 20 and 70 percent of height
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
		if( (auxv.size() == auxt.size()) && auxt.size()!=0 )
		{
			if(auxv[0]>auxv[auxv.size()-1]) //check if first point is greater then last point (means wvf goes down)
				timedwns.push_back(auxt);
			else if (auxv[0]<auxv[auxv.size()-1])  //check if first point is less then last point (means wvf goes up)
				timeups.push_back(auxt);			
			else cout<<"error analyze"<<endl;
		}//TODO add a check to find the ones that jumps from under threshold to above threshold and therefore are not saved		
	}
	ROOT::VecOps::RVec<double> tup_center, tdwn_center, tup_cerr, tdwn_cerr; 
	for( unsigned int i = 0; i<timedwns.size(); i++ )
	{
		auto dwnstCE = Center_and_Err(timedwns[i]);
		tdwn_center.push_back( dwnstCE.first ), tup_cerr.push_back( dwnstCE.second );
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
	
	auto t0 = std::chrono::high_resolution_clock::now(); //to evaluate execution time
	if ( argc != 3 ){
		std::cout<<"\nINPUT ERROR"<<endl
		<<"needs:\noutputDir \nrootInputFile\n"<<endl
		<<"as an example \n./executable path/to/output/folder/ path/to/file/fileName.root"<<endl;
		return 1;
	}
  
	//set I/O file names
	string 	outDir  = argv[1],
		rootIn  = argv[2],
		rootOut = rootIn; //update existing TTree with new branches (columns)
	
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
			df_wvf.Snapshot("amulet_aux",rootOut.c_str(), {"idx","measN","ch0_wvf_amp","ch0_wvf_time","ch1_wvf_amp","ch1_wvf_time"}, opt);
		}else
			throw std::runtime_error(e.what());
	}

	//reopen new tree
	const char* name = (alreadyExecuted) ? "amulet_aux" : "amulet";
	ROOT::RDataFrame df_dec(name, rootIn);
	//define new columns with the decoded informations and create new tree
	df_dec	.Define("ch0Decoded"    , [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return analyze(ADC_ch, sampfreq);     }, {"ch0_wvf_amp"} )
		.Define("ch1Decoded"    , [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return analyze(ADC_ch, sampfreq);     }, {"ch1_wvf_amp"} )
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
