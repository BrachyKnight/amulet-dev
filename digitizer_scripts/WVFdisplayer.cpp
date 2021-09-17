/* 
 * Author: Massimo Girola
 * Date: June 2021
 * Developed using ROOT 6.24/00
 * c++ -Wall -o analogicWVFdisplayer analogicWVFdisplayer.cpp `root-config --cflags --glibs`
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
	if ( argc != 2 ){
		std::cout<<"\nINPUT ERROR"<<endl
		<<"needs:\noutputDir \nrootInputFile\n"<<endl
		<<"as an example \n./executable path/to/file/fileName.root"<<endl;
		return 1;
	}
  
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
			.Define("ch0_wvf_time"	, [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return temporalize(ADC_ch, sampfreq); }, {"ch0_wvf_amp"} ) //add time point to wvf ch0
			.Define("ch1_wvf_time"	, [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return temporalize(ADC_ch, sampfreq); }, {"ch1_wvf_amp"} ) //add time point to wvf ch1
			.Define("ch2_wvf_time"	, [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return temporalize(ADC_ch, sampfreq); }, {"ch2_wvf_amp"} ) //add time point to wvf ch2
			.Define("ch3_wvf_time"	, [sampfreq]( ROOT::VecOps::RVec<int> ADC_ch ){ return temporalize(ADC_ch, sampfreq); }, {"ch3_wvf_amp"} ) //add time point to wvf ch3
			.Snapshot("amulet_wvf",rootOut.c_str(), {"idx","ch0_wvf_amp","ch0_wvf_time","ch1_wvf_amp","ch1_wvf_time","ch2_wvf_amp","ch2_wvf_time","ch3_wvf_amp","ch3_wvf_time"}, opt);
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
			df_wvf.Snapshot("amulet_aux",rootOut.c_str(),  {"idx","ch0_wvf_amp","ch0_wvf_time","ch1_wvf_amp","ch1_wvf_time","ch2_wvf_amp","ch2_wvf_time","ch3_wvf_amp","ch3_wvf_time"}, opt);
		}else
			throw std::runtime_error(e.what());
	}

	TFile f(rootOut.c_str(),"UPDATE");
	auto *tmain = f.Get<TTree>("amulet");
	tmain->Write("",TObject::kOverwrite);
	if(alreadyExecuted) gDirectory->Delete("amulet_aux;*");
	f.Close();
	
	cout<<endl<<"\nDATAFRAME UPDATED IN FILE "<<rootOut.c_str()<<endl<<endl;
	
	auto t3 = std::chrono::high_resolution_clock::now();
	auto total_duration = std::chrono::duration_cast<std::chrono::seconds>( t3 - t0 ).count();
	std::cout<<"Total execution time "<<total_duration<<" seconds"<<endl<<endl;
	return 0;
}
