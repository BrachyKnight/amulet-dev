/* 
 * Author: Massimo Girola
 * Date: June 2021
 * Developed using ROOT 6.24/00
 * c++ -Wall -o xmltoTTreeRDF xmltoTTreeRDF.cpp `root-config --cflags --glibs`
 */
 
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <chrono>

#include <ROOT/RDataFrame.hxx>

#include <TFile.h>

using namespace std;

//#define KEEP_TXT //uncomment this line if you want to keep the txt files at the end of the execution

//function to get infos from auxiliary txt file
void getInfos (std::string infoFile, unsigned int &Nevts, unsigned int &SamplesInOneEvent, unsigned int &Nchannels)
{
	ifstream inFile;
	inFile.open(infoFile);
	if( !inFile.is_open() )
		throw std::runtime_error( "info file not found" );
	inFile>>Nevts>>SamplesInOneEvent>>Nchannels;
	inFile.close();
}

//function to execute python script that will generate auxiliary txt file with the xml data
void execute_python_script(string script_dir, std::vector<string> python_inputs ){
	auto t1 = std::chrono::high_resolution_clock::now();
	string command = "python3 " + script_dir + " ";
	for (unsigned int i = 0; i < python_inputs.size(); i++)
		command = command + python_inputs[i] + " ";
	std::cout<<endl<<"executing python script:"<<endl<<command<<endl;
	system(command.c_str());
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>( t2 - t1 ).count();
	std::cout<<"execution python script terminated in "<<duration<<" minutes"<<endl<<endl;
}

int main(int argc, char** argv)
{
	if ( argc != 3 ){
		std::cout<<"\nINPUT ERROR"<<endl
		<<"needs:\noutputDir \nrootInputFile\n"<<endl
		<<"as an example \n./executable path/to/output/folder/ path/to/file/fileName.xml"<<endl;
		return 1;
	}
	
	//setting I/O file names:
	string 	outDir      = argv[1],
	       	xmlIn       = argv[2],
	       	txtOut      = xmlIn  ,
		rootOut     = xmlIn  ,
		txtOutName  = txtOut.substr(txtOut.find_last_of("/")+1),
		rootOutName = rootOut.substr(rootOut.find_last_of("/")+1);
	if( outDir.begin() + outDir.find_last_of("/") != outDir.end() )
		outDir += "/";	
	txtOut      = outDir + txtOutName ;
	txtOut.replace(txtOut.end()-4, txtOut.end(), ".txt"); 
	rootOut     = outDir + rootOutName;
	rootOut.replace(rootOut.end()-4, rootOut.end(), "_amulet.root");
	//just an auxiliary info file produced by the python script from which to take important data to make the RDataFrame
	string txtFileInfo = outDir + txtOutName.replace(txtOutName.end()-4, txtOutName.end(), "_infos.txt"), 
	txtIn = txtOut;
	
	//print some infos:
	std::cout<<endl<<"\nReading RAW data (waveforms) from: "<<xmlIn
		 <<endl<<"Saving data in TTree in root file: "<<rootOut<<endl;
	
	//run python code to write data in external auziliary txt file:
	string pythonScriptDir = "../digitizer_scripts/python_modules/xmltotxt.py";
	execute_python_script(pythonScriptDir, {xmlIn, txtOut});

	unsigned int nEvents = 0, sampInOneEvent = 0, nExistingChs = 0;
	getInfos( txtFileInfo, nEvents, sampInOneEvent, nExistingChs );

	//open input file
	ifstream inFile(txtIn);
	if( !inFile.is_open() )
		throw std::runtime_error( "impossible to open file" );
	
	//prepare variables to store data
	int voltCh0 = 0, voltCh1 = 0, voltCh2 = 0, voltCh3 = 0, aux=0;
	ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> WVFofThisEventChannels;
	ROOT::VecOps::RVec<int> WVFofThisEventch0, WVFofThisEventch1, WVFofThisEventch2, WVFofThisEventch3;
	
	//prepare a cpp lambda function to get the actual WVF
	//returns indented list containing each WVF in each channel
	auto getWVFs = [&](){
		WVFofThisEventch0.clear();
		WVFofThisEventch1.clear();
		WVFofThisEventch2.clear();
		WVFofThisEventch3.clear();
		for(unsigned int i = 0; i<sampInOneEvent; i++ )
		{
			if( !inFile.is_open() )
				throw std::runtime_error( "file is not open" );
			if (nExistingChs == 1){
				inFile>>voltCh0;
				WVFofThisEventch0.push_back(voltCh0);
			}else if (nExistingChs == 2){
				inFile>>voltCh0>>voltCh1;
				WVFofThisEventch0.push_back(voltCh0);
				WVFofThisEventch1.push_back(voltCh1);
			}else if (nExistingChs == 3){
				inFile>>voltCh0>>voltCh1>>voltCh2;
				WVFofThisEventch0.push_back(voltCh0);
				WVFofThisEventch1.push_back(voltCh1);
				WVFofThisEventch2.push_back(voltCh2);
			}else if (nExistingChs == 4){
				inFile>>voltCh0>>voltCh1>>voltCh2>>voltCh3;
				WVFofThisEventch0.push_back(voltCh0);
				WVFofThisEventch1.push_back(voltCh1);
				WVFofThisEventch2.push_back(voltCh2);
				WVFofThisEventch3.push_back(voltCh3);
			}else
				throw runtime_error( "nExistingChs is not in the possible values" );
		}
		WVFofThisEventChannels = {WVFofThisEventch0, WVFofThisEventch1, WVFofThisEventch2, WVFofThisEventch3};
		if( inFile.eof() )
			throw std::runtime_error( "reached end of file while reading data" );
		return WVFofThisEventChannels;
	};

	//tell snapshot to update root file
	ROOT::RDF::RSnapshotOptions opt;
	opt.fMode="RECREATE";

	//here is where the computation takes place: RDataFrame will call the lambda function
    	auto t1 = std::chrono::high_resolution_clock::now();
	ROOT::RDataFrame df_wvf(nEvents); //create empty df
	std::cout<<"Creating and saving TTree (using RDataFrame)..."<<endl;
	//define columns (TTree branches)
	df_wvf	.Define("allChannels", getWVFs   ) //get list of lists containing ALL the four digitizer channels (even empty ones, then you decide wich to keep)
		.Define("ch0_wvf_amp",[=](ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> allChs){ return allChs[0]; }, {"allChannels"} ) //get channel 0
		.Define("ch1_wvf_amp",[=](ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> allChs){ return allChs[1]; }, {"allChannels"} ) //get channel 1
		.Snapshot("amulet"   , rootOut.c_str(), {"ch0_wvf_amp","ch1_wvf_amp"}, opt); //save TTree in external file
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>( t2 - t1 ).count();
	std::cout<<"Snapshot has been saved in "<<rootOut.c_str()<<endl
		 <<"in "<<duration<<" minutes"<<endl;

	//check that file has been completely read
	inFile>>aux;
	if( inFile.is_open() && !inFile.eof() )
		std::cout<<"\n!!!!!!!!!WARNING: you probably haven't read all the data!!!!!!!!!!!!\n"<<endl;
	inFile.close();
	
	//define KEEP_TXT to keep auxiliary txt filese
	#ifndef KEEP_TXT
	if (!(remove(txtOut.c_str()) == 0 && remove(txtFileInfo.c_str()) == 0))
		std::cout<<"problems in eliminating auxiliary txt files"<<std::endl;
	#endif
	
	return 0;
}
