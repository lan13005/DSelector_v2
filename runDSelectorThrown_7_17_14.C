// macro to process analysis TTree with DSelector
// We cannot just run this macro, the library doesnt load properly. We can run the following two lines of code
//.x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C
//.x runDSelector.C

#include <iostream> 
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"

R__LOAD_LIBRARY(libDSelector.so) 
   
void runDSelectorThrown_7_17_14(bool proof = 1, string path = "") 
{
	// Load DSelector library
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	int proof_Nthreads = 36;
	//int proof_Nthreads = 50;

	// open ROOT files and TTree
	//TChain *chain = new TChain("ksks_Tree");
	//Check the name of the tree in the root files
	TString nameOfTree = "Thrown_Tree"; // pi0eta__B4_Tree is the old one
	TChain *chain = new TChain(nameOfTree);

        // Dont delete this space since background sampels need to be in a specific place
        // Dont delete this space since background sampels need to be in a specific place
        // Dont delete this space since background sampels need to be in a specific place
        // Dont delete this space since background sampels need to be in a specific place
        // Dont delete this space since background sampels need to be in a specific place
        // Dont delete this space since background sampels need to be in a specific place
        // Dont delete this space since background sampels need to be in a specific place
        // Dont delete this space since background sampels need to be in a specific place

	// **********************************************************************************	
	// ************************** ------ BACKGROUND SAMPLES ---------**************************	
	// **********************************************************************************	
//        chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/a2pi/tree_thrown*");
//        chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/b1/tree_thrown*");
//        chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/etap_to_etapipi/tree_thrown*");
//        chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/eta_to_3pi/tree_thrown*");
//        chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/f1_to_etapipi/tree_thrown*");
//        chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/omega_pi0g_2018_8/thrown/tree_thrown*.root");
//        chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/pi0pi0/tree_thrown*");

	//a0
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0_091519/tree_thrown.root");

	// a0a2
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2_noPlugin_Geant4_30730/tree_thrown.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a2_10M_071620/tree_thrown.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0_10M_071620/tree_thrown.root");

	// Flat
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_noPlugin_Geant4_30730_8to9GeV/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_8GeVPlus_lustre_upTo3GeVResMass/tree_thrown.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_2.1t/tree_thrown.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2017/tree_thrown.root");
        // newer flat narrow Mpi0eta range
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2017_Mpi0eta16to3/tree_thrown.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2018_8_Mpi0eta16to3/tree_thrown.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2018_1_Mpi0eta16to3/tree_thrown.root");
        // newer flat full Mpi0eta range
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_8to9GeV/thrown/tree_thrown_gen_amp_*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_8to9GeV_2018_8/thrown/tree_thrown_gen_amp_*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_tslope_matching/tslope5/thrown/tree_thrown_gen_amp_*");

        // Much larger flat sample
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/120921/2018_1_400M/merged/tree_thrown*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_130M/merged/tree_thrown*");
        
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2_posRefD0D2_negRefD1_posM/thrown/tree*");
        
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2a2pi1_polarized_largerPi1/thrown/tree_thrown_gen_amp*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/pi1_polarized/thrown/tree_thrown_gen_amp_*");

	// BA studies. Hddm filtered to make phi have some cos dependence
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_21t_hddmFiltered_8288_1628_uniquePolar/tree_thrown.root");

        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/nonres_eff_tests/etapi0_zlm_d2/thrown/tree_*"); 
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/nonres_eff_tests/matchingFlat2018_8/etapi0_zlm_d2/thrown/*"); 

	// test
	//chain->Add("/d/grid13/ln16/MC/pi0eta_flat_2.3t/hddm/tree_thrown.root");

        // Justin VPS
        chain->Add("/d/grid17/ln16/b1_vps/degALL_thrown_as_etapi_trees_DSelector.root");
        
	string degAngle="degALL";
        //string tag="_bkgndSample_gen";
        string tag="_b1vps_as_4g_gen";
        //string tag="_flat_gen_2018_1";
        //string tag="_flat_gen_2018_8";
        //string tag="_flat_gen_2018_8_tslope1";
        //string tag="_a0a2_posRefD0D2_negRefD1_posM_gen";
        //string tag="_nonres_eff_test_zlm_d2_matchingFlat2018_8_gen";

	// The following two files work separately but not together...

	//string options = sample;
	string options = "";
	if(proof) { // add TTree to chain and use PROOFLiteManager
		//string outputHistFileName = Form("flatUpTo3GeVResMass_2_gen_hists_DSelector_pi0eta.root");//_GEANT4.root");
                //string outputTreeFileName = Form("flatUpTo3GeVResMass_2_gen_trees_DSelector_pi0eta.root");//_GEANT4.root");
		string outputHistFileName = degAngle+tag+"_hists_DSelector.root";//_GEANT4.root");
                string outputTreeFileName = degAngle+tag+"_trees_DSelector.root";//_GEANT4.root");
		DPROOFLiteManager::Process_Chain(chain, "DSelector_thrown.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
		//DPROOFLiteManager::Process_Chain(chain, "DSelector_newTest.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		//TFile *f = TFile::Open(fileName);
		//TTree *tree = (TTree*)f->Get("omega_skim_Tree");
		chain->Process("DSelector_thrown.C+", options.data());
		//chain->Process("DSelector_newTest.C+", options.data());
	}
	
	return;
}
