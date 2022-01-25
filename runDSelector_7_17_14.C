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


//R__LOAD_LIBRARY(/group/halld/Software/builds/Linux_CentOS7-x86_64-gcc4.8.5/gluex_root_analysis/gluex_root_analysis-0.5/Linux_CentOS7-x86_64-gcc4.8.5/lib/libDSelector.so)
//R__LOAD_LIBRARY(/d/home/ln16/gluex_top/gluex_root_analysis/gluex_root_analysis_1.7.0/Linux_CentOS7-x86_64-gcc4.8.5/lib/libDSelector.so) 
R__LOAD_LIBRARY(libDSelector.so)
   
void runDSelector_7_17_14(bool useproof = 1, string path = "") 
{
	cout << "Loaded using R__LOAD_LIBRARY" << endl;
	// Load DSelector library
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	// change the directory that proof saves the data to
	//gEnv->SetValue("ProofLite.Sandbox", "/d/grid15/ln16/.proof");
	int proof_Nthreads = 48;
	//int proof_Nthreads = 50;

	// open ROOT files and TTree
	TString nameOfTree = "pi0eta__B4_M17_M7_Tree"; 
        //TString nameOfTree = "pi0eta__B4_M7_M17_Tree";
        //TString nameOfTree = "pi0eta__B4_Tree";
        //TString nameOfTree = "pi0eta__B4_M7_Tree";
        //TString nameOfTree = "pi0eta__B4_M7_Tree";
	TChain *chain = new TChain(nameOfTree);
	
	// **********************************************************************************	
	// ************************** ------ PI0ETA BELOW ---------**************************	
	// **********************************************************************************	
	// omega -> gamma pi0 as pi0eta
	//chain->Add("/d/grid15/ln16/rootFiles/omega_gammaPi0/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/pi0eta/testing/omegaPi0/rootFiles/fullRecon_20M/0303_0304/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/pi0eta/testing/omegaPi0/rootFiles/fullRecon_20M/0305_0306/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/pi0eta/testing/omegaPi0/rootFiles/fullRecon_20M/0307_0308/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/pi0eta/testing/omegaPi0/rootFiles/fullRecon_20M/0309_0310/tree_pi0eta__B4_M17_M7.root");
	//
	// pi0pi0 as pi0eta
	//chain->Add("/d/grid15/ln16/rootFiles/pi0pi0/noResonances/tree_pi0eta__B4_M17_M7.root");
	
        // MC flat pi0pi0 reconstructed as pi0eta
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_20M_2017_8GeVp/tree_pi0eta__B4_M17_M7.root");
        
        // MC flat pi0eta
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2017/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2017_Mpi0eta16to3/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2018_8_Mpi0eta16to3/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/etapi0_flat_8GeVp_2018_1_Mpi0eta16to3/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat/trees/tree_pi0eta__B4_M17_M7_gen_amp*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_8to9GeV/trees/tree_pi0eta__B4_M17_M7_gen_amp*");
        //
        //Much larger flat MC
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_260M/merged/tree_pi0eta__B4_M17_M7_merged*");
        chain->Add("/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_130M/merged/tree_pi0eta__B4_M17_M7_merged10.root");
        chain->Add("/d/grid17/ln16/rootFiles/pi0eta/120921/2018_8_130M/merged/tree_pi0eta__B4_M17_M7_merged11.root");

        // Testing tslope matching to data
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_tslope_matching/tslope1/trees/tree_pi0eta__B4_M17_M7_gen_amp*");

        //
        // Mass (un)constrained flat etapi0
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_massConstrained/trees/tree_*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_constrainedVsNonConstrained/unconstrained/trees/tree*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_constrainedVsNonConstrained/constrained/trees/tree*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_2018_8_kinFitStudies/etapi_flat_pi0MassConstrainedOnly/trees/tree*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/etapi_flat_2018_8_kinFitStudies/etapi_flat_etaMassConstrainedOnly/trees/tree*");
        
	// BA studies. Hddm filtered to make phi have some cos dependence
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/flat_21t_hddmFiltered_8288_1628_uniquePolar/tree_pi0eta__B4_M17_M7.root");
        
        // a0a2a2pi1 polarized amplitudes
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2a2pi1_polarized_largerPi1/trees/tree_pi0eta__B4_M17_M7_gen_amp_*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/malte_MC_pol035/trees/tree_pi0eta__B4_M17_M7_gen_amp*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2a2pi1_polarized_D2neg/trees/tree_pi0eta__B4_M17_M7_gen_amp_03*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2a2pi1_polarized_moreD2p_lessD10p_moreS/trees/tree_pi0eta__B4_M7_M17_gen_amp_03*");
        // a0a2 polarized amplitudes
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2_posRefD0D2_negRefD1_posM/trees/tree*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a2nonres_2017/trees/tree_*");
        
        // seans bkg MC
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/rootTrees/30M/b1/tree_pi0eta__B4_M17_M7_gen_amp*");
	
	// a0a2 recon_2017
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0a2_a2pi1/a0a2a2pi1_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0_091519/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a2_10M/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a2_10M_071620/tree_pi0eta__B4_M17_M7.root"); // newer a2
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/a0_10M_072120/tree_pi0eta__B4_M17_M7.root"); // newer a0
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2_D0p_moreD2p_D1m/trees/tree_pi0eta__B4_M17_M7_gen_amp_03*");
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/a0a2_D0p_moreD2p_D1p/trees/tree_pi0eta__B4_M17_M7_gen_amp_03*");
	
	// vincent
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30461/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30730/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30461_withP/tree_pi0eta__B4_M17_M7.root");
	//chain->Add("/d/grid15/ln16/rootFiles/pi0eta/vincent_noPlugin_Geant4_30730_withP/tree_pi0eta__B4_M17_M7.root");
	//
	// sean resolution test
	//chain->Add("/d/home/sdobbs/grid13/MCtest/pi0eta_a2/hddm/tree_pi0eta__B4_F1_M7_M17_00*");
	// sean new resolution test
	//chain->Add("/d/home/sdobbs/grid13/MCtest/71001_0/tree_pi0eta__B4_F1_M7_M17.root");
	//chain->Add("/d/home/sdobbs/grid13/MCtest/71001_3.save/tree_pi0eta__B4_F1_M7_M17.root");

	// 2017 DATA
	// allData
	// This one contains the showerQuality variables
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/analysis-ver27/tree_pi0eta__B4_M17_M7/merged/*");
        // ver52
	//chain->Add("/d/grid17/ln16/rootFiles/pi0eta/RunPeriod-2017-01/analysis-ver52/tree_pi0eta__B4_M17_M7/merged/*");
	// 2018 DATA
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2018-01/analysis-ver02/tree_pi0eta__B4_M17_M7/merged/*");
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2018-08/analysis-ver02/tree_pi0eta__B4_M17_M7/merged/tree_pi0eta__B4_M17_M7_*");

        // ********************* USING THIS FOR Q-values so far *********		
	//chain->Add("/d/grid15/ln16/pi0eta/092419/boolSpectPi0EtaSelectedFiles/pi0eta_data_tree_DSelector.root");


        // ************* MALTE KMATRIX ********
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/malte_kmatrix/tree_pi0eta__B4_M7_M17*");

        // *********** From mergeEvents.C, merging etapi and b1 events that prefer the 4g hypothesis
        //chain->Add("/d/grid17/ln16/mergingTrees/etapi_and_b1_as_4g.root");
        //chain->Add("/d/grid17/ln16/mergingTrees/reject5gHypoStudy/etapi_primary_etapi_as_b1_secondary.root");
        //chain->Add("/d/grid17/ln16/mergingTrees/etapi/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid17/ln16/mergingTrees/reject5gHypoStudy/b1_as_etapi_primary_b1_secondary.root");
        //chain->Add("/d/grid17/ln16/mergingTrees/b1_as_etapi/tree_pi0eta__B4_M17_M7.root");


        // ********************** v2 base cuts - using Justin's vertex selections - looser photon theta - ******************** 
        //chain->Add("./zBaseCuts/baseCuts_with_looseChiUE_v2/degALL_data_2018_1_baseCuts_looseChiUE_tree_DSelector.root");
        // ********************** STD PROTON/PHOTON/EXCLUSIVITY LOOSE CHI CUTS ********************
        //chain->Add("./zThesisCuts/degALL_data_2017_stdProtonPhotonExclusivity_looseChi_for_thesis_tree_DSelector.root");
        // ********************** NEW BASE CUTS APPLIED ONLY LOOSE CHI UE - Think only difference is with UE cut compared with above******************** 
        //chain->Add("./zBaseCuts/baseCuts_with_looseChiUE/degALL_data_2017_baseCuts_looseChiUE_tree_DSelector.root");
        // ********************** NEW BASE CUTS APPLIED ONLY TIGHT CHI UE ********************
        //chain->Add("./zBaseCuts/baseCuts_with_tightChiUE/degALL_data_2017_baseCuts_tree_DSelector.root");
        //
        // ********************** OLD  BASE CUTS APPLIED ONLY ********************
	//// 2017
        //chain->Add("/d/grid15/ln16/pi0eta/092419/zSelectedBaseCuts/pi0eta_data_tree_DSelector.root");
	//// 2017 - loose ChiSq and UE cut
	//chain->Add("/d/grid15/ln16/pi0eta/092419/zSelectedLooseChiSqUE/pi0eta_looseCuts_tree_DSelector.root");
	// 2018_8
	//chain->Add("/d/grid15/ln16/pi0eta/092419/zSelectedBaseCuts_2018_8/pi0eta_baseCuts_2018_tree_DSelector.root");
	// 2018_1
	//chain->Add("/d/grid15/ln16/pi0eta/092419/zSelectedBaseCuts_2018_1/pi0eta_baseCuts_2018_tree_DSelector.root");
	//
	// BA for Double regge MC studies
	// Rejection sampled the uniform phi distribution. 
	//chain->Add("/d/grid15/ln16/hddm_filtered_root/zFiltered/fullRecon/tree_pi0eta__B4_M17_M7.root");

        // quick run over b1 and flat etapi trees 
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/runAll/b1to5g/degALL_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_DSelector.root");
        //chain->Add("/d/grid15/ln16/rootFiles/pi0eta/seansBkgMC/runAll/flat_etapi/degALL_mEllipse_8288_chi13_tpLT05_pipicut_omegacut_tree_DSelector.root");

        // vector pseudoscalar b1 simulations from justin
        //chain->Add("/d/grid17/ln16/rootFiles/vectorps_b1_justin/b1_vectorps_justin_as_etapi/tree_pi0eta__B4_M17_M7.root");
        //chain->Add("/d/grid17/ln16/b1_vps/degALL_recon_as_etapi_tree_DSelector.root");
        
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/nonres_eff_tests/etapi0_zlm_/trees/*"); 
        //chain->Add("/d/grid17/ln16/rootFiles/pi0eta/010820/nonres_eff_tests/matchingFlat2018_8/etapi0_zlm_d2/trees/*");
        
	// **********************************************************************************	
	// ************************** ------ PI0PI0 BELOW ---------**************************	
	// **********************************************************************************	
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/analysis-ver20/tree_pi0pi0__B3_F1_M7/merged/tree_pi0pi0__B3_F1_M7_*");
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/tree_pi0pi0__B3_F1_U1_M7/merged/tree_pi0pi0__B3_F1_U1_M7_03*");
	//chain->Add("/d/home/sdobbs/GlueX/gluex_data/RunPeriod-2017-01/tree_pi0pi0__B3_F1_U1_M7/merged/tree_pi0pi0__B3_F1_U1_M7_03033*");
	//
	// MC with f2 resonance
	//chain->Add("/d/grid15/ln16/rootFiles/pi0pi0/Aug27Sep09/tree_pi0pi0__B3_F1_U1_M7.root");
	//chain->Add("/d/grid13/ln16/MC/pi0pi0_f2_interactive/hddm/combined/tree_pi0pi0__B3_F1_U1_M7.root");
        //
        
        //string tag="_data_2017_BA";
        //string tag="_flat_2018_1_mEllipse_mt_mDelta_mOmega_mEbeam";
        //string tag="_b1_vectorps_as_etapi_BA";
        //string tag="_b1_vps_as_etapi_selectOmega_mEllipse";
        //string tag="_b1_vps_as_etapi_stdProtonPhotonExclusivity_chi13_for_thesis";
        //string tag="_data_2018_1_mEllipse_mMandelstam_t_mEbeam_mDelta_chi13";
        //string tag="_flat_2018_8_mEllipse_tslope1";
        //string tag="_flat_2018_8_mEllipse";
        //string tag="_data_2017_BA_loosePi0EtaCut";
        //string tag="_flat_2018_8_mEllipse_mMandelstam_t_8288_chi13_omegacut";
        //string tag="_flat_2018_1_BA";
        //string tag="_nonres_eff_test_d2_mEllipse";
        //string tag="_nonres_eff_test_zlm_d2_matchingFlat2018_8_mEllipse";
        
        //string tag="_flat_2017_noSelections_for_thesis";
        //string tag="_flat_2017_stdProtonPhoton_looseChi_for_thesis";
        //string tag="_flat_2017_stdProtonPhotonExclusivity_chi13_for_thesis";

	// should change the name below from data to reco when running over MC
	string degAngle="degALL";
        //string tag="_b1vps_as_4g_mEllipse_8288_tLT1_chi13";
        //string tag="_kmatrix_mEllipse_8288_tLT1_chi13";

        //string tag="_flat_2017_mEllipse_mt_mDelta_mOmega_mEbeam";
        //string tag="_flat_2017_mEllipse_mEbeam_mt";
        //string tag="_data_2018_1_mEllipse_mEbeam_mt_mDelta_mOmega";
        string tag="_flat_2018_8_mEllipse_8288_tLT1";
        //string tag="_data_2017_BA";
        //string tag="_data_2017_baseCuts_looseChiUE";
        //string tag="_data_2017_mEllipse_8288_chi13_tLT1_pipicut_omegacut";
        //string tag="_malte_kmatrix_2018_8_mEllipse_8288_tLT1_chi13_omegacut";//_mPhotonETheta_mDelta";
        //string tag="_flat_2017_mEllipse_8288_tpLT1_chi13_omegacut_rfpeak";
        //
        //string tag="_b1_as_4g_all_mEllipse_8288_tLT1_chi13_omegacut";
        //string tag="_etapi_as_4g_mEllipse_8288_tLT1_chi13_omegacut";
        //string tag="_b1_as_4g_chiSqUEOnly";
        //
        //string tag="_flat_M7_2018_8_mEllipse_8288_tpLT1_chi13_omegacut_rfpeak";
        //string tag="_malte_kmatrix_2018_8_mEllipse_8288_tLT1_chi13_omegacut";//_mPhotonETheta_mDelta";

	//  ===== Section is for pulling in data by polarization for asymmetry ===== /////
	// dividing allData into gorups by polarization. Actually the files we import have
	// the root files separated by polarizations. 
	//ifstream myReadFile;
	//TString degFile = "separateRuns/"+degAngle+".txt";
	//myReadFile.open(degFile);
	//char output[120];
	//if (myReadFile.is_open()) {
	//	// read the degAngle files until the end of it... and Add the file to the chain
	//	 while (!myReadFile.eof()) {
	//		myReadFile >> output;
	//		chain->Add(output);
	//		cout << output << endl;
	//	}
	//}
	/// ========================================================================///////


	//char oldHome[80]="HOME=/d/home/ln16";
	//char newHome[80]="HOME=/d/grid15/ln16";
	//putenv(newHome);
	//system("echo New Home At:");
	//system("echo $HOME");

	//string options = sample;
	string options = "";
	if(useproof) { // add TTree to chain and use PROOFLiteManager
		string outputHistFileName;
                string outputTreeFileName;
		outputHistFileName = degAngle+tag+"_hists_DSelector.root";
		outputTreeFileName = degAngle+tag+"_tree_DSelector.root"; 
	
		DPROOFLiteManager::Process_Chain(chain, "DSelector_ver20.C+",  proof_Nthreads, outputHistFileName, outputTreeFileName, options);
	}
	else { // get TTree and use standard TTree::Process
		//TFile *f = TFile::Open(fileName);
		//TTree *tree = (TTree*)f->Get("omega_skim_Tree");
		chain->Process("DSelector_ver20.C+", options.data());
		
	}
	
	//putenv(oldHome);
	//system("echo Returning to  Home At:");
	//system("echo $HOME");
	cout << "\033[1;31mMAKE SURE IF YOU ARE RUNNING ON PI0PI0 DATA WE USE mergePi0.root WITH THE CORRECT FILE!\n1)Change file name in mergePi0.root\n2)root -l -b -q mergePi0.root \033[0m\n";
	cout << "\033[1;31mRemember that Mpi0, Meta, Mpi0eta fundamental branches have a cut applied on them! This is for SPlotting to define a better region \033[0m\n";

	return;
}
