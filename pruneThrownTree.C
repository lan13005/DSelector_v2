// Load the library at macro parsing time: we need this to use its content in the code
//R__LOAD_LIBRARY($ROOTSYS/test/libEvent.so)
 
////////////
// The goal of this program is to remove a bunch of useless branches. 
//  This is partly done to allow the use of uproot since it doesn't like 4-vectors here for w.e. reason
//  I think flat_to_amptools.C also doesnt like the full thrown trees


void pruneThrownTree()
{
   //vector<string> tags={"d2","d1","d0","p0","p1","s0"};
   //string baseFolder="/d/grid17/ln16/myDSelector/nonres_eff_test/attempt7_matching/";
   //vector<string> tags={"2018_8"};
   //string baseFolder="/d/grid17/ln16/myDSelector/nonres_eff_test/attempt5_matching/";
   //vector<string> tags={"2018_8"};
   //string baseFolder="/d/grid17/ln16/myDSelector/nonres_eff_test/attempt7_matching/";
   //vector<string> tags={"2018_1"};
   //string baseFolder="/d/grid17/ln16/myDSelector/ztSlopeMatching/";
   vector<string> tags={""};
   string baseFolder="/d/grid17/ln16/myDSelector/";
   bool includeP4=false;
   for (auto tag : tags){
        //string oldFileName="selected_degALL_nonres_eff_test_zlm_"+tag+"_matchingFlat2018_8_gen_trees_DSelector.root";
        //string oldFileName="selected_degALL_flat_gen_"+tag+"_tslope1_trees_DSelector.root";
        //string oldFileName="selected_deg000_flat_gen_2018_8_trees_DSelector.root";
        //string oldFileName="selected_degALL_flat_gen_2018_1_trees_DSelector.root";
        string oldFileName="degALL_b1vps_as_4g_gen_trees_DSelector.root";
        cout << "running over: " << baseFolder+oldFileName << endl;
        string filename=baseFolder+oldFileName;
        TFile oldfile(filename.c_str());
        TTree *oldtree;
        oldfile.GetObject("Thrown_Tree", oldtree);
        
        TLorentzVector* g1_p4_kin=0;
        TLorentzVector* g2_p4_kin=0;
        TLorentzVector* g3_p4_kin=0;
        TLorentzVector* g4_p4_kin=0;
        TLorentzVector etaP4;
        TLorentzVector pi0P4;
        double Mpi0eta;
        double Mpi0p;
        double Metap;
        double mandelstam_t;
        double mandelstam_tp;
        double Ebeam;
        double cosTheta_eta_gj;
        double cosTheta_eta_hel;
        double phi_eta_gj;
        double phi_eta_hel;
        cout << "-----\nPRE\n-----" << endl;
        if (includeP4){
            oldtree->SetBranchAddress("g1_p4_kin",&g1_p4_kin);
            oldtree->SetBranchAddress("g2_p4_kin",&g2_p4_kin);
            oldtree->SetBranchAddress("g3_p4_kin",&g3_p4_kin);
            oldtree->SetBranchAddress("g4_p4_kin",&g4_p4_kin);
            for (int i=0; i<10; i++){
                oldtree->GetEntry(i);
                etaP4=*g3_p4_kin+*g4_p4_kin;
                pi0P4=*g1_p4_kin+*g2_p4_kin;
                cout << etaP4.M() << ", " << pi0P4.M() << endl;
            }
        }
        oldtree->SetBranchAddress("Mpi0eta",&Mpi0eta);
        oldtree->SetBranchAddress("Mpi0p",&Mpi0p);
        oldtree->SetBranchAddress("Metap",&Metap);
        oldtree->SetBranchAddress("mandelstam_tp",&mandelstam_tp);
        oldtree->SetBranchAddress("mandelstam_t",&mandelstam_t);
        oldtree->SetBranchAddress("Ebeam",&Ebeam);
        oldtree->SetBranchAddress("cosTheta_eta_gj",&cosTheta_eta_gj);
        oldtree->SetBranchAddress("cosTheta_eta_hel",&cosTheta_eta_hel);
        oldtree->SetBranchAddress("phi_eta_gj",&phi_eta_gj);
        oldtree->SetBranchAddress("phi_eta_hel",&phi_eta_hel);

        // Deactivate all branches
        oldtree->SetBranchStatus("*", 0);
        
        // Activate some branches
        //vector<string> activeBranches = {"beam_p4_kin", "p_p4_kin", "g1_p4_kin", "g2_p4_kin", "g3_p4_kin", "g4_p4_kin", 
        //                                    "Mpi0eta", "mandelstam_t", "mandelstam_tp", "Ebeam"};
        vector<string> activeBranches = {"Mpi0eta", "Mpi0p", "Metap", "mandelstam_t", "mandelstam_tp", "Ebeam",
            "cosTheta_eta_gj", "cosTheta_eta_hel", "phi_eta_gj", "phi_eta_hel"
        };
        for (auto activeBranchName : activeBranches)
           oldtree->SetBranchStatus(activeBranchName.c_str(), 1);
        
        // Create a new file + a clone of old tree in new file
        TFile newfile((baseFolder+"pruned_"+oldFileName).c_str(), "recreate");
        auto newtree = oldtree->CloneTree();

        cout << "-----\nPOST\n-----" << endl;
        if (includeP4){
            newtree->SetBranchAddress("g1_p4_kin",&g1_p4_kin);
            newtree->SetBranchAddress("g2_p4_kin",&g2_p4_kin);
            newtree->SetBranchAddress("g3_p4_kin",&g3_p4_kin);
            newtree->SetBranchAddress("g4_p4_kin",&g4_p4_kin);
            for (int i=0; i<10; i++){
                newtree->GetEntry(i);
                etaP4=*g3_p4_kin+*g4_p4_kin;
                pi0P4=*g1_p4_kin+*g2_p4_kin;
                cout << etaP4.M() << ", " << pi0P4.M() << endl;
            }
        }
        //newtree->Print();
        newfile.Write();
   }
}
