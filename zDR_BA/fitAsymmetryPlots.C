// This one is for Asymmetry vs t_recoil which is equal to u3 in vincent/colin language.
//#include "/d/grid13/gluex/gluex_top/gluex_style.C"

double degToRad=TMath::Pi()/180;
// par[3] is used to shift phase by the para or perp orientation, either 0 for para or 90 for perp. 0/-45 is para and 45/90 is perp. 
int numDOFsig_sc = 3;
Double_t shiftedCos(Double_t *x, Double_t *par){
	return par[0]*(1.0 - par[1]*TMath::Cos(2*degToRad*(x[0]-par[2])));
}
int numDOFsig_flat = 1;
Double_t flat(Double_t *x, Double_t *par){
	return par[0];
}

// par0 and par1 = polarization fraction of perp and para components respectively 
int numDOFsig_asym=4;
Double_t asymmetry(Double_t *x, Double_t *par){
	return ((par[0]+par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3]))/(2+(par[0]-par[1])*par[2]*TMath::Cos(2*degToRad*(x[0]-par[3]))));
}

// Number of t1 bins to consider, have to manually enter the names
int nt1bins=5; 
string mint1s[5] = {"01", "03", "05", "07", "09"};
//int nt1bins=1;
//string mint1s[1] = {"05"}; 

// This tag will be appended to the resulting csv file and all folders, can use to denote specific studies like differences in weighting
//string tag="_sig_AS"; // do accidentals only so we need to select on signal region also. Either through the rectangular selection of insideEllipse bool
//string weightVar="AccWeight";
//string tag="_ASBS"; // do accidentals and sideband subtraction 
//string weightVar="weightASBS";
string tag="_etaRight_ASBS"; // select left of eta peak and sideband subtract that side also
string weightVar="weightASBS";


string fileType="png";
string folder="/d/grid17/ln16/myDSelector/zDR_BA/fitAsymmetry_trees_loosePi0EtaCut_mEllipse_v2/";
static const int nDataSets = 3;
string dataSetTag[nDataSets] = { "2017" , "2018_1", "2018_8" };
string dataSetTag2[nDataSets+1] = { "2017" , "2018_1", "2018_8", "total" };
double fluxRatios_90_0[nDataSets] = {  4.346818e+12/4.188001e+12, 0.965429, 0.918503 };
double fluxRatios_45_135[nDataSets] = {  4.076065e+12/4.095013e+12 , 1.02261, 1.03254 };
string polnames[7] = {"phi000","phi045","phi090","phi135","phi000_090","phi045_135","phiAMO"};
// We dont actually have a pol offset for the total dataset so we will just choose the largest dataset, 2018_1 
//    These numbers are from Alex's rho SDME analysis
float Phi0_0_90[nDataSets+1] = {3.1,4.2,3.1,4.2}; 
float Phi0_45_135[nDataSets+1] = {-45+3.2,-45+2.8,-45+3.1,-45+2.8};
// Last 2 elements is or {0,90} and {45/-45} summed which should be ~flat. Tegan/Will uses the asymmetry offsets. These -999 will be updated in the dataset loop
float poloffset[6] = {0, 45, 90, -45, -999, -999}; 
int minEntriesInPhiHists=30;

// WL for weighted log-likelihood but probably doenst work for asymmetry histogram since the entries are not counts
// M for improved MINUIT results
// E for error estimation - Minos
string fitOption = "E S Q"; 

// From mike dugger's script
vector<vector<float>>  polfractions = { //2017, 2018_1, 2018_8, total
    {0.3537,0.3484,0.3472,0.3512}, // 0, 45, 90, -45 
    {0.3458,0.3416,0.3448,0.3597},
    {0.3563,0.3403,0.3430,0.3523},
    {0.3458,0.3416,0.3448,0.3597} // copy the 2018_1 fractions for the totals
};
vector<vector<float>>  polfractions_err = {
    {0.0100,0.0102,0.0098,0.0102},  
    {0.0055,0.0056,0.0055,0.0057},
    {0.0070,0.0075,0.0074,0.0077},
    {0.0055,0.0056,0.0055,0.0057} // copy the 2018_1 fractions for the totals
};

struct criteria{
    string variable;
    float minval;
    float maxval;
};

void constructAndFit(
    map<string, float*> &variable_map, 
    map<string, vector<float>> &array_variable_map,  
    map<string, vector<int>> &array_BeamAngles_map, 
    vector<criteria> selections,
    string fitSaveLoc,
    ofstream* saveCsv
    ){
    string phis_orientation[5] =  {"phi000", "phi045", "phi090", "phi135", "phiAMO"};
    float t1HalfWidth=1.0/nt1bins/2;

    gStyle->SetOptFit(111);
    gStyle->SetStatY(1);
    gStyle->SetStatX(1);
    gStyle->SetStatW(0.18);
    gStyle->SetStatH(0.18);
    gStyle->SetFitFormat("3.2g");

    // *****************************
    // INITALIZE THE HISTOGRAMS
    // *****************************
    TH1F *phi000_eta[nDataSets+1][nt1bins]; // +1 is for the total phase 1 dataset
    TH1F *phi045_eta[nDataSets+1][nt1bins]; 
    TH1F *phi090_eta[nDataSets+1][nt1bins]; 
    TH1F *phi135_eta[nDataSets+1][nt1bins]; 
    TH1F *phiAMO_eta[nDataSets+1][nt1bins]; 
    TH1F *phi000_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi045_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi090_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi135_pi0[nDataSets+1][nt1bins]; 
    TH1F *phiAMO_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi000_090_eta[nDataSets+1][nt1bins]; 
    TH1F *phi045_135_eta[nDataSets+1][nt1bins]; 
    TH1F *phi000_090_pi0[nDataSets+1][nt1bins]; 
    TH1F *phi045_135_pi0[nDataSets+1][nt1bins]; 
    
    float minphi=-180;
    float maxphi=180;
    float nphi=30;
    int binphi=9;
    for (int iData=0; iData<nDataSets; ++iData){
        for (int it1bin=0; it1bin<nt1bins; ++it1bin){
            string yaxis="Entries / "+to_string(binphi)+" degrees";
            phi000_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi000_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi045_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi045_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi); 
            phi090_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi090_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi135_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi135_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phiAMO_eta[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phiAMO_eta"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi000_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi000_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi045_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi045_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi090_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi090_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phi135_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phi135_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
            phiAMO_pi0[iData][it1bin] = new TH1F((dataSetTag[iData]+"_phiAMO_pi0"+to_string(it1bin)).c_str(),(";#phi; "+yaxis).c_str(), nphi, minphi, maxphi);
        }
    }
    
    ////////////////////////////////////
    // FILL THE HISTOGRAMS
    ////////////////////////////////////
    float weight;
    float selectionVariable;
    float teta;
    float tpi0;
    bool selected;
    int BeamAngle;
    int total_selected=0;
    int iteta;
    int itpi0;
    float p0; // parameter for shiftedCosine 
    int counter;
    for (int iData=0; iData<nDataSets; ++iData){
        Long64_t nentries=(Long64_t)array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"].size();
        for (Long64_t ientry=0; ientry<nentries; ++ientry){
            selected=true;
            for (auto selection: selections){
                selectionVariable=array_variable_map[dataSetTag[iData]+"_"+selection.variable][ientry]; 
                if(selectionVariable<=0.0){
                    cout << "selectionVar negative!" << endl;
                    exit;
                }
                selected *= (selectionVariable>selection.minval)*(selectionVariable<selection.maxval);
            } 
            if(selected){
                weight=array_variable_map[dataSetTag[iData]+"_"+weightVar][ientry];
                BeamAngle=array_BeamAngles_map[dataSetTag[iData]][ientry];
                teta=array_variable_map[dataSetTag[iData]+"_mandelstam_teta"][ientry];
                tpi0=array_variable_map[dataSetTag[iData]+"_mandelstam_tpi0"][ientry];
                //cout << BeamAngle << " " << teta << " " << tpi0 << endl;
             
                iteta=int(ceil(teta/(1.0/nt1bins))-1);
                itpi0=int(ceil(tpi0/(1.0/nt1bins))-1);   

                if(BeamAngle==0){
                    if (iteta<nt1bins)
                        phi000_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if (itpi0<nt1bins)
                        phi000_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==45){
                    if (iteta<nt1bins)
                        phi045_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if (itpi0<nt1bins)
                        phi045_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==90){
                    if (iteta<nt1bins)
                        phi090_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if (itpi0<nt1bins)
                        phi090_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==135){
                    if (iteta<nt1bins)
                        phi135_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if (itpi0<nt1bins)
                        phi135_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                else if(BeamAngle==-1){
                    if (iteta<nt1bins)
                        phiAMO_eta[iData][iteta]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_eta_lab"][ientry],weight); 
                    if (itpi0<nt1bins)
                        phiAMO_pi0[iData][itpi0]->Fill(array_variable_map[dataSetTag[iData]+"_"+"phi_pi0_lab"][ientry],weight); 
                }
                ++total_selected;
            }
        }
    }
    cout << "Total passed selections: " << total_selected << endl;
    
    TF1 *fit_asym = new TF1("fit_asym",asymmetry,-180,180,numDOFsig_asym); 
    TF1 *fit_sc = new TF1("fit_flat",shiftedCos,-180,180,numDOFsig_sc); 
    fit_asym->SetLineColor(kRed);
    fit_sc->SetLineColor(kRed);
    //////////////////////////////////////
    // SCALE THE HISTOGRAMS BY FLUX RATIO
    //////////////////////////////////////
    for (int it1bin=0; it1bin<nt1bins; ++it1bin){
        phi000_eta[3][it1bin] = new TH1F(("total_phi000_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phi045_eta[3][it1bin] = new TH1F(("total_phi045_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180); 
        phi090_eta[3][it1bin] = new TH1F(("total_phi090_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phi135_eta[3][it1bin] = new TH1F(("total_phi135_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phiAMO_eta[3][it1bin] = new TH1F(("total_phiAMO_eta"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phi000_pi0[3][it1bin] = new TH1F(("total_phi000_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phi045_pi0[3][it1bin] = new TH1F(("total_phi045_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phi090_pi0[3][it1bin] = new TH1F(("total_phi090_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phi135_pi0[3][it1bin] = new TH1F(("total_phi135_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);
        phiAMO_pi0[3][it1bin] = new TH1F(("total_phiAMO_pi0"+to_string(it1bin)).c_str(),";#phi; Entries / 9 degrees", 30, -180, 180);

        for (int iData=0; iData<nDataSets; ++iData){
            // First we scale all the para orientations with the flux ratios
            cout << "(" << dataSetTag[iData] << " 000 eta) entries/maximum before scaling: " << phi000_eta[iData][it1bin]->GetEntries() << "/" 
                                            << phi000_eta[iData][it1bin]->GetMaximum() << endl;
            phi000_eta[iData][it1bin]->Scale(fluxRatios_90_0[iData]);
            phi000_pi0[iData][it1bin]->Scale(fluxRatios_90_0[iData]);
            phi135_eta[iData][it1bin]->Scale(fluxRatios_45_135[iData]);
            phi135_pi0[iData][it1bin]->Scale(fluxRatios_45_135[iData]);
            cout << "(" << dataSetTag[iData] << " 000 eta) entries/maximum after scaling: " << phi000_eta[iData][it1bin]->GetEntries() << "/" 
                                            << phi000_eta[iData][it1bin]->GetMaximum() << endl;
        
            // Add up all the datasets for the total
            phi000_eta[3][it1bin]->Add(phi000_eta[iData][it1bin]);
            phi045_eta[3][it1bin]->Add(phi045_eta[iData][it1bin]); 
            phi090_eta[3][it1bin]->Add(phi090_eta[iData][it1bin]);
            phi135_eta[3][it1bin]->Add(phi135_eta[iData][it1bin]);
            phiAMO_eta[3][it1bin]->Add(phiAMO_eta[iData][it1bin]);
            phi000_pi0[3][it1bin]->Add(phi000_pi0[iData][it1bin]);
            phi045_pi0[3][it1bin]->Add(phi045_pi0[iData][it1bin]);
            phi090_pi0[3][it1bin]->Add(phi090_pi0[iData][it1bin]);
            phi135_pi0[3][it1bin]->Add(phi135_pi0[iData][it1bin]);
            phiAMO_pi0[3][it1bin]->Add(phiAMO_pi0[iData][it1bin]);
        }
        cout << "Phase 1 000 eta entries " << phi000_eta[3][it1bin]->GetEntries() << endl;
        
        //////////////////////////////////
        // FIT PRELIMINARIES
        //////////////////////////////////
        TCanvas *allCanvases = new TCanvas("","",1440,900);
        Int_t fitStatus;
        
        TFitResultPtr fitPointer;
        float asymmetries_000_eta;
        float asymmetries_000_eta_err;
        float asymmetries_045_eta;
        float asymmetries_045_eta_err;
        float asymmetries_000_pi0;
        float asymmetries_000_pi0_err;
        float asymmetries_045_pi0;
        float asymmetries_045_pi0_err;
        float pol_frac_000; 
        float pol_frac_090; 
        float pol_frac_045; 
        float pol_frac_135; 
        float asymmetry000_090_eta_chiPerDOF;
        float asymmetry045_135_eta_chiPerDOF;
        float asymmetry000_090_pi0_chiPerDOF;
        float asymmetry045_135_pi0_chiPerDOF;
        for(int iData=0; iData<4; ++iData){
            pol_frac_000 = polfractions[iData][0];//+-polfractions_err[iData][0]; 
            pol_frac_045 = polfractions[iData][1];//+-polfractions_err[iData][1]; 
            pol_frac_090 = polfractions[iData][2];//+-polfractions_err[iData][2]; 
            pol_frac_135 = polfractions[iData][3];//+-polfractions_err[iData][3]; 
            poloffset[4]=Phi0_0_90[iData];
            poloffset[5]=Phi0_45_135[iData];
            // *****************************
            // GETTING ASYMMETRIES
            // *****************************
            TH1* asymmetry000_090_eta = phi090_eta[iData][it1bin]->GetAsymmetry(phi000_eta[iData][it1bin]);
            TH1* asymmetry045_135_eta = phi045_eta[iData][it1bin]->GetAsymmetry(phi135_eta[iData][it1bin]);
            asymmetry000_090_eta->SetTitle("0/90 Asymmetry");
            asymmetry045_135_eta->SetTitle("45/135 Asymmetry");
            asymmetry000_090_eta->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            asymmetry045_135_eta->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            TH1* asymmetry000_090_pi0 = phi090_pi0[iData][it1bin]->GetAsymmetry(phi000_pi0[iData][it1bin]);
            TH1* asymmetry045_135_pi0 = phi045_pi0[iData][it1bin]->GetAsymmetry(phi135_pi0[iData][it1bin]);
            asymmetry000_090_pi0->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            asymmetry045_135_pi0->GetYaxis()->SetTitle(("Yield Asymmetry / "+to_string(binphi)+" degrees").c_str());
            asymmetry000_090_pi0->SetTitle("0/90 Asymmetry");
            asymmetry045_135_pi0->SetTitle("45/135 Asymmetry");
            // Constructing sums of paired polarizations for systematics
            phi000_090_eta[iData][it1bin]=(TH1F*)phi000_eta[iData][it1bin]->Clone("eta_000_090");
            phi045_135_eta[iData][it1bin]=(TH1F*)phi045_eta[iData][it1bin]->Clone("eta_045_135");
            phi000_090_pi0[iData][it1bin]=(TH1F*)phi000_pi0[iData][it1bin]->Clone("pi0_000_090");
            phi045_135_pi0[iData][it1bin]=(TH1F*)phi045_pi0[iData][it1bin]->Clone("pi0_045_135");
            phi000_090_eta[iData][it1bin]->Add(phi090_eta[iData][it1bin]);
            phi045_135_eta[iData][it1bin]->Add(phi135_eta[iData][it1bin]);
            phi000_090_pi0[iData][it1bin]->Add(phi090_pi0[iData][it1bin]);
            phi045_135_pi0[iData][it1bin]->Add(phi135_pi0[iData][it1bin]);
            // Loading phi distributions into a vector
            TH1F *phis_eta[7] =  {phi000_eta[iData][it1bin], phi045_eta[iData][it1bin], phi090_eta[iData][it1bin], 
                                    phi135_eta[iData][it1bin], phi000_090_eta[iData][it1bin], phi045_135_eta[iData][it1bin], phiAMO_eta[iData][it1bin]} ;
            TH1F *phis_pi0[7] =  {phi000_pi0[iData][it1bin], phi045_pi0[iData][it1bin], phi090_pi0[iData][it1bin], 
                                    phi135_pi0[iData][it1bin], phi000_090_pi0[iData][it1bin], phi045_135_pi0[iData][it1bin], phiAMO_pi0[iData][it1bin]};
            vector<float> eta_psig;
            vector<float> eta_psig_err;
            vector<float> pi0_psig;
            vector<float> pi0_psig_err;

        
            // *****************************
            // Fitting asymmetry for eta
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(2,1);
            allCanvases->cd(1);
            fit_asym->SetParameters(pol_frac_090,pol_frac_000,0,Phi0_0_90[iData]);
            fit_asym->FixParameter(0,pol_frac_090);
            fit_asym->FixParameter(1,pol_frac_000);
            fit_asym->SetParLimits(2,-1,1),
            fit_asym->FixParameter(3,Phi0_0_90[iData]);

            if ((phi090_eta[iData][it1bin]->GetEntries()>minEntriesInPhiHists)*(phi000_eta[iData][it1bin]->GetEntries()>minEntriesInPhiHists)){
                fitPointer = asymmetry000_090_eta->Fit(fit_asym,fitOption.c_str());
                if ( (int)fitPointer!=0 )
                    cout << "A fit was not successful! exiting..." << endl;
                asymmetries_000_eta = fit_asym->GetParameter(2);
                asymmetries_000_eta_err = fit_asym->GetParError(2);
                asymmetry000_090_eta_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_000_eta = -999;
                asymmetries_000_eta_err = -999;
                asymmetry000_090_eta_chiPerDOF = -999;
            }
            asymmetry000_090_eta->Draw("SAME");
            asymmetry000_090_eta->SetAxisRange(-1,1,"Y");
        
            allCanvases->cd(2);
            fit_asym->SetParameters(pol_frac_135,pol_frac_045,0.5,Phi0_45_135[iData]);
            fit_asym->FixParameter(0,pol_frac_135);
            fit_asym->FixParameter(1,pol_frac_045);
            fit_asym->FixParameter(3,Phi0_45_135[iData]);
        
            if ((phi135_eta[iData][it1bin]->GetEntries()>minEntriesInPhiHists)*(phi045_eta[iData][it1bin]->GetEntries()>minEntriesInPhiHists)){
                fitPointer = asymmetry045_135_eta->Fit(fit_asym,fitOption.c_str());
                if ( (int)fitPointer!=0 )
                    cout << "A fit was not successful! exiting..." << endl;
                asymmetries_045_eta = fit_asym->GetParameter(2);
                asymmetries_045_eta_err = fit_asym->GetParError(2);
                asymmetry045_135_eta_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_045_eta = -999;
                asymmetries_045_eta_err = -999;
                asymmetry045_135_eta_chiPerDOF = -999;

            }
            asymmetry045_135_eta->Draw("SAME");
            asymmetry045_135_eta->SetAxisRange(-1,1,"Y");
            allCanvases->SaveAs((fitSaveLoc+"_eta"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"."+fileType).c_str()); 

            // *****************************
            // Fitting shifted cos - phi for eta
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(4,2);
            counter=0;
            for (auto phi: phis_eta){
            	allCanvases->cd(counter+1);
            	phi->SetTitle((polnames[counter]).c_str());
            	p0 = phi->GetEntries()/phi->GetNbinsX();
                if (counter<6){
            	    fit_sc->SetParameters(p0,0.1,poloffset[counter]);
                    fit_sc->FixParameter(2,poloffset[counter]);
                }
                else{ // index 6 corresponds to AMO dataset
                    fit_sc->SetParameters(p0,0.1,0);
                    fit_sc->SetParLimits(2,-1,1);
                }
            	fitStatus = phi->Fit(fit_sc,fitOption.c_str());
                if ( fitStatus!=0 ){
                    cout << "A fit was not successful! exiting..." << endl;
                    exit;
                }
            	eta_psig.push_back(fit_sc->GetParameter(1));
            	eta_psig_err.push_back(fit_sc->GetParError(1));
            	phi->Draw("SAME");
            	++counter;
            }
            allCanvases->SaveAs((fitSaveLoc+"_eta"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"_phis."+fileType).c_str()); 
            
            // *****************************
            // Fitting asymmetry for pi0
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(2,1);
            allCanvases->cd(1);
            fit_asym->SetParameters(pol_frac_090,pol_frac_000,0.5,Phi0_0_90[iData]);
            fit_asym->FixParameter(0,pol_frac_090);
            fit_asym->FixParameter(1,pol_frac_000);
            fit_asym->FixParameter(3,Phi0_0_90[iData]);
            if ((phi090_pi0[iData][it1bin]->GetEntries()>minEntriesInPhiHists)*(phi000_pi0[iData][it1bin]->GetEntries()>minEntriesInPhiHists)){
                fitPointer = asymmetry000_090_pi0->Fit(fit_asym,fitOption.c_str());
                if ( (int)fitPointer!=0 )
                    cout << "A fit was not successful! exiting..." << endl;
                asymmetries_000_pi0 = fit_asym->GetParameter(2);
                asymmetries_000_pi0_err = fit_asym->GetParError(2);
                asymmetry000_090_pi0_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_000_pi0 = -999;
                asymmetries_000_pi0_err = -999;
                asymmetry000_090_pi0_chiPerDOF = -999;
            }
            asymmetry000_090_pi0->Draw("SAME");
            asymmetry000_090_pi0->SetAxisRange(-1,1,"Y");
            
            allCanvases->cd(2);
            fit_asym->SetParameters(pol_frac_135,pol_frac_045,0.5,Phi0_45_135[iData]);
            fit_asym->FixParameter(0,pol_frac_135);
            fit_asym->FixParameter(1,pol_frac_045);
            fit_asym->FixParameter(3,Phi0_45_135[iData]);
            if ((phi135_pi0[iData][it1bin]->GetEntries()>minEntriesInPhiHists)*(phi045_pi0[iData][it1bin]->GetEntries()>minEntriesInPhiHists)){
                fitPointer = asymmetry045_135_pi0->Fit(fit_asym,fitOption.c_str());
                if ( (int)fitPointer!=0 )
                    cout << "A fit was not successful! exiting..." << endl;
                asymmetries_045_pi0 = fit_asym->GetParameter(2);
                asymmetries_045_pi0_err = fit_asym->GetParError(2);
                asymmetry045_135_pi0_chiPerDOF = fit_asym->GetChisquare()/fit_asym->GetNDF();
            }
            else{
                asymmetries_045_pi0 = -999; 
                asymmetries_045_pi0_err = -999; 
                asymmetry045_135_pi0_chiPerDOF = -999;
            }
            asymmetry045_135_pi0->Draw("SAME");
            asymmetry045_135_pi0->SetAxisRange(-1,1,"Y");
            allCanvases->SaveAs((fitSaveLoc+"_pi0"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"."+fileType).c_str()); 

            // *****************************
            // Fitting shifted cos - phi for pi0
            // *****************************
            allCanvases->Clear();
            allCanvases->Divide(4,2);
            counter=0;
            for (auto phi: phis_pi0){
            	allCanvases->cd(counter+1);
            	phi->SetTitle((polnames[counter]).c_str());
            	p0 = phi->GetEntries()/phi->GetNbinsX();
                if (counter<6){
            	    fit_sc->SetParameters(p0,0.1,poloffset[counter]);
                    fit_sc->FixParameter(2,poloffset[counter]);
                }
                else{
                    fit_sc->SetParameters(p0,0.1,0);
                    fit_sc->SetParLimits(2,-1,1);
                }
            	fitStatus = phi->Fit(fit_sc,fitOption.c_str());
                if ( fitStatus!=0 ){
                    cout << "A fit was not successful! exiting..." << endl;
                    exit;
                }
            	pi0_psig.push_back(fit_sc->GetParameter(1));
            	pi0_psig_err.push_back(fit_sc->GetParError(1));
            	phi->Draw("SAME");
            	++counter;
            }
            allCanvases->SaveAs((fitSaveLoc+"_pi0"+mint1s[it1bin]+"_"+dataSetTag2[iData]+"_phis."+fileType).c_str()); 
        
            // ********************************
            // Writing out the results
            // ********************************
            //float midt1=0.2*it1bin+0.1;
            float midt1=t1HalfWidth*2*(it1bin+0.5);
            *saveCsv << selections[0].variable << " " << dataSetTag2[iData] << " "
                     << midt1 << " " << t1HalfWidth << " "
                     << asymmetries_000_eta << " " << asymmetries_000_eta_err << " " 
                     << asymmetries_045_eta << " " << asymmetries_045_eta_err << " "
                     << asymmetries_000_pi0 << " " << asymmetries_000_pi0_err << " " 
                     << asymmetries_045_pi0 << " " << asymmetries_045_pi0_err << " "
                     << asymmetry000_090_eta_chiPerDOF << " " << asymmetry045_135_eta_chiPerDOF << " "
                     << asymmetry000_090_pi0_chiPerDOF << " " << asymmetry045_135_pi0_chiPerDOF << " "
                     << phi000_eta[iData][it1bin]->GetEntries() << " " << phi045_eta[iData][it1bin]->GetEntries() << " "
                     << phi090_eta[iData][it1bin]->GetEntries() << " " << phi135_eta[iData][it1bin]->GetEntries() << " "
                     << phi000_pi0[iData][it1bin]->GetEntries() << " " << phi045_pi0[iData][it1bin]->GetEntries() << " "
                     << phi090_pi0[iData][it1bin]->GetEntries() << " " << phi135_pi0[iData][it1bin]->GetEntries();
            for (int i=0; i<(int)eta_psig.size(); ++i){
                *saveCsv << " " << eta_psig[i] << " " << eta_psig_err[i] << " " << pi0_psig[i] << " " << pi0_psig_err[i];
            }
            for (auto selection : selections)
                *saveCsv << " " << selection.variable << " " << selection.minval << " " << selection.maxval;
            *saveCsv << endl;

        } // end data loop
    } // end t loop
    
    ///////////////////////////////////
    // CLEAN UP THE POINTERS FOR NEXT ITER
    ///////////////////////////////////
    for(int iData=0; iData<4; ++iData){
        for (int it1bin=0; it1bin<nt1bins; ++it1bin){
            delete phi000_eta[iData][it1bin];
            delete phi045_eta[iData][it1bin]; 
            delete phi090_eta[iData][it1bin];
            delete phi135_eta[iData][it1bin];
            delete phiAMO_eta[iData][it1bin];
            delete phi000_pi0[iData][it1bin];
            delete phi045_pi0[iData][it1bin];
            delete phi090_pi0[iData][it1bin];
            delete phi135_pi0[iData][it1bin];
            delete phiAMO_pi0[iData][it1bin];
        }
    }
}


void fitAsymmetryPlots(){
    ///////////////////////////////////
    // Storage variables to set branch addresses to 
    ///////////////////////////////////
    float phi_eta_lab;
    float phi_pi0_lab;
    float Mpi0;
    float Meta;
    float Mpi0p;
    float Metap;
    float mandelstam_t;
    float mandelstam_tp;
    float mandelstam_teta;
    float mandelstam_tpi0;
    float Mpi0eta;
    float AccWeight;
    float weightASBS;
    int BeamAngle;
    bool insideEllipse;

    map<string, float*> variable_map={
        {"phi_eta_lab",&phi_eta_lab},
        {"phi_pi0_lab",&phi_pi0_lab},
        {"Mpi0",&Mpi0},
        {"Meta",&Meta},
        {"Mpi0p",&Mpi0p},
        {"Metap",&Metap},
        {"mandelstam_t",&mandelstam_t},
        {"mandelstam_tp",&mandelstam_tp},
        {"mandelstam_teta",&mandelstam_teta},
        {"mandelstam_tpi0",&mandelstam_tpi0},
        {"Mpi0eta",&Mpi0eta},
        {"AccWeight",&AccWeight},
        {"weightASBS",&weightASBS},
    };

    vector<string> variables={"phi_eta_lab","phi_pi0_lab","Mpi0","Meta","Mpi0p","Metap",
        "mandelstam_t","mandelstam_tp","Mpi0eta","AccWeight","weightASBS","mandelstam_tpi0","mandelstam_teta"};
    map<string, vector<float>> array_variable_map;
    map<string, vector<int>> array_BeamAngles_map;

    // **************************************************
    // Define eta signal and sideband regions
    //   Basically non-existent pi0 sidebands, just cut signal out
    // **************************************************
    float eta_peak=0.548;
    float eta_std=0.0123;
    float etasbRL=eta_peak+5*eta_std;
    float etasbRR=eta_peak+11*eta_std;
    float etasbLL=eta_peak-11*eta_std;
    float etasbLR=eta_peak-5*eta_std;
    float etasigL=eta_peak-3*eta_std;
    float etasigR=eta_peak+3*eta_std;
    float pi0_peak=0.136;
    float pi0_std=0.006425;
    float pi0sigL=pi0_peak-3*pi0_std;
    float pi0sigR=pi0_peak+3*pi0_std;

    TH1F* dHist_Mpi0 = new TH1F("dHist_Mpi0","dHist_Mpi0;Mpi0",200,0.05,0.25);
    TH1F* dHist_Meta = new TH1F("dHist_Meta","dHist_Meta;Meta",300,0.25,0.85);
    TCanvas* c1 = new TCanvas("","",1440,900);
    TBox* box = new TBox();

    // *************************************
    // Loading the data in from the histogram 
    // *************************************
    TTree* tree;
    for (int iData=0; iData <3; ++iData){
        string dataFileName=folder+"degALL_data_"+dataSetTag[iData]+"_BA_treeFlat_DSelector.root";
    	cout << "LOADING ROOT FILE: " << dataFileName << endl; 
    	TFile *dataFile = new TFile(dataFileName.c_str());
        dataFile->GetObject("tree_4g_flat",tree);
        
        for(auto variable: variables){
            array_variable_map[dataSetTag[iData]+"_"+variable]=vector<float>{};
            tree->SetBranchAddress(variable.c_str(),variable_map[variable]);
        }
        tree->SetBranchAddress("BeamAngle",&BeamAngle);
        tree->SetBranchAddress("insideEllipse",&insideEllipse);
    
        Long64_t nentries=tree->GetEntries();
        //nentries=3000;
        cout << dataSetTag[iData] << " nentries: " << nentries << endl;
        bool selectPi0;
        bool selectEta;
        for (Long64_t ientry=0; ientry<nentries; ++ientry){
            tree->GetEntry(ientry);
            dHist_Mpi0->Fill(*variable_map["Mpi0"],*variable_map[weightVar]);
            dHist_Meta->Fill(*variable_map["Meta"],*variable_map[weightVar]);

            //selectPi0=(*variable_map["Mpi0"]>pi0sigL)*(*variable_map["Mpi0"]<pi0sigR);

            // ******* Select Eta Sidebands
            //selectEta=(*variable_map["Meta"]>etasbLL)*(*variable_map["Meta"]<etasbLR)||(*variable_map["Meta"]>etasbRL)*(*variable_map["Meta"]<etasbRR);
            // ******* Select Eta Right Sideband
            //selectEta=(*variable_map["Meta"]>etasbRL)*(*variable_map["Meta"]<etasbRR);
            // ******* Select Eta Left Sideband
            //selectEta=(*variable_map["Meta"]>etasbLL)*(*variable_map["Meta"]<etasbLR);
            // ******* Select Eta Signal
            //selectEta=(*variable_map["Meta"]>etasigL)*(*variable_map["Meta"]<etasigR);

            // ******* Select Eta Right Half Peak
            selectEta=(*variable_map["Meta"]>eta_peak);
            // ******* Select Eta Left Half Peak
            //selectEta=(*variable_map["Meta"]<eta_peak);
            
            // Do not select on Mpi0 nor Meta if we are going to do sideband subtraction also
            selectPi0=true; 
            //selectEta=true;

            if(selectPi0*selectEta){
            //if(insideEllipse){
                for(auto variable: variables){
                    array_variable_map[dataSetTag[iData]+"_"+variable].push_back(*variable_map[variable]); 
                }
                array_BeamAngles_map[dataSetTag[iData]].push_back(BeamAngle);
            }
        }
    }

    //gSystem->Exec("rm -rf fitAsymmetryPlots_results");
    gSystem->Exec("mkdir -p fitAsymmetryPlots_results");

    // Draw a diagnostic plot to denote the Meta and Mpi0 distributions
    dHist_Mpi0->SetMinimum(0);
    dHist_Meta->SetMinimum(0);
    c1->Divide(2,1);
    c1->cd(1);
    dHist_Mpi0->Draw();
    box->SetFillColorAlpha(kGreen+2,0.3);
    box->DrawBox(pi0sigL,0,pi0sigR,dHist_Mpi0->GetMaximum()*1.15);
    c1->cd(2);
    dHist_Meta->Draw();
    box->SetFillColorAlpha(kGreen+2,0.3);
    box->DrawBox(etasigL,0,etasigR,dHist_Meta->GetMaximum()*1.15);
    box->SetFillColorAlpha(kRed+1,0.3);
    box->DrawBox(etasbLL,0,etasbLR,dHist_Meta->GetMaximum()*1.15);
    box->DrawBox(etasbRL,0,etasbRR,dHist_Meta->GetMaximum()*1.15);
    c1->SaveAs(("fitAsymmetryPlots_results/mass_plots"+tag+".pdf").c_str());

    ofstream saveCsv;
    saveCsv.open(("fitAsymmetryPlots_results/results"+tag+".csv").c_str());
    saveCsv << "binVar data midt1 t1_err asym_000_eta asym_000_eta_err asym_045_eta asym_045_eta_err ";
    saveCsv << "asym_000_pi0 asym_000_pi0_err asym_045_pi0 asym_045_pi0_err ";
    saveCsv << "000090_eta_chiPerDOF 045135_eta_chiPerDOF 000090_pi0_chiPerDOF 045135_pi0_chiPerDOF ";
    saveCsv << "entries_000_eta entries_045_eta entries_090_eta entries_135_eta ";
    saveCsv << "entries_000_pi0 entries_045_pi0 entries_090_pi0 entries_135_pi0";
    for (int i=0; i<(int)(sizeof(polnames)/sizeof(polnames[0])); ++i){
        saveCsv << " " << polnames[i] << "_eta " << polnames[i] << "_eta_err " << polnames[i] << "_pi0 " << polnames[i] << "_pi0_err";
    }
    for (int i=1; i<3; ++i){
        saveCsv << " selection" << i << " minval" << i << " maxval" << i;
    }
    saveCsv << endl;

    ////////////////////////////////////////////////////////////
    //   BINNED IN U3
    ////////////////////////////////////////////////////////////
    string fitFolder="fitAsymmetryPlots_results/binned_u3"+tag+"/";
    gSystem->Exec(("mkdir -p "+fitFolder).c_str());
    
    const int ntbins=4;
    float mints[ntbins]={0.0,0.5,1.000,0.000};
    float maxts[ntbins]={0.5,1.0,100.0,100.0};
    //float mints[ntbins]={0.000};
    //float maxts[ntbins]={100.0};
    for(int it=0; it<ntbins; ++it){
        vector<criteria> selections={
            {"mandelstam_t",mints[it],maxts[it]},
            {"Mpi0p",1.4,10}
        };
        string fitSaveLoc=fitFolder+"asym_tB"+to_string(it); 
        constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,fitSaveLoc,&saveCsv);
    }

    ////////////////////////////////////////////////////////////
    //   BINNED IN PI0P
    ////////////////////////////////////////////////////////////
    fitFolder="fitAsymmetryPlots_results/binned_spi0p"+tag+"/";
    gSystem->Exec(("mkdir -p "+fitFolder).c_str());
    //float minPi0P[5] = {1.15, 1.4, 1.6, 1.75, 1.95};
    //float maxPi0P[5] = {1.4, 1.6, 1.75, 1.95, 2.4};
    float minPi0P[3] = {1.15, 1.4, 2.20};
    float maxPi0P[3] = {1.40, 2.2, 10.0};
    for(int it=0; it<3; ++it){
        vector<criteria> selections={
            {"Mpi0p",minPi0P[it],maxPi0P[it]},
            {"Mpi0p",0,10}
        };
        string fitSaveLoc=fitFolder+"asym_Mpi0pB"+to_string(it); 
        constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,fitSaveLoc,&saveCsv);
    }

    ////////////////////////////////////////////////////////////
    //   BINNED IN METAP
    ////////////////////////////////////////////////////////////
    fitFolder="fitAsymmetryPlots_results/binned_setap"+tag+"/";
    gSystem->Exec(("mkdir -p "+fitFolder).c_str());
    //float minEtaP[5] = {1.5, 1.65, 1.9, 2.2, 2.4};
    //float maxEtaP[5] = {1.65, 1.9, 2.2, 2.4, 2.8};
    float minEtaP[3] = {0.0, 2.1, 2.60};
    float maxEtaP[3] = {2.1, 2.6, 10.0};
    for(int it=0; it<3; ++it){
        vector<criteria> selections={
            {"Metap",minEtaP[it],maxEtaP[it]},
            {"Mpi0p",1.4,10}
        };
        string fitSaveLoc=fitFolder+"asym_Mpi0pB"+to_string(it); 
        constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,fitSaveLoc,&saveCsv);
    }

    ////////////////////////////////////////////////////////////
    //   BINNED IN MPI0ETA
    ////////////////////////////////////////////////////////////
    fitFolder="fitAsymmetryPlots_results/binned_s12"+tag+"/";
    gSystem->Exec(("mkdir -p "+fitFolder).c_str());
    //float minMpi0eta[5] = {1.65, 1.9, 2.15, 2.4, 2.65};
    //float maxMpi0eta[5] = {1.9, 2.15, 2.4, 2.65, 2.9};
    float minMpi0eta[7] = {1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8};
    float maxMpi0eta[7] = {1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
    for(int it=0; it<7; ++it){
        vector<criteria> selections={
            {"Mpi0eta",minMpi0eta[it],maxMpi0eta[it]},
            {"Mpi0p",1.4,10}
        };
        string fitSaveLoc=fitFolder+"asym_Mpi0etaB"+to_string(it); 
        constructAndFit(variable_map,array_variable_map,array_BeamAngles_map,selections,fitSaveLoc,&saveCsv);
    }

    saveCsv.close();
}













