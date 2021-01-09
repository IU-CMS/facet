void table() {
    #include "TPad.h"
  

    //TString infile = "pyt8_10MEvt_13TeVpp_J_Psy_4Vec_Psy_pTHat_0_10.root";
    TString infile = "pyt8_10MEvt_13TeVpp_Upsilon_4Vec_Ups_pTHat_0_10.root";
    
    
    gSystem->Load("libPhysics.so");

    TFile f(infile);
    f.SetCompressionLevel(1);
    //auto file = new TFile(infile);

    TTree* tree =  (TTree*)f.Get("tree");
    TLorentzVector *particle_p4    = new TLorentzVector(0., 0., 0., 0.);
    
    //tree->SetBranchAddress("J_Psy_p4", &particle_p4);
    tree->SetBranchAddress("Upsilon_p4", &particle_p4);

    int nevents = tree->GetEntriesFast();
    cout << "nevents: " << nevents << endl; 

    
    double nXbin = 50;
    double nYbin = nXbin;


    double theta_min = 0.000001;
    double theta_max = 0.004; // 0,010 
    double delta_theta = (theta_max - theta_min)/nXbin;
    
    vector<double> v_theta;
    //fill v_theta
    for (double theta = theta_min; theta < theta_max + delta_theta; theta = theta + delta_theta){
        v_theta.push_back(theta);
    }

    
    double pt_min=0;
    double pt_max=10;
    double delta_pt = ( pt_max - pt_min ) / nYbin;

    double delta_P = ( pt_max/sin(theta_max) - pt_min/sin(theta_min) ) / nYbin;
    
    vector<double> v_pt;
    //fill v_pt
    for (double pt = pt_min; pt < pt_max + delta_pt; pt = pt + delta_pt){
        v_pt.push_back(pt);
    }


    
    //TH2D * hist = new TH2D ("J/ #Psi", "J/ #Psi distribution in #theta,pt plane", nXbin, &v_theta[0], nYbin,  &v_pt[0]);
    TH2D * hist = new TH2D ("#Upsilon", "#Upsilon distribution in #theta,pt plane", nXbin, &v_theta[0], nYbin,  &v_pt[0]);
    double count(0);

    //for over all events
    for (int i = 0; i < nevents; ++i){
        tree->GetEntry(i);
        
        if(particle_p4->Theta() < theta_min) continue;
        if(particle_p4->Pt() < pt_min)       continue;
        if(particle_p4->Theta() > theta_max) continue;
        if(particle_p4->Pt() > pt_max)       continue;

        hist->Fill(particle_p4->Theta(), particle_p4->Pt());
        count +=1;  
    }
    cout << "# of passed: " << count << endl;

    TCanvas *c1 = new TCanvas("c1","c1",200,10,700,500);

    hist->SetXTitle("#theta");
    hist->SetYTitle("Pt");
    hist->Draw("colz");
    
    ofstream file("Psi_table.dat");
    

    file<< "J/ #Psi theta-P cross section (pb) table" << endl;

    file <<"#theta_(rad) "<<"theta+Delta_theta "<<"P "<<"P+deltaP "<<"cross_section_(pb)"<< endl;
    
    for (double i = 0; i < nXbin; ++i){
        for (double j = 0; j < nYbin; ++j){

             double P = v_pt[j]/sin(v_theta[i]);
             double crx = hist->GetBinContent(i+1, j+1)*k;
             
             file << v_theta[i]                <<" "<<
                     v_theta[i] + delta_theta  <<" "<<
                     P                         <<" "<<
                     P + delta_P               <<" "<<
                     crx                            << endl;
        }
    }
    


    file.close();
    //cout<< "" <<hist->GetYaxis()->GetBinLowEdge(3)<<endl;
    

    c1->Update();
    TPaveStats * ps1 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");

    ps1->SetY1NDC(0.70); ps1->SetY2NDC(0.90);
    ps1->SetX1NDC(0.75); ps1->SetX2NDC(0.90);
    c1->Modified();
    
    c1->Print("Psi_theta_pt.pdf");

}