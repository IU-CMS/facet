

void filter(int pid, TString infile, TString outfile) {

	if (pid == 443) string par="J_Psi";
    if (pid == 553) string par="Upsilon";
	//int pid = 553; //eta 221, pi0 111, J/Psi(1S) 443, e 11, mu 13, Upsilon 553

    auto file = new TFile(infile);
    TTreeReader myReader("LHEF", file);

    gSystem->Load("libPhysics.so");
    TFile f(outfile, "recreate");
    f.SetCompressionLevel(1); 
    TTree tdat("tree", "tree");

    TLorentzVector *J_Psy_p4    = new TLorentzVector(0., 0., 0., 0.);

    tdat.Branch("J_Psy_p4", "TLorentzVector", &J_Psy_p4);

    TTreeReaderValue<int>    rv_Particle_size(myReader, "Particle_size");
    TTreeReaderArray<int>    ra_Particle_PID(myReader, "Particle.PID");
    TTreeReaderArray<int>    ra_Particle_Mother1(myReader, "Particle.Mother1");
    TTreeReaderArray<int>    ra_Particle_Mother2(myReader, "Particle.Mother2");

    TTreeReaderArray<double> ra_Particle_Px(myReader, "Particle.Px");
    TTreeReaderArray<double> ra_Particle_Py(myReader, "Particle.Py");
    TTreeReaderArray<double> ra_Particle_Pz(myReader, "Particle.Pz");
    TTreeReaderArray<double> ra_Particle_E(myReader, "Particle.E");
    TTreeReaderArray<double> ra_Particle_pT(myReader, "Particle.PT");
    TTreeReaderArray<double> ra_Particle_Eta(myReader, "Particle.Eta");
    TTreeReaderArray<double> ra_Particle_Phi(myReader, "Particle.Phi");
    TTreeReaderArray<double> ra_Particle_Mass(myReader, "Particle.M");

    int J_Psy(0);
    int n_events;
   
    n_events = myReader.GetEntries(1);


    for (int i_event = 0; i_event < n_events; ++i_event) {

        myReader.Next();

        for (int i_particle = 0; i_particle < *rv_Particle_size; ++i_particle) {

            //if the particle is a *pid
            if (abs(ra_Particle_PID.At(i_particle)) == pid) {

            //    if(ra_Particle_Eta.At(i_particle) > 6.17){
                J_Psy+=1;

                    J_Psy_p4->SetPxPyPzE(ra_Particle_Px.At(i_particle),
                                      ra_Particle_Py.At(i_particle),
                                      ra_Particle_Pz.At(i_particle),
                                      ra_Particle_E.At(i_particle));

                tdat.Fill();    
            //    }
            }
        }
    }

    tdat.Write();
    cout<<"# of J_Psy: " <<J_Psy<<endl;
        

}

int main() {
 return 0;
}

