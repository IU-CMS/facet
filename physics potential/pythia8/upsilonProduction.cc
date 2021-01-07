// Bu macro LHC 13 TeV de upsilon(bbar) ve J/Psi(ccbar) üretimi için yazıldı. 
// Üretilen ağır mezonlar daha sonra ağır dark photon üretimi için kullanılacak.
// Bu dark photonlar forward bölgede araştırılacak(theta ~1-4mrad).

using namespace std;

#include "Pythia8/Pythia.h"
using namespace Pythia8;

int main(int argc, char** argv) 
{


string lhe_file_name  = argv[1];
string eCM            = argv[2];
string pTHatMin       = argv[3];
string pTHatMax       = argv[4];
string seed           = argv[5];
int event_number = stoi(argv[6]);


Pythia pythia;

// Switch off generation of steps subsequent to the process level one.
// (These will not be stored anyway, so only steal time.)
pythia.readString("PartonLevel:all = off");

// LHC 13 TeV initialization.
pythia.readString("Beams:eCM =" + eCM); 


pythia.readString("Bottomonium:all = on");
pythia.readString("Charmonium:all = on");

//pythia.readString("phaseSpace:pTHatMin =" + pTHatMin);
pythia.readString("phaseSpace:pTHatMax =" + pTHatMax);


pythia.readString("Random:setSeed = on");
pythia.readString("Random:seed =" + seed);

//CommonSettings//
pythia.readString("Tune:preferLHAPDF = 2");
pythia.readString("Main:timesAllowErrors = 10000");
pythia.readString("Check:epTolErr = 0.01");
pythia.readString("Beams:setProductionScalesFromLHEF = off");
pythia.readString("SLHA:minMassSM = 1000.");
pythia.readString("ParticleDecays:limitTau0 = on");
pythia.readString("ParticleDecays:tau0Max = 10");
pythia.readString("ParticleDecays:allowPhotonRadiation = on");
//CommonSettings//

//CUEP8M1Settings//
pythia.readString("Tune:pp 14");
pythia.readString("Tune:ee 7");
pythia.readString("MultipartonInteractions:pT0Ref=2.4024");
pythia.readString("MultipartonInteractions:ecmPow=0.25208");
pythia.readString("MultipartonInteractions:expPow=1.6");
//CUEP8M1Settings//


pythia.init();



// Create an LHAup object that can access relevant information in pythia.
LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);

// Open a file on which LHEF events should be stored, and write header.
myLHA.openLHEF(lhe_file_name);

// Store initialization info in the LHAup object.
myLHA.setInit();

// Write out this initialization info on the file.
myLHA.initLHEF();

  // Loop over events.
  for (int i = 0; i < event_number; ++i) {

      // Generate an event.
      pythia.next();

      // Store event info in the LHAup object.
      myLHA.setEvent();
      myLHA.eventLHEF();

  }

  // Statistics: full printout.
  pythia.stat();

  // Update the cross section info based on Monte Carlo integration during run.
  myLHA.updateSigma();

  // Write endtag. Overwrite initialization info with new cross sections.
  myLHA.closeLHEF(true);

  cout<< "\ncross section: " << pythia.info.sigmaGen() * 1e+9 << " (pb)" << endl;

  return 0;
}
