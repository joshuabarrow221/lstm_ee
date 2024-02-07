#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Var.h"

#include "TestBeamAna/Cuts/ToFCuts.h"
#include "TestBeamAna/Cuts/Cuts.h"
#include "TestBeamAna/Vars/Vars.h"

#include "NDAna/numucc_inc/NumuCCInCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLegendEntry.h"
#include "TMarker.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace ana;

void new_kin_vars_plots()
{
  //SpectrumLoader loaderMC("p3.tbmc.caf.root");
  SpectrumLoader loaderMC("/pnfs/nova/persistent/production/concat/R20-11-25-prod5.1reco.a/nd/sumdecaf/ndphysics_contain/genie/prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_145_of_600.root");

  /*const Var kN3dpngPass([](const caf::SRProxy* sr)
   {
     float ret = 0;
     for(int i = 0; i < sr->vtx.elastic.fuzzyk.npng; i++)
       {
	 if(sr->vtx.elastic.fuzzyk.png[i].dir.z>0.9 && sr->vtx.elastic.fuzzyk.png[i].nhit<23)
	   {
	     ret++;
	   }
       }
     return ret;
     });*/

  /*const Cut kTruthExists([](const caf::SRProxy* sr)
   {
     bool ret = false;
     for(int i = 0; i < sr->vtx.elastic.fuzzyk.npng; i++)
       {
	 if(sr->vtx.elastic.fuzzyk.png[i].dir.z>0.9 && sr->vtx.elastic.fuzzyk.png[i].nhit<23)
	   {
	     if (sr->vtx.elastic.fuzzyk.png[i].truth.p.E>0&&sr->vtx.elastic.fuzzyk.png[i].truth.p.E<5000)
	       {
		 ret=true;
	       }
	     break;
	   }
       }
     return ret;
     });*/
  
  //const HistAxis kTBNPrngPassAxis("Num. 3D prongs passing cuts", Binning::Simple(10, 0, 10), kN3dpngPass);

  ///////////////////////////////////////////////////////
  //Variables
  //////////////////////////////////////////////////////
  const Var var_kNuEnergyLSTM([](const caf::SRProxy* sr)
   {
     float ret = 0;
     for(int i = 0; i < sr->vtx.elastic.fuzzyk.npng; i++)
       {
	 if(sr->vtx.elastic.fuzzyk.png[i].dir.z>0.9 && sr->vtx.elastic.fuzzyk.png[i].nhit<23)
	   {
	     ret++;
	   }
       }
     return ret;
     
   });

  ///////////////////////////////////////////////////////
  //Special Cuts
  //////////////////////////////////////////////////////
  
  const Cut kSlcTimeMC([](const caf::SRProxy* sr)
   {
     return (sr->slc.meantime > 215000 && sr->slc.meantime < 230000);
   });

  ///////////////////////////////////////////////////////
  //Axes and Spectra
  //////////////////////////////////////////////////////

  const HistAxis axis_kNuEnergyLSTM("LSTM Neutrino Energy (GeV)", Binning::Simple(100, 0., 0.), var_kNuEnergyLSTM);
  
  Spectrum spec_kNuEnergyLSTM(loaderMC, axis_kNuEnergyLSTM, kTruthExists&&kSlcTimeMC);
  Spectrum spec_kNuEnergyLSTM1(loaderMC, axis_kNuEnergyLSTM, /*kTruthExists&&*/kSlcTimeMC);
  
  loaderMC.Go();
  

  ///////////////////////////////////////////////////////
  //Histograms and ROOT
  //////////////////////////////////////////////////////
  
  TH1* hist_kNuEnergyLSTM = spec_kNuEnergyLSTM.ToTH1(spec_kNuEnergyLSTM.POT());
  TH1* hist_kNuEnergyLSTM1 = spec_kNuEnergyLSTM.ToTH1(spec_kNuEnergyLSTM1.POT());
  
  TFile fout("plots.root", "recreate");
  
  hnpngPass->Write("hist_kNuEnergyLSTM");
  hnpngPass1->Write("hist_kNuEnergyLSTM1");
  
  fout.Close();
}

