#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Var.h"

#include "TestBeamAna/Cuts/ToFCuts.h"
#include "TestBeamAna/Cuts/Cuts.h"
#include "TestBeamAna/Vars/Vars.h"

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

void make_tb_plots()
{
  //SpectrumLoader loaderTBMC("p3.tbmc.caf.root");
  SpectrumLoader loaderTBMC("/pnfs/nova/persistent/production/concat/R20-11-25-prod5.1reco.a/nd/sumdecaf/ndphysics_contain/genie/prod_sumdecaf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_145_of_600.root");

  const Var kN3dpngPass([](const caf::SRProxy* sr) {
      float ret = 0;
      for(int i = 0; i < sr->vtx.tbvtx[0].fuzzyk.npng; i++) {
	if(sr->vtx.tbvtx[0].fuzzyk.png[i].dir.z>0.9 && sr->vtx.tbvtx[0].fuzzyk.png[i].nhit<23) {
	  ret++;
	}
      }
      return ret;
    });

  const Cut kSlcTimeMC([](const caf::SRProxy* sr) {
      return (sr->slc.meantime > 84020 && sr->slc.meantime < 84170);
    });

  const Cut kTruthExists([](const caf::SRProxy* sr) {
      bool ret = false;
      for(int i = 0; i < sr->vtx.tbvtx[0].fuzzyk.npng; i++) {
	if(sr->vtx.tbvtx[0].fuzzyk.png[i].dir.z>0.9 && sr->vtx.tbvtx[0].fuzzyk.png[i].nhit<23) {
	  if (sr->vtx.tbvtx[0].fuzzyk.png[i].truth.p.E>0&&sr->vtx.tbvtx[0].fuzzyk.png[i].truth.p.E<5000) {
	    ret=true;
	  }
	  break;
	}
      }
      return ret;
    });

  const HistAxis kTBNPrngPassAxis("Num. 3D prongs passing cuts", Binning::Simple(10, 0, 10), kN3dpngPass);

  Spectrum npngPass(loaderTBMC, kTBNPrngPassAxis, kTruthExists&&kSlcTimeMC);

  loaderTBMC.Go();

  TH1* hnpngPass = npngPass.ToTH1(npngPass.POT());

  TFile fout("tbprotons_plots.root", "recreate");

  hnpngPass->Write("npngPass");

  fout.Close();
}

