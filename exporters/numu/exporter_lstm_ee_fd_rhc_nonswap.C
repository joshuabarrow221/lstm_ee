#include <string>
#include <iostream>

#include "CAFAna/Analysis/CSVMaker.h"

#include "3FlavorAna/Cuts/NumuCuts2018.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "CAFAna/Weights/GenieWeights.h"
#include "3FlavorAna/Vars/LSTMEHelperVars.h"
#include "3FlavorAna/Vars/NueEnergy2020.h"
#include "3FlavorAna/Vars/NumuEFxs.h"
#include "3FlavorAna/Vars/NumuVars.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"

#include "StandardRecord/Proxy/SRProxy.h"

using namespace ana;

const std::string DATA =
    "dataset_def_name_newest_snapshot "
    " prod_caf_R19-11-18-prod5reco.f_fd_genie_N1810j0211a_nonswap_rhc_nova_v08_period4_v1"
    ;

const Cut kTrueEbelow7GeV = kTrueE < 7.0;

const Cut SanityCut(
    [] (const caf::SRProxy *sr)
    {
        return (sr->mc.nnu > 0) && (! sr->mc.nu[0].prim.empty());
    }
);

const Cut kNumuLoosePID(
    [] (const caf::SRProxy* sr)
    {
        return (
               (sr->sel.remid.pid > 0.5)
            && (sr->sel.cvnloosepreselptp.numuid > 0.5)
        );
    }
);

const Cut cut    =
       kIsNumuCC
    && (
           kNumuBasicQuality
        && kNumuContainFD2017
        && kNumuLoosePID
    )
    && kTrueEbelow7GeV
    && SanityCut;

const std::vector<std::pair<std::string, Var>> TRUTH_VAR_DEFS({
    { "mode",     SIMPLEVAR(mc.nu[0].mode) },
    { "trueE",    SIMPLEVAR(mc.nu[0].E) },
    { "trueLepE", SIMPLEVAR(mc.nu[0].prim[0].p.E) },
    { "trueHadE",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { return sr->mc.nu[0].E - sr->mc.nu[0].prim[0].p.E; }
      )
    }
});

const std::vector<std::pair<std::string, Var>> RECO_VAR_DEFS({
    { "numuRecoMuonE",    kNumuMuE2018   },
    { "numuRecoHadE",     kNumuHadE2018  },
    { "numuRecoE",        kNumuE2018     },
    { "nueRecoLepE",      kRecoEME       },
    { "nueRecoHadE",      kRecoHADE      },
    { "nueRecoE",         kNueEnergy2018 },
});

const std::vector<std::pair<std::string, Var>> EXTRA_VAR_DEFS({
    { "trkLen",      kTrkLength },
    { "remID",       kRemID     },
    { "cvn.numuid",  kCVNm      },
    { "cvn.nueid",   kCVNe      },
    { "cvn.nutauid", kCVNt      },
    { "cvn.ncid",    kCVNnc     },
    { "trkNHit",     kTrkNhits  },
    { "run",         SIMPLEVAR(hdr.run) },
    { "hadCalE",     SIMPLEVAR(energy.numu.hadcalE) },
    { "hadTrkE",     SIMPLEVAR(energy.numu.hadtrkE) },
});

const Weight weight = kPPFXFluxCVWgt * kXSecCVWgt2020;

void exporter_lstm_ee_fd_rhc_nonswap()
{
    CSVMaker maker(DATA, "dataset_lstm_ee_fd_rhc_nonswap.csv");
    maker.setPrecision(6);

    maker.addVars(kSliceVarDefs);

    maker.addVars(TRUTH_VAR_DEFS);
    maker.addVars(RECO_VAR_DEFS);
    maker.addVars(EXTRA_VAR_DEFS);

    maker.addMultiVars(kPng2dVarDefs);
    maker.addMultiVars(kPng3dVarDefs);

    maker.addVar("weight", VarFromWeight(weight));

    maker.SetSpillCut(kStandardSpillCuts);
    maker.setCut(cut);

    maker.Go();
}
