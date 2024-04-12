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
    " prod_caf_R19-11-18-prod5reco.f_fd_genie_N1810j0211a_nonswap_fhc_nova_v08_period3_v1"
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
    /*   ORIGINAL VARIALBES   */
    //Original variable--should describe the scattering mode of the neutrino (coded as QE, RES, DIS, COH, etc.)
    { "mode",     SIMPLEVAR(mc.nu[0].mode) },
    //Original variable--should describe the interaction type of the neutrino (coded as CC, NC, etc.)
    //This may be deprecated...may now be nu.inttype, but I'm not sure...the range of values seems off
    //{ "inttype",  SIMPLEVAR(mc.nu[0].inttype) },
    { "interaction", SIMPLEVAR(mc.nu[0].interaction) },
    //Original variable--should describe the neutrino flavor (coded as numu, nue, nutau, etc.)
    { "flavor",   SIMPLEVAR(mc.nu[0].flavor) },
    //Original variable--should describe the incoming neutrino energy
    //We are good here as this is always the zeroth particle in the list from the GHEP record
    { "trueE",    SIMPLEVAR(mc.nu[0].E) },
    //Original variable--should describe the outgoing lepton momentum
    //We should be good here as this is always the zeroth particle in the list of final state particle from the GHEP record
    { "trueLepMom", SIMPLEVAR(mc.nu[0].prim[0].p) },
    //Original variable--should describe the outgoing lepton energy
    { "trueLepE", SIMPLEVAR(mc.nu[0].prim[0].p.E) },
    //Original variable--should describe the outgoing lepton angle
    { "trueLepCos", SIMPLEVAR(mc.nu[0].prim[0].p.costh) },
    //Original variable--should describe the outgoing hadronic system energy by conservation between neutrino and lepton
    { "trueHadE",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { return sr->mc.nu[0].E - sr->mc.nu[0].prim[0].p.E;}
      ) },
      /*   ORIGINAL VARIALBES   */
      
      
      /*     NEW  VARIALBES     */
      //This is the total energy of the hadronic and lepton system, overtly excluding the nuclear remnants
      //As I understand it, the cafs don't actually have anything directly contain anything in the
      //mc.nu.prim branches which are remnants, so this should be fine and actually redundant with
      //the above definition of trueHadE
      //I now think that this lambda function is not necessary and even incorrect given that the final
      //entry in the list should not actually be a nuclear remnant, as it supposedly doesn't exist in
      //cafs as a final state particle
      /*{"trueHadE_NoRemnants",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { return sr->mc.nu[0].E - sr->mc.nu[0].prim[0].p.E - sr->mc.nu[0].prim[sr->mc.nu[0].prim.size()-1].p.E; }
      ) },*/
    
    
    //////////////////////////////////////////////
    //Lets consider only the hadronic system first
    //////////////////////////////////////////////
    //This is the total x-direction momentum px of the hadronic system
    { "trueHadTotMomX",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { return sr->mc.nu[0].p.px - sr->mc.nu[0].prim[0].p.px;}
      ) },
    //This is the total y-direction momentum py of the hadronic system
    { "trueHadTotMomY",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { return sr->mc.nu[0].p.py - sr->mc.nu[0].prim[0].p.py;}
      ) },
    //This is the total z-direction momentum pz of the hadronic system
    { "trueHadTotMomZ",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { return sr->mc.nu[0].p.pz - sr->mc.nu[0].prim[0].p.pz;}
      ) },
    //This is the total momentum p_Tot of the hadronic system
    { "trueHadTotMom",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { return sqrt(  pow(sr->mc.nu[0].p.px - sr->mc.nu[0].prim[0].p.px,2) 
                        + pow(sr->mc.nu[0].p.py - sr->mc.nu[0].prim[0].p.py,2)
                        + pow(sr->mc.nu[0].p.pz - sr->mc.nu[0].prim[0].p.pz,2));}
      ) },

    ////////////////////////////////////////////////////
    //Now lets consider all particles in the final state
    ////////////////////////////////////////////////////
    //This is the total x-direction momentum px of the hadronic and lepton system, overtly including all final state particles
    { "trueTotMomX_all",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { double sum_momentum_x = 0.0;
            for (const auto& particle_i : sr->mc.nu.@prim.size())
                {
                    //This statement should always be true, but I'm including it just in case
                    if (abs(mc.nu[0].prim[particle_i].pdg) < 1000000000 /*this should cut out any weird stragglers which shouldn't be there*/)
                        {
                          sum_momentum_x += mc.nu[0].prim[particle_i].p.px;
                        }
                }
            return sum_momentum_x;
          }
      ) },
    //This is the total y-direction momentum py of the hadronic and lepton system, overtly including all final state particles
    { "trueTotMomY_all",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { double sum_momentum_y = 0.0;
            for (const auto& particle_i : sr->mc.nu.@prim.size())
                {
                    if (abs(mc.nu[0].prim[particle_i].pdg) < 1000000000 /*this should cut out any weird stragglers which shouldn't be there*/)
                        {
                          sum_momentum_y += mc.nu[0].prim[particle_i].p.py;
                        }
                }
            return sum_momentum_y;
          }
      ) },
    //This is the total z-direction momentum pz of the hadronic and lepton system, overtly including all final state particles
    { "trueTotMomZ_all",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { double sum_momentum_z = 0.0;
            for (const auto& particle_i : sr->mc.nu.@prim.size())
                {
                    if (abs(mc.nu[0].prim[particle_i].pdg) < 1000000000 /*this should cut out any weird stragglers which shouldn't be there*/)
                        {
                          sum_momentum_z += mc.nu[0].prim[particle_i].p.pz;
                        }
                }
            return sum_momentum_z;
          }
      ) },
    //This is the total momentum p_Tot of the hadronic and lepton system, overtly including all final state particles
    { "trueTotMom_all",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { double sum_momentum_tot = 0.0;
            double sum_momentum_x   = 0.0;
            double sum_momentum_y   = 0.0;
            double sum_momentum_z   = 0.0;
            for (const auto& particle_i : sr->mc.nu.@prim.size())
                {
                    if (abs(mc.nu[0].prim[particle_i].pdg) < 1000000000 /*this should cut out any weird stragglers which shouldn't be there*/)
                        {
                          sum_momentum_x += mc.nu[0].prim[particle_i].p.px;
                          sum_momentum_y += mc.nu[0].prim[particle_i].p.py;
                          sum_momentum_z += mc.nu[0].prim[particle_i].p.pz;
                          sum_momentum_tot = sqrt(pow(sum_momentum_x,2) + pow(sum_momentum_y,2) + pow(sum_momentum_z,2));
                        }
                }
            return sum_momentum_tot;
          }
      ) },

    /////////////////////////////////////////////////////////////////////////////////////
    //Now lets do the same thing, but exclude the neutrons from the final state particles
    /////////////////////////////////////////////////////////////////////////////////////
    //This is the total x-direction momentum px of the hadronic and lepton system, overtly excluding neutrons but keeping everything else
    { "trueTotMomX_no_neutrons",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { double sum_momentum_x_no_neutrons = 0.0;
            for (const auto& particle_i : sr->mc.nu.@prim.size())
                {
                    if (mc.nu[0].prim[particle_i].pdg != 2112 /*neutron*/
                        && abs(mc.nu[0].prim[particle_i].pdg) < 1000000000 /*this should cut out any weird stragglers which shouldn't be there*/)
                        {
                          sum_momentum_x_no_neutrons += mc.nu[0].prim[particle_i].p.px;
                        }   
                }
            return sum_momentum_x_no_neutrons;
          }
      ) },
    //This is the total y-direction momentum py of the hadronic and lepton system, overtly excluding neutrons but keeping everything else
    { "trueTotMomY_no_neutrons",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { double sum_momentum_y_no_neutrons = 0.0;
            for (const auto& particle_i : sr->mc.nu.@prim.size())
                {
                    if (mc.nu[0].prim[particle_i].pdg != 2112 /*neutron*/
                        && abs(mc.nu[0].prim[particle_i].pdg) < 1000000000 /*this should cut out any weird stragglers which shouldn't be there*/)
                        {
                          sum_momentum_y_no_neutrons += mc.nu[0].prim[particle_i].p.py;
                        }   
                }
            return sum_momentum_y_no_neutrons;
          }
      ) },
    { "trueTotMomZ_no_neutrons",
      Var(
          [] (const caf::SRProxy *sr) -> double
          { double sum_momentum_z_no_neutrons = 0.0;
            for (const auto& particle_i : sr->mc.nu.@prim.size())
                {
                    if (mc.nu[0].prim[particle_i].pdg != 2112 /*neutron*/
                        && abs(mc.nu[0].prim[particle_i].pdg) < 1000000000 /*this should cut out any weird stragglers which shouldn't be there*/)
                        {
                          sum_momentum_z_no_neutrons += mc.nu[0].prim[particle_i].p.pz;
                        }   
                }
            return sum_momentum_z_no_neutrons;
          }
      ) },
});
      /*       NEW  VARIALBES       */

const std::vector<std::pair<std::string, Var>> RECO_VAR_DEFS({
    /* ORIGINAL VARIABLES */
    { "numuRecoMuonE",    kNumuMuE2018   },
    { "numuRecoHadE",     kNumuHadE2018  },
    { "numuRecoE",        kNumuE2018     },
    { "nueRecoLepE",      kRecoEME       },
    { "nueRecoHadE",      kRecoHADE      },
    { "nueRecoE",         kNueEnergy2018 },
    /* ORIGINAL VARIABLES */


    /*   NEW  VARIABLES   */
    { "SlcVisE", SIMPLEVAR(rec.slc.calE) },

    /*   NEW  VARIABLES   */
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

void exporter_lstm_ee_fd_fhc_nonswap()
{
    CSVMaker maker(DATA, "dataset_lstm_ee_fd_fhc_nonswap.csv");
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
      /*     ORIGINAL VARIALBES     */