// $Id: ntupleFilter.cc,v 1.10 2010/09/21 08:13:48 dlevans Exp $

#include <assert.h>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"
#include "../../CORE/Tools/goodrun.cc"
#include "TMath.h"
#include "../../MT2Analysis/MT2CORE/mt2tree.cc"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;


#include "Rtypes.h"
typedef ULong64_t uint64;

float getPx(float pt, float phi) {return pt*TMath::Cos(phi);}
float getPy(float pt, float phi) {return pt*TMath::Sin(phi);}
float getPz(float pt, float eta) {return pt*TMath::SinH(eta);}
float getE(float pt, float eta, float m) {return TMath::Sqrt(pt*pt + getPz(pt,eta)*getPz(pt,eta) + m*m);}

bool jetLepOverlap(float jeta, float jphi, float leta, float lphi);

void ntupleFilter (TChain *chain, const std::string &outfile, bool printPass=false, bool isData = true)  
{
    // set good run list
    set_goodrun_file("/home/users/fgolf/run2/ambulance_chasing/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_snt.txt");

    // output file and tree
    TFile *output =TFile::Open(outfile.c_str(), "RECREATE");
    assert(output != 0);
    TTree *newtree = 0;

    const long long max_tree_size = 20000000000000000LL;
    TTree::SetMaxTreeSize(max_tree_size);

    FILE *log = 0; //for keeping any output desired on selection
    if( printPass ) {
        size_t pos = outfile.find(".root");
        assert( pos != string::npos );
        std::string outcpy = outfile;
        log = fopen( outcpy.replace(pos, 5, "_run_lumi_event").c_str(), "w" );
    }

    TObjArray *listOfFiles = chain->GetListOfFiles();
    const uint64 nEventsChain = chain->GetEntries();
    uint64 nEventsTotal = 0;
    uint64 nEventsSelected = 0;

    mt2tree t;
    
    // file loop
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    bool first = true;
    int i_permille_old = 0;
    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile f(currentFile->GetTitle());
        TTree *tree = (TTree*)f.Get("mt2");
        // for the first file, clone the tree
        if ( first ) {
            newtree = chain->CloneTree(0);
            newtree->SetDirectory(output);
            first = false;

        }

        // init
        t.Init(newtree);
        t.Init(tree);

        // Event Loop
        const unsigned int nEvents = tree->GetEntries();
        for (unsigned int event = 0; event < nEvents; ++event, ++nEventsTotal) {
            int i_permille = (int)floor(10000 * nEventsTotal / float(nEventsChain));
            if (i_permille != i_permille_old) {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%5.2f%%"
                            "\033[0m\033[32m <---\033[0m\015", i_permille/100.);
                    fflush(stdout);
                }
                i_permille_old = i_permille;
            }

            t.GetEntry(event);
            //set condition to skip event
            if (t.nlep < 2) continue;
            
            //
            // start by finding all of Z candidates
            //
            std::vector<LorentzVector> Zp4s;
            std::vector<std::pair<int, int> > Zindices;
            for (int idx = 0; idx < t.nlep; idx++)
            {
                if (t.lep_pt[idx] < 20.) continue;
                if (fabs(t.lep_eta[idx]) > 2.4) continue;
                if (t.lep_relIso03[idx] > 0.4) continue;
                LorentzVector candLep_p4;
                candLep_p4.SetPx(getPx(t.lep_pt[idx],t.lep_phi[idx]));
                candLep_p4.SetPy(getPy(t.lep_pt[idx],t.lep_phi[idx]));
                candLep_p4.SetPx(getPz(t.lep_pt[idx],t.lep_eta[idx]));
                candLep_p4.SetE(getE(t.lep_pt[idx],t.lep_eta[idx],t.lep_mass[idx]));
          
                for (int sidx = idx+1; sidx < t.nlep; sidx++)
                {
                    if (t.lep_pt[sidx] < 20.) continue;
                    if (fabs(t.lep_eta[sidx]) > 2.4) continue;
                    if (t.lep_relIso03[sidx] > 0.4) continue;              
                    if (t.lep_charge[idx] * t.lep_charge[sidx] > 0) continue;
                    if (abs(t.lep_pdgId[idx]) != abs(t.lep_pdgId[sidx])) continue;
              
                    LorentzVector candLep2_p4;
                    candLep2_p4.SetPx(getPx(t.lep_pt[sidx],t.lep_phi[sidx]));
                    candLep2_p4.SetPy(getPy(t.lep_pt[sidx],t.lep_phi[sidx]));
                    candLep2_p4.SetPz(getPz(t.lep_pt[sidx],t.lep_eta[sidx]));
                    candLep2_p4.SetE(getE(t.lep_pt[sidx],t.lep_eta[sidx],t.lep_mass[sidx]));
                    
                    LorentzVector tmpZ_p4 = candLep_p4 + candLep2_p4;
                    if (tmpZ_p4.mass() < 70 || tmpZ_p4.mass() > 110) continue;

                    Zp4s.push_back(tmpZ_p4);
                    Zindices.push_back(std::make_pair(idx, sidx));
                }
            } // end loop over leptons

            //
            // if no Z candidates found, skip event
            //
            if (Zp4s.size() == 0) continue;           

            bool has_two_jets = false;
            for (int idx = 0; idx < Zp4s.size(); idx++)
            {
                int idx1 = Zindices.at(idx).first;
                int idx2 = Zindices.at(idx).second;

                std::vector<int> jet_indices;
                for (int jidx = 0; jidx < t.njet; jidx++)
                {
                    if (t.jet_pt[jidx] < 30) continue;
                    if (fabs(t.jet_eta[jidx]) > 2.5) continue;                        
                    if (jetLepOverlap(t.jet_eta[jidx], t.jet_phi[jidx], t.lep_eta[idx1], t.lep_phi[idx1])) continue;
                    if (jetLepOverlap(t.jet_eta[jidx], t.jet_phi[jidx], t.lep_eta[idx2], t.lep_phi[idx2])) continue;                    
                    jet_indices.push_back(jidx);
                }
                    
                if (jet_indices.size() > 1)
                {
                    has_two_jets = true;
                    break;
                }
            }

            if (t.nlep < 4 || t.met_pt < 50 || !has_two_jets) continue;
            
            ++nEventsSelected;
            if( printPass ) {
                //cout << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl;
                fprintf(log, "%i %i %i\n", t.run, t.lumi, t.evt );
            }

            // cms2.LoadAllBranches();

            // fill the new tree
            newtree->Fill();
        }
    }

    if( printPass ) {
        fprintf(log, "\nTotal events run on: %llu\n", nEventsTotal);
        fprintf(log, "Num events selected: %llu\n", nEventsSelected ); //need two fprintf statements bc of some gcc bug
        //cout << endl
        //	  << "Total events run on: " << nEventsTotal << endl
        //	  << "Num events selected: " << nEventsSelected << endl;
        //<< "Copy finished. Closing Files" << endl;
    }

    output->cd();
    newtree->Write();
    delete output;
}

bool jetLepOverlap(float jet_eta, float jet_phi, float lep_eta, float lep_phi)
{
    float deta = jet_eta - lep_eta;
    float dphi = TMath::ACos(TMath::Cos(jet_phi - lep_phi));
    float dr = TMath::Sqrt(std::pow(deta,2) + std::pow(dphi,2));
    return (dr < 0.4);
}
