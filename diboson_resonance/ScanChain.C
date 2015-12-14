// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"
#include "Math/LorentzVector.h"
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

// CMS3
#include "../MT2Analysis/MT2CORE/mt2tree.cc"
#include "../CORE/Tools/goodrun.cc"
#include "../CORE/Tools/dorky/dorky.cc"

using namespace std;
// using namespace tas;
using namespace duplicate_removal;

float getPx(float pt, float phi) {return pt*TMath::Cos(phi);}
float getPy(float pt, float phi) {return pt*TMath::Sin(phi);}
float getPz(float pt, float eta) {return pt*TMath::SinH(eta);}
float getE(float pt, float eta, float m) {return TMath::Sqrt(pt*pt + getPz(pt,eta)*getPz(pt,eta) + m*m);}

bool jetLepOverlap(float jeta, float jphi, float leta, float lphi);

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

    //
    // set good run list
    //
    set_goodrun_file("/home/users/fgolf/run2/AmbulanceChasing/json/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_snt.txt");
    
    
    // Benchmark
    TBenchmark *bmark = new TBenchmark();
    bmark->Start("benchmark");

    // Example Histograms
    // TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    TFile outfile("histos.root","RECREATE");
    outfile.cd();

    //
    // event level histograms
    //
    TH1F* h_nleps   = new TH1F("h_nleps", "h_nleps", 6,-0.5,5.5);    
    TH1F* h_nZs     = new TH1F("h_nZs", "h_nZs", 4, -0.5, 3.5);
    
    TH1F* h_4l_mass    = new TH1F("h_4l_mass", "h_4l_mass", 200,0,1000);
    TH1F* h_4l_2z_mass    = new TH1F("h_4l_2z_mass", "h_4l_2z_mass", 200,0,1000);
    TH1F* h_4l_z1_mass = new TH1F("h_4l_z1_mass", "h_4l_z1_mass", 40,50,130);
    TH1F* h_4l_z2_mass = new TH1F("h_4l_z2_mass", "h_4l_z2_mass", 100,0,200);  
    TH1F* h_4l_met     = new TH1F("h_4l_met", "h_4l_met", 20,0,400);
    TH1F* h_4l_njets   = new TH1F("h_4l_njets", "h_4l_njets", 6,-0.5,5.5);
    TH1F* h_4l_nbtags  = new TH1F("h_4l_nbtags", "h_4l_nbtags", 4,-0.5,3.5);    
    
    TH1F* h_Zmet_mass    = new TH1F("h_Zmet_mass", "h_Zmet_mass", 200,0,1000);
    TH1F* h_Zmet_z_mass = new TH1F("h_Zmet_z1_mass", "h_Zmet_z1_mass", 40,50,130);
    TH1F* h_Zmet_met     = new TH1F("h_Zmet_met", "h_Zmet_met", 20,0,400);
    TH1F* h_Zmet_Zpt     = new TH1F("h_Zmet_Zpt", "h_Zmet_Zpt", 20,0,400);    
    TH1F* h_Zmet_njets   = new TH1F("h_Zmet_njets", "h_Zmet_njets", 6,-0.5,5.5);
    TH1F* h_Zmet_nbtags  = new TH1F("h_Zmet_nbtags", "h_Zmet_nbtags", 4,-0.5,3.5);    

    TH1F* h_2L2b_mass    = new TH1F("h_2L2b_mass", "h_2L2b_mass", 200,0,1000);
    TH1F* h_2L2b_2z_mass = new TH1F("h_2L2b_2z_mass", "h_2L2b_2z_mass", 200,0,1000);    
    TH1F* h_2L2b_zll_mass = new TH1F("h_2L2b_zll_mass", "h_2L2b_zll_mass", 40,50,130);
    TH1F* h_2L2b_zbb_mass = new TH1F("h_2L2b_zbb_mass", "h_2L2b_zbb_mass", 40,0,400);
    TH1F* h_2L2b_met     = new TH1F("h_2L2b_met", "h_2L2b_met", 20,0,400);
    TH1F* h_2L2b_njets   = new TH1F("h_2L2b_njets", "h_2L2b_njets", 6,-0.5,5.5);
    TH1F* h_2L2b_nbtags  = new TH1F("h_2L2b_nbtags", "h_2L2b_nbtags", 5,-0.5,4.5);    

    TH1F* h_2L2j_mass    = new TH1F("h_2L2j_mass", "h_2L2j_mass", 200,0,1000);
    TH1F* h_2L2j_2z_mass = new TH1F("h_2L2j_2z_mass", "h_2L2j_2z_mass", 200,0,1000);    
    TH1F* h_2L2j_zll_mass = new TH1F("h_2L2j_zll_mass", "h_2L2j_zll_mass", 40,50,130);
    TH1F* h_2L2j_zjj_mass = new TH1F("h_2L2j_zjj_mass", "h_2L2j_zjj_mass", 40,0,400);
    TH1F* h_2L2j_met     = new TH1F("h_2L2j_met", "h_2L2j_met", 20,0,400);
    TH1F* h_2L2j_njets   = new TH1F("h_2L2j_njets", "h_2L2j_njets", 6,-0.5,5.5);
    TH1F* h_2L2j_nbtags  = new TH1F("h_2L2j_nbtags", "h_2L2j_nbtags", 4,-0.5,3.5);    

    TH1F* h_2L2j_all_mass    = new TH1F("h_2L2j_all_mass", "h_2L2j_all_mass", 200,0,1000);
    TH1F* h_2L2j_all_2z_mass = new TH1F("h_2L2j_all_2z_mass", "h_2L2j_all_2z_mass", 200,0,1000);    
    TH1F* h_2L2j_all_zll_mass = new TH1F("h_2L2j_all_zll_mass", "h_2L2j_all_zll_mass", 40,50,130);
    TH1F* h_2L2j_all_zjj_mass = new TH1F("h_2L2j_all_zjj_mass", "h_2L2j_all_zjj_mass", 40,0,400);
    TH1F* h_2L2j_all_met     = new TH1F("h_2L2j_all_met", "h_2L2j_all_met", 20,0,400);
    TH1F* h_2L2j_all_njets   = new TH1F("h_2L2j_all_njets", "h_2L2j_all_njets", 6,-0.5,5.5);
    TH1F* h_2L2j_all_nbtags  = new TH1F("h_2L2j_all_nbtags", "h_2L2j_all_nbtags", 4,-0.5,3.5);    
    
    // samplehisto->SetDirectory(rootdir);
  
    // Loop over events to Analyze
    unsigned int nEventsTotal = 0;
    unsigned int nEventsChain = chain->GetEntries();
    unsigned int nEventsPassedZ = 0;    
    unsigned int nEventsPassed4L = 0;
    unsigned int nEventsPassed2Lmet = 0;
    unsigned int nEventsPassed2L2b = 0;
    unsigned int nEventsPassed2L2j = 0;        
    if( nEvents >= 0 ) nEventsChain = nEvents;
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    mt2tree t;

    bool print = true;
    
    // File Loop
    while ( (currentFile = (TFile*)fileIter.Next()) ) {

        // Get File Content
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("mt2");
        if(fast) TTreeCache::SetLearnEntries(10);
        if(fast) tree->SetCacheSize(128*1024*1024);
        t.Init(tree);
        
        // Loop over Events in current file
        if( nEventsTotal >= nEventsChain ) continue;

        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {
            
            // Get Event Content
            if( nEventsTotal >= nEventsChain ) continue;
            if(fast) tree->LoadTree(event);
            t.GetEntry(event);
            ++nEventsTotal;
            
            if (!goodrun(t.run, t.lumi)) continue;
            
            if(t.isData)
            {
                DorkyEventIdentifier id(t.run, t.evt, t.lumi);
                if (is_duplicate(id)) continue;
            }
            
            // Progress
            mt2tree::progress( nEventsTotal, nEventsChain );
            
            //
            // start with the 4 lepton search
            //
            h_nleps->Fill(std::min(t.nlep,5));
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
                    nEventsPassedZ++;
                }
            } // end loop over leptons

            //
            // if no Z candidates found, skip event
            //
            h_nZs->Fill(Zp4s.size());
            if (Zp4s.size() == 0) continue;           

            //
            // get good jets without overlap removal
            //
            std::vector<int> all_btag_indices;
            std::vector<int> all_lf_indices;
            for (int jidx = 0; jidx < t.njet; jidx++)
            {
                if (t.jet_pt[jidx] < 30) continue;
                if (fabs(t.jet_eta[jidx]) > 2.5) continue;
                if (t.jet_btagCSV[jidx] > 0.89)
                    all_btag_indices.push_back(jidx);
                else
                    all_lf_indices.push_back(jidx);
            }
            
            //
            // now find the rest of the di-lepton pairs
            //
            std::vector<std::vector<LorentzVector> > extraZ_p4s;
            std::vector<std::vector<std::pair<int, int> > > extraZ_indices;
            for (int zidx = 0; zidx < Zp4s.size(); zidx++)
            {                
                std::vector<LorentzVector> tmp_p4s;
                std::vector<std::pair<int, int> > tmp_indices;
                for (int idx = 0; idx < t.nlep; idx++)
                {
                    if (idx == Zindices.at(zidx).first || idx == Zindices.at(zidx).second) continue;

                    if (t.lep_pt[idx] < 20.) continue;
                    if (fabs(t.lep_eta[idx]) > 2.4) continue;
                    if (t.lep_relIso03[idx] > 0.4) continue;                        
                    LorentzVector candLep_p4;
                    candLep_p4.SetPx(getPx(t.lep_pt[idx],t.lep_phi[idx]));
                    candLep_p4.SetPy(getPy(t.lep_pt[idx],t.lep_phi[idx]));
                    candLep_p4.SetPz(getPz(t.lep_pt[idx],t.lep_eta[idx]));
                    candLep_p4.SetE(getE(t.lep_pt[idx],t.lep_eta[idx],t.lep_mass[idx]));
          
                    for (int sidx = idx+1; sidx < t.nlep; sidx++)
                    {
                        if (sidx == Zindices.at(zidx).first || sidx == Zindices.at(zidx).second) continue;
              
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
                        tmp_p4s.push_back(tmpZ_p4);
                        tmp_indices.push_back(std::make_pair(idx, sidx));
                    }                    
                }

                extraZ_p4s.push_back(tmp_p4s);
                extraZ_indices.push_back(tmp_indices);
            } // end loop over Z candidates

            //
            // ok, let's look for 4L events
            //
            for (int idx = 0; idx < Zp4s.size(); idx++)
            {
                if (extraZ_p4s.at(idx).size() != 0)
                {
                    LorentzVector secondZ(0,0,0,0);
                    std::pair<int, int> secondZ_indices = std::make_pair(-1,-1);
                    for (int zidx = 0; zidx < extraZ_p4s.at(idx).size(); zidx++)
                    {
                        secondZ = extraZ_p4s.at(idx).at(0);
                        secondZ_indices = extraZ_indices.at(idx).at(0);

                        LorentzVector new_p4 = Zp4s.at(idx) + secondZ;
                
                        //
                        // fill histograms
                        //
                        nEventsPassed4L++;            
                        h_4l_mass->Fill(std::min(new_p4.mass(),(float)999.9));
                        h_4l_z1_mass->Fill(Zp4s.at(idx).mass());
                        h_4l_z2_mass->Fill(std::min(secondZ.mass(),(float)199.9));
                        h_4l_met->Fill(std::min(t.met_pt,(float)399.9));
                        h_4l_njets->Fill(std::min(t.nJet30,5));
                        h_4l_nbtags->Fill(std::min(t.nBJet20,3));
                        
                        if (secondZ.mass() > 70 && secondZ.mass() < 110)
                            h_4l_2z_mass->Fill(std::min(new_p4.mass(),(float)999.9));
                    }
                }
            } // end 4L section
                
            //
            // beginning of 2L MET search section
            //
            if (t.met_pt > 75)
            {
                for (int idx = 0; idx < Zp4s.size(); idx++)
                {
                    nEventsPassed2Lmet++;
                    
                    LorentzVector Zp4 = Zp4s.at(idx);
                    
                    h_Zmet_z_mass ->Fill(Zp4.mass());
                    h_Zmet_Zpt    ->Fill(std::min(Zp4.pt(), (float)399.9));
                    h_Zmet_met    ->Fill(std::min(t.met_pt,(float)399.9));
                    h_Zmet_njets  ->Fill(std::min(t.nJet30,5));
                    h_Zmet_nbtags ->Fill(std::min(t.nBJet20,3));                

                    LorentzVector met_p4(getPx(t.met_pt,t.met_phi),getPy(t.met_pt,t.met_phi),0,getE(t.met_pt,0.,91.2));
                    LorentzVector new_p4 = met_p4 + Zp4;
                                                  
                    h_Zmet_mass->Fill(std::min(new_p4.mass(),(float)999.9));
                }
            } // end 2L+MET section
            
            //
            // beginning of 2L + 2b search section
            //
            if (all_btag_indices.size() > 1)
            {
                for (int idx = 0; idx < Zp4s.size(); idx++)
                {
                    int idx1 = Zindices.at(idx).first;
                    int idx2 = Zindices.at(idx).second;

                    std::vector<int> btag_indices;
                    for (int jidx = 0; jidx < all_btag_indices.size(); jidx++)
                    {                        
                        int bidx1 = all_btag_indices.at(jidx);
                        if (jetLepOverlap(t.jet_eta[bidx1], t.jet_phi[bidx1], t.lep_eta[idx1], t.lep_phi[idx1])) continue;
                        if (jetLepOverlap(t.jet_eta[bidx1], t.jet_phi[bidx1], t.lep_eta[idx2], t.lep_phi[idx2])) continue;                    
                        btag_indices.push_back(bidx1);
                    }
                    
                    if (btag_indices.size() < 2) continue;
                    for (int bidx = 0; bidx < btag_indices.size(); bidx++)
                    {
                        int bidx1 = btag_indices.at(bidx);
                        LorentzVector bjet1_p4(getPx(t.jet_pt[bidx1], t.jet_phi[bidx1]),
                                              getPy(t.jet_pt[bidx1], t.jet_phi[bidx1]),
                                              getPz(t.jet_pt[bidx1], t.jet_eta[bidx1]),
                                              getE(t.jet_pt[bidx1], t.jet_eta[bidx1], t.jet_mass[bidx1]));
                        
                        for (int jidx = bidx+1; jidx < btag_indices.size(); jidx++)
                        {
                            int bidx2 = btag_indices.at(jidx);
                            LorentzVector bjet2_p4(getPx(t.jet_pt[bidx2], t.jet_phi[bidx2]),
                                                  getPy(t.jet_pt[bidx2], t.jet_phi[bidx2]),
                                                  getPz(t.jet_pt[bidx2], t.jet_eta[bidx2]),
                                                  getE(t.jet_pt[bidx2], t.jet_eta[bidx2], t.jet_mass[bidx2]));
                            LorentzVector tmp_p4 = bjet1_p4 + bjet2_p4;
                            LorentzVector new_p4 = tmp_p4 + Zp4s.at(idx);

                            nEventsPassed2L2b++;
                            h_2L2b_mass->Fill(std::min(new_p4.mass(), (float)999.9));
                            h_2L2b_zll_mass->Fill(Zp4s.at(idx).mass());
                            h_2L2b_zbb_mass->Fill(std::min(tmp_p4.mass(),(float)399.9));
                            h_2L2b_met    ->Fill(std::min(t.met_pt,(float)399.9));
                            h_2L2b_njets  ->Fill(std::min(t.nJet30,5));
                            h_2L2b_nbtags ->Fill(std::min((int)btag_indices.size(),3));                

                            if (tmp_p4.mass() > 60 && tmp_p4.mass() < 120)
                                h_2L2b_2z_mass->Fill(std::min(new_p4.mass(),(float)999.9));
                        }
                    }
                }
            } // end 2L2b section

            //
            // beginning of 2L + 2j (not b-tagged) search section
            //
            if (all_lf_indices.size() > 1)
            {
                for (int idx = 0; idx < Zp4s.size(); idx++)
                {
                    int idx1 = Zindices.at(idx).first;
                    int idx2 = Zindices.at(idx).second;

                    std::vector<int> jet_indices;
                    for (int jidx = 0; jidx < all_lf_indices.size(); jidx++)
                    {
                        int lidx = all_lf_indices.at(jidx);                        
                        if (jetLepOverlap(t.jet_eta[lidx], t.jet_phi[lidx], t.lep_eta[idx1], t.lep_phi[idx1])) continue;
                        if (jetLepOverlap(t.jet_eta[lidx], t.jet_phi[lidx], t.lep_eta[idx2], t.lep_phi[idx2])) continue;                    
                        jet_indices.push_back(lidx);
                    }
                    
                    if (jet_indices.size() < 2) continue;
                    for (int jidx = 0; jidx < jet_indices.size(); jidx++)
                    {
                        int jidx1 = jet_indices.at(jidx);
                        LorentzVector jet1_p4(getPx(t.jet_pt[jidx1], t.jet_phi[jidx1]),
                                              getPy(t.jet_pt[jidx1], t.jet_phi[jidx1]),
                                              getPz(t.jet_pt[jidx1], t.jet_eta[jidx1]),
                                              getE(t.jet_pt[jidx1], t.jet_eta[jidx1], t.jet_mass[jidx1]));
                        
                        for (int sidx = jidx+1; sidx < jet_indices.size(); sidx++)
                        {
                            int jidx2 = jet_indices.at(sidx);
                            LorentzVector jet2_p4(getPx(t.jet_pt[jidx2], t.jet_phi[jidx2]),
                                                  getPy(t.jet_pt[jidx2], t.jet_phi[jidx2]),
                                                  getPz(t.jet_pt[jidx2], t.jet_eta[jidx2]),
                                                  getE(t.jet_pt[jidx2], t.jet_eta[jidx2], t.jet_mass[jidx2]));
                            LorentzVector tmp_p4 = jet1_p4 + jet2_p4;
                            LorentzVector new_p4 = tmp_p4 + Zp4s.at(idx);

                            nEventsPassed2L2j++;
                            h_2L2j_mass->Fill(std::min(new_p4.mass(), (float)999.9));
                            h_2L2j_zll_mass->Fill(Zp4s.at(idx).mass());
                            h_2L2j_zjj_mass->Fill(std::min(tmp_p4.mass(), (float)399.9));
                            h_2L2j_met    ->Fill(std::min(t.met_pt,(float)399.9));
                            h_2L2j_njets  ->Fill(std::min((int)jet_indices.size(),5));
                            h_2L2j_nbtags ->Fill(std::min(t.nBJet20,3));

                            if (tmp_p4.mass() > 60 && tmp_p4.mass() < 120)
                                h_2L2j_2z_mass->Fill(std::min(new_p4.mass(),(float)999.9));                            
                        }
                    }
                }
            } // end 2L2j section
            
            //
            // beginning of 2L + 2j (all) search section
            //
            if (t.njet > 1)
            {
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
                    
                    if (jet_indices.size() < 2) continue;
                    for (int jidx = 0; jidx < jet_indices.size(); jidx++)
                    {
                        int jidx1 = jet_indices.at(jidx);
                        LorentzVector jet1_p4(getPx(t.jet_pt[jidx1], t.jet_phi[jidx1]),
                                              getPy(t.jet_pt[jidx1], t.jet_phi[jidx1]),
                                              getPz(t.jet_pt[jidx1], t.jet_eta[jidx1]),
                                              getE(t.jet_pt[jidx1], t.jet_eta[jidx1], t.jet_mass[jidx1]));
                        
                        for (int sidx = jidx+1; sidx < jet_indices.size(); sidx++)
                        {
                            int jidx2 = jet_indices.at(sidx);
                            LorentzVector jet2_p4(getPx(t.jet_pt[jidx2], t.jet_phi[jidx2]),
                                                  getPy(t.jet_pt[jidx2], t.jet_phi[jidx2]),
                                                  getPz(t.jet_pt[jidx2], t.jet_eta[jidx2]),
                                                  getE(t.jet_pt[jidx2], t.jet_eta[jidx2], t.jet_mass[jidx2]));
                            LorentzVector tmp_p4 = jet1_p4 + jet2_p4;
                            LorentzVector new_p4 = tmp_p4 + Zp4s.at(idx);

                            nEventsPassed2L2j++;
                            h_2L2j_all_mass->Fill(std::min(new_p4.mass(), (float)999.9));
                            h_2L2j_all_zll_mass->Fill(Zp4s.at(idx).mass());
                            h_2L2j_all_zjj_mass->Fill(std::min(tmp_p4.mass(), (float)399.9));
                            h_2L2j_all_met    ->Fill(std::min(t.met_pt,(float)399.9));
                            h_2L2j_all_njets  ->Fill(std::min((int)jet_indices.size(),5));
                            h_2L2j_all_nbtags ->Fill(std::min(t.nBJet20,3));

                            if (tmp_p4.mass() > 60 && tmp_p4.mass() < 120)
                                h_2L2j_all_2z_mass->Fill(std::min(new_p4.mass(),(float)999.9));                            
                        }
                    }
                }
            } // end 2L2j section
            
            
        } // End loop over events
  
        // Clean Up
        delete tree;
        file->Close();
        delete file;
    }
    if ( nEventsChain != nEventsTotal ) {
        cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
    }
  
    // Example Histograms
    // samplehisto->Draw();
    outfile.Write();
    outfile.Close();

    // return
    bmark->Stop("benchmark");
    cout << endl;
    cout << nEventsTotal << " Events Processed" << endl;
    // cout << nEventsPassedGoodRunList << " Events passed good run list" << endl;
    // cout << nEventsPassedDuplicate << " Events passed duplicate removal" << endl;    
    cout << nEventsPassedZ << " Z candidates found" << endl;    
    cout << nEventsPassed4L << " candidates passing 4L selection" << endl;
    cout << nEventsPassed2Lmet << " candidates passing 2L+MET selection" << endl;
    cout << nEventsPassed2L2b << " candidates passing 2L+2b selection" << endl;
    cout << nEventsPassed2L2j << " candidates passing 2L+2j selection" << endl;    
    cout << "------------------------------" << endl;
    cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
    cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
    cout << endl;
    delete bmark;
    return 0;
}


bool jetLepOverlap(float jet_eta, float jet_phi, float lep_eta, float lep_phi)
{
    float deta = jet_eta - lep_eta;
    float dphi = TMath::ACos(TMath::Cos(jet_phi - lep_phi));
    float dr = TMath::Sqrt(std::pow(deta,2) + std::pow(dphi,2));
    return (dr < 0.4);
}
