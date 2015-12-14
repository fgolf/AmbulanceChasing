#include <iostream>
#include <string>
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"

void dumpPlots (std::string fname)
{
    TCanvas c1("c1","c1",600,400);
    gStyle->SetOptStat("emrou");

    TFile file(fname.c_str());
    file.cd();
    TList *flist = file.GetListOfKeys();
    for (auto item : *flist)
    {
        TObject *dobj = file.Get(item->GetName());
        if (dobj->InheritsFrom(TH1::Class()))
        {
            TString hname(dobj->GetName());            
            if (hname.Contains("2L2b_2z_mass") || hname.Contains("2L2j_2z_mass") || hname.Contains("2L2j_all_2z_mass"))
            {
                TH1F* htemp = (TH1F*)dobj;
                TH1F* hobj = (TH1F*)htemp->Rebin(4);
                hobj->Draw();
                c1.SetLogy(0);
                c1.Print(Form("plots/%s.pdf", dobj->GetName()));
                c1.Print(Form("plots/%s.png", dobj->GetName()));
                c1.SetLogy(1);
                c1.Print(Form("plots/%s_log.pdf", dobj->GetName()));
                c1.Print(Form("plots/%s_log.png", dobj->GetName()));
            }
            else
            {
                dobj->Draw();
                c1.SetLogy(0);
                c1.Print(Form("plots/%s.pdf", dobj->GetName()));
                c1.Print(Form("plots/%s.png", dobj->GetName()));
                c1.SetLogy(1);
                c1.Print(Form("plots/%s_log.pdf", dobj->GetName()));
                c1.Print(Form("plots/%s_log.png", dobj->GetName()));
            }
        }
    }
}
