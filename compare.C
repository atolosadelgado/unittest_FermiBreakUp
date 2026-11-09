
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>
#include <vector>
#include <string>

void compare(const char* file1Name,
             const char* file2Name,
             int rebin = 1)
{
    if (rebin < 1) {
        std::cerr << "Error: rebin must be >= 1" << std::endl;
        return;
    }

    TFile file1(file1Name, "READ");
    TFile file2(file2Name, "READ");

    if (file1.IsZombie() || file2.IsZombie()) {
        std::cerr << "Error opening input files." << std::endl;
        return;
    }

    const std::vector<std::string> histNames = {
        "hEexcPreco",
        "hEexcDeex",
        "hAEexcPreco",
        "hAEexcDeex",
        "hAEexcDeexNopreco",
        "hApreco",
        "hZpreco",
        "hAdeex",
        "hZdeex"
    };

    gStyle->SetOptStat(0);

    for (const auto& name : histNames) {

        TH1D* h1 = dynamic_cast<TH1D*>(file1.Get(name.c_str()));
        TH1D* h2 = dynamic_cast<TH1D*>(file2.Get(name.c_str()));

        if (!h1 || !h2) {
            std::cerr << "Skipping " << name
                      << ": histogram missing in one of the files."
                      << std::endl;
            continue;
        }

        // Clone histograms to avoid modifying the input files.
        TH1D* hist1 = static_cast<TH1D*>(h1->Clone(
            (name + "_file1").c_str()));
        TH1D* hist2 = static_cast<TH1D*>(h2->Clone(
            (name + "_file2").c_str()));

        hist1->SetDirectory(nullptr);
        hist2->SetDirectory(nullptr);

        if (rebin > 1) {
            hist1->Rebin(rebin);
            hist2->Rebin(rebin);
        }

        // Canvas and pads.
        auto* canvas = new TCanvas(
            ("c_" + name).c_str(),
            name.c_str(), 800, 800);

        auto* pad1 = new TPad("pad1", "Histograms",
                              0.0, 0.30, 1.0, 1.0);
        auto* pad2 = new TPad("pad2", "Ratio",
                              0.0, 0.0, 1.0, 0.30);

        pad1->SetBottomMargin(0.02);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.30);

        pad1->Draw();
        pad2->Draw();

        // Upper pad.
        pad1->cd();

        hist1->SetLineColor(kBlue + 1);
        hist2->SetLineColor(kRed + 1);

        hist1->SetLineWidth(2);
        hist2->SetLineWidth(2);

        hist1->SetTitle(name.c_str());
        hist1->GetYaxis()->SetTitle("Events");
        hist1->GetXaxis()->SetLabelSize(0);

        // Normalize for shape comparison.
        //if (hist1->Integral() > 0)
        //    hist1->Scale(1.0 / hist1->Integral());

        //if (hist2->Integral() > 0)
        //    hist2->Scale(1.0 / hist2->Integral());

        hist1->SetMaximum(
            1.15 * std::max(hist1->GetMaximum(),
                            hist2->GetMaximum()));

        hist1->Draw("HIST");
        hist2->Draw("HIST SAME");

        auto* legend = new TLegend(0.50, 0.5, 0.90, 0.90);
        legend->AddEntry(hist1, file1Name, "l");
        legend->AddEntry(hist2, file2Name, "l");
        legend->Draw();

        // Lower pad: ratio.
				pad2->cd();

				auto* ratio = static_cast<TH1D*>(
						hist1->Clone((name + "_ratio").c_str()));

				ratio->SetDirectory(nullptr);
				ratio->Divide(hist2);

				ratio->SetTitle("");
				ratio->SetStats(0);
				ratio->SetLineColor(kBlack);
				ratio->SetMarkerColor(kBlack);
				ratio->SetMarkerStyle(20);
				ratio->SetMarkerSize(0.5);

				ratio->GetYaxis()->SetTitle("File 1 / File 2");
				ratio->GetYaxis()->SetNdivisions(505);
				ratio->GetYaxis()->SetTitleSize(0.10);
				ratio->GetYaxis()->SetTitleOffset(0.45);
				ratio->GetYaxis()->SetLabelSize(0.09);

				ratio->GetXaxis()->SetTitle(
						hist1->GetXaxis()->GetTitle());
				ratio->GetXaxis()->SetTitleSize(0.11);
				ratio->GetXaxis()->SetLabelSize(0.09);

				ratio->SetMinimum(0.0);
				ratio->SetMaximum(2.0);

				// Draw explicitly with axes and error bars.
				ratio->Draw("E1");

				auto* line = new TLine(
						ratio->GetXaxis()->GetXmin(), 1.0,
						ratio->GetXaxis()->GetXmax(), 1.0);

				line->SetLineStyle(2);
				line->Draw("SAME");

				pad2->Modified();
				pad2->Update();
				canvas->Modified();
				canvas->Update();

        canvas->Update();
    }
}

