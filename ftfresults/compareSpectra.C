
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

enum ParticleType {
    PROTON,
    NEUTRON,
    PI0,
    PIPLUS,
    PIMINUS,
    EM,
    FRAGMENT,
    OTHER,
    NPARTICLES
};

const std::vector<std::string> particleNames = {
    "Proton",
    "Neutron",
    "Pi0",
    "Pi+",
    "Pi-",
    "e^{#pm} + #gamma",
    "Fragments",
    "Other"
};

const std::vector<int> colors = {
    kRed + 1,
    kBlue + 1,
    kGreen + 2,
    kMagenta + 1,
    kOrange + 7,
    kCyan + 2,
    kBlack,
    kGray + 2
};

ParticleType classify(int pdg)
{
    switch (pdg) {
        case 2212:  return PROTON;
        case 2112:  return NEUTRON;
        case 111:   return PI0;
        case 211:   return PIPLUS;
        case -211:  return PIMINUS;
        case 11:
        case -11:
        case 22:    return EM;
        default:    break;
    }

    // Nuclear PDG codes: 10LZZZAAAI
    if (std::abs(static_cast<long long>(pdg)) >= 1000000000LL)
        return FRAGMENT;

    return OTHER;
}

int getA(int pdg)
{
    return (std::abs(pdg) / 10) % 1000;
}

int getZ(int pdg)
{
    return (std::abs(pdg) / 10000) % 1000;
}

// ------------------------------------------------------------
// Calculate energy flow from histogram
//
double getEnergyFlow(TH1D* h)
{
    double eflow = 0.;

    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        eflow += h->GetBinCenter(bin) * h->GetBinContent(bin);
    }

    return eflow;
}

// ------------------------------------------------------------
// Dedicated function to create canvas of A and Z distribution
//
void drawFragmentDistribution(
    TH1D* h1,
    TH1D* h2,
    const char* label1,
    const char* label2,
    const char* observable,
    const char* output,
    double xmax, 
    double nEvents
    )
{
    TCanvas* c = new TCanvas(
        Form("c_%s", observable),
        observable, 900, 800);

    // Upper pad
    TPad* upper = new TPad("upperFrag", "", 0., 0.30, 1., 1.);
    upper->SetBottomMargin(0.02);
    upper->SetLeftMargin(0.13);
    upper->SetRightMargin(0.05);
    upper->SetLogy();
    upper->Draw();
    upper->cd();

    h1->SetLineColor(kBlue + 1);
    h2->SetLineColor(kRed + 1);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);

    h1->SetTitle(
        Form(";Fragment %s;Counts", observable));

    h1->SetMinimum(1);

    h1->SetMaximum(
        5 * std::max(h1->GetMaximum(), h2->GetMaximum()));

    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.65, 0.78, 0.90, 0.91);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(h1, label1, "l");
    legend->AddEntry(h2, label2, "l");
    legend->Draw();
    // Fragment statistics
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.030);
    latex.SetTextFont(42);
    latex.SetTextAlign(12);

		double max_fragment_property = 5; // max Z or A for "light fragment"
		double nLight1 = h1->Integral(1, h1->FindFixBin(max_fragment_property));
		double nLight2 = h2->Integral(1, h2->FindFixBin(max_fragment_property));
		double nHeavy1 = h1->Integral(h1->FindFixBin(max_fragment_property), -1);
		double nHeavy2 = h2->Integral(h2->FindFixBin(max_fragment_property), -1);

		latex.SetTextColor(kBlue + 1);
    latex.DrawLatex(
        0.16, 0.85,
        Form("%s:  <N_{frag}> = %.1f,  <N_{%s>5}> = %.1f",
             label1, double(h1->GetEntries())/nEvents, observable, nHeavy1/nEvents));


		latex.SetTextColor(kRed + 1);
    latex.DrawLatex(
        0.16, 0.80,
        Form("%s:  <N_{frag}> = %.1f,  <N_{%s>5}> = %.1f",
             label2, double(h2->GetEntries())/nEvents, observable, nHeavy2/nEvents));
             
    // Lower ratio pad
    c->cd();

    TPad* lower = new TPad("lowerFrag", "", 0., 0., 1., 0.30);
    lower->SetTopMargin(0.03);
    lower->SetBottomMargin(0.30);
    lower->SetLeftMargin(0.13);
    lower->SetRightMargin(0.05);
    lower->Draw();
    lower->cd();

    TH1D* ratio = (TH1D*)h2->Clone(
        Form("ratio_%s", observable));

    ratio->SetDirectory(nullptr);
    ratio->Divide(h1);

    ratio->SetTitle("");
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(0.45);

    ratio->GetYaxis()->SetTitle(
        Form("%s / %s", label2, label1));

    ratio->GetXaxis()->SetTitle(
        Form("Fragment %s", observable));

    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.10);
    ratio->GetYaxis()->SetLabelSize(0.085);
    ratio->GetYaxis()->SetTitleOffset(0.55);

    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetXaxis()->SetTitleOffset(1.0);

    ratio->SetMinimum(0.);
    ratio->SetMaximum(4.);

    ratio->Draw("E1");

    TLine* unity = new TLine(0., 1., xmax, 1.);
    unity->SetLineColor(kRed);
    unity->SetLineStyle(2);
    unity->Draw("SAME");

    c->Print(output);

    delete legend;
    delete unity;
    delete ratio;
    delete upper;
    delete lower;
    delete c;
}

// ------------------------------------------------------------
// Read one file and fill the particle spectra
// ------------------------------------------------------------

Long64_t readFile(
    const char* filename,
    const char* label,
    std::vector<TH1D*>& histograms,
    TH1D* hA,
    TH1D* hZ)
{
    TFile* file = TFile::Open(filename, "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open " << filename << std::endl;
        return -1;
    }

    TTree* tree = nullptr;
    file->GetObject("events", tree);

    if (!tree) {
        std::cerr << "Cannot find tree 'events' in "
                  << filename << std::endl;
        file->Close();
        return -1;
    }

    TTreeReader reader(tree);

    TTreeReaderValue<std::vector<int>> pdg(
        reader, "ftfp_secPDG");

    TTreeReaderValue<std::vector<double>> kinetic(
        reader, "ftfp_secKinetic_i_MeV");

    Long64_t nEvents = 0;

    while (reader.Next()) {

        ++nEvents;

        const auto& pdgs = *pdg;
        const auto& energies = *kinetic;

        if (pdgs.size() != energies.size()) {
            std::cerr << "Inconsistent vector sizes at event "
                      << nEvents << std::endl;
            continue;
        }

        for (size_t i = 0; i < pdgs.size(); ++i) {

            const ParticleType type = classify(pdgs[i]);

            const double E = energies[i];

            if (!std::isfinite(E) || E < 0.)
                continue;

            histograms[type]->Fill(E);
						if (type == FRAGMENT) {

								const int A = getA(pdgs[i]);
								const int Z = getZ(pdgs[i]);

								hA->Fill(A);
								hZ->Fill(Z);
						}
        }
    }

    std::cout << label << ": processed "
              << nEvents << " events" << std::endl;

    file->Close();
    return nEvents;
}

// ------------------------------------------------------------
// Main comparison macro
// ------------------------------------------------------------

void compareSpectra(
    const char* file1 = "run1.root",
    const char* file2 = "run2.root",
    const char* label1 = "Run 1",
    const char* label2 = "Run 2",
    const char* output = "spectra_comparison.pdf")
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.045);
		gStyle->SetTickLength(-0.03, "Y");
		gStyle->SetLabelOffset(0.03,"Y");
		
    const int nBins = 100;
    const double Emin = 0.;
    const double Emax = 10000.;

    std::vector<TH1D*> h1(NPARTICLES);
    std::vector<TH1D*> h2(NPARTICLES);

    for (int i = 0; i < NPARTICLES; ++i) {

        const std::string name1 = "h1_" + std::to_string(i);
        const std::string name2 = "h2_" + std::to_string(i);

        h1[i] = new TH1D(
            name1.c_str(), "", nBins, Emin, Emax);

        h2[i] = new TH1D(
            name2.c_str(), "", nBins, Emin, Emax);

        h1[i]->SetDirectory(nullptr);
        h2[i]->SetDirectory(nullptr);

        h1[i]->Sumw2();
        h2[i]->Sumw2();
    }
    
		TH1D* hA1 = new TH1D("hA1", ";Mass number A;Fragments",
				                 50, 0., 250.);
		TH1D* hA2 = new TH1D("hA2", ";Mass number A;Fragments",
				                 50, 0., 250.);

		TH1D* hZ1 = new TH1D("hZ1", ";Atomic number Z;Fragments",
				                 20, 0., 100.);
		TH1D* hZ2 = new TH1D("hZ2", ";Atomic number Z;Fragments",
				                 20, 0., 100.);
                     
		Long64_t nEvents1 = readFile(
				file1, label1, h1, hA1, hZ1);

		Long64_t nEvents2 = readFile(
				file2, label2, h2, hA2, hZ2);

    // --------------------------------------------------------
    // PDF with one page per particle category
    // --------------------------------------------------------

    TCanvas* canvas = new TCanvas(
        "canvas", "Hadronic spectra", 900, 800);

    canvas->Print(Form("%s[", output));

    for (int i = 0; i < NPARTICLES; ++i) {

        canvas->Clear();

        // Upper spectrum pad
        TPad* upper = new TPad(
            "upper", "upper", 0., 0.30, 1., 1.);

        upper->SetBottomMargin(0.02);
        upper->SetLeftMargin(0.13);
        upper->SetRightMargin(0.05);
        upper->SetLogy();

        upper->Draw();
        upper->cd();

				// Print total particle multiplicities
				// including overflow and underflow!
				const double integral1 = h1[i]->Integral(0, nBins + 1);
				const double integral2 = h2[i]->Integral(0, nBins + 1);

				std::cout << particleNames[i] << std::endl;
				std::cout << "  " << label1 << ": " << integral1
								  << " particles ("
								  << integral1 / nEvents1 << " per event)" << std::endl;

				std::cout << "  " << label2 << ": " << integral2
								  << " particles ("
								  << integral2 / nEvents2 << " per event)" << std::endl;
        //if (integral1 > 0.)
        //    h1[i]->Scale(1. / integral1);

        //if (integral2 > 0.)
        //    h2[i]->Scale(1. / integral2);

        h1[i]->SetLineColor(kBlue + 1);
        h2[i]->SetLineColor(kRed + 1);

        h1[i]->SetLineWidth(2);
        h2[i]->SetLineWidth(2);

        h1[i]->SetTitle(
            Form("%s;Kinetic energy [MeV];"
                 "Counts",
                 particleNames[i].c_str()));

        h1[i]->SetMinimum(1.0);

        if (h2[i]->GetMaximum() > 0.)
            h1[i]->SetMaximum(
                5. * std::max(
                    h1[i]->GetMaximum(),
                    h2[i]->GetMaximum()));

        h1[i]->Draw("HIST");
        h2[i]->Draw("HIST SAME");

        TLegend* legend = new TLegend(
            0.65, 0.78, 0.90, 0.91);

        legend->SetBorderSize(0);
        legend->SetFillStyle(0);

        legend->AddEntry(h1[i], label1, "l");
        legend->AddEntry(h2[i], label2, "l");

        legend->Draw();
        
				// --------------------------------------------------------
				// Particle multiplicity and energy flow
				// --------------------------------------------------------

				const double N1 = h1[i]->Integral(1, nBins);
				const double N2 = h2[i]->Integral(1, nBins);

				const double Eflow1 = getEnergyFlow(h1[i]);
				const double Eflow2 = getEnergyFlow(h2[i]);

				TLatex latex;
				latex.SetNDC();
				latex.SetTextSize(0.030);
				latex.SetTextAlign(12);

				latex.SetTextColor(kBlue + 1);

				latex.DrawLatex(
						0.15, 0.87,
						Form("%s: <N> = %.1f, <E_{flow}> = %.3g MeV",
								 label1, N1/nEvents1, Eflow1/nEvents1));

				latex.SetTextColor(kRed + 1);

				latex.DrawLatex(
						0.15, 0.82,
						Form("%s: <N> = %.1f, <E_{flow}> = %.3g MeV",
								 label2, N2/nEvents2, Eflow2/nEvents2));

        // ----------------------------------------------------
        // Lower ratio pad
        // ----------------------------------------------------

        canvas->cd();

        TPad* lower = new TPad(
            "lower", "lower", 0., 0., 1., 0.30);

        lower->SetTopMargin(0.03);
        lower->SetBottomMargin(0.30);
        lower->SetLeftMargin(0.13);
        lower->SetRightMargin(0.05);

        lower->Draw();
        lower->cd();

        TH1D* ratio = (TH1D*)h2[i]->Clone(
            Form("ratio_%d", i));

        ratio->SetDirectory(nullptr);
        ratio->Divide(h1[i]);

        ratio->SetTitle("");
        ratio->SetLineColor(kBlack);
        ratio->SetMarkerColor(kBlack);
        ratio->SetMarkerStyle(20);
        ratio->SetMarkerSize(0.45);

        ratio->GetYaxis()->SetTitle( Form("%s / %s", label2, label1) );
        ratio->GetXaxis()->SetTitle(
            "Kinetic energy [MeV]");

        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetYaxis()->SetTitleSize(0.10);
        ratio->GetYaxis()->SetLabelSize(0.085);
        ratio->GetYaxis()->SetTitleOffset(0.55);

        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetLabelSize(0.10);
        ratio->GetXaxis()->SetTitleOffset(1.0);

        ratio->SetMinimum(0.);
        ratio->SetMaximum(2.);

        ratio->Draw("E1");

        TLine* unity = new TLine(
            Emin, 1., Emax, 1.);

        unity->SetLineColor(kRed);
        unity->SetLineStyle(2);
        unity->Draw("SAME");

        canvas->cd();
        canvas->Print(output);

        delete legend;
        delete unity;
        delete ratio;
        delete upper;
        delete lower;
    }
    // --------------------------------------------------------
		// Fragment mass and charge distributions
		// --------------------------------------------------------

		drawFragmentDistribution(
				hA1, hA2,
				label1, label2,
				"A",
				output,
				250.,
				nEvents1);

		drawFragmentDistribution(
				hZ1, hZ2,
				label1, label2,
				"Z",
				output,
				100.,
				nEvents2);

    canvas->Print(Form("%s]", output));

    std::cout << "Created " << output << std::endl;

    delete canvas;

    for (int i = 0; i < NPARTICLES; ++i) {
        delete h1[i];
        delete h2[i];
    }
}

