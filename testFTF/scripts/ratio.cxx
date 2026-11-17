
#include <iostream>
#include <random>

#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

int main()
{
    constexpr int nEvents = 1'000'000;
    constexpr int maxX1 = 100;

    // Random number generator
    std::mt19937 generator(12345);

    // Histograms
    TH1D* hRatio = new TH1D(
        "hRatio",
        "Ratio x_{2}/x_{1};x_{2}/x_{1};Events",
        200, 0.0, 1.0
    );

    TH1D* hX1 = new TH1D(
        "hX1",
        "Initial integer x_{1};x_{1};Events",
        maxX1, 0.5, maxX1 + 0.5
    );

    TH1D* hX2 = new TH1D(
        "hX2",
        "Final integer x_{2};x_{2};Events",
        maxX1, 0.5, maxX1 + 0.5
    );

    for (int i = 0; i < nEvents; ++i)
    {
        // x1 uniformly distributed between 1 and 100
        int x1 = 1 + generator() % maxX1;

        // x2 uniformly distributed between 0 and x1
        int x2 = generator() % (x1 + 1);

        double ratio = static_cast<double>(x2) / x1;

        hX1->Fill(x1);
        hX2->Fill(x2);
        hRatio->Fill(ratio);
    }

    // Plot
    TCanvas* canvas = new TCanvas("canvas", "Integer ratios", 1200, 800);
    canvas->Divide(2, 2);

    canvas->cd(1);
    hX1->Draw();

    canvas->cd(2);
    hX2->Draw();

    canvas->cd(3);
    gPad->SetLogy();
    hRatio->Draw();

    canvas->SaveAs("integer_ratios.png");

    std::cout << "Generated " << nEvents << " events." << std::endl;

    return 0;
}

