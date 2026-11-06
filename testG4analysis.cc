

#include "G4AnalysisManager.hh"

int main()
{
    auto ana = G4AnalysisManager::Instance();

    ana->OpenFile("debug.root");

    ana->CreateH1(
        "hEexcPreco",
        "Excitation energy of remnant (Before PreEq)",
        200, 0., 20.
    );
    ana->FillH1(0, 5);
    ana->Write();
    ana->CloseFile();


}
