#ifndef G4DEBUGPRECOMPOUNDMODEL_HH
#define G4DEBUGPRECOMPOUNDMODEL_HH

#include "G4PreCompoundModel.hh"
#include "G4SystemOfUnits.hh"


#include "G4AnalysisManager.hh"
#include "G4VPreCompoundTransitions.hh"
#include "G4PreCompoundEmission.hh"

#include "HadronicEventTuple.hh"

class G4DebugPreCompoundModel : public G4PreCompoundModel
{
public:

    explicit G4DebugPreCompoundModel(
        G4ExcitationHandler* ptr = nullptr)
      : G4PreCompoundModel(ptr)
    {
        auto ana = G4AnalysisManager::Instance();

        hEexcPrecoID = ana->CreateH1(
            "hEexcPreco", "Excitation energy of remnant (Before PreEq);E* (GeV);Events",
            200, 0., 20);

        hEexcDeexID = ana->CreateH1(
            "hEexcDeex", "Excitation energy of remnant (Before Equilibrium);E* (GeV);Events",
            200, 0., 20);

        hAEexcPrecoID = ana->CreateH1(
            "hAEexcPreco", "Excitation energy of remnant (Before PreEq);E*/A (MeV);Events",
            200, 0., 50);

        hAEexcDeexID = ana->CreateH1(
            "hAEexcDeex", "Excitation energy of remnant (Before Equilibrium);E*/A (MeV);Events",
            200, 0., 50);

        hAEexcDeexNoprecoID = ana->CreateH1(
            "hAEexcDeexNopreco", "Excitation energy of remnant (Before Equilibrium AND skipped preco);E*/A (MeV);Events",
            200, 0., 50);

        hAprecoID = ana->CreateH1(
            "hApreco", "Residual atomic mass A (Before PreEq);   A;Events",
            250, 0., 250);

        hZprecoID = ana->CreateH1(
            "hZpreco", "Residual element number Z (Before PreEq);   Z;Events",
            100, 0., 100);

        hAdeexID = ana->CreateH1(
            "hAdeex", "Residual atomic mass A (Before Equilibrium);   A;Events",
            250, 0., 250);

        hZdeexID = ana->CreateH1(
            "hZdeex", "Residual element number Z (Before Equilibrium);   Z;Events",
            100, 0., 100);

        tuple = std::make_unique<HadronicEventTuple>();
    }

    G4ReactionProductVector*
    DeExcite(G4Fragment& fragment) override;
    void MyDoIt(G4ReactionProductVector* result, G4Fragment& aFragment);

protected:
    void PerformEquilibriumEmission(const G4Fragment& aFragment,
                                           G4ReactionProductVector* result) const override;

public:
    G4int hEexcPrecoID;
    G4int hEexcDeexID;

    G4int hAEexcPrecoID;
    G4int hAEexcDeexID;
    G4int hAEexcDeexNoprecoID;

    G4int hAprecoID;
    G4int hZprecoID;
    G4int hAdeexID;
    G4int hZdeexID;

    std::unique_ptr<HadronicEventTuple> tuple;
    bool MyDebug=false;


};

#endif
