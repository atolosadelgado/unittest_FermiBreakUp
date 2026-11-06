
#include "G4PreCompoundModel.hh"
#include "G4SystemOfUnits.hh"


#include "TH1D.h"
#include "TFile.h"

TH1D* hEexcFTFremnant = new TH1D("hEexcFTFremnant", "FTF wounded nucleus (after dummy INC) excitation energy;E* (GeV);Events",
                        200, 0., 20);
TH1D* hEexcEquilibrium = new TH1D("hEexcEquilibrium", "Equilibrium excitation energy;E* (GeV);Events",
                        200, 0., 20);

TH1D* hAEexcFTFremnant = new TH1D("hAEexcFTFremnant", "FTF wounded nucleus (after dummy INC) excitation energy;E*/A (MeV);Events",
                        200, 0., 50);
TH1D* hAEexcEquilibrium = new TH1D("hAEexcEquilibrium", "Equilibrium excitation energy;E*/A (MeV);Events",
                        200, 0., 50);
TH1D* hAEexcEquilibriumOutsideLimit = new TH1D("hAEexcEquilibriumULimit", "Equilibrium excitation energy outside PreEquilibrium limits;E*/A (MeV);Events",
                        200, 0., 50);

class G4DebugPreCompoundModel : public G4PreCompoundModel
{
public:

    explicit G4DebugPreCompoundModel(
        G4ExcitationHandler* ptr = nullptr)
      : G4PreCompoundModel(ptr)
    {}

    G4ReactionProductVector*
    DeExcite(G4Fragment& fragment) override;
protected:
    void PerformEquilibriumEmission(const G4Fragment& aFragment,
                                           G4ReactionProductVector* result) const override;
};

#include "G4PreCompoundInterface.hh"
#include "G4VMultiFragmentation.hh"

G4ReactionProductVector*
G4DebugPreCompoundModel::DeExcite(G4Fragment& aFragment)
{
    if (!isInitialised)
    {
        InitialiseModel();
    }

    if (usePrecoInterface)
    {
        return fInterface->DeExcite(aFragment);
    }

    G4ReactionProductVector* result =
        new G4ReactionProductVector();

    G4int Z = aFragment.GetZ_asInt();
    G4int A = aFragment.GetA_asInt();
    G4double U = aFragment.GetExcitationEnergy();

    // FTF remnant entering precompound
    // G4cout << "\n===== FTF REMNANT =====\n"
    //        << "A = " << A << "\n"
    //        << "Z = " << Z << "\n"
    //        << "E* = " << U / MeV << " MeV\n"
    //        << "E*/A = " << U / A / MeV << " MeV/nucleon\n"
    //        << "=======================\n";
    hEexcFTFremnant->Fill(U / GeV);
    hAEexcFTFremnant->Fill(U / A / MeV);


    // if (!isActive)
    //     G4cout << "Skipping: inactive" << G4endl;
    //
    // if (Z < minZ && A < minA)
    //     G4cout << "Skipping: nucleus too light" << G4endl;
    //
    // if (U < fLowLimitExc * A)
    //     G4cout << "Skipping: excitation " << U/MeV << " below low limit. Limit per nucleon (MeV) = " << fLowLimitExc/MeV
    //           << G4endl;
    //
    // if (U > fHighLimitExc * A)
    //     G4cout << "Skipping: excitation " << U/MeV << " above high limit. Limit per nucleon (MeV) = " << fHighLimitExc/MeV
    //           << G4endl;
    //
    // if (aFragment.GetNumberOfLambdas() > 0)
    //     G4cout << "Skipping: hypernucleus" << G4endl;
    //
    // if (Z == 0 || Z >= A)
    //     G4cout << "Skipping: invalid Z/A" << G4endl;

    // Conditions to skip pre-compound and perform equilibrium emission
    if (!isActive || (Z < minZ && A < minA)
        || U < fLowLimitExc * A
        || U > fHighLimitExc * A
        || 0 < aFragment.GetNumberOfLambdas()
        || 0 == Z || Z >= A)
    {
        // G4cout << "===== EXCITATION HANDLER INPUT =====\n"
        //        << "A = " << A << "\n"
        //        << "Z = " << Z << "\n"
        //        << "E* = " << U / MeV << " MeV\n"
        //        << "E*/A = " << U / A / MeV << " MeV/nucleon\n"
        //        << "====================================\n";
        if (U < fLowLimitExc * A || U > fHighLimitExc * A)
          hAEexcEquilibriumOutsideLimit->Fill(U / A / MeV);
        PerformEquilibriumEmission(aFragment, result);
        return result;
    }

    // Apply multi-fragmentation above threshold
    G4FragmentVector* theTempResult;

    if (U > A * fMinEForMultiFrag)
    {
        G4cout << "===== MULTIFRAGMENTATION INPUT =====\n"
               << "A = " << A << "\n"
               << "Z = " << Z << "\n"
               << "E* = " << U / MeV << " MeV\n"
               << "E*/A = " << U / A / MeV << " MeV/nucleon\n"
               << "====================================\n";
        theTempResult = theMultiFrag->BreakItUp(aFragment);
    }
    else
    {
        theTempResult = new G4FragmentVector();
        theTempResult->push_back(new G4Fragment(aFragment));
    }

    // Apply pre-compound on the list of fragments
    for (auto& ptr : *theTempResult)
    {
        DoIt(result, *ptr);
        delete ptr;
    }

    delete theTempResult;

    return result;
}

void G4DebugPreCompoundModel::PerformEquilibriumEmission(
    const G4Fragment& aFragment,
    G4ReactionProductVector* result) const
{
    // G4cout << "\n===== EQUILIBRIUM EMISSION =====\n";
    // G4cout << "A = " << aFragment.GetA_asInt() << G4endl;
    // G4cout << "Z = " << aFragment.GetZ_asInt() << G4endl;
    // G4cout << "E* = "
    //        << aFragment.GetExcitationEnergy()/MeV
    //        << " MeV" << G4endl;
    G4double U_GeV = aFragment.GetExcitationEnergy()/GeV;
    G4int A = aFragment.GetA_asInt();
    hEexcEquilibrium->Fill(U_GeV);
    hAEexcEquilibrium->Fill(1000*U_GeV/A);
    G4PreCompoundModel::PerformEquilibriumEmission(
        aFragment, result);
}
