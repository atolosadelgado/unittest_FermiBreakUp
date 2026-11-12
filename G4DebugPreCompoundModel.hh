
#include "G4PreCompoundModel.hh"
#include "G4SystemOfUnits.hh"


#include "G4AnalysisManager.hh"

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
    auto ana = G4AnalysisManager::Instance();
    ana->FillH1(hEexcPrecoID, U / GeV);
    ana->FillH1(hAEexcPrecoID, U / A / MeV);
    ana->FillH1(hAprecoID, A);
    ana->FillH1(hZprecoID, Z);

    ana->FillNtupleIColumn( tuple->ids.AremPreco, A );
    ana->FillNtupleIColumn( tuple->ids.ZremPreco, Z );
    ana->FillNtupleDColumn( tuple->ids.EexcPreco, U/MeV );
    ana->FillNtupleIColumn( tuple->ids.NexcPreco, aFragment.GetNumberOfExcitons() );
    ana->FillNtupleIColumn( tuple->ids.NparPreco, aFragment.GetNumberOfParticles() );
    ana->FillNtupleIColumn( tuple->ids.NholPreco, aFragment.GetNumberOfHoles() );

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
        ana->FillH1(hAEexcDeexNoprecoID, U / A / MeV);
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
    G4double U = aFragment.GetExcitationEnergy();
    G4int A = aFragment.GetA_asInt();
    G4int Z = aFragment.GetZ_asInt();
    auto ana = G4AnalysisManager::Instance();
    ana->FillH1(hEexcDeexID, U/GeV);
    ana->FillH1(hAEexcDeexID, U/MeV/A);
    ana->FillH1(hAdeexID, A);
    ana->FillH1(hZdeexID, Z);


    ana->FillNtupleIColumn( tuple->ids.AremDeex, A );
    ana->FillNtupleIColumn( tuple->ids.ZremDeex, Z );
    ana->FillNtupleDColumn( tuple->ids.EexcDeex, U/MeV );
    ana->FillNtupleIColumn( tuple->ids.NexcDeex, aFragment.GetNumberOfExcitons() );
    ana->FillNtupleIColumn( tuple->ids.NparDeex, aFragment.GetNumberOfParticles() );
    ana->FillNtupleIColumn( tuple->ids.NholDeex, aFragment.GetNumberOfHoles() );

    G4PreCompoundModel::PerformEquilibriumEmission(
        aFragment, result);
}
