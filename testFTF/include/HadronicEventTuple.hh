
#ifndef HADRONIC_EVENT_TUPLE_HH
#define HADRONIC_EVENT_TUPLE_HH

#include <vector>
#include <numeric>
#include <string>

#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4AnalysisManager.hh"
#include "G4DynamicParticle.hh"
#include "G4Exception.hh"

class HadronicEventTuple
{
public:

    // Information about the nuclear remnant
    struct RemnantInfo {
        G4int A = -1;
        G4int Z = -1;
        G4double Eexc = -1.0;

        G4int Nexc = -1;
        G4int Npar = -1;
        G4int Nhol = -1;
    };

    HadronicEventTuple()
    {
        auto* man = G4AnalysisManager::Instance();

        // Create ntuple
        ntupleID = man->CreateNtuple(
            "events", "Hadronic events");

        // --------------------------------------------------
        // Remnant before precompound
        // --------------------------------------------------

        AremBeforePreco = man->CreateNtupleIColumn(
            ntupleID, "AremBeforePreco");

        ZremBeforePreco = man->CreateNtupleIColumn(
            ntupleID, "ZremBeforePreco");

        EexcBeforePreco = man->CreateNtupleDColumn(
            ntupleID, "EexcBeforePreco_MeV");

        NexcBeforePreco = man->CreateNtupleIColumn(
            ntupleID, "NexcBeforePreco");

        NparBeforePreco = man->CreateNtupleIColumn(
            ntupleID, "NparBeforePreco");

        NholBeforePreco = man->CreateNtupleIColumn(
            ntupleID, "NholBeforePreco");


        // --------------------------------------------------
        // Remnant before deexcitation
        // --------------------------------------------------

        AremBeforeDeex = man->CreateNtupleIColumn(
            ntupleID, "AremBeforeDeex");

        ZremBeforeDeex = man->CreateNtupleIColumn(
            ntupleID, "ZremBeforeDeex");

        EexcBeforeDeex = man->CreateNtupleDColumn(
            ntupleID, "EexcBeforeDeex_MeV");

        NexcBeforeDeex = man->CreateNtupleIColumn(
            ntupleID, "NexcBeforeDeex");

        NparBeforeDeex = man->CreateNtupleIColumn(
            ntupleID, "NparBeforeDeex");

        NholBeforeDeex = man->CreateNtupleIColumn(
            ntupleID, "NholBeforeDeex");


        // --------------------------------------------------
        // Secondary collections
        // --------------------------------------------------

        CreateSecondaryColumns(man, "ftf_", ftf);
        CreateSecondaryColumns(man, "ftfp_", ftfp);

        man->FinishNtuple(ntupleID);
    }


    // ======================================================
    // Public interface
    // ======================================================

    void SetRemnantBeforePrecompound(const RemnantInfo& info)
    {
        remnantBeforePreco = info;
    }

    void SetRemnantBeforeDeexcitation(const RemnantInfo& info)
    {
        remnantBeforeDeex = info;
    }


    void AddFTFSecondary(const G4DynamicParticle* particle)
    {
        AddSecondary(ftf, particle);
    }

    void AddHadronicSecondary(const G4DynamicParticle* particle)
    {
        AddSecondary(ftfp, particle);
    }


    void FillEvent()
    {
        auto* man = G4AnalysisManager::Instance();

        if (!CheckConsistency(ftf) || !CheckConsistency(ftfp)) {
            G4Exception(
                "HadronicEventTuple::FillEvent",
                "HadronicTuple001",
                JustWarning,
                "Inconsistent secondary vector sizes.");
            return;
        }

        // --------------------------------------------------
        // Remnant before precompound
        // --------------------------------------------------

        man->FillNtupleIColumn(
            ntupleID, AremBeforePreco, remnantBeforePreco.A);

        man->FillNtupleIColumn(
            ntupleID, ZremBeforePreco, remnantBeforePreco.Z);

        man->FillNtupleDColumn(
            ntupleID, EexcBeforePreco,
            remnantBeforePreco.Eexc / MeV);

        man->FillNtupleIColumn(
            ntupleID, NexcBeforePreco, remnantBeforePreco.Nexc);

        man->FillNtupleIColumn(
            ntupleID, NparBeforePreco, remnantBeforePreco.Npar);

        man->FillNtupleIColumn(
            ntupleID, NholBeforePreco, remnantBeforePreco.Nhol);


        // --------------------------------------------------
        // Remnant before deexcitation
        // --------------------------------------------------

        man->FillNtupleIColumn(
            ntupleID, AremBeforeDeex, remnantBeforeDeex.A);

        man->FillNtupleIColumn(
            ntupleID, ZremBeforeDeex, remnantBeforeDeex.Z);

        man->FillNtupleDColumn(
            ntupleID, EexcBeforeDeex,
            remnantBeforeDeex.Eexc / MeV);

        man->FillNtupleIColumn(
            ntupleID, NexcBeforeDeex, remnantBeforeDeex.Nexc);

        man->FillNtupleIColumn(
            ntupleID, NparBeforeDeex, remnantBeforeDeex.Npar);

        man->FillNtupleIColumn(
            ntupleID, NholBeforeDeex, remnantBeforeDeex.Nhol);


        // --------------------------------------------------
        // Secondary summary information
        // --------------------------------------------------

        FillSecondarySummary(man, ftf);
        FillSecondarySummary(man, ftfp);


        // Write one row per event
        man->AddNtupleRow(ntupleID);

        Clear();
    }


    void Clear()
    {
        remnantBeforePreco = RemnantInfo{};
        remnantBeforeDeex = RemnantInfo{};

        ClearSecondaryData(ftf);
        ClearSecondaryData(ftfp);
    }


private:

    // ======================================================
    // Internal data structures
    // ======================================================

    struct SecondaryData {

        // Scalar column IDs
        G4int nSecondaries;
        G4int sumKinetic;
        G4int sumMass;
        G4int sumTotal;

        // Vector column IDs
        G4int PDG;
        G4int charge;
        G4int baryonNumber;

        G4int kinetic;
        G4int mass;
        G4int total;

        G4int px;
        G4int py;
        G4int pz;

        // Vector buffers
        std::vector<G4int> secPDG;
        std::vector<G4int> secBaryonNumber;

        std::vector<G4double> secCharge;

        std::vector<G4double> secKinetic;
        std::vector<G4double> secMass;
        std::vector<G4double> secTotal;

        std::vector<G4double> secPx;
        std::vector<G4double> secPy;
        std::vector<G4double> secPz;
    };


    // ======================================================
    // Ntuple column IDs
    // ======================================================

    G4int ntupleID;

    // Remnant before precompound
    G4int AremBeforePreco;
    G4int ZremBeforePreco;
    G4int EexcBeforePreco;

    G4int NexcBeforePreco;
    G4int NparBeforePreco;
    G4int NholBeforePreco;

    // Remnant before deexcitation
    G4int AremBeforeDeex;
    G4int ZremBeforeDeex;
    G4int EexcBeforeDeex;

    G4int NexcBeforeDeex;
    G4int NparBeforeDeex;
    G4int NholBeforeDeex;


    // ======================================================
    // Event data
    // ======================================================

    RemnantInfo remnantBeforePreco;
    RemnantInfo remnantBeforeDeex;

    SecondaryData ftf;
    SecondaryData ftfp;


    // ======================================================
    // Internal methods
    // ======================================================

    void CreateSecondaryColumns(
        G4AnalysisManager* man,
        const std::string& prefix,
        SecondaryData& data)
    {
        // Scalar columns
        data.nSecondaries = man->CreateNtupleIColumn(
            ntupleID, prefix + "nSecondaries");

        data.sumKinetic = man->CreateNtupleDColumn(
            ntupleID, prefix + "sumKinetic_MeV");

        data.sumMass = man->CreateNtupleDColumn(
            ntupleID, prefix + "sumMass_MeV");

        data.sumTotal = man->CreateNtupleDColumn(
            ntupleID, prefix + "sumTotal_MeV");


        // Vector columns
        data.PDG = man->CreateNtupleIColumn(
            ntupleID, prefix + "secPDG", data.secPDG);

        data.charge = man->CreateNtupleDColumn(
            ntupleID, prefix + "secCharge", data.secCharge);

        data.baryonNumber = man->CreateNtupleIColumn(
            ntupleID, prefix + "secBaryonNumber",
            data.secBaryonNumber);

        data.kinetic = man->CreateNtupleDColumn(
            ntupleID, prefix + "secKinetic_i_MeV",
            data.secKinetic);

        data.mass = man->CreateNtupleDColumn(
            ntupleID, prefix + "secMass_i_MeV",
            data.secMass);

        data.total = man->CreateNtupleDColumn(
            ntupleID, prefix + "secTotal_i_MeV",
            data.secTotal);

        data.px = man->CreateNtupleDColumn(
            ntupleID, prefix + "secPx_i_MeV",
            data.secPx);

        data.py = man->CreateNtupleDColumn(
            ntupleID, prefix + "secPy_i_MeV",
            data.secPy);

        data.pz = man->CreateNtupleDColumn(
            ntupleID, prefix + "secPz_i_MeV",
            data.secPz);
    }


    void AddSecondary(
        SecondaryData& data,
        const G4DynamicParticle* particle)
    {
        if (!particle) return;

        const auto* def = particle->GetDefinition();
        const auto& p = particle->GetMomentum();

        const G4double ekin = particle->GetKineticEnergy();
        const G4double mass = particle->GetMass();

        data.secPDG.push_back(def->GetPDGEncoding());
        data.secCharge.push_back(def->GetPDGCharge());
        data.secBaryonNumber.push_back(def->GetBaryonNumber());

        data.secKinetic.push_back(ekin / MeV);
        data.secMass.push_back(mass / MeV);
        data.secTotal.push_back(particle->GetTotalEnergy() / MeV);

        data.secPx.push_back(p.x() / MeV);
        data.secPy.push_back(p.y() / MeV);
        data.secPz.push_back(p.z() / MeV);
    }


    void FillSecondarySummary(
        G4AnalysisManager* man,
        const SecondaryData& data)
    {
        const auto sum = [](const std::vector<G4double>& v) {
            return std::accumulate(
                v.begin(), v.end(), G4double{0.0});
        };

        man->FillNtupleIColumn(
            ntupleID, data.nSecondaries, data.secPDG.size());

        man->FillNtupleDColumn(
            ntupleID, data.sumKinetic, sum(data.secKinetic));

        man->FillNtupleDColumn(
            ntupleID, data.sumMass, sum(data.secMass));

        man->FillNtupleDColumn(
            ntupleID, data.sumTotal, sum(data.secTotal));
    }


    G4bool CheckConsistency(const SecondaryData& data) const
    {
        const auto n = data.secPDG.size();

        return data.secCharge.size() == n &&
               data.secBaryonNumber.size() == n &&
               data.secKinetic.size() == n &&
               data.secMass.size() == n &&
               data.secTotal.size() == n &&
               data.secPx.size() == n &&
               data.secPy.size() == n &&
               data.secPz.size() == n;
    }


    void ClearSecondaryData(SecondaryData& data)
    {
        data.secPDG.clear();
        data.secCharge.clear();
        data.secBaryonNumber.clear();

        data.secKinetic.clear();
        data.secMass.clear();
        data.secTotal.clear();

        data.secPx.clear();
        data.secPy.clear();
        data.secPz.clear();
    }

};

#endif
