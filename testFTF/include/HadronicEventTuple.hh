
#include <vector>
#include <memory>
#include <algorithm>
#include <utility>

#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4AnalysisManager.hh"
#include "G4DynamicParticle.hh"

class HadronicEventTuple
{
public:

    struct IDs {
        G4int ntuple;

        // Event-level
        G4int AremPreco;
        G4int ZremPreco;
        G4int EexcPreco;
        G4int AremDeex;
        G4int ZremDeex;
        G4int EexcDeex;

        G4int NexcPreco;
        G4int NparPreco;
        G4int NholPreco;
        G4int NexcDeex;
        G4int NparDeex;
        G4int NholDeex;

        G4int nSecondaries;
        G4int secKinetic;
        G4int secMass;
        G4int secTotal;

        // Secondary-level vectors
        G4int secPDG;
        G4int secCharge;
        G4int secBaryonNumber;

        G4int secKineticVector;
        G4int secMassVector;
        G4int secTotalVector;

        G4int secPx;
        G4int secPy;
        G4int secPz;
    };

    HadronicEventTuple()
    {
        auto* man = G4AnalysisManager::Instance();

        ids.ntuple = man->CreateNtuple("events", "Hadronic events");

        ids.AremPreco = man->CreateNtupleIColumn(ids.ntuple, "AremPreco");
        ids.ZremPreco = man->CreateNtupleIColumn(ids.ntuple, "ZremPreco");
        ids.EexcPreco = man->CreateNtupleDColumn(ids.ntuple, "EexcPreco_MeV");

        ids.AremDeex = man->CreateNtupleIColumn(ids.ntuple, "AremDeex");
        ids.ZremDeex = man->CreateNtupleIColumn(ids.ntuple, "ZremDeex");
        ids.EexcDeex = man->CreateNtupleDColumn(ids.ntuple, "EexcDeex_MeV");

        ids.NexcPreco = man->CreateNtupleIColumn(ids.ntuple, "NexcPreco");
        ids.NparPreco = man->CreateNtupleIColumn(ids.ntuple, "NparPreco");
        ids.NholPreco = man->CreateNtupleIColumn(ids.ntuple, "NholPreco");

        ids.NexcDeex = man->CreateNtupleIColumn(ids.ntuple, "NexcDeex");
        ids.NparDeex = man->CreateNtupleIColumn(ids.ntuple, "NparDeex");
        ids.NholDeex = man->CreateNtupleIColumn(ids.ntuple, "NholDeex");

        ids.nSecondaries = man->CreateNtupleIColumn(
            ids.ntuple, "nSecondaries");

        ids.secKinetic = man->CreateNtupleDColumn(
            ids.ntuple, "secKinetic_MeV");

        ids.secMass = man->CreateNtupleDColumn(
            ids.ntuple, "secMass_MeV");

        ids.secTotal = man->CreateNtupleDColumn(
            ids.ntuple, "secTotal_MeV");

        // Vector columns
        ids.secPDG = man->CreateNtupleIColumn(
            ids.ntuple, "secPDG", secPDG);

        ids.secCharge = man->CreateNtupleDColumn(
            ids.ntuple, "secCharge", secCharge);

        ids.secBaryonNumber = man->CreateNtupleIColumn(
            ids.ntuple, "secBaryonNumber", secBaryonNumber);

        ids.secKineticVector = man->CreateNtupleDColumn(
            ids.ntuple, "secKinetic_i_MeV", secKinetic_i);

        ids.secMassVector = man->CreateNtupleDColumn(
            ids.ntuple, "secMass_i_MeV", secMass_i);

        ids.secTotalVector = man->CreateNtupleDColumn(
            ids.ntuple, "secTotal_i_MeV", secTotal_i);

        ids.secPx = man->CreateNtupleDColumn(
            ids.ntuple, "secPx_MeV", secPx);

        ids.secPy = man->CreateNtupleDColumn(
            ids.ntuple, "secPy_MeV", secPy);

        ids.secPz = man->CreateNtupleDColumn(
            ids.ntuple, "secPz_MeV", secPz);

        man->FinishNtuple(ids.ntuple);
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    IDs ids;

    // Vector buffers
    std::vector<G4int> secPDG;
    std::vector<G4int> secBaryonNumber;

    std::vector<G4double> secCharge;
    std::vector<G4double> secKinetic_i;
    std::vector<G4double> secMass_i;
    std::vector<G4double> secTotal_i;

    std::vector<G4double> secPx;
    std::vector<G4double> secPy;
    std::vector<G4double> secPz;

    void AddSecondary(const G4DynamicParticle* particle)
    {
        const auto* def = particle->GetDefinition();
        const auto& p = particle->GetMomentum();

        secPDG.push_back(def->GetPDGEncoding());
        secCharge.push_back(def->GetPDGCharge());
        secBaryonNumber.push_back(def->GetBaryonNumber());

        secKinetic_i.push_back(particle->GetKineticEnergy() / MeV);
        secMass_i.push_back(particle->GetMass() / MeV);
        secTotal_i.push_back(particle->GetTotalEnergy() / MeV);

        secPx.push_back(p.x() / MeV);
        secPy.push_back(p.y() / MeV);
        secPz.push_back(p.z() / MeV);

    }

    void FillEvent()
    {
        auto* man = G4AnalysisManager::Instance();

        // Assuming scalar columns are filled already...
        auto sum_ff = [](std::vector<G4double> & vector){
            G4double sum=0;
            for( auto number : vector) sum+=number;
            return sum;
        };

        G4double secKinetic_MeV = sum_ff(secKinetic_i);
        G4double secMass_MeV = sum_ff(secMass_i);
        G4double secTotal_MeV = sum_ff(secTotal_i);

        man->FillNtupleDColumn( ids.secKinetic, secKinetic_MeV);
        man->FillNtupleDColumn( ids.secMass, secMass_MeV);
        man->FillNtupleDColumn( ids.secTotal, secTotal_MeV);

        man->FillNtupleIColumn(ids.ntuple, ids.nSecondaries,
                            secPDG.size());

        man->AddNtupleRow(ids.ntuple);

        Clear();
    }

    void Clear()
    {
        secPDG.clear();
        secCharge.clear();
        secBaryonNumber.clear();

        secKinetic_i.clear();
        secMass_i.clear();
        secTotal_i.clear();

        secPx.clear();
        secPy.clear();
        secPz.clear();
    }

};
