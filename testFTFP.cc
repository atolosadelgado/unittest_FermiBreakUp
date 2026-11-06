#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4HadronicParameters.hh"

#include <CLHEP/Random/RanecuEngine.h>

#include <iostream>
#include <memory>

//----
#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
G4TheoFSGenerator* BuildFTFP();

// Headers to initialize particles
#include "G4DecayPhysics.hh"
#include "G4ParticleTable.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4IonTable.hh"

// custom precompound
#include "G4DebugPreCompoundModel.hh"

#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"


void InitializeParticles();

int main()
{
    // Random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4Random::setTheSeed(123456);

    InitializeParticles();

    // Hadronic parameters
    G4HadronicParameters::Instance()->SetVerboseLevel(0);


    G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
    G4cout << "G4 precompound type: " << param->GetPreCompoundType() << G4endl;
    param->SetPreCompoundType(G4PreCompoundType::eDefault);

    G4cout << "My precompound type: " << param->GetPreCompoundType() << G4endl;

    // Initialize projectile: proton, 100 GeV
    G4DynamicParticle dynamicParticle(
        G4Proton::Proton(),
        G4ThreeVector(0., 0., 1.),
        100. * GeV
    );

    G4HadProjectile projectile(dynamicParticle);

    // Initialize target: carbon-12
    G4Nucleus target;
    target.SetParameters(208, 82);

    // Initialize FTFP
    std::unique_ptr<G4TheoFSGenerator> ftfp(BuildFTFP());

    for(int i =0; i<20000; ++i)
    {    // Generate interaction
        G4HadFinalState* finalState =
        ftfp->ApplyYourself(projectile, target);

        // G4cout << "End iteration = " << i << G4endl;
        // G4cerr << std::flush;

    }
    // auto nucleons = ftfp->GetHighEnergyGenerator()->GetWoundedNucleus()->GetNucleons();
    // std::cout << "number of nucleons " << nucleons.size() << std::endl;
    // auto * ftfmodel = static_cast<const G4FTFModel*>(ftfp->GetHighEnergyGenerator());
    // std::cout << "GetNumberOfTargetSpectatorNucleons " << ftfmodel->GetNumberOfTargetSpectatorNucleons() << std::endl;
    // std::cout << "GetNumberOfNNcollisions " << ftfmodel->GetNumberOfNNcollisions() << std::endl;

    // auto * preco = static_cast<const G4GeneratorPrecompoundInterface*>(ftfp->GetTransport());
    // preco->Ge

    // Print secondaries
    // std::cout << "\nFinal state:\n";
    //
    // for (G4int i = 0;
    //      i < finalState->GetNumberOfSecondaries(); ++i)
    // {
    //     auto* particle = finalState->GetSecondary(i)->GetParticle();
    //
    //     std::cout
    //         << particle->GetDefinition()->GetParticleName()
    //         << "  "
    //         << particle->GetKineticEnergy() / GeV
    //         << " GeV\n";
    // }
    //
    // std::cout << "\nNumber of secondaries: "
    //           << finalState->GetNumberOfSecondaries()
    //           << std::endl;

    TFile* file = new TFile("output.root", "RECREATE");
    hEexcFTFremnant->Write();
    hEexcEquilibrium->Write();
    hAEexcFTFremnant->Write();
    hAEexcEquilibrium->Write();
    hAEexcEquilibriumOutsideLimit->Write();
    file->Close();

    delete file;
    return 0;
}


void InitializeParticles()
{


// Define generic ion
auto* gion = G4GenericIon::Definition();
gion->SetProcessManager(new G4ProcessManager(gion));

// Define standard particles
auto* decays = new G4DecayPhysics;
decays->ConstructParticle();

// Mark particle table as ready
auto* particleTable = G4ParticleTable::GetParticleTable();
particleTable->SetReadiness();

// Initialize ions
auto* ions = particleTable->GetIonTable();
ions->CreateAllIon();
ions->CreateAllIsomer();

}

G4TheoFSGenerator* BuildFTFP()
{
    auto* excitationHandler = new G4ExcitationHandler;

    auto* precompound =
        // new G4PreCompoundModel(excitationHandler);
        new G4DebugPreCompoundModel(excitationHandler);

    auto* transport =
        new G4GeneratorPrecompoundInterface;

    transport->SetDeExcitation(precompound);

    auto* fragmentation =
        new G4LundStringFragmentation;

    auto* stringDecay =
        new G4ExcitedStringDecay(fragmentation);

    auto* ftf = new G4FTFModel;

    ftf->SetFragmentationModel(stringDecay);

    auto* ftfp = new G4TheoFSGenerator("FTFP");

    ftfp->SetHighEnergyGenerator(ftf);
    ftfp->SetTransport(transport);

    return ftfp;
}
