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

// custom precompound
#include "G4DebugPreCompoundModel.hh"

void BuildFTFP(std::unique_ptr<G4TheoFSGenerator>       & ftfp,
               std::unique_ptr<G4DebugPreCompoundModel> & precompound);

// Headers to initialize particles
#include "G4DecayPhysics.hh"
#include "G4ParticleTable.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4IonTable.hh"


#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"


void InitializeParticles();

int main()
{

    int nEvents = 1000;
    int nEvents4printout = nEvents/10;
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
    std::unique_ptr<G4TheoFSGenerator> ftfp;
    std::unique_ptr<G4DebugPreCompoundModel> precompound;
    BuildFTFP(ftfp,precompound);

    for(int i =0; i<nEvents; ++i)
    {    // Generate interaction
        G4HadFinalState* finalState =
        ftfp->ApplyYourself(projectile, target);

        if( i%nEvents4printout ==0 )
            G4cout << "Progress... " << 100*i/nEvents << "%" << G4endl;
        // G4cerr << std::flush;
    }

    TFile* file = new TFile("output2.root", "RECREATE");
    precompound->hEexcFTFremnant->Write();
    precompound->hEexcEquilibrium->Write();
    precompound->hAEexcFTFremnant->Write();
    precompound->hAEexcEquilibrium->Write();
    precompound->hAEexcEquilibriumOutsideLimit->Write();
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


void BuildFTFP(std::unique_ptr<G4TheoFSGenerator>& ftfp,
               std::unique_ptr<G4DebugPreCompoundModel>& precompound)
{
    auto* excitationHandler = new G4ExcitationHandler;

    precompound = std::make_unique<G4DebugPreCompoundModel>(excitationHandler);

    auto* transport = new G4GeneratorPrecompoundInterface;
    transport->SetDeExcitation(precompound.get());

    auto* fragmentation = new G4LundStringFragmentation;
    auto* stringDecay = new G4ExcitedStringDecay(fragmentation);

    auto* ftf = new G4FTFModel;
    ftf->SetFragmentationModel(stringDecay);

    ftfp = std::make_unique<G4TheoFSGenerator>("FTFP");

    ftfp->SetHighEnergyGenerator(ftf);
    ftfp->SetTransport(transport);
}
