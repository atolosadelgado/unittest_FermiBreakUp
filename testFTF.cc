#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4HadronicParameters.hh"

#include <CLHEP/Random/RanecuEngine.h>

#include <iostream>
#include <memory>
#include <string>
#include <cstdlib>

//----
#include "G4TheoFSGenerator.hh"
#include "G4FTFModel.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4LundStringFragmentation.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
// custom precompound
#include "G4DebugPreCompoundModel.hh"

void BuildFTF(std::unique_ptr<G4TheoFSGenerator>       & theogen,
              std::unique_ptr<G4DebugPreCompoundModel> & precompound,
              std::string theogenName);

// Headers to initialize particles
#include "G4DecayPhysics.hh"
#include "G4ParticleTable.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"
#include "G4IonTable.hh"

void InitializeParticles();

#include "G4AnalysisManager.hh"

struct InputParameters
{
    int targetMassNumberA = 208;
    int targetAtomicNumberZ = 82;
    double projectileEnergyGeV = 2.0;
    std::string projectileName = "proton";
    std::string theoGeneratorName = "FTFB";
    int nEvents = 1000;

    InputParameters(int argc, char* argv[])
    {
        if (argc > 1) targetMassNumberA = std::atoi(argv[1]);
        if (argc > 2) targetAtomicNumberZ = std::atoi(argv[2]);
        if (argc > 3) projectileEnergyGeV = std::atof(argv[3]);
        if (argc > 4) projectileName = argv[4];
        if (argc > 5) theoGeneratorName = argv[5];
        if (argc > 6) nEvents = std::atoi(argv[6]);
    }
    void Show(){
        std::cout << "---------------------------------" << std::endl;
        std::cout << "Input parameters:\n";
        std::cout << "\tTarget A = "   << targetMassNumberA << '\n';
        std::cout << "\tTarget Z = "   << targetAtomicNumberZ << '\n';
        std::cout << "\tEnergy = "     << projectileEnergyGeV << " GeV\n";
        std::cout << "\tProjectile = " << projectileName << '\n';
        std::cout << "\tModel = "      << theoGeneratorName << '\n';
        std::cout << "\tEvents = "     << nEvents << '\n';
        std::cout << "---------------------------------" << std::endl;
    }
};


int main(int argc, char* argv[])
{

    InputParameters input(argc, argv);
    input.Show();

    int targetMassNumberA = input.targetMassNumberA;
    int targetAtomicNumberZ = input.targetAtomicNumberZ;
    double projectileEnergyGeV = input.projectileEnergyGeV;
    std::string projectileName = input.projectileName;
    std::string theoGeneratorName = input.theoGeneratorName;
    int nEvents = input.nEvents;

    std::string output_filename("test_");
    output_filename += theoGeneratorName;
    output_filename += "_";
    output_filename += (std::to_string(targetMassNumberA) + "A");
    output_filename += "_";
    output_filename += (std::to_string(targetAtomicNumberZ) + "Z");
    output_filename += "_";
    output_filename += (projectileName + std::to_string(int(projectileEnergyGeV)) + "GeV");
    output_filename += "_";
    output_filename += (std::to_string(nEvents) + "evts");
    output_filename += ".root";

    auto ana = G4AnalysisManager::Instance();
    ana->OpenFile(output_filename);

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


   G4ParticleDefinition * projectileDef = G4ParticleTable::GetParticleTable()->FindParticle( projectileName );

    G4DynamicParticle dynamicParticle(
        projectileDef,
        G4ThreeVector(0., 0., 1.),
        projectileEnergyGeV * GeV
    );

    G4HadProjectile projectile(dynamicParticle);

    // Initialize target: carbon-12
    G4Nucleus target;
    target.SetParameters(targetMassNumberA, targetAtomicNumberZ);

    // Initialize FTFP
    std::unique_ptr<G4TheoFSGenerator> ftfp;
    std::unique_ptr<G4DebugPreCompoundModel> precompound;
    BuildFTF(ftfp,precompound, theoGeneratorName);

    for(int i =0; i<nEvents; ++i)
    {
        // Generate interaction
        G4HadFinalState* finalState =
        ftfp->ApplyYourself(projectile, target);

        if( i%nEvents4printout ==0 )
            G4cout << "Progress... " << 100*i/nEvents << "%" << G4endl;
    }

    ana->Write();
    ana->CloseFile();
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


void BuildFTF(std::unique_ptr<G4TheoFSGenerator>& theogen,
               std::unique_ptr<G4DebugPreCompoundModel>& precompound,
               std::string theogenName)
{
    auto* excitationHandler = new G4ExcitationHandler;

    precompound = std::make_unique<G4DebugPreCompoundModel>(excitationHandler);

    G4VIntraNuclearTransportModel *transport = nullptr;
    if("FTFP" == theogenName)
        transport = new G4GeneratorPrecompoundInterface;
    else if ("FTFB" == theogenName)
        transport = new G4BinaryCascade;
    else{
        throw std::runtime_error("Theo generator name <" + theogenName + ">not known");
    }

    transport->SetDeExcitation(precompound.get());

    auto* fragmentation = new G4LundStringFragmentation;
    auto* stringDecay = new G4ExcitedStringDecay(fragmentation);

    auto* ftf = new G4FTFModel;
    ftf->SetFragmentationModel(stringDecay);

    theogen = std::make_unique<G4TheoFSGenerator>(theogenName);

    theogen->SetHighEnergyGenerator(ftf);
    theogen->SetTransport(transport);
}

