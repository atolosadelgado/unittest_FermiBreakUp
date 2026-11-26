#include "G4HadronInelasticFTFB_BIC.hh"

#include "G4BGGNucleonInelasticXS.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4BinaryCascade.hh"
#include "G4BuilderType.hh"
#include "G4CascadeInterface.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4HadParticles.hh"
#include "G4HadProcesses.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronicBuilder.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronGeneralProcess.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronRadCapture.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleInelasticXS.hh"
#include "G4PreCompoundModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4TheoFSGenerator.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronInelasticFTFB_BIC);

G4HadronInelasticFTFB_BIC::G4HadronInelasticFTFB_BIC(G4int ver) : G4VHadronPhysics("hInelasticFTFB_BIC")
{
  SetPhysicsType(bHadronInelastic);
  auto param = G4HadronicParameters::Instance();
  // param->SetEnableBCParticles(true);
  param->SetEnableNeutronGeneralProcess(false);
  // param->SetUseRFilesForXS(true);
  param->SetVerboseLevel(ver);
}

void G4HadronInelasticFTFB_BIC::ConstructProcess()
{
  G4HadronicParameters* param = G4HadronicParameters::Instance();
  G4bool useFactorXS = param->ApplyFactorXS();
  G4PhysicsListHelper* physicsHelper = G4PhysicsListHelper::GetPhysicsListHelper();

  // configure models
  const G4double eminFtf = param->GetMinEnergyTransitionFTF_Cascade();
  const G4double emaxBic = param->GetMaxEnergyTransitionFTF_Cascade();
  const G4double emax = param->GetMaxEnergy();

  if (G4Threading::IsMasterThread() && param->GetVerboseLevel() > 0)
  {
    G4cout << __PRETTY_FUNCTION__ << std::endl;
    G4cout << "\tEmin(FTFB)= " << eminFtf / CLHEP::GeV << " GeV" << std::endl;
    G4cout << "\tEmax(FTFB)= " << emax    / CLHEP::GeV << " GeV" << std::endl;
    G4cout << "\rEmax(BIC) = " << emaxBic / CLHEP::GeV << " GeV" << std::endl;
  }

  auto ff_builder = [&](G4ParticleDefinition * particle)->void {

    // Create hadronic process for hadronic inelastic interactions
    auto procname = particle->GetParticleName() + "Inelastic";
    G4HadronicProcess* hp = new G4HadronInelasticProcess(procname, particle);

    // Register process in the physics list G4PhysicsListHelper
    physicsHelper->RegisterProcess(hp, particle);

    // Start filling the process object:

    // 1. Attach BIC and FTFB models
    // 1.1 FTFB model
    G4FTFModel * theStringModel = new G4FTFModel();
    theStringModel->SetFragmentationModel(new G4ExcitedStringDecay());
    G4TheoFSGenerator * theFTFB = new G4TheoFSGenerator("FTFB");
    theFTFB->SetHighEnergyGenerator(theStringModel);
    auto * transport_for_FTF = new G4BinaryCascade(/*thePreCompound*/);
    theFTFB->SetTransport(transport_for_FTF);
    theFTFB->SetMinEnergy(eminFtf);
    theFTFB->SetMaxEnergy(emax);

    // 1.2 BIC model
    auto theBIC = new G4BinaryCascade(/*thePreCompound*/);
    theBIC->SetMaxEnergy(emaxBic);

    // 1.3 Register models into the process
    hp->RegisterMe(theFTFB);
    hp->RegisterMe(theBIC);

    // 2. Attach dataset of hadronic inelastic to process
    G4VCrossSectionDataSet * cross_section_inel = nullptr;
    if(G4Neutron::Neutron() == particle ) cross_section_inel = new G4NeutronInelasticXS();
    else if(G4PionMinus::PionMinus() == particle) cross_section_inel = new G4BGGPionInelasticXS(particle);
    else if(G4PionPlus::PionPlus() == particle)   cross_section_inel = new G4BGGPionInelasticXS(particle);
    else                                  cross_section_inel = new G4ParticleInelasticXS(particle);

    hp->AddDataSet(cross_section_inel);

    // 3. Customization per particle
    if(G4Proton::Proton() == particle)
    {
        if (useFactorXS) hp->MultiplyCrossSectionBy(param->XSFactorNucleonInelastic());
    }

    if(G4Neutron::Neutron() == particle)
    {
        if (useFactorXS) hp->MultiplyCrossSectionBy(param->XSFactorNucleonInelastic());

        // Neutrons require capture as extra process
        // (ignoring neutron general process mechanism)
        auto nCap = new G4NeutronCaptureProcess("nCapture");
        nCap->RegisterMe(new G4NeutronRadCapture());
        physicsHelper->RegisterProcess(nCap, particle);
    }
    if(G4PionPlus::PionPlus() == particle)
    {
         hp->MultiplyCrossSectionBy(param->XSFactorPionInelastic());
    }
    if(G4PionMinus::PionMinus() == particle)
    {
         hp->MultiplyCrossSectionBy(param->XSFactorPionInelastic());
    }

};

    ff_builder(G4Proton::Proton());
    ff_builder(G4Neutron::Neutron());
    ff_builder(G4PionPlus::PionPlus());
    ff_builder(G4PionMinus::PionMinus());

    // keep FTFP_BERT builders for particles different from pions and nucleons
    // kaons
    G4HadronicBuilder::BuildKaonsFTFP_BERT();

    // high energy particles
    if (emax > param->EnergyThresholdForHeavyHadrons())
    {
        // pbar, nbar, anti light ions
        G4HadronicBuilder::BuildAntiLightIonsFTFP();

        // hyperons
        G4HadronicBuilder::BuildHyperonsFTFP_BERT();

        // b-, c- baryons and mesons
        if (param->EnableBCParticles())
        {
        G4HadronicBuilder::BuildBCHadronsFTFP_BERT();
        }
    }
}
