
#include "G4ExcitationHandler.hh"
#include "G4VFermiBreakUp.hh"
#include "G4Fragment.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NucleiProperties.hh"
#include "G4LorentzVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include "InputParser.hh"

int main(int argc, char** argv)
{
    InputParser input(argc, argv);

    const int massNumberA = input.massNumberA;
    const int atomicNumberZ = input.atomicNumberZ;
    const double excitationEnergyMeV = input.excitationEnergyMeV;
    const std::string fermiBreakUpModelName = input.fermiBreakUpModel;
    const int nEvents = input.nEvents;
    const int fermiBreakUpModelVerbosity = 0;
    const int testVerbosity=0;

    std::cout << "A = " << massNumberA << '\n'
                << "Z = " << atomicNumberZ << '\n'
                << "Excitation energy = " << excitationEnergyMeV << " MeV\n"
                << "Fermi BreakUp model = " << fermiBreakUpModelName << '\n'
                << "Number of events = " << nEvents << '\n';

    using Isotope = std::pair<G4int, G4int>;
    using Multiplicity = G4int;

    std::map<Isotope, Multiplicity> fragmentCounts;
		
		
    // ----------------------------------------------------
    // Configure de-excitation, choose Fermi breakup model
    // ----------------------------------------------------

    G4DeexPrecoParameters* parameters =
        G4NuclearLevelData::GetInstance()->GetParameters();

    // Default implementation: G4FermiBreakUpVI
    if( "VI" == fermiBreakUpModelName )
    {
        parameters->SetFermiBreakUpType(bModelVI);
    }
    else if( "AN" == fermiBreakUpModelName )
    {
        parameters->SetFermiBreakUpType(bModelAN);
    }
    else
    {
        std::cerr << "Fermi Breakup model <" << fermiBreakUpModelName << "> not known" << '\n';
        return 1;
    }

    G4cout << "GetFBUEnergyLimit = " << parameters->GetFBUEnergyLimit()/MeV << G4endl;
    G4cout << "GetFermiEnergy = " << parameters->GetFermiEnergy()/MeV << G4endl;

    // --------------------------------------------------
    // Initialise the excitation handler
    // --------------------------------------------------

    G4ExcitationHandler handler;
    handler.Initialise();

    // Get the Fermi Breakup model
    auto* fermi = handler.GetFermiModel();

    if (!fermi) {
        G4cerr << "ERROR: Fermi Breakup model not available"
               << G4endl;
        return 1;
    }
    fermi->SetVerbose(fermiBreakUpModelVerbosity);

    // --------------------------------------------------
    // Test applicability of Fermi Breakup
    // --------------------------------------------------
    const G4double groundMass =
        G4NucleiProperties::GetNuclearMass(massNumberA, atomicNumberZ);
    {
        // Fragment at rest, with excitation energy
        G4LorentzVector momentum(
            0.0, 0.0, 0.0, groundMass + excitationEnergyMeV*MeV
        );

        G4Fragment* initial =
            new G4Fragment(massNumberA, atomicNumberZ, momentum);

        if (!fermi->IsApplicable(
                atomicNumberZ, massNumberA, initial->GetExcitationEnergy())) {

            G4cerr << "[localtest] Fermi Breakup is not applicable!"
                << G4endl;

            delete initial;
            return 1;
        }
    }
    
    // --------------------------------------------------
    // Wrapper on the call of the model
    // --------------------------------------------------

    auto unittest_break = [&fermi, &fragmentCounts, massNumberA, atomicNumberZ, excitationEnergyMeV, groundMass]() {

        G4LorentzVector momentum(
            0.0, 0.0, 0.0, groundMass + excitationEnergyMeV
        );

        // New fragment for every call
        auto* initial = new G4Fragment(massNumberA, atomicNumberZ, momentum);

        G4FragmentVector products;

        // initial is consumed during the call
        fermi->BreakFragment(&products, initial);
        if(testVerbosity>0){
            G4cout << "\nFinal state:" << G4endl;
            G4cout << "Number of fragments: "
                    << products.size() << G4endl;
        }

        for (const auto* fragment : products) {

                if(testVerbosity>0){
                    G4cout << "  A = " << fragment->GetA_asInt()
                                        << ", Z = " << fragment->GetZ_asInt()
                                        << ", E* = "
                                        << fragment->GetExcitationEnergy() / MeV
                                        << " MeV"
                                        << ", E = "
                                        << fragment->GetMomentum().e() / MeV
                                        << " MeV"
                                        << G4endl;
                }
                Isotope isotope(
                        fragment->GetA_asInt(),
                        fragment->GetZ_asInt()
                );

                ++fragmentCounts[isotope];
            delete fragment;
        }
    }; // end call of Fermi Breakup model

    for (int i = 0; i < nEvents; ++i) {
            unittest_break();
    }

    
    G4cout << "\nFragment multiplicities:" << G4endl;

		for (const auto& [isotope, multiplicity] : fragmentCounts) {

				G4cout << "  A = " << isotope.first
				       << ", Z = " << isotope.second
				       << ", multiplicity = " << multiplicity
				       << G4endl;
		}

    return 0;
}

