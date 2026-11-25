
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4IonisParamMat.hh"
#include "G4Version.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <iomanip>
#include <vector>

int main()
{
    auto* nist = G4NistManager::Instance();
    nist->SetVerbose(0);

    const auto& names = nist->GetNistMaterialNames();

    std::cout << std::setprecision(17);

    std::cout << "# Geant4 version: "
              << G4VERSION_NUMBER << "\n";

    std::cout << "# Number of NIST materials: "
              << names.size() << "\n\n";

    for (const auto& name : names)
    {
        auto* mat = nist->FindOrBuildMaterial(name);

        if (!mat)
        {
            std::cerr << "WARNING: Cannot build "
                      << name << "\n";
            continue;
        }

        std::cout << "================================================\n";
        std::cout << "MATERIAL " << name << " \n";

        // General properties
        std::cout << "Density       " << mat->GetDensity() / (CLHEP::g/CLHEP::cm3)<< "\n";
        std::cout << "State         " << mat->GetState() << "\n";
        std::cout << "Temperature   " << mat->GetTemperature() / (CLHEP::kelvin) << "\n";
        std::cout << "Pressure      " << mat->GetPressure() / (CLHEP::pascal)<< "\n";

        // Radiation and interaction lengths
        std::cout << "RadLength     "
          << mat->GetRadlen() / CLHEP::mm << " mm\n";

        std::cout << "IntLength     "
          << mat->GetNuclearInterLength() / CLHEP::mm << " mm\n";
        // Atomic properties
        std::cout << "ElectronDensity "
          << mat->GetElectronDensity() * CLHEP::cm3
          << " electrons/cm3\n";

        std::cout << "AtomsPerVolume "
          << mat->GetTotNbOfAtomsPerVolume() * CLHEP::cm3
          << " atoms/cm3\n";
	
	// Mean excitation energy
        auto* ion = mat->GetIonisation();

        if (ion)
        {
          std::cout << "MeanExcitationEnergy "
            << ion->GetMeanExcitationEnergy() / CLHEP::eV
            << " eV\n";
	  // Density effect parameters
	  std::cout << "PlasmaEnergy "
            << ion->GetPlasmaEnergy() / CLHEP::eV
            << " eV\n";

	  std::cout << "AdjustmentFactor "
            << ion->GetAdjustmentFactor() << "\n";

	  std::cout << "Cdensity " << ion->GetCdensity() << "\n";
	  std::cout << "Mdensity " << ion->GetMdensity() << "\n";
	  std::cout << "Adensity " << ion->GetAdensity() << "\n";
	  std::cout << "X0density " << ion->GetX0density() << "\n";
	  std::cout << "X1density " << ion->GetX1density() << "\n";
	  std::cout << "D0density " << ion->GetD0density() << "\n";
	  for (G4double x = -1.0; x <= 4.0; x += 0.1)
	  {
	    std::cout << "DensityCorrection "
        	      << x << " "
	              << ion->DensityCorrection(x) << "\n";
	  }
	  std::cout << "LogMeanExcEnergy "
          << ion->GetLogMeanExcEnergy() << "\n";

	  std::cout << "Taul "
            << ion->GetTaul() << "\n";

	  const auto* shell = ion->GetShellCorrectionVector();

	  for (G4int i = 0; i < 3; ++i)
	  {
	    std::cout << "ShellCorrection[" << i << "] "
        	      << shell[i] << "\n";
	  }
	  std::cout << "F1fluct " << ion->GetF1fluct() << "\n";
	  std::cout << "F2fluct " << ion->GetF2fluct() << "\n";

	  std::cout << "Energy1fluct "
          	    << ion->GetEnergy1fluct() / CLHEP::eV << " eV\n";

	  std::cout << "Energy2fluct "
          	    << ion->GetEnergy2fluct() / CLHEP::eV << " eV\n";

	  std::cout << "Energy0fluct "
          	    << ion->GetEnergy0fluct() / CLHEP::eV << " eV\n";

	  std::cout << "Rateionexcfluct "
                    << ion->GetRateionexcfluct() << "\n";
	  std::cout << "Zeffective " << ion->GetZeffective() << "\n";
	  std::cout << "FermiEnergy "
          	    << ion->GetFermiEnergy() / CLHEP::eV << " eV\n";
	  std::cout << "LFactor " << ion->GetLFactor() << "\n";
	  std::cout << "InvA23 " << ion->GetInvA23() << "\n";
	  std::cout << "BirksConstant "
	            << ion->GetBirksConstant() / (CLHEP::mm / CLHEP::MeV)
        	    << " mm/MeV\n";

	  std::cout << "MeanEnergyPerIonPair "
          	    << ion->GetMeanEnergyPerIonPair() / CLHEP::eV
	            << " eV\n";	  
	}
        // Composition
        const G4int nElements = mat->GetNumberOfElements();

        std::cout << "NumberOfElements " << nElements << "\n";

        const auto* elements = mat->GetElementVector();
        const auto* fractions = mat->GetFractionVector();
        const auto* atoms = mat->GetAtomsVector();

        for (G4int i = 0; i < nElements; ++i)
        {
            const auto* el = (*elements)[i];

            std::cout << "ELEMENT "
                      << el->GetName()
                      << " Z=" << el->GetZ()
                      << " A=" << el->GetA() / (CLHEP::g / CLHEP::mole)
                      << " MassFraction=" << fractions[i];

            if (atoms)
                std::cout << " Atoms=" << atoms[i];

            std::cout << "\n";
        }

        std::cout << "\n";
    }

    return 0;
}

