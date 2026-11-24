
#include "G4DebugPreCompoundModel.hh"

#include "G4NuclearLevelData.hh"
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

    // Ntuple
    tuple->SetRemnantBeforePrecompound({
        A,
        Z,
        U,
        aFragment.GetNumberOfExcitons(),
        aFragment.GetNumberOfParticles(),
        aFragment.GetNumberOfHoles()
    });

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
        // MyDoIt(result, *ptr);
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


    tuple->SetRemnantBeforeDeexcitation({
        A,
        Z,
        U,
        aFragment.GetNumberOfExcitons(),
        aFragment.GetNumberOfParticles(),
        aFragment.GetNumberOfHoles()
    });

    G4PreCompoundModel::PerformEquilibriumEmission(
        aFragment, result);
}

void G4DebugPreCompoundModel::MyDoIt(G4ReactionProductVector* result, G4Fragment& aFragment)
{
  // main loop
  G4int count = 0;
  const G4double ldfact = 12.0 / CLHEP::pi2;
  const G4int countmax = 300;
  for (;;)
  {
    if(MyDebug)G4cout << "\n-------------\nNEW PRECOMPOUND LOOP" << G4endl;
    G4double U = aFragment.GetExcitationEnergy();
    G4int Z = aFragment.GetZ_asInt();
    G4int A = aFragment.GetA_asInt();
    G4int nExcitons = aFragment.GetNumberOfExcitons();

    // check conditions for pre-compound model for updated fragment
    if (nExcitons <= 0 || (Z < minZ && A < minA) || U > fHighLimitExc * A || U <= fLowLimitExc * A)
    {
      PerformEquilibriumEmission(aFragment, result);
      return;
    }
    G4int eqExcitonNumber = G4lrint(std::sqrt(ldfact * U * fNuclData->GetLevelDensity(Z, A, U)));
    if (2 < fVerbose)
    {
      G4cout << "   1st loop " << count << ". Z=" << Z << " A=" << A << " U(MeV)=" << U
             << " eqExcitationNumber=" << eqExcitonNumber << G4endl;
    }
    // J. M. Quesada (Jan. 08)  equilibrium hole number could be used as preeq.
    // evap. delimiter (IAEA report)

    // Loop for transitions, it is performed while there are
    // preequilibrium transitions.
    G4bool isTransition = true;
    do
    {
      ++count;

      // soft cutoff criterium as an "ad-hoc" solution to force
      // increase in  evaporation
      nExcitons = aFragment.GetNumberOfExcitons();
      if (nExcitons <= 0)
      {
        PerformEquilibriumEmission(aFragment, result);
        return;
      }

      // J. M. Quesada (Apr. 2008): soft-cutoff switched off by default
      G4bool go_ahead = (nExcitons < eqExcitonNumber);
      G4double x = (G4double)(nExcitons - eqExcitonNumber) / (G4double)eqExcitonNumber;
      if (useSCO && go_ahead)
      {

        if (G4UniformRand() < 1.0 - G4Exp(-x * x / 0.32))
        {
          go_ahead = false;
        }
      }
      if (!go_ahead)
      {
        if(MyDebug)G4cout << "STOP: Nex >= Neq" << G4endl;
        PerformEquilibriumEmission(aFragment, result);
        return;
      }
      // compute transition probability
      G4double transProbability = theTransition->CalculateProbability(aFragment);
      G4double P1 = theTransition->GetTransitionProb1();
      G4double P2 = theTransition->GetTransitionProb2();
      G4double P3 = theTransition->GetTransitionProb3();
      if(MyDebug)G4cout << "\n--- PRECOMPOUND TRANSITION ---" << G4endl;
      if(MyDebug)G4cout << "A = " << A
                        << " Z = " << Z
                        << " E* = " << U / CLHEP::MeV << " MeV"
                        << " Np = " << aFragment.GetNumberOfParticles()
                        << " Nh = " << aFragment.GetNumberOfHoles()
                        << " Nex = " << nExcitons
                        << " Neq = " << eqExcitonNumber
                        << " Nex-Neq = " << nExcitons-eqExcitonNumber << G4endl
                        << "P1 (n -> n+2) = " << P1 << G4endl
                        << "P2 (n -> n-2)   = " << P2 << G4endl
                        << "P3 (n -> n) = " << P3 << G4endl
                        << G4endl;

      if (2 < fVerbose)
      {
        G4cout << "   2nd loop " << count << ". Nex=" << nExcitons << " P1=" << P1 << " P2=" << P2
               << " P3=" << P3 << G4endl;
        G4cout << aFragment << G4endl;
      }
      // J.M. Quesada (May 2008) Physical criterium (lamdas)  PREVAILS over
      //                         approximation (critical exciton number)
      if (P1 <= P2 + P3)
      {
        if(MyDebug) G4cout << "STOP: P1 <= P2 + P3" << G4endl;
        PerformEquilibriumEmission(aFragment, result);
        return;
      }
      // compute emission probability
      G4double emissionProbability = theEmission->GetTotalProbability(aFragment);
      G4double totProbability = emissionProbability + transProbability;
      if (2 < fVerbose)
      {
        G4cout << "             prob_em=" << emissionProbability << " prob_tot=" << totProbability
               << G4endl;
      }

        // G4cout << "Ptotal        = " << transProbability << G4endl;
      // Select subprocess
      if (totProbability * G4UniformRand() > emissionProbability)
      {
        if(MyDebug) G4cout << "PRECOMPOUND INTERNAL TRANSITION" << G4endl;
        // It will be transition to state with a new number of excitons
        theTransition->PerformTransition(aFragment);
      }
      else
      {
        if(MyDebug) G4cout << "PRECOMPOUND PARTICLE EMISSION" << G4endl;
        // pre-compound fragment emission, primary fragment is changed
        isTransition = false;
        result->push_back(theEmission->PerformEmission(aFragment));
      }
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (isTransition);  // end of do loop with emission

    // stop if too many iterations
    if (count >= countmax)
    {
      G4ExceptionDescription ed;
      ed << "G4PreCompoundModel loop over " << countmax << " iterations; "
         << "current G4Fragment: \n"
         << aFragment;
      G4Exception("G4PreCompoundModel::DeExcite()", "had0034", JustWarning, ed, "");
      break;
    }
  }  // end of for (;;) loop
  PerformEquilibriumEmission(aFragment, result);
}
