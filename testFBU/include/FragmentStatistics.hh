
#ifndef FRAGMENT_STATISTICS_HH
#define FRAGMENT_STATISTICS_HH

#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <ostream>
#include <utility>

class FragmentStatistics
{
public:

    struct Moments {
        double sum = 0.0;
        double sum2 = 0.0;
        unsigned long long n = 0;

        void Add(double x) {
            sum += x;
            sum2 += x * x;
            ++n;
        }

        double Mean() const {
            return n ? sum / n : 0.0;
        }

        double RMS() const {
            if (!n) return 0.0;
            const double mean = Mean();
            return std::sqrt(std::max(0.0, sum2 / n - mean * mean));
        }
    };

    struct IsotopeStatistics {
        unsigned long long totalProduction = 0;
        Moments multiplicityPerEvent;
    };

public:
    FragmentStatistics(int A, int Z)
        : initialA(A), initialZ(Z)
    {}
    // Global statistics
    unsigned long long nEvents = 0;
    unsigned long long totalFragments = 0;

    Moments multiplicity;
    Moments kineticEnergy;
    Moments excitationEnergy;

    // Isotope statistics: key = (A, Z)
    std::map<std::pair<int,int>, IsotopeStatistics> isotopes;

    // Analyze one event
    void Analyze(const G4FragmentVector& products)
    {
        ++nEvents;

        multiplicity.Add(products.size());
        totalFragments += products.size();

        double eventKineticEnergy = 0.0;
        double eventExcitationEnergy = 0.0;
        int eventA = 0;
        int eventZ = 0;

        // Temporary event isotope multiplicities
        std::map<std::pair<int,int>, unsigned long long> eventIsotopes;

        for (const auto* fragment : products) {

            const int A = fragment->GetA_asInt();
            const int Z = fragment->GetZ_asInt();
            eventA += A;
            eventZ += Z;
            const auto isotope = std::make_pair(A, Z);

            isotopes[isotope].totalProduction++;
            eventIsotopes[isotope]++;
            const auto & p4 = fragment->GetMomentum();
            const double FragmentKineticEnergy = p4.e() - p4.m();
            const double FragmentExcitatEnergy = fragment->GetExcitationEnergy();
            const double FragmentGroundMass    = fragment->GetGroundStateMass();
            const double FragmentInvariantMass = p4.m();
            constexpr double massTolerance = 0.05;
            const double massDifference = FragmentInvariantMass - FragmentExcitatEnergy - FragmentGroundMass;
            if( std::fabs(massDifference) > massTolerance*FragmentInvariantMass )
            {
                G4cout  << fragment
                        << "Ground-state mass: " << FragmentGroundMass / MeV    << " MeV"
                        << " | Invariant mass: " << FragmentInvariantMass / MeV << " MeV"
                        << " | Excitation: "     << FragmentExcitatEnergy / MeV << " MeV"
                        << " | Difference: "     << massDifference / keV        << " keV"
                        << G4endl;
            }

            eventKineticEnergy += FragmentKineticEnergy;
            eventExcitationEnergy += FragmentExcitatEnergy;
        }

        // Check conservation event by event
        const int deltaA = eventA - initialA;
        const int deltaZ = eventZ - initialZ;

        if (deltaA == 0) ++conservedA;
        if (deltaZ == 0) ++conservedZ;
        if (deltaA == 0 && deltaZ == 0) ++conservedAZ;

        maxDeltaA = std::max(maxDeltaA, std::abs(deltaA));
        maxDeltaZ = std::max(maxDeltaZ, std::abs(deltaZ));

        kineticEnergy.Add(eventKineticEnergy);
        excitationEnergy.Add(eventExcitationEnergy);

        // Record multiplicity for every known isotope,
        // including zero production in this event.
        for (auto& [isotope, stats] : isotopes) {
            const auto it = eventIsotopes.find(isotope);

            const double count = (it != eventIsotopes.end())
                               ? static_cast<double>(it->second)
                               : 0.0;

            stats.multiplicityPerEvent.Add(count);
        }
    }

    // Print accumulated statistics to any output stream
    void PrintOutput(std::ostream& out) const
    {
        out << std::fixed << std::setprecision(4);

        out << "========================================\n";
        out << "       FRAGMENTATION STATISTICS\n";
        out << "========================================\n\n";

        out << "Number of events: " << nEvents << '\n';
        out << "Total number of fragments: "
            << totalFragments << "\n\n";

        out << "--- Fragment multiplicity ---\n";
        out << "Mean: " << multiplicity.Mean() << '\n';
        out << "RMS:  " << multiplicity.RMS() << "\n\n";

        out << "--- Total kinetic energy [MeV] ---\n";
        out << "Mean: " << kineticEnergy.Mean() / MeV << '\n';
        out << "RMS:  " << kineticEnergy.RMS() / MeV << "\n\n";

        out << "--- Total excitation energy [MeV] ---\n";
        out << "Mean: " << excitationEnergy.Mean() / MeV << '\n';
        out << "RMS:  " << excitationEnergy.RMS() / MeV << "\n\n";

        out << "--- Isotope production ---\n";
        out << " A    Z    Total production"
            << "    Mean multiplicity/event"
            << "    RMS/event\n";

        for (const auto& [isotope, stats] : isotopes) {

            out << std::setw(3) << isotope.first << "  "
                << std::setw(3) << isotope.second << "  "
                << std::setw(15) << stats.totalProduction << "  "
                << std::setw(22)
                << stats.multiplicityPerEvent.Mean() << "  "
                << std::setw(10)
                << stats.multiplicityPerEvent.RMS()
                << '\n';
        }

        out << "\n========================================\n";

        out << "--- Conservation checks ---\n";

        out << "Mass number A: "
            << conservedA << " / " << nEvents
            << " events ("
            << (nEvents ? 100.0 * conservedA / nEvents : 0.0)
            << "%)\n";

        out << "Charge Z:      "
            << conservedZ << " / " << nEvents
            << " events ("
            << (nEvents ? 100.0 * conservedZ / nEvents : 0.0)
            << "%)\n";

        out << "Both A and Z:  "
            << conservedAZ << " / " << nEvents
            << " events ("
            << (nEvents ? 100.0 * conservedAZ / nEvents : 0.0)
            << "%)\n";

        out << "Maximum |delta A|: " << maxDeltaA << '\n';
        out << "Maximum |delta Z|: " << maxDeltaZ << "\n\n";

        out << "\n========================================\n";
    }
    private:
        const int initialA;
        const int initialZ;
        // Conservation statistics
        unsigned long long conservedA = 0;
        unsigned long long conservedZ = 0;
        unsigned long long conservedAZ = 0;

        int maxDeltaA = 0;
        int maxDeltaZ = 0;
};

#endif
