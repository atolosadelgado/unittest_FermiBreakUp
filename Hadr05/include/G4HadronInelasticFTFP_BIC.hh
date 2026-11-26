#ifndef G4HADRONINELASTICFTFPBIC_HH
#define G4HADRONINELASTICFTFPBIC_HH

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4HadronInelasticFTFP_BIC : public G4VHadronPhysics
{
  public:

    // constructor
    G4HadronInelasticFTFP_BIC(G4int ver = 1);

    virtual ~G4HadronInelasticFTFP_BIC() = default;

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructProcess() override;

    // copy constructor and hide assignment operator
    G4HadronInelasticFTFP_BIC(G4HadronInelasticFTFP_BIC&) = delete;
    G4HadronInelasticFTFP_BIC& operator=(const G4HadronInelasticFTFP_BIC& right) = delete;
};

#endif
