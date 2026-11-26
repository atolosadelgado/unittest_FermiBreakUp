#ifndef G4HADRONINELASTICFTFBBIC_HH
#define G4HADRONINELASTICFTFBBIC_HH

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4HadronInelasticFTFB_BIC : public G4VHadronPhysics
{
  public:

    // constructor
    G4HadronInelasticFTFB_BIC(G4int ver = 1);

    virtual ~G4HadronInelasticFTFB_BIC() = default;

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    void ConstructProcess() override;

    // copy constructor and hide assignment operator
    G4HadronInelasticFTFB_BIC(G4HadronInelasticFTFB_BIC&) = delete;
    G4HadronInelasticFTFB_BIC& operator=(const G4HadronInelasticFTFB_BIC& right) = delete;
};

#endif
