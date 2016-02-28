
#include <iostream>
#include <fstream>
#include <vector>

#include "generaltools.cc"
#include "structures.cc"
#include "librarytools.cc"

using namespace std;

int main(int argc, char* argv[]) {

    cout << "--New simulations starts--" << endl;

    /// Read fuel info and libraries
    // Generate necessary libraries and path strings
    LibInfo lwrlib, moxlib;
    string lwrpath = "libs/extLWR", moxpath = "libs/PWRMOX";

    // Populate the libraries
    LibraryReader("lwrlib", lwrpath, lwrlib);
    LibraryReader("moxlib", moxpath, moxlib);

    // Build reactor object and fuel isos
    TypeInfo lwrtype;

    RegionInfo regionlwr;
    regionlwr.mass_ = 180;
    regionlwr.fractions[922350] = 0.0530;
    regionlwr.fractions[922380] = 0.9470;
    regionlwr.BuildIso(lwrlib);

    lwrtype.batch.assign(4, regionlwr);

    TypeInfo moxtype;

    RegionInfo regionmox;
    regionmox.mass_ = 180;
    regionmox.fractions[922350] = 0.00654120;
    regionmox.fractions[922380] = 0.91342476;
    regionmox.fractions[942380] = 0.00216000;
    regionmox.fractions[942390] = 0.04480000;
    regionmox.fractions[942400] = 0.02072000;
    regionmox.fractions[942410] = 0.00592000;
    regionmox.fractions[942420] = 0.00584000;
    regionmox.fractions[952440] = 0.00056000;
    regionmox.BuildIso(lwrlib);

    moxtype.batch.assign(4, regionmox);

    ReactorXInfo reactor;
    reactor.type.push_back(lwrtype);
    reactor.type.push_back(moxtype);
    reactor.core_mass_ = regionlwr.mass_ + regionmox.mass_;
    reactor.abs_flux_tol_ = 0.001;
    reactor.base_flux_ = 3E14;
    reactor.DA_mode_ = 0;
    reactor.fluence_timestep_ = 50;
    reactor.flux_mode_ = 1;
    reactor.SS_tol_ = 0.001;


    /// Run steady state calculations

    /// Output results

    cout << "--Simulation reached the end--" << endl;
    return 0;
}








































