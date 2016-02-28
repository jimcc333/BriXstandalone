
#include <iostream>
#include <fstream>
#include <vector>

#include "generaltools.cc"
#include "structures.cc"
#include "librarytools.cc"

using namespace std;

int main(int argc, char* argv[]) {

    cout << "New simulations starts" << endl;

    /// Read fuel info and libraries
    // Generate necessary libraries and path strings
    LibInfo lwrlib, moxlib;
    string lwrpath = "libs/extLWR", moxpath = "libs/PWRMOX";

    // Populate the libraries
    LibraryReader("lwrlib", lwrpath, lwrlib);
    LibraryReader("moxlib", moxpath, moxlib);

    // Build reactor object
    RegionInfo regionlwr;
    regionlwr.mass_ = 180;
    regionlwr.fractions[922350] = 0.0530;
    regionlwr.fractions[922380] = 0.9470;

    ReactorXInfo reactor;
    reactor.region.push_back(regionlwr);

    // Fuel compositions and build fuel isos

    // Run steady state calculations

    // Output results

    cout << "--Simulation reached the end--" << endl;
    return 0;
}
