/**
by: Cem Bagdatlioglu

A simple code to test BriXsuite cases for simulation in Cyclus

Also will be used for reactor library testing

Assumptions:
- Batches of a type of fuel cannot have different masses
- Removed structural material calculations from kcalc
- Different libraries have the same fluence steps (equal power flux)

**/
#include <iostream>
#include <fstream>
#include <vector>

#include "generaltools.cc"
#include "structures.cc"
#include "librarytools.cc"

using namespace std;

void LWRsolution(ReactorXInfo &reactor);

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
    lwrtype.mass = 180;

    RegionInfo regionlwr;
    regionlwr.fractions[922350] = 0.0530;
    regionlwr.fractions[922380] = 0.9470;
    regionlwr.BuildIso(lwrlib);

    lwrtype.batch.assign(4, regionlwr);

    TypeInfo moxtype;
    moxtype.mass = 84; // watchout

    RegionInfo regionmox;
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
    reactor.type.push_back(lwrtype);
    reactor.pnl = 0.8215;
    reactor.thermal_pow_ = 41;
    reactor.core_mass_ = lwrtype.mass + moxtype.mass;
    reactor.abs_flux_tol_ = 0.01;
    reactor.base_flux_ = 3E18;
    reactor.DA_mode_ = 0;
    reactor.fluence_timestep_ = 50;
    reactor.flux_mode_ = 1;
    reactor.SS_tol_ = 0.001;
    reactor.batches = 0;
    for(int type_i = 0; type_i < reactor.type.size(); type_i++) {
        for(int batch_i = 0; batch_i < reactor.type[type_i].batch.size(); batch_i++) {
            reactor.batches++;
        }
    }


    /// Run steady state calculations
/*
    for(int cycle = 0; cycle < 3; cycle++) {
        reactor.PrintFluences();
        CriticalityBurn(reactor);
        //reactor.PrintFluences();
        for(int type_i = 0; type_i < reactor.type.size(); type_i++) {
            for(int batch_i = 0; batch_i < reactor.type[type_i].batch.size() - 1; batch_i++) {
                reactor.type[type_i].batch[batch_i].fluence_ = reactor.type[type_i].batch[batch_i+1].fluence_;
            }
            reactor.type[type_i].batch.back().fluence_ = 0;
        }
    }
*/

    reactor.type[1].batch[3] = regionmox;
    reactor.type[1].mass = 84;

    for(int cycle = 0; cycle < 4; cycle++) {
        cout << "BURN" << endl;
        CriticalityBurn(reactor);

        reactor.type[0].batch[0].fluence_ = reactor.type[0].batch[1].fluence_;
        reactor.type[0].batch[1].fluence_ = reactor.type[0].batch[2].fluence_;
        reactor.type[0].batch[2].fluence_ = reactor.type[0].batch[3].fluence_;
        reactor.type[0].batch[3].fluence_ = 0;

        reactor.type[1].batch[0] = reactor.type[1].batch[1];
        reactor.type[1].batch[1] = reactor.type[1].batch[2];
        reactor.type[1].batch[2] = reactor.type[1].batch[3];
        reactor.type[1].batch[3] = regionmox;
    }


    cout << " UO2 burnup: " << reactor.type[0].batch[0].CalcBU() << endl;
    cout << " MOX burnup: " << reactor.type[1].batch[0].CalcBU() << endl;
    cout << " core burnup: " << reactor.CalcBU() << endl;

    float assembly_BU = (reactor.type[0].batch[0].CalcBU() * reactor.type[0].mass + reactor.type[1].batch[0].CalcBU() * reactor.type[1].mass)
                        / reactor.core_mass_;

    cout << " assembly burnup: " << reactor.AssemblyBU(0) << endl;

    /// Output compositions
    reactor.type[0].batch[0].UpdateComp();
    reactor.type[1].batch[0].UpdateComp();

    cout << "UO2 composition: " << endl;
    cout << " U234  " << reactor.type[0].batch[0].comp["U234"] << endl;
    cout << " U235  " << reactor.type[0].batch[0].comp["U235"] << endl;
    cout << " U236  " << reactor.type[0].batch[0].comp["U236"] << endl;
    cout << " U238  " << reactor.type[0].batch[0].comp["U238"] << endl;
    cout << " PU238 " << reactor.type[0].batch[0].comp["PU238"] << endl;
    cout << " PU239 " << reactor.type[0].batch[0].comp["PU239"] << endl;
    cout << " PU240 " << reactor.type[0].batch[0].comp["PU240"] << endl;
    cout << " PU241 " << reactor.type[0].batch[0].comp["PU241"] << endl;
    cout << " PU242 " << reactor.type[0].batch[0].comp["PU242"] << endl;
    cout << " CS137 " << reactor.type[0].batch[0].comp["CS137"] << endl;

    cout << "MOX composition: " << endl;
    cout << " U234  " << reactor.type[1].batch[0].comp["U234"] << endl;
    cout << " U235  " << reactor.type[1].batch[0].comp["U235"] << endl;
    cout << " U236  " << reactor.type[1].batch[0].comp["U236"] << endl;
    cout << " U238  " << reactor.type[1].batch[0].comp["U238"] << endl;
    cout << " PU238 " << reactor.type[1].batch[0].comp["PU238"] << endl;
    cout << " PU239 " << reactor.type[1].batch[0].comp["PU239"] << endl;
    cout << " PU240 " << reactor.type[1].batch[0].comp["PU240"] << endl;
    cout << " PU241 " << reactor.type[1].batch[0].comp["PU241"] << endl;
    cout << " PU242 " << reactor.type[1].batch[0].comp["PU242"] << endl;
    cout << " CS137 " << reactor.type[1].batch[0].comp["CS137"] << endl;
    cout << "--Simulation reached the end--" << endl;
    return 0;
}



void LWRsolution(ReactorXInfo &reactor) {
    /// A about 42 MWd/kg burnup LWR reactor example
    /// Yes, the author recognizes this is a TERRIBLE way to save "inputs"

    // Generate necessary libraries and path strings
    LibInfo lwrlib;
    string lwrpath = "libs/extLWR";

    // Populate the libraries
    LibraryReader("lwrlib", lwrpath, lwrlib);

    // Build reactor object and fuel isos
    TypeInfo lwrtype;
    lwrtype.mass = 112000;

    RegionInfo regionlwr;
    regionlwr.fractions[922350] = 0.036;
    regionlwr.fractions[922380] = 0.964;
    regionlwr.BuildIso(lwrlib);

    lwrtype.batch.assign(4, regionlwr);

    reactor.type.push_back(lwrtype);
    reactor.pnl = 0.97;
    reactor.thermal_pow_ = 4170;
    reactor.core_mass_ = lwrtype.mass ;
    reactor.abs_flux_tol_ = 0.01;
    reactor.base_flux_ = 3E18;
    reactor.DA_mode_ = 0;
    reactor.fluence_timestep_ = 50;
    reactor.flux_mode_ = 1;
    reactor.SS_tol_ = 0.001;
    reactor.batches = 0;
    for(int type_i = 0; type_i < reactor.type.size(); type_i++) {
        for(int batch_i = 0; batch_i < reactor.type[type_i].batch.size(); batch_i++) {
            reactor.batches++;
        }
    }


    for(int cycle = 0; cycle < 20; cycle++) {

        CriticalityBurn(reactor);

        for(int type_i = 0; type_i < reactor.type.size(); type_i++) {
            for(int batch_i = 0; batch_i < reactor.type[type_i].batch.size() - 1; batch_i++) {
                reactor.type[type_i].batch[batch_i].fluence_ = reactor.type[type_i].batch[batch_i+1].fluence_;
            }
            reactor.type[type_i].batch.back().fluence_ = 0;
        }
    }

    cout << " discharge burnup: " << reactor.type[0].batch[0].CalcBU() << endl;
    cout << " core burnup: " << reactor.CalcBU() << endl;


}




































