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

float SSCalc(ReactorXInfo &reactor);
void LWRsolution(ReactorXInfo &reactor);

int main(int argc, char* argv[]) {

    cout << "--New simulations starts--" << endl;

    /// Read fuel info and libraries

    // Generate necessary libraries and path strings
    LibInfo lwrlib, moxlib;
    string lwrpath = "libs/E5_50", moxpath = "libs/PWRMOX";

    // Populate the libraries
    LibraryReader("lwrlib", lwrpath, lwrlib);
    LibraryReader("moxlib", moxpath, moxlib);

    // Build reactor object and fuel isos
    TypeInfo lwrtype;
    lwrtype.mass = 180;

    RegionInfo regionlwr;
    float enrichment = 0.0530;
    regionlwr.fractions[922350] = enrichment;
    regionlwr.fractions[922380] = 1-enrichment;
    regionlwr.BuildIso(lwrlib);

    lwrtype.batch.assign(4, regionlwr);

    TypeInfo moxtype;
    moxtype.mass = 84;

    RegionInfo regionmox;
    regionmox.fractions[922350] = 0.00654120;
    regionmox.fractions[922380] = 0.91342476;
    regionmox.fractions[942380] = 0.00216000;
    regionmox.fractions[942390] = 0.04480000;
    regionmox.fractions[942400] = 0.02072000;
    regionmox.fractions[942410] = 0.00592000;
    regionmox.fractions[942420] = 0.00584000;
    regionmox.fractions[952410] = 0.00056000;
    regionmox.BuildIso(moxlib);

    moxtype.batch.assign(4, regionmox);

    ReactorXInfo reactor;
    reactor.type.push_back(lwrtype);
    reactor.type.push_back(moxtype);
    reactor.pnl = 0.95;
    reactor.thermal_pow_ = 401;
    reactor.core_mass_ = lwrtype.mass + moxtype.mass;
    reactor.abs_flux_tol_ = 0.01;
    reactor.base_flux_ = 3E19;
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
    RegionInfo nextmox;


    /// Run steady state calculations
    // Run 1
    float target = 60;


    while(SSCalc(reactor) > target) {
        reactor.pnl -= 0.001;
    }
    cout << "Final pnl: " << reactor.pnl << " SS BU:" << SSCalc(reactor) << endl;


    for(int i = 0; i < 6; i++) {
        // Update compositions
        reactor.type[0].batch[0].UpdateComp();
        reactor.type[1].batch[0].UpdateComp();

        float total = reactor.type[1].batch[0].comp["U235"] + reactor.type[1].batch[0].comp["U238"]  + reactor.type[1].batch[0].comp["PU238"] +
        reactor.type[1].batch[0].comp["PU239"] + reactor.type[1].batch[0].comp["PU240"] + reactor.type[1].batch[0].comp["PU241"] +
        reactor.type[1].batch[0].comp["PU242"] + reactor.type[1].batch[0].comp["AM241"];

        cout << "total:" << total << endl;

        nextmox.fractions[922350] = reactor.type[1].batch[0].comp["U235"] + (1-total)*regionmox.fractions[922350];
        nextmox.fractions[922380] = reactor.type[1].batch[0].comp["U238"] + (1-total)*regionmox.fractions[922380];
        nextmox.fractions[942380] = reactor.type[1].batch[0].comp["PU238"] + (1-total)*regionmox.fractions[942380];
        nextmox.fractions[942390] = reactor.type[1].batch[0].comp["PU239"] + (1-total)*regionmox.fractions[942390];
        nextmox.fractions[942400] = reactor.type[1].batch[0].comp["PU240"] + (1-total)*regionmox.fractions[942400];
        nextmox.fractions[942410] = reactor.type[1].batch[0].comp["PU241"] + (1-total)*regionmox.fractions[942410];
        nextmox.fractions[942420] = reactor.type[1].batch[0].comp["PU242"] + (1-total)*regionmox.fractions[942420];
        nextmox.fractions[952410] = reactor.type[1].batch[0].comp["AM241"] + (1-total)*regionmox.fractions[952410];

        nextmox.iso.Clear();
        nextmox.BuildIso(moxlib);

        reactor.type[1].batch[0] = nextmox;
        reactor.type[1].batch[1] = nextmox;
        reactor.type[1].batch[2] = nextmox;
        reactor.type[1].batch[3] = nextmox;

        float newtot = nextmox.fractions[922350] + nextmox.fractions[922380] + nextmox.fractions[942380] + nextmox.fractions[942390]
        + nextmox.fractions[942400]+ nextmox.fractions[942410] + nextmox.fractions[942420]+    nextmox.fractions[952410];

        cout << "new total: " << newtot << endl << endl;
        /*
        cout << regionmox.fractions[922350] << "  " << nextmox.fractions[922350] << endl;
        cout << regionmox.fractions[922380] << "  " << nextmox.fractions[922380] << endl;
        cout << regionmox.fractions[942380] << "  " << nextmox.fractions[942380] << endl;
        cout << regionmox.fractions[942390] << "  " << nextmox.fractions[942390] << endl;
        cout << regionmox.fractions[942400] << "  " << nextmox.fractions[942400] << endl;
        cout << regionmox.fractions[942410] << "  " << nextmox.fractions[942410] << endl;
        cout << regionmox.fractions[942420] << "  " << nextmox.fractions[942420] << endl;
        cout << regionmox.fractions[952410] << "  " << nextmox.fractions[952410] << endl;*/

        while(SSCalc(reactor) < target) {
            enrichment += 0.0001;

            regionlwr.iso.Clear();
            regionlwr.fractions[922350] = enrichment;
            regionlwr.fractions[922380] = 1-enrichment;

            regionlwr.BuildIso(lwrlib);

            reactor.type[0].batch[0] = regionlwr;
            reactor.type[0].batch[1] = regionlwr;
            reactor.type[0].batch[2] = regionlwr;
            reactor.type[0].batch[3] = regionlwr;
        }

        cout << "Cycle " << i+2 << " enrichment: " << enrichment << endl;

        cout << nextmox.fractions[942380] << endl;
        cout << nextmox.fractions[942390] << endl;
        cout << nextmox.fractions[942400] << endl;
        cout << nextmox.fractions[942410] << endl;
        cout << nextmox.fractions[942420] << endl;
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



/*
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
    cout << " CS137 " << reactor.type[1].batch[0].comp["CS137"] << endl;*/
    cout << "--Simulation reached the end--" << endl;
    return 0;
}


float SSCalc(ReactorXInfo &reactor) {
    /// uses the burnup of the oldest assembly for this calculation!
    // all fuel isos are correct when this is called
    float bu_prev = reactor.AssemblyBU(0);
    float burnup = 5;
    int iter = 0;

    while(abs(bu_prev - burnup)/burnup > 0.004) {
        iter++;
        if(iter > 60) {
            cout << "!Too many iterations for SS calc!" << endl;
            return burnup;
        }

        bu_prev = reactor.AssemblyBU(0);

        // New cycle
        reactor.type[0].batch[0].fluence_ = reactor.type[0].batch[1].fluence_;
        reactor.type[0].batch[1].fluence_ = reactor.type[0].batch[2].fluence_;
        reactor.type[0].batch[2].fluence_ = reactor.type[0].batch[3].fluence_;
        reactor.type[0].batch[3].fluence_ = 0;

        reactor.type[1].batch[0].fluence_ = reactor.type[1].batch[1].fluence_;
        reactor.type[1].batch[1].fluence_ = reactor.type[1].batch[2].fluence_;
        reactor.type[1].batch[2].fluence_ = reactor.type[1].batch[3].fluence_;
        reactor.type[1].batch[3].fluence_ = 0;

        // Burn the fuel (update fluences)
        CriticalityBurn(reactor);
        burnup = reactor.AssemblyBU(0);

    }
    if(iter < 4) {
        cout << "  ----Warning! SS burnup found after " << iter+1 << " iterations" << endl;
    }

    return burnup;
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




































