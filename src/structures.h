#ifndef STRUCTURES_H_INCLUDED
#define STRUCTURES_H_INCLUDED

#include <utility>
#include <map>
#include <string>
#include <cmath>
#include <sstream>


struct Daughter {
    unsigned int name;
    string string_name;
    std::vector<float> mass;
    float fraction;
};

class IsoInfo {
public:
    unsigned int name;          // Nucid if uncollapsed
    float fraction;             // Fraction of this isotope in the inherited class
    std::vector<float> fluence;
    std::vector<float> neutron_prod;
    std::vector<float> neutron_dest;
    std::vector<float> BU;      //total BU
    std::vector<float> fission_products; // Fission product lump mass
    std::vector<Daughter> iso_vector;

    void Print(int times = 5);
};

struct nonActinide {
    int name;
    float sng;
    float scattering;
    float sn2n;
    float snp;
    float sngx;
    float sn2nx;
    float yyn;
    float total_prod;
    float total_dest;
};

// Thermal disadvantage factor parameters for a given region
struct DisadvParams {
    float a;            // Fuel radius [cm]
    float b;            // Moderator radius [cm]
    float mod_Sig_a;     // Moderator Sig a
    float mod_Sig_s;     // Mod Sig s
    float fuel_Sig_s;    // Fuel Sig s
};

// Spatial flux calculation parameters for a given region
struct SpatialParamsX {
    float region_area;  // Total area of the region in [cm2]
    float Sig_a;        // Macroscopic absorption  cross section [cm-1]
    float Sig_tr;       // Macroscopic transport cross section [cm-1]
    float nuSig_f;      // Macroscopic fission cross section times nu [cm-1]
};

// Spatial flux calculation parameters for a given region
struct SpatialParamsLite {
    float delta;                // [cm]
    float fuel_area;            // Total area of the region in [cm2]
    float spatial_mod_thickness;// [cm]

    // All cross-sections in [cm]
    float spatial_mod_Sig_a;
    float spatial_mod_Sig_tr;
    float spatial_mod_Sig_f;
    float spatial_fuel_Sig_tr;
};

struct LibInfo {
    std::string name;               // Matches library in folder

    // Parameters used in library database generation if available
    float lib_flux;        // Flux in [n/cm2]
    float lib_res_t;       // Fuel residence time [day]
    float lib_power;       // Power in [MW(th)]
    float lib_mass;        // Mass in [kg]
    float lib_BU;          // Burnup in [MWd/kgIHM]
    float lib_CR;          // Conversion ratio

    std::vector<IsoInfo> all_iso;   // All iso's in the library go here
    float fraction;        // Fraction of this library in core region
};

// Information about the fuel region
class RegionInfo {
public:
    float mass_;             // Mass of region

    IsoInfo iso;             // Collapsed, isoinfo for region
    std::map<int,float> fractions; // Name and fraction of each isotope for iso building
    std::map<int,float> comp; // Composition of the region

    unsigned int location_;  // Radial location of region, 1:center

    float fluence_ = 0;      // Fluence of this region [n/cm2]
    float rflux_ = 1;        // Relative flux of region
    ///TODO DA missing underscore! ffs
    float DA = 1;            // Disadvantage factor

    void Print();            // Displays info on the region on to terminal
    void BuildIso(LibInfo library);
    void UpdateComp();              // Updates the composition of isotopes at region fluence
    float CalcBU();                 // Returns the burnup/mass at region fluence
    float CalcBU(float fluence);    // Returns the burnup/mass at given fluence
    float CalcProd();               // Returns the neutron production at region fluence
    float CalcProd(float fluence);  // Returns the neutron production at given fluence
    float CalcDest();               // Return the neutron destruction at region fluence
    float CalcDest(float fluence);  // Return the neutron destruction at given fluence
    float CalcNuSigf();             // Return the macroscopic fission cross section [cm-1] times nu
    float CalcSiga();               // Return the macroscopic absorption cross section [cm-1]
};

class TypeInfo {
public:
    std::vector<RegionInfo> batch; // number of batches for this type of fuel
};

class ReactorLiteInfo {
public:
    // Initialized during startup
    unsigned int regions_;  // Total number of regions/batches
    float cycles_ = 0;          // Number of cycles fuel resides in core
    float thermal_pow_;     // Reactor thermal power [MWth]
    float core_mass_;       // Total mass of all fuel in [kg]
    float target_BU_;       // Target burnup in [MWd/kgIHM]
    float target_CR_;       // Target conversion ratio
    float pnl;              // Nonleakage probability
    float fluence_timestep_;// Fluence propagation time step [second]
    float base_flux_;       // Library base flux or last cycle flux
    std::vector<int> CR_fissile_; // List of fissile isotopes for CR calc

    float abs_flux_tol_;    // Convergence tolerance for absolute flux calculation
    float SS_tol_;          // Convergence tolerance for steady state fluence calculation

    unsigned int flux_mode_;// Flux calculation mode:
    // 0: Equal Power Share, 1:Uniform, 2:Inv.Neut.Prod, 3:Spatial
    unsigned int DA_mode_;  // Disadvantage factor. 0:OFF, 1:ON

    SpatialParamsLite spatial_; // Spatial flux calculation parameters
    DisadvParams DA_;           // DA calculation parameters

    float struct_prod_ = 0;  // Non-fuel material neutron prod
    float struct_dest_ = 0;  // Non-fuel material neutron dest

    std::vector<std::string> libraries_;
    LibInfo library_;
    float CR_;

    // Regions are populated based on reactor parameters
    std::vector<RegionInfo> region;  // region[0] is oldest

    bool ConsistencyCheck();    // Makes sure there are no contradicting inputs
    void PrintRegionIsos();
    void PrintFluences();
    void BuildRegionIsos();
    void Reorder();             // Reorders regions from lowest k to highest
    void UpdateComp();          // Updates the composition of isotopes in all regions
    float CalcBU();             // Calculates the burnup of the core
    float CalcBU(float flux);   // Calculates the burnup if reactor at given flux
};

// A fuel type is based on one library
class FuelType {
    unsigned int f_groups_; // Fuel groups of this fuel type
    LibInfo library_;       // The library of the fuel group


};

class ReactorXInfo {

public:
    // Initialized during startup
    float thermal_pow_;     // Reactor thermal power [MWth]
    //float target_BU_;       // Target burnup in [MWd/kgIHM]
    //float target_CR_;       // Target conversion ratio
    //float pnl;              // Nonleakage probability
    float fluence_timestep_;// Fluence propagation time step [second]
    float base_flux_;       // Library base flux or last cycle flux
    std::vector<int> CR_fissile_; // List of fissile isotopes for CR calc

    float abs_flux_tol_;    // Convergence tolerance for absolute flux calculation
    float SS_tol_;          // Convergence tolerance for steady state fluence calculation

    unsigned int flux_mode_;// Flux calculation mode:
    // 0: Equal Power Share, 1:Uniform, 2:Inv.Neut.Prod, 3:Spatial
    unsigned int DA_mode_;  // Disadvantage factor. 0:OFF, 1:ON

    std::vector<std::string> libraries_;
    LibInfo library_;
    float CR_;

    // Regions are populated based on reactor parameters
    unsigned int regions_;  // calculated total number of regions
    float core_mass_;       // Total mass of all fuel in [kg]
    std::vector<TypeInfo> type;  // fuel types, such as UO2 and MOX

    void PrintRegionIsos();
    void PrintFluences();
    void BuildRegionIsos();
    void Reorder();             // Reorders regions from lowest k to highest
    void UpdateComp();          // Updates the composition of isotopes in all regions
    float CalcBU();             // Calculates the burnup of the core
    float CalcBU(float flux);   // Calculates the burnup if reactor at given flux

};

extern std::map<std::string, LibInfo> global_libs;
#endif // STRUCTURES_H_INCLUDED







