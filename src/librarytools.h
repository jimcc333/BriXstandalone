#ifndef LIBRARYTOOLS_H
#define LIBRARYTOOLS_H

#include <string>
#include <sstream>
#include <cmath>


// Used for DACalc
#include <boost/math/special_functions/bessel.hpp>

// Used for spatial flux calc
//#include <eigen3/Eigen/Core>
//#include <eigen3/Eigen/Eigenvalues>

using namespace std;

void LibraryReader(string library_name, string library_path, LibInfo &library);
void IsoBuilder(string library_path, IsoInfo &iso);
float FluxFinder(string library_path);
void StructReader(string library_path, float &struct_prod, float &struct_dest);
/*
void DACalc(ReactorLiteInfo &reactor_core);
void BurnFuel(ReactorLiteInfo &reactor_core);
void CriticalityBurn(ReactorLiteInfo &reactor_core);
void FluxCalc(ReactorLiteInfo &reactor_core);
float kCalc(ReactorLiteInfo &reactor_core);
void EqPowPhi(ReactorLiteInfo &reactor_core);
void InvProdPhi(ReactorLiteInfo &core);
float RegionCRCalc(ReactorLiteInfo &core, unsigned const int reg_i);
float CoreCRCalc(ReactorLiteInfo &core);
float AbsFluxCalc(ReactorLiteInfo &core, float abs_flux);
void SpatialPhi(ReactorLiteInfo &core);
float SteadyStateFluence(ReactorLiteInfo core, IsoInfo iso);
*/
#endif // LIBRARYTOOLS_H
