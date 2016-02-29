#include "structures.h"


std::map<std::string, LibInfo> global_libs;
// Prints the isotope data to terminal
void IsoInfo::Print(int times) {
    std::cout << std::endl << "Isotope lib name: " << name << " fraction: "
              << fraction << " k[0]: " << neutron_prod[0]/neutron_dest[0] << std::endl;

    if (fluence.size() < 1) {
        std::cout << "  !Region not built." << std::endl;
        return;
    }

    if (times > fluence.size()) {times = fluence.size();}

    std::cout << "  Fluence: ";
    for (int i = 0; i < times; i++) {
        std::cout << fluence[i] << " ";
    } std::cout << std::endl;

    std::cout << "  NProd: ";
    for (int i = 0; i < times; i++) {
        std::cout << neutron_prod[i] << " ";
    } std::cout << std::endl;

    std::cout << "  NDest: ";
    for (int i = 0; i < times; i++) {
        std::cout << neutron_dest[i] << " ";
    } std::cout << std::endl;

    std::cout << "  k: ";
    for (int i = 0; i < times; i++) {
        std::cout << neutron_prod[i]/neutron_dest[i] << " ";
    } std::cout << std::endl;

    std::cout << "  BU: ";
    for (int i = 0; i < times; i++) {
        std::cout << BU[i] << " ";
    } std::cout << std::endl;

    std::cout << "  ------------" << std::endl;

}

void RegionInfo::Print() {
    std::cout << "Fluence: " << fluence_ << " rFlux: " << rflux_ << " mass: " << mass_;
    iso.Print();
}

// Updates the iso in that region with the fraction
void RegionInfo::BuildIso(LibInfo library) {

    iso.name = 0; // Zero represents a collapsed iso and not a specific isotope
    //For every available isotope in the library
    for (unsigned int lib_i = 0; lib_i < library.all_iso.size(); lib_i++) {
        // Save the name of the current library isotope
        unsigned const int iso_name = library.all_iso[lib_i].name;
        //
        if (library.all_iso[lib_i].iso_vector.size() > 0) {
            if (iso.fluence.size() < 1) {
                for (int i = 0; i < library.all_iso[lib_i].fluence.size(); i++) {
                    iso.fluence.push_back(library.all_iso[lib_i].fluence[i]);
                }
                for (int i = 0; i < library.all_iso[lib_i].neutron_prod.size(); i++) {
                    iso.neutron_prod.push_back(fractions[iso_name] * library.all_iso[lib_i].neutron_prod[i]);
                }
                for (int i = 0; i < library.all_iso[lib_i].neutron_dest.size(); i++) {
                    iso.neutron_dest.push_back(fractions[iso_name] * library.all_iso[lib_i].neutron_dest[i]);
                }
                for (int i = 0; i < library.all_iso[lib_i].BU.size(); i++) {
                    iso.BU.push_back(fractions[iso_name] * library.all_iso[lib_i].BU[i]);
                }
                for (int i = 0; i < library.all_iso[lib_i].iso_vector.size(); i++) {
                    iso.iso_vector.push_back(library.all_iso[lib_i].iso_vector[i]);
                    for(int k = 0; k < iso.iso_vector[i].mass.size(); k++) {
                        iso.iso_vector[i].mass[k] = fractions[iso_name] * iso.iso_vector[i].mass[k];
                    }
                }
            } else {

                for (int i = 0; i < library.all_iso[lib_i].neutron_prod.size(); i++) {
                    iso.neutron_prod[i] += fractions[iso_name] * library.all_iso[lib_i].neutron_prod[i];
                }
                for (int i = 0; i < library.all_iso[lib_i].neutron_dest.size(); i++) {
                    iso.neutron_dest[i] += fractions[iso_name] * library.all_iso[lib_i].neutron_dest[i];
                }
                for (int i = 0; i < library.all_iso[lib_i].BU.size(); i++) {
                    iso.BU[i] += fractions[iso_name] * library.all_iso[lib_i].BU[i];
                }
                for (int i = 0; i < library.all_iso[lib_i].iso_vector.size(); i++) {
                    bool iso_check = true;
                    for (int j = 0; j < iso.iso_vector.size(); j++) {
                        if (library.all_iso[lib_i].iso_vector[i].name == iso.iso_vector[j].name) {
                            for (int k = 0; k < iso.iso_vector[j].mass.size(); k++) {
                                for (int ii = 0; ii < library.all_iso[lib_i].iso_vector[i].mass.size(); ii ++) {
                                    if ( k ==ii ) {
                                        iso.iso_vector[j].mass[k] += fractions[iso_name] *library.all_iso[lib_i].iso_vector[i].mass[ii];
                                    }
                                }
                            }
                            iso_check = false;
                        }
                    }
                    if (iso_check == true) {
                        iso.iso_vector.push_back(library.all_iso[lib_i].iso_vector[i]);
                        for(int k = 0; k < iso.iso_vector[iso.iso_vector.size()-1].mass.size()-1; k++) {
                            iso.iso_vector[iso.iso_vector.size()-1].mass[k] = iso.iso_vector[iso.iso_vector.size()-1].mass[k]*fractions[iso_name] ;
                        }
                    }
                }
            }
        } else {
            if (iso.neutron_prod.size() > 1) {
                for (int i = 0; i < library.all_iso[lib_i].neutron_prod.size(); i++) {
                    iso.neutron_prod[i] = iso.neutron_prod[i] + fractions[iso_name] *library.all_iso[lib_i].neutron_prod[i];
                }
                for (int i = 0; i < library.all_iso[lib_i].neutron_dest.size(); i++) {
                    iso.neutron_dest[i] = iso.neutron_dest[i] + fractions[iso_name] *library.all_iso[lib_i].neutron_dest[i];
                }
            } else {
                for (int i = 0; i < library.all_iso[lib_i].neutron_prod.size(); i++) {
                    iso.neutron_prod.push_back(fractions[iso_name] *library.all_iso[lib_i].neutron_prod[i]+(fractions[iso_name] )*library.all_iso[1].neutron_prod[i]);
                }
                for (int i = 0; i < library.all_iso[lib_i].neutron_dest.size(); i++) {
                    iso.neutron_dest.push_back(fractions[iso_name] *library.all_iso[lib_i].neutron_dest[i]+(fractions[iso_name] )*library.all_iso[1].neutron_dest[i]);
                }
            }
        }
    }
}

// Determines the composition of the given isotope at region fluence
void RegionInfo::UpdateComp() {
    if(fluence_ <= 0) {return;}
    if(fluence_ >= iso.fluence.back()) {
        std::cout << " Warning! Region fluence exceeds library range in composition calculation!"
                << " Composition will not be updated." << std::endl;
        return;
    }

    unsigned int ii;

    for(ii = 1; iso.fluence[ii] <= fluence_; ii++){}

    // Finds the interpolation slope to use for faster interpolation
    const double slope = (fluence_ - iso.fluence[ii-1]) / (iso.fluence[ii] - iso.fluence[ii-1]);

    for(int iso_i = 0; iso_i < iso.iso_vector.size(); iso_i++) {
        comp[iso.iso_vector[iso_i].name] = (iso.iso_vector[iso_i].mass[ii-1] +
                (iso.iso_vector[iso_i].mass[ii] - iso.iso_vector[iso_i].mass[ii-1]) * slope) / 1000;
    }
}

// Determines the burnup/mass of the region based on the region fluence
float RegionInfo::CalcBU() {
    return CalcBU(fluence_);
}

// Determines the burnup/mass of the region based on the given fluence
float RegionInfo::CalcBU(float fluence) {
    if(fluence <= 0) {return 0;}
    if(fluence >= iso.fluence.back()) {return iso.BU.back();}

    unsigned int ii;
    float burnup;

    for(ii = 1; iso.fluence[ii] <= fluence; ii++){}

    burnup = Interpolate(iso.BU[ii-1], iso.BU[ii], iso.fluence[ii-1], iso.fluence[ii], fluence);
    if(burnup < 0) {return 0;}
    return burnup;
}

// Determines the neutron production of the region based on region fluence
float RegionInfo::CalcProd() {
    return CalcProd(fluence_);
}

// Determines the neutron production of the region based on given fluence
float RegionInfo::CalcProd(float fluence) {
    if(fluence < 0){return 0;}
    if(fluence >= iso.fluence.back()){return iso.neutron_prod.back();}

    unsigned int ii;

    for(ii = 1; iso.fluence[ii] <= fluence; ii++){}

    return Interpolate(iso.neutron_prod[ii-1], iso.neutron_prod[ii],
                       iso.fluence[ii-1], iso.fluence[ii], fluence);
}

// Determines the neutron production of the region based on region fluence
float RegionInfo::CalcDest() {
    return CalcDest(fluence_);
}

// Determines the neutron destruction of the region based on given fluence
float RegionInfo::CalcDest(float fluence) {
    if(fluence < 0){return 0;}
    if(fluence >= iso.fluence.back()){return iso.neutron_dest.back();}

    unsigned int ii;

    for(ii = 1; iso.fluence[ii] <= fluence; ii++){}

    return Interpolate(iso.neutron_dest[ii-1], iso.neutron_dest[ii],
                       iso.fluence[ii-1], iso.fluence[ii], fluence);
}

// Determines nu Sig_f [cm-1] from neutron production rate
float RegionInfo::CalcNuSigf() {
    return CalcProd() * 0.01097;
}

// Determines Sig_a [cm-1] from neutron destruction rate
float RegionInfo::CalcSiga() {
    return CalcDest() * 0.01097;
}

// Checks variables for consistency
bool ReactorLiteInfo::ConsistencyCheck() {
    unsigned int regions_;  // Total number of regions/batches
    float thermal_pow_;     // Reactor thermal power [MWth]
    float core_mass_;       // Total mass of all fuel in [kg]
    float target_BU_;       // Target burnup in [MWd/kgIHM]
    float target_CR_;       // Target conversion ratio
    float pnl;              // Nonleakage probability
    float fluence_timestep_;// Fluence propagation time step [second]
    float base_flux_;       // Library base flux or last cycle flux
    std::vector<int> CR_fissile_; // List of fissile isotopes for CR calc

    float SS_tol_;          // Convergence tolerance for steady state fluence calculation

    if(flux_mode_ < 0 || flux_mode_ > 3) {
        std::cout << "Unrecognized flux mode " << flux_mode_ << std:: endl;
        return false;
    }

    if(DA_mode_ != 0 || DA_mode_ != 1) {
        std::cout << "Unrecognized DA mode " << DA_mode_ << std::endl;
        return false;
    }

    if(abs_flux_tol_ < 0.00001 || SS_tol_ < 0.00001) {
        return false;
    }

    if(abs_flux_tol_ > 0.4 || SS_tol_ > 0.4) {
        return false;
    }


}

// Prints the iso of each region
void ReactorLiteInfo::PrintRegionIsos() {
    std::cout << global_libs[libraries_[0]].name << " region isos:" << std::endl;
    for(unsigned int reg_i = 0; reg_i < region.size(); reg_i++) {
        region[reg_i].Print();;
    }
    std::cout << "-------------------------------" << std::endl;
}

// Prints the fluence of each region to terminal
void ReactorLiteInfo::PrintFluences() {
    std::cout << global_libs[libraries_[0]].name << " fluences: ";
    for(unsigned int reg_i = 0; reg_i < region.size(); reg_i++) {
        std::cout << region[reg_i].fluence_ << "  ";
    }
    std::cout << std::endl;
}

// Goes through all regions and builds their isos
void  ReactorLiteInfo::BuildRegionIsos() {
    for (unsigned int reg_i = 0; reg_i < region.size(); reg_i++) {
        // If a region has zero fluence the iso must not be built yet
        if (region[reg_i].fluence_ == 0) {
            region[reg_i].BuildIso(global_libs[libraries_[0]]);
        }
    }
}

// Reorders regions so that lowest k is at entry zero
void ReactorLiteInfo::Reorder() {
    std::vector<RegionInfo> temp_region = region;
    region.clear();
    double k0, k1;

    while(temp_region.size() != 0) {
        unsigned int lowest = 0;
        for(unsigned int remain_i = 1; remain_i < temp_region.size(); remain_i++) {
            // The first two discrete points are used as opposed to the very first
            k0 = temp_region[lowest].iso.neutron_prod[0]
                 / temp_region[lowest].iso.neutron_dest[0];
            k1 = (temp_region[remain_i].iso.neutron_prod[0])
                 / (temp_region[remain_i].iso.neutron_dest[0]);

            if(k0 > k1) {
                lowest = remain_i;
            }
        }
        region.push_back(temp_region[lowest]);
        temp_region.erase(temp_region.begin() + lowest);
    }
}

// Updates the composition of isotopes in all regions
void ReactorLiteInfo::UpdateComp() {
    for(unsigned int reg_i = 0; reg_i < region.size(); reg_i++) {
        region[reg_i].UpdateComp();
    }
}

// Returns the total burnup/mass of the core by going through each batch
float ReactorLiteInfo::CalcBU() {
    unsigned const int regions = region.size();
    float total_BU = 0;
    for(unsigned int reg_i = 0; reg_i < regions; reg_i++) {
        total_BU += region[reg_i].CalcBU() * region[reg_i].mass_;
    }
    return total_BU / core_mass_;
}

// Returns the burnup if the reactor at given flux; using rflux_ and fluence_timestep_
float ReactorLiteInfo::CalcBU(float flux) {
    unsigned const int regions = region.size();
    float burnup = 0;
    float fluence;

    for(int reg_i = 0; reg_i < regions; reg_i++) {
        fluence = region[reg_i].fluence_ + (region[reg_i].rflux_ * flux * fluence_timestep_);
        burnup += region[reg_i].CalcBU(fluence) * region[reg_i].mass_;
    }

    return burnup / core_mass_;
}


// Returns the total burnup/mass of the core by going through each batch
float ReactorXInfo::CalcBU() {

    return CalcBU(0);
}

float ReactorXInfo::CalcBU(float flux) {
    float total_BU = 0;
    float fluence;

    for(int type_i = 0; type_i < type.size(); type_i++) {
        for(int batch_i = 0; batch_i < type[type_i].batch.size(); batch_i++) {
            fluence = type[type_i].batch[batch_i].fluence_ +
                      type[type_i].batch[batch_i].rflux_ * flux * fluence_timestep_;
            total_BU += type[type_i].batch[batch_i].CalcBU(fluence);
            //cout << "BUcalc; fluence:" << fluence << " calcBU(fluence):" << type[type_i].batch[batch_i].CalcBU(fluence) << endl;
            //cout << "   batchfluence:" << type[type_i].batch[batch_i].fluence_ << " rflux:" << type[type_i].batch[batch_i].rflux_ << " fluencetimestep:" << fluence_timestep_ << " flux:" << flux << endl;
        }
        total_BU = total_BU * type[type_i].mass / core_mass_ / type[type_i].batch.size();
    }
    return total_BU;
}











