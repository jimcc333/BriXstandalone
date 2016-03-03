#include "librarytools.h"


// Reads from library in library_path to build all_iso
void LibraryReader(string library_name, string library_path, LibInfo &library) {
    library.name = library_name;

    ifstream inf (library_path + "/manifest.txt"); // Opens manifest file

    if(!inf){
        cout << "Error! LibraryReader failed to read from " << library_path << endl;
    }

    string line;
    string iso_name;

    // Put the name of every available isotope in the library in all_iso
    while (getline(inf, line)) {
        IsoInfo iso; // Temporary iso to push

        istringstream iss(line);

        iss >> iso_name;
        iso.name = atoi(iso_name.c_str()); // Adds name to iso

        iso.fraction = 0; // Default fraction

        IsoBuilder(library_path, iso);

        library.all_iso.push_back(iso);
    }

}

// Reads the isotope information and puts it in iso
void IsoBuilder(string library_path, IsoInfo &iso) {
    //type = cyclus::Env::GetInstallPath() + "/share/brightlite/" + type;
    const float flux_value = FluxFinder(library_path);

    ifstream inf(library_path + "/" + to_string(iso.name) + ".txt");
    if(!inf){
        cout << "Failed to read file for " + library_path + " " +  to_string(iso.name) << endl;
        return;
    }

    int i = 0;
    float value;
    string buffer;
    string line;

    while(getline(inf, line)){
        istringstream iss(line);
        iss >> buffer;
        if (i >= 4){
            Daughter daughter;
            ///changed here without knowing what it actually will do, from this:
            //daughter.name = pyne::nucname::zzaaam(buffer);
            /// to this:
            daughter.string_name = buffer; // (also had to change add daughter.string_name)
            while (iss >> value){
                daughter.mass.push_back(value);
            }
            iso.iso_vector.push_back(daughter);
        } else {
            while (iss >> value){
                if (i == 0){
                    iso.fluence.push_back(value*flux_value*84600);
                } else if (i == 1){
                    iso.neutron_prod.push_back(value);
                } else if (i == 2){
                    iso.neutron_dest.push_back(value);
                } else if (i == 3){
                    iso.BU.push_back(value);
                }
            }
            i++;
        }
    }
    inf.close();

    // Erases zero beginnings since these zero-zero points are nonphysical
    if(iso.neutron_prod[0] == 0 || iso.neutron_dest[0] == 0){
        iso.neutron_prod.erase(iso.neutron_prod.begin());
        iso.neutron_dest.erase(iso.neutron_dest.begin());
        iso.BU.erase(iso.BU.begin());
        iso.fluence.erase(iso.fluence.begin());
        for(int j = 0; j < iso.iso_vector.size(); j++){
            iso.iso_vector[j].mass.erase(iso.iso_vector[j].mass.begin());
        }
    }

    // Assumes derivative if the burnup vector is not increasing
    if(DecreaseChecker(iso.BU)) {
        CumulativeAdd(iso.BU); // Turns it into cumulative vector
    }
}



// Returns the library base flux value
float FluxFinder(string library_path) {

    // Open the params file where flux info is stored
    ifstream inf(library_path + "/params.txt");

    if(!inf){
        cout << "Error attempting to open " << library_path << "/params.txt" << endl;
        return 0;
    }
    string buffer;
    float value;
    string line;
    while(getline(inf, line)){

        istringstream iss(line);
        iss >> buffer >> value;
        if (buffer == "FLUX"){
            return value;
        }
    }
    cout << "Error in finding flux in " << library_path << "/params.txt" << endl;
    return 0;
}

// Finds the structural material contribution
void StructReader(string library_path, float &struct_prod, float &struct_dest) {
// Uses the compositions in structural.txt to combine cross section info in TAPE9.INP

    ifstream struct_file(library_path + "/structural.txt");
    if (!struct_file) {
        cout << "Error! StructReader can't open " << library_path << "/structural.txt  "
        << "Structural neutron production and destruction info will not be used." << endl;
        return;
    }

    ifstream tape9_file(library_path + "/TAPE9.INP");
    if (!tape9_file) {
        cout << "Error! StructReader can't open " << library_path << "/TAPE9.INP  "
        << "Structural neutron production and destruction info will not be used." << endl;
        return;
    }

    unsigned int nucid;
    float fraction;
    float tot_dest = 0;
    float tot_prod = 0;
    string line, library, iso;
    nonActinide x_sections;                     // This stores tape9 cross-sections
    vector<nonActinide> nonActinides;           // This stores isotope compositions
    float sng, sn2n, sna, snp, sngx, sn2nx, yyn;

    // Read the TAPE9 file for cross sections
    while (getline(tape9_file, line)) {
        istringstream iss(line);
        iss >> iso >> iso >> iso;

        if (iso == "MATERIAL"){
            while (getline(tape9_file, line)) {
                istringstream iss1(line);
                iss1 >> library;
                if (library == "-1"){
                    tape9_file.close();
                    for (int i = 0; i < nonActinides.size(); i++) {
                        nonActinides[i].total_prod = 2*nonActinides[i].sn2n + 2*nonActinides[i].sn2nx;
                        nonActinides[i].total_dest = nonActinides[i].snp + nonActinides[i].sng +
                        nonActinides[i].sngx + nonActinides[i].sn2n + nonActinides[i].sn2nx;
                    }
                    goto next;
                }
                iss1 >> iso >> sng >> sn2n >> snp >> sngx >> sn2nx >> yyn;
                x_sections.name = atoi(iso.c_str());
                x_sections.sng = sng;
                x_sections.sn2n = sn2n;
                x_sections.snp = snp;
                x_sections.sngx = sngx;
                x_sections.sn2nx = sn2nx;
                x_sections.yyn = yyn;
                nonActinides.push_back(x_sections);
            }
        }
    }

    next: // Breaks out of the previous nested loop

    // Combine total cross sections with mass compositions
	while (getline(struct_file, line)) {
        istringstream iss(line);
        iss >> nucid >> fraction;

        for(int i = 0; i < nonActinides.size(); i++){
            if(nonActinides[i].name == nucid){
                nucid = nucid % 10000;
                nucid = nucid / 10;

                tot_prod += nonActinides[i].total_prod * fraction * 1000*0.602/nucid;
                tot_dest += nonActinides[i].total_dest * fraction * 1000*0.602/nucid;
            }
        }
    }

    struct_prod = tot_prod;
    struct_dest = tot_dest;
}

/*
// Calculates fuel disadvantage factor
void DACalc(ReactorLiteInfo &core){
// DA = phi_thermal_Mod / phi_thermal_Fuel

    const float a = core.DA_.a;               // radius of the fuel rod
    const float b = core.DA_.b;               // radius of the equivalent cell
    const float Sig_sF = core.DA_.fuel_Sig_s; // macroscopic scatter CS of fuel
    const float Sig_sM = core.DA_.mod_Sig_s;  //macroscopic scatter CS of moderator
    const float Sig_aM = core.DA_.mod_Sig_a;  // macroscopic abs. CS of moderator

    float Sig_trF;          // macroscopic transport CS of fuel
    float Sig_tF;           // macroscopic total CS of fuel
    float D_F;              // diffusion coef. of fuel
    const float A_F = 235;  // A number of fuel
    float L_F;              // diffusion length of fuel
    float Sig_aF;           // macroscopic abs. CS of fuel
    const float V_F = pow(a,2)*3.141592; // Fuel volume

    const float Sig_tM = Sig_aM + Sig_sM;       // Macroscopic total CS of moderator
    const float A_M = 18;                       // A of moderator
    const float Sig_trM = Sig_tM - 2/3/A_M*Sig_sM;  // macroscopic transport CS of moderator
    const float D_M = 1 / (3 * Sig_trM);        // diffusion coef. of moderator
    const float L_M = sqrt(D_M/Sig_aM);         // diffusion length of moderator
    const float V_M = pow(b,2)*3.141592 - V_F;  // volume of moderator

    // calculated equivalent dimensions
    float x;
    const float y = a/L_M;
    const float z = b/L_M;
    float F, E;     // lattice functions
    float f;        // flux of fuel divided by total flux(fuel+moderator)

    for (int i = 0; i < core.region.size(); i++) {
        Sig_aF = core.region[i].CalcSiga();

        Sig_tF = Sig_aF+Sig_sF;
        Sig_trF = Sig_tF - 2/3/A_F*Sig_sF;
        D_F = 1 / (3 * Sig_trF);
        L_F = sqrt(D_F/Sig_aF);
        x = a/L_F;

        F = x * boost::math::cyl_bessel_i(0,x) / (2 * boost::math::cyl_bessel_i(1, x));

        E = (z*z - y*y) / (2 * y) * ( (boost::math::cyl_bessel_i(0, y) *
                boost::math::cyl_bessel_k(1, z)+ boost::math::cyl_bessel_k(0, y) *
                boost::math::cyl_bessel_i(1, z)) / (boost::math::cyl_bessel_i(1, z) *
                boost::math::cyl_bessel_k(1, y) - boost::math::cyl_bessel_k(1, z) *
                                                    boost::math::cyl_bessel_i(1, y)));

        f = pow((((Sig_aM * V_M)/(Sig_aF * V_F)) * F + E), (-1.));

        core.region[i].DA = (Sig_aF*V_F - f*Sig_aF*V_F)/(f*Sig_aM*V_M);
    }
}

// Determines the operating mode of reactor and burns fuel accordingly
void BurnFuel(ReactorLiteInfo &core) {

    if(core.target_CR_ < 0) {
        // Reactor is in stop at k=1 mode
        CriticalityBurn(core);
    } else {
        ///TODO CR burn
    }

}
*/

// Burns the fuel by advancing fluence based only on criticality
void CriticalityBurn(ReactorXInfo &core) {
    float kcore = 1.5;
    float kcore_prev;
    float abs_flux = core.base_flux_;
    unsigned const int batches = core.batches;

    core.base_flux_ = AbsFluxCalc(core, abs_flux, batches);

    while(kcore > 1) {
        kcore_prev = kcore; // Save previous k for final interpolation
        //FluxCalc(core); // Update relative flux of regions
        EqPowPhi(core);

        for(int i = 0; i < core.type[0].batch.size(); i++) {
            ///TODO this will change!!! just to see if it works, for now
            core.type[1].batch[i].rflux_ = core.type[0].batch[i].rflux_ * 0.95;
        }

        // Calculate DA
        //    DACalc(core);

        abs_flux = AbsFluxCalc(core, abs_flux, batches);
        core.base_flux_ = abs_flux;
        //std::cout << "Absolute flux: " << abs_flux << " k: " << kcore << std::endl;

        // Update fluences
        for(int type_i = 0; type_i < core.type.size(); type_i++) {
            for(int batch_i = 0; batch_i < core.type[type_i].batch.size(); batch_i++) {
                core.type[type_i].batch[batch_i].fluence_ += core.type[type_i].batch[batch_i].rflux_
                                                             * core.fluence_timestep_ * abs_flux;
            }
        }

        // Calculate k
        kcore = kCalc(core);
        //std::cout << "k: " << kcore << " BU: " << core.region[0].CalcBU() << " CR: " << CoreCRCalc(core) << std::endl;
    }

    // Update base_flux_ for next time
    core.base_flux_ = abs_flux;

    // Find the discharge fluences
    for(int type_i = 0; type_i < core.type.size(); type_i++) {
        for(int batch_i = 0; batch_i < core.type[type_i].batch.size(); batch_i++) {
        ///The subtraction here is meh
            core.type[type_i].batch[batch_i].fluence_ = Interpolate(
                                                                core.type[type_i].batch[batch_i].fluence_ -
                                                                core.type[type_i].batch[batch_i].rflux_
                                                                * core.fluence_timestep_ * abs_flux,
                                                                core.type[type_i].batch[batch_i].fluence_,
                                                                kcore_prev, kcore, 1);
        }
    }
}
/*
// Determines the flux calculation method and calls flux function accordingly
void FluxCalc(ReactorXInfo &core) {
    unsigned int mode = core.flux_mode_;

    // Override to mode=1 if any region exceeds library limit
    for(unsigned int reg_i = 0; reg_i < core.region.size(); reg_i++) {
        if(core.region[reg_i].fluence_ > core.region[reg_i].iso.fluence.back()) {
            mode = 1;
        }
    }

    if(mode == 1) {
        // Simplest mode, all fluxes are 1
        for(unsigned int reg_i = 0; reg_i < core.region.size(); reg_i++){
            core.region[reg_i].rflux_ = 1;
            return;
        }

    }
    else if(mode == 0) {EqPowPhi(core);}
    else if(mode == 2) {InvProdPhi(core);}
    else if(mode == 3) {SpatialPhi(core);}
    else {
        std::cout << "  Error in flux mode input for ReactorLite" << std::endl;
    }
}
*/

// Calculates the k-value of the core
float kCalc(ReactorXInfo &core) {
    unsigned const int batches = core.batches;
    float prod_tot = 0;
    float dest_tot = 0;

    for(int type_i = 0; type_i < core.type.size(); type_i++) {
        for(int batch_i = 0; batch_i < core.type[type_i].batch.size(); batch_i++) {
            prod_tot += core.type[type_i].batch[batch_i].CalcProd()
                        * core.type[type_i].batch[batch_i].rflux_ * core.type[type_i].mass;

            dest_tot += core.type[type_i].batch[batch_i].CalcDest()
                        * core.type[type_i].batch[batch_i].rflux_ * core.type[type_i].mass;
        }
    }


    if(dest_tot <= 0) {return 0;}
    return core.pnl * prod_tot / dest_tot;
}


// Calculates relative fluxes based on the equal power sharing assumption (1)
void EqPowPhi(ReactorXInfo &core) {
    /// works for two fuel types with equal number of batches
    ///   -> its assuming each assembly has two types of fuel
    ///   - saves the calculated fluxes to fuel type [0]
    // Operates on: Regions of equal power have equal burnup
    // this function updates relative fluxes instead of directly calculating burnup

    float max_fluence = core.fluence_timestep_ * core.base_flux_;
    float bu_old, bu_next, delta_bu, batch_bu;
    float batch_fluence;
    unsigned const int N = core.type[0].batch.size();
    float max_flux = -1;
    float min_flux = 10;
    unsigned int jk;
    core.type[0].batch[0].rflux_ = 1;
    core.type[1].batch[0].rflux_ = 1;

    // Find the current burnup of the oldest assembly
    bu_old = core.AssemblyBU(0);

    // Find the burnup for next step
    // This assumes oldest batch will have the least burnup for a given change in fluence
    bu_next = core.AssemblyBU(0, core.type[0].batch[0].fluence_ + max_fluence);
    delta_bu = bu_next - bu_old;

    for(int unsigned i = 0; i < N; i++) {
        batch_bu = core.AssemblyBU(i);

        // find the discrete points before and after batch bu
        for(jk = 0; core.type[0].batch[i].iso.BU[jk] < batch_bu + delta_bu; jk++){
        }

        batch_fluence = Interpolate(core.type[0].batch[i].iso.fluence[jk-1],
                                    core.type[0].batch[i].iso.fluence[jk],
                                    core.AssemblyBU(i, core.type[0].batch[i].iso.fluence[jk-1]),
                                    core.AssemblyBU(i, core.type[0].batch[i].iso.fluence[jk]),
                                    batch_bu + delta_bu);

        core.type[0].batch[i].rflux_ = (batch_fluence - core.type[0].batch[i].fluence_)/(max_fluence);

        if(core.type[0].batch[i].rflux_ > max_flux){max_flux = core.type[0].batch[i].rflux_;}
        if(core.type[0].batch[i].rflux_ < 0){core.type[0].batch[i].rflux_ = 0;}
        if(core.type[0].batch[i].rflux_ < min_flux){min_flux = core.type[0].batch[i].rflux_;}
    }

    for(int unsigned i = 0; i < N; i++) {
        if(core.type[0].batch[i].rflux_ == 0) {
            core.type[0].batch[i].rflux_ = min_flux/max_flux;
        } else {
            core.type[0].batch[i].rflux_ = core.type[0].batch[i].rflux_ / max_flux;
        }
    }

    //cout << "EqPow fluxes: " << endl;
    for(int unsigned i = 0; i < N; i++) {
        //cout << "               " << core.type[0].batch[i].rflux_ << endl;
    }

}
/*
// Calculates relative fluxes based on the inverse of neutron production assumption (2)
void InvProdPhi(ReactorLiteInfo &core) {
// Updates the rflux of each region in core.region
// Assumes the flux of a region is proportional to the inverse neutron prod rate

    double maxphi = 0;

    // Find the inverse of neutron production
    for(unsigned int reg_i = 0; reg_i < core.region.size(); reg_i++){

        core.region[reg_i].rflux_ = 1. / core.region[reg_i].CalcProd();

        if(maxphi < core.region[reg_i].rflux_){
            maxphi = core.region[reg_i].rflux_;
        }
    }

    // Normalize all the flux values
    for(unsigned int reg_i = 0; reg_i < core.region.size(); reg_i++){
        core.region[reg_i].rflux_ /= maxphi;
    }
}

// Calculates the CR for reg_i
float RegionCRCalc(ReactorLiteInfo &core, unsigned const int reg_i) {
    // reg_i is the region number, starting from zero

    float FP = 0, FP0 = 0, FP1 = 0;
    float fissile = 0, fissile0 = 0, fissile1 = 0;
    float ini_fissile = 0;
    float CR;
    unsigned int ii, ZZ;

    const unsigned int CR_upper = 160, CR_lower = 70;

    // Find the discrete point index for region fluence
    if(core.region[reg_i].fluence_ > core.region[reg_i].iso.fluence.back()) {
        ii = core.region[reg_i].iso.fluence.size()-1;
    } else {
        for(ii = 1; core.region[reg_i].iso.fluence[ii] < core.region[reg_i].fluence_; ii++){}
    }

    for(int iso_j = 0; iso_j < core.region[reg_i].iso.iso_vector.size(); iso_j++) {
        // Convert name to mass number
        ZZ = core.region[reg_i].iso.iso_vector[iso_j].name;
        ZZ = ZZ % 10000;
        ZZ /= 10;

        // Add up the FP
        if(ZZ < CR_upper && ZZ > CR_lower) {
            // Interpolation will be done at the end
            FP0 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii-1];
            FP1 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii];
        }

        // Add up fissiles
        for(int fis = 0; fis < core.CR_fissile_.size(); fis++){
            if(core.region[reg_i].iso.iso_vector[iso_j].name == core.CR_fissile_[fis]){
                fissile0 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii-1];
                fissile1 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii];

                ini_fissile += core.region[reg_i].iso.iso_vector[iso_j].mass[0];
            }
        }
    }

    // recycling variable FP0 here to check greater than zero
    FP0 = Interpolate(FP0, FP1, core.region[reg_i].iso.fluence[ii-1],
                      core.region[reg_i].iso.fluence[ii], core.region[reg_i].fluence_);
    if(FP0 > 0) {FP += FP0;}

    fissile = Interpolate(fissile0, fissile1, core.region[reg_i].iso.fluence[ii-1],
                          core.region[reg_i].iso.fluence[ii], core.region[reg_i].fluence_);

    if(FP > 0){
        CR = (FP+fissile-ini_fissile)/FP;
    } else {
        CR  = 0;
    }

    if(CR < 0){CR = 0;}

    return CR;
}

// Calculates the CR for the core
float CoreCRCalc(ReactorLiteInfo &core) {
    float FP = 0, FP0 = 0, FP1 = 0;
    float fissile = 0, fissile0 = 0, fissile1 = 0;
    float ini_fissile = 0;
    float CR;
    unsigned int ii, ZZ;

    const unsigned int CR_upper = 160, CR_lower = 70;

    for(unsigned int reg_i = 0; reg_i < core.region.size(); reg_i++) {
        // Find the discrete point index for region fluence
        if(core.region[reg_i].fluence_ > core.region[reg_i].iso.fluence.back()) {
            ii = core.region[reg_i].iso.fluence.size()-1;
        } else {
            for(ii = 1; core.region[reg_i].iso.fluence[ii] < core.region[reg_i].fluence_; ii++){}
        }

        for(int iso_j = 0; iso_j < core.region[reg_i].iso.iso_vector.size(); iso_j++) {
            // Convert name to mass number
            ZZ = core.region[reg_i].iso.iso_vector[iso_j].name;
            ZZ = ZZ % 10000;
            ZZ /= 10;

            // Add up the FP
            if(ZZ < CR_upper && ZZ > CR_lower) {
                // Interpolation will be done at the end
                FP0 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii-1];
                FP1 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii];
            }

            // Add up fissiles
            for(int fis = 0; fis < core.CR_fissile_.size(); fis++){
                if(core.region[reg_i].iso.iso_vector[iso_j].name == core.CR_fissile_[fis]){
                    fissile0 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii-1];
                    fissile1 += core.region[reg_i].iso.iso_vector[iso_j].mass[ii];

                    ini_fissile += core.region[reg_i].iso.iso_vector[iso_j].mass[0];
                }
            }
        }

        // recycling variable FP0 here to check greater than zero
        FP0 = Interpolate(FP0, FP1, core.region[reg_i].iso.fluence[ii-1],
                          core.region[reg_i].iso.fluence[ii], core.region[reg_i].fluence_);
        if(FP0 > 0) {FP += FP0;}

        fissile += Interpolate(fissile0, fissile1, core.region[reg_i].iso.fluence[ii-1],
                              core.region[reg_i].iso.fluence[ii], core.region[reg_i].fluence_);

        FP0 = 0;
        FP1 = 0;
        fissile0 = 0;
        fissile1 = 0;
    }

    if(FP > 0){
        CR = (FP+fissile-ini_fissile)/FP;
    } else {
        CR  = 0;
    }

    if(CR < 0){CR = 0;}

    return CR;
}
*/
// Determines the absolute flux (correct flux units) for the timestep
float AbsFluxCalc(ReactorXInfo &core, float abs_flux, int regions) {
    //std::cout << "ABSFLUXCALC BEGIN" << std::endl;

    const float power = core.thermal_pow_;          // [MWth]
    const float mass = core.core_mass_;             // [kg]
    const float timestep = core.fluence_timestep_;  // [day]
    const float BU_prev = core.CalcBU();            // [MWd/kgIHM]
    float BU_next, delta_BU;
    float fluence;
    float step_power1 = 0, step_power2 = 0;
    float abs_flux1 = abs_flux, abs_flux2 = abs_flux*2, temp_flux;
    unsigned int times = 0;

    delta_BU = core.CalcBU(abs_flux1) - BU_prev;
    step_power1 = delta_BU * mass / timestep;

    //cout << "delta BU: " << delta_BU << " BU1: " << core.CalcBU(abs_flux1) << " BU2: " << BU_prev << endl;
    //cout << "abs flux1: " << abs_flux1 << endl;

    // Quickly scale up flux if abs_flux was too far off
    while(step_power1 < power && times < 20) {
        abs_flux1 *= 1.5;

        delta_BU = core.CalcBU(abs_flux1) - BU_prev;
        step_power1 = delta_BU * mass / timestep;

        times++;
    }

    times = 0;
    abs_flux2 = abs_flux1 * 0.9;


    while(times < 20) {
        // Determine the timestep power implied by the current abs_flux
    //cout << "abs flux2: " << abs_flux2 << endl;
        delta_BU = core.CalcBU(abs_flux2) - BU_prev;
        //cout << " delta BU: " << delta_BU << " BUprev: " << BU_prev << endl;
        step_power2 = delta_BU * mass / timestep;   // [MWd/kgIHM] * [kgIHM] / [day]

        if(std::abs(step_power2 - power)/power < core.abs_flux_tol_) {
            return abs_flux2;
        }

        //cout << "f2:" << abs_flux2 << " steppow2:" << step_power2 << " f1:" << abs_flux1 << " steppow1:" << step_power1 << endl;
        temp_flux = abs_flux2 + (power - step_power2)*(abs_flux2-abs_flux1)/(step_power2-step_power1);
        //cout << "temp flux:" << temp_flux << endl;
        if(temp_flux < 0) {temp_flux = 0;}

        step_power1 = step_power2;
        abs_flux1 = abs_flux2;
        abs_flux2 = temp_flux;

        times++;
    }

    std::cout << "Warning! Absolute flux calculation reached max iterations!" << std::endl;
    return abs_flux2;
}




/*
// Returns the steady state fluence of the core by refueling with IsoInfo iso
float SteadyStateFluence(ReactorLiteInfo core, IsoInfo iso) {
    // Note that the ReactorLiteInfo object is a copy here
    // Function needs the original ReactorLiteInfo object from tick() as first argument
    // Assumes region[0] is oldest

    float fluence1 = 0, fluence2 = 0;
    const unsigned int region = core.region.size();

    // Assign iso to all regions
    for(int reg_i = 0; reg_i < region; reg_i++) {
        core.region[reg_i].iso = iso;

        ///TODO think about this more, ran in to this optimization frontier before
        core.region[reg_i].fluence_ = 0;
    }

    // Burn fuel once before entering loop
    BurnFuel(core);
    fluence2 = core.region[0].fluence_;

    while(std::abs(fluence1 - fluence2)/fluence2 > core.SS_tol_) {
        fluence1 = fluence2;

        // "reload" fuel by adjusting fluences
        for(int reg_i = 0; reg_i < region-1; reg_i++) {
            core.region[reg_i].fluence_ = core.region[reg_i+1].fluence_;
        }
        core.region[region-1].fluence_ = 0;

        BurnFuel(core);
        fluence2 = core.region[0].fluence_;
    }

    return fluence2;
}
*/





































