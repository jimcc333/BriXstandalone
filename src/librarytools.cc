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
            daughter.name = pyne::nucname::zzaaam(buffer);
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

// Burns the fuel by advancing fluence based only on criticality
void CriticalityBurn(ReactorLiteInfo &core) {
// yes this IS old burnupcalc
    float kcore = 1.5;
    float kcore_prev;
    unsigned const int regions = core.region.size();
    float abs_flux = core.base_flux_;

    while(kcore > 1) {
        kcore_prev = kcore; // Save previous k for final interpolation
        FluxCalc(core); // Update relative flux of regions

        // Calculate DA
        if(core.DA_mode_ == 1) {
            DACalc(core);
        }

        abs_flux = AbsFluxCalc(core, abs_flux);
        //std::cout << "Absolute flux: " << abs_flux << std::endl;

        // Update fluences
        for(unsigned int reg_i = 0; reg_i < regions; reg_i++) {
            core.region[reg_i].fluence_ += core.region[reg_i].rflux_
                    * core.fluence_timestep_ * abs_flux;
        }

        // Calculate k
        kcore = kCalc(core);
        //std::cout << "k: " << kcore << " BU: " << core.region[0].CalcBU() << " CR: " << CoreCRCalc(core) << std::endl;
    }

    // Update base_flux_ for next time
    core.base_flux_ = abs_flux;

    // Find the discharge fluences
    for(unsigned int reg_i = 0; reg_i < regions; reg_i++) {
        ///The subtraction here is meh
        core.region[reg_i].fluence_ = Interpolate((core.region[reg_i].fluence_ -
                   core.region[reg_i].rflux_ * core.fluence_timestep_ * abs_flux),
                   core.region[reg_i].fluence_, kcore_prev, kcore, 1);
    }

}

// Determines the flux calculation method and calls flux function accordingly
void FluxCalc(ReactorLiteInfo &core) {
    unsigned int mode = core.flux_mode_;

    // Override to mode=zero if any region exceeds library limit
    for(unsigned int reg_i = 0; reg_i < core.region.size(); reg_i++) {
        if(core.region[reg_i].fluence_ > core.region[reg_i].iso.fluence.back()) {
            mode = 0;
        }
    }

    if(mode == 0) {
        // Simplest mode, all fluxes are 1
        for(unsigned int reg_i = 0; reg_i < core.region.size(); reg_i++){
            core.region[reg_i].rflux_ = 1;
            return;
        }

    }
    else if(mode == 1) {EqPowPhi(core);}
    else if(mode == 2) {InvProdPhi(core);}
    else if(mode == 3) {SpatialPhi(core);}
    else {
        std::cout << "  Error in flux mode input for ReactorLite" << std::endl;
    }
}

// Calculates the k-value of the core
float kCalc(ReactorLiteInfo &core) {
    unsigned const int regions = core.region.size();
    float prod_tot = 0;
    float dest_tot = 0;

    for(unsigned int reg_i = 0; reg_i < regions; reg_i++) {
        prod_tot += ( (core.region[reg_i].CalcProd() + core.struct_prod_ * core.region[reg_i].DA)
                    * core.region[reg_i].rflux_ * core.region[reg_i].mass_);

        dest_tot += ( (core.region[reg_i].CalcDest() +
                       core.struct_dest_ * core.region[reg_i].DA)
                    * core.region[reg_i].rflux_) * core.region[reg_i].mass_;
    }

    if(dest_tot <= 0) {return 0;}
    return core.pnl * prod_tot / dest_tot;
}

// Calculates relative fluxes based on the equal power sharing assumption (1)
void EqPowPhi(ReactorLiteInfo &core) {
    // Operates on: Regions of equal power have equal burnup
    // this function updates relative fluxes instead of directly calculating burnup

    float max_fluence = core.fluence_timestep_ * core.base_flux_;
    float bu_old, bu_next, delta_bu, batch_bu;
    float batch_fluence;
    unsigned const int N = core.region.size();
    float max_flux = -1;
    float min_flux = 10;
    unsigned int jk;
    core.region[0].rflux_ = 1;

    // Find the current burnup of the oldest batch
    bu_old = core.region[0].CalcBU();

    // Find the burnup for next step
    // This assumes oldest batch will have the least burnup for a given change in fluence
    bu_next = core.region[0].CalcBU(core.region[0].fluence_ + max_fluence);
    delta_bu = bu_next - bu_old;

    for(int unsigned i = 0; i < N; i++) {
        batch_bu = core.region[i].CalcBU();
        // find the discrete points before and after batch bu
        for(jk = 0; core.region[i].iso.BU[jk] < batch_bu + delta_bu; jk++){
        }

        batch_fluence = Interpolate(core.region[i].iso.fluence[jk-1], core.region[i].iso.fluence[jk],
                                    core.region[i].iso.BU[jk-1], core.region[i].iso.BU[jk], batch_bu + delta_bu);

        core.region[i].rflux_ = (batch_fluence - core.region[i].fluence_)/(max_fluence);
        if(core.region[i].rflux_ > max_flux){max_flux = core.region[i].rflux_;}
        if(core.region[i].rflux_ < 0){core.region[i].rflux_ = 0;}
        if(core.region[i].rflux_ < min_flux){min_flux = core.region[i].rflux_;}
    }

    for(int unsigned i = 0; i < N; i++) {
        if(core.region[i].rflux_ == 0) {
            core.region[i].rflux_ = min_flux/max_flux;
        } else {
            core.region[i].rflux_ = core.region[i].rflux_ / max_flux;
        }
    }
}

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

// Determines the absolute flux (correct flux units) for the timestep
float AbsFluxCalc(ReactorLiteInfo &core, float abs_flux) {
    //std::cout << "ABSFLUXCALC BEGIN" << std::endl;
    const int regions = core.regions_;
    const float power = core.thermal_pow_;          // [MWth]
    const float mass = core.core_mass_;
    const float timestep = core.fluence_timestep_/86400.;  // [day]
    const float BU_prev = core.CalcBU();            // [MWd/kgIHM]
    float BU_next, delta_BU;
    float fluence;
    float step_power1 = 0, step_power2 = 0;
    float abs_flux1 = abs_flux, abs_flux2 = abs_flux*2, temp_flux;
    unsigned int times = 0;

    delta_BU = core.CalcBU(abs_flux1) - BU_prev;
    step_power1 = delta_BU * mass / timestep;

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
        delta_BU = core.CalcBU(abs_flux2) - BU_prev;
        step_power2 = delta_BU * mass / timestep;   // [MWd/kgIHM] * [kgIHM] / [day]

        if(std::abs(step_power2 - power)/power < core.abs_flux_tol_) {
            return abs_flux2;
        }

        temp_flux = abs_flux2 + (power - step_power2)*(abs_flux2-abs_flux1)/(step_power2-step_power1);
        if(temp_flux < 0) {temp_flux = 0;}

        step_power1 = step_power2;
        abs_flux1 = abs_flux2;
        abs_flux2 = temp_flux;

        times++;
    }

    std::cout << "Warning! Absolute flux calculation reached max iterations!" << std::endl;
    return abs_flux2;
}

// Performs the spatial flux calculation
void SpatialPhi(ReactorLiteInfo &core) {
// Assumes an outermost water region

    unsigned const int region = core.region.size();
    float delta = core.spatial_.delta;
    float R[region+1];          //radial thickness of each region
    int N[region+1];            //number of mesh points in each region
    int NC[region+1];           //cumulative N
    int NTotal;                 //total number of mesh points
    int iter;
    float dd2[region+1];        // D/delta^2
    float Sigma_a[region+1];    //mac. abs. cs of each region
    float NuSigma_f[region+1];  //nu sigma f of each region
    float Sigma_tr[region+1];   //mac. transport cs of each region
    float D[region+1];          //diff coef. for each region
    float LSquared[region+1];
    float k = 1;
    float k_prev = 0.9;
    float prod, prod_prev;
    float flux[region+1], maxflux=0;
    float sum = 0;

    // Set the radial thickness of each region
    R[0] = std::sqrt(core.spatial_.fuel_area/region/3.141592);

    for(unsigned int reg_i = 1; reg_i < region; reg_i++) {
        R[reg_i] = sqrt(core.spatial_.fuel_area/region/3.141592*(reg_i+1));
    }
    R[region] = R[region-1] + core.spatial_.spatial_mod_thickness;

    // Assign fuel cross sections
    for(int reg_i = 0; reg_i < region; reg_i++) {

        NuSigma_f[reg_i] = core.region[reg_i].CalcNuSigf();
        Sigma_a[reg_i] = core.region[reg_i].CalcSiga();

    /*    // comment out!
        Sigma_a[0] = 0.0230;
        Sigma_a[1] = 0.0246;
        Sigma_a[2] = 0.0324;
        NuSigma_f[0] = 0.0184;
        NuSigma_f[1] = 0.0217;
        NuSigma_f[2] = 0.0382;
        /// till here!  */

        Sigma_tr[reg_i] = core.spatial_.spatial_fuel_Sig_tr;
        D[reg_i] = 1/(Sigma_tr[reg_i]*3.);
        LSquared[reg_i] = D[reg_i]/Sigma_a[reg_i];
    }

    // Assign moderator cross sections
    Sigma_a[region] = core.spatial_.spatial_mod_Sig_a;
    Sigma_tr[region] = core.spatial_.spatial_mod_Sig_tr;
    D[region] = 1/(Sigma_tr[region]*3.);
    LSquared[region] = D[region]/Sigma_a[region];
    NuSigma_f[region] = core.spatial_.spatial_mod_Sig_f;

    // Populate dd2
    for(int reg_i = 0; reg_i < region+1; reg_i++) {
        dd2[reg_i] = D[reg_i]/delta/delta;
    }

    if((R[region-1]-R[region-2])/delta < 2.99) {
        std::cout << "  Warning, too few discrete points in spatial flux calc." << std::endl;
        unsigned int last = region;
        if(region == 1){last = 2;}
        delta = (R[last-1]-R[last-2])/3;
        core.spatial_.delta = delta;
        std::cout << "    Delta changed to " << delta << " [cm]." << std::endl;
    }

    // Populate N, number of mesh points in each region
    N[0] = round(R[0]/delta);
    NC[0] = N[0];
    NTotal = N[0];

    for(unsigned int reg_i = 1; reg_i < region+1; reg_i++) {
        N[reg_i] = std::round((R[reg_i] - R[reg_i-1]) / delta);
        if(N[reg_i] < 3) {
            std::cout << "  Warning! Region " << reg_i+1 << " has too few discrete points ("
                    << N[reg_i] << ").  - Increase fuel_area or mod_thickness." << std::endl;
            EqPowPhi(core);
            return;
        }
        NC[reg_i] = NC[reg_i-1] + N[reg_i];
        NTotal += N[reg_i];
    }
    NC[region] += 1;
    NTotal += 1;
/*
    Eigen::MatrixXf A(NTotal, NTotal);
    Eigen::MatrixXf F(NTotal, 1);
    Eigen::MatrixXf phi(NTotal, 1);
    Eigen::MatrixXf phi_prev(NTotal, 1);
    Eigen::MatrixXf S(NTotal, 1);
    Eigen::MatrixXf S_prev(NTotal, 1);

    A.setZero();
    F.setZero();

    int jprev = 0;
    int j;

    int r = 0; //region index
    for(int i = 1; i < NTotal-1; i++) {
        A(i, i-1) = (-1.)*dd2[r]*(2*i-1)/(2*i);
        A(i,i) = dd2[r]*2. + Sigma_a[r];
        A(i, i+1) = (-1.)*dd2[r]*(2*i+1)/(2*i);
        if(i == NC[r]){
            r += 1;
        }
        if(core.region[r].rflux_ < 1.1 && core.region[r].rflux_ > 0) {
            phi_prev(i) = 1;//core.region[r].rflux; //uses last runs results if available
        } else {
            phi_prev(i) = 1;
        }

    }
    A(0,0) = 1;
    A(0,1) = -1;
    A(NTotal-1,NTotal-1) = 1;
    phi_prev(NTotal-1) = 0;

    // Boundary conditions
    for(r = 0; r < region; r++) {
        A(NC[r],NC[r]-1) = D[r];
        A(NC[r],NC[r]) = -D[r]-D[r+1];
        A(NC[r],NC[r]+1) = D[r+1];
    }

    r = 0;
    for(int i = 1; i < NTotal; i++) {
        if(i != NC[r]) {
            F(i) = NuSigma_f[r];
        }
        if(i == NC[r]) {r += 1;}
        S_prev(i) = F(i)*phi_prev(i);
    }

    for(iter = 0; iter < 100; iter++) {
        phi = A.colPivHouseholderQr().solve(S_prev)/k_prev;

        if(!phi.allFinite()) {phi = phi_prev;}

        S = F.array() * phi.array();

        prod = (0.25)*3.141592*NuSigma_f[0]*phi(0)*delta*delta;
        prod_prev = (0.25)*3.141592*NuSigma_f[0]*phi_prev(0)*delta*delta;
        r = 0;

        for(int i = 0; i < NTotal; i++) {
            if(i == NC[r]){
                prod += 3.141592*NuSigma_f[r]*phi(i)*(i-0.25)*delta*delta;
                prod_prev += 3.141592*NuSigma_f[r]*phi_prev(i)*(i-0.25)*delta*delta;
                r += 1;
                prod += 3.141592*NuSigma_f[r]*phi(i)*(i-0.25)*delta*delta;
                prod_prev += 3.141592*NuSigma_f[r]*phi_prev(i)*(i-0.25)*delta*delta;
            } else {
                prod += 2.*3.141592*NuSigma_f[r]*phi(i)*i*delta*delta;
                prod_prev += 2.*3.141592*NuSigma_f[r]*phi_prev(i)*i*delta*delta;
            }
        }
        ///TODO check this
        if(abs((k_prev-k)/k) < 0.001 && iter > 3) {break;}
        //cout << "prod: " << prod << "  prod_prev: " << prod_prev << "  k: " << prod/prod_prev << endl;
        k = prod/prod_prev*k_prev;
        phi_prev = phi;
        k_prev = k;
        S_prev = S;
    }

    // Find area weighted average phi per batch
    r = 0;
    flux[0] = 0;
    for(int i = 0; i < NTotal; i++) {
        if(phi(i) < 0) {phi(i) = 0;}

        flux[r] += phi(i)*(2*(i+1)-1);
        sum += (2*(i+1)-1);

        // uses some trickery to switch between regions
        if(i == NC[r] || i == NTotal-1) {
            // divide by total area
            flux[r] /= sum;

            // reset sum and flux of next region
            sum = 0;
            r += 1;
            flux[r] = 0;
        }
    }
    for(r = 0; r < region+1; r++) {
        if(flux[r] > maxflux){
            maxflux = flux[r];
        }
    }

    //cout << " Iterations:" << iter << " k:" << k << endl;
    //cout << "--- A ---" << endl << A << endl << " --- ----" << endl;
    //cout << "---phi---" << endl<< phi << endl << "--------" << endl;

    //find area weighted average phi per batch
    r = 0;
    flux[0] = 0;
    for(int i = 0; i < NTotal; i++) {
        flux[r] += phi(i)*(2*(i+1)-1);
        sum += (2*(i+1)-1);

        if(i == NC[r] || i == NTotal-1) {
            flux[r] /= sum;
            sum = 0;
            r += 1;
            flux[r] = 0;
        }
    }
    for(r = 0; r < region+1; r++) {
        if(flux[r] > maxflux){
            maxflux = flux[r];
        }
    }

    // Normalize the fluxes
    for(r = 0; r < region+1; r++) {
        flux[r] /= maxflux;
    }

    for(int reg_i = 0; reg_i < region; reg_i++) {
        core.region[reg_i].rflux_ = flux[reg_i];
    }*/
}

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






































