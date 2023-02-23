# Import python libs:
import numpy as np
import os

# Import auxiliary libraries:
import juliet
import exoctk

# Import libraries for this script:
import utils

# Exposure time of 2-min cadence data, in days:
exp_time = (2./60.)/24.

# Load target list. This will only load the planet name and TIC IDs:
target_list = utils.read_data('TESS_ID.txt')
density_list = utils.read_data_density('TESS_Densities.txt')

# Define fitting method. If blank (""), fit GP and transit lightcurve together. 
# If set to "fit_out", fit out-of-transit lightcurve first with a GP, then use posteriors of that fit to 
# fit the in-transit lightcurve. The factor variable defines what is "in-transit" in units of the transit duration:
method = "fit_out"
factor = 5.

# Run multi-sector fit? This will run a global joint fit of all sectors:
run_multisector = True

# Do we want to run catwoman fits?
fit_catwoman = True

# Number of threads for multi-sector fits:
#nthreads = 4
nthreads = 10

# Iterate through all planets:
for planet in target_list.keys():

    print('Working on ',planet)
    # Load all data for this particular planet; if data exists, go forward. If not, skip 
    # target:
    with open('completed_planets.txt') as f:
        if str(planet) in f.read():
            print("Already complete:",planet)
            continue
    with open('problem_planets.txt') as f:
        if str(planet) in f.read():
            print("Already complete:",planet)
            continue
    try:

        # If Kepler planet, don't analyze it. Too much blending:
        if 'Kepler' not in planet:
        #if 'WASP' in planet:
            planet_data, url = exoctk.utils.get_target_data(planet)

            # Extract useful data:
            tdur = planet_data['transit_duration']
            tdepth = planet_data['transit_depth']
            period = planet_data['orbital_period']
            period_err = (planet_data['orbital_period_upper'] + planet_data['orbital_period_lower'])*0.5
            t0 = planet_data['transit_time'] + 2400000.5
            t0_err = (planet_data['transit_time_upper'] + planet_data['transit_time_lower'])*0.5
            
            density=density_list[planet]['density']
            dens_err_up=density_list[planet]['dens_err_up']
            dens_err_low=density_list[planet]['dens_err_low']
            
            if 'N/A' not in density and 'N/A' not in dens_err_low and 'N/A' not in dens_err_up:
                density=float(density)*1000
                dens_err_low=float(dens_err_low)*1000
                dens_err_up=float(dens_err_up)*1000
            
            if (type(density) is float) and (type(dens_err_up) and float) and (type(dens_err_low) is float):
                print('Fit for {} will include stellar density as a prior'.format(planet))
                rho = density
                print('rho is {}'.format(rho))
                rho_sig = np.mean([dens_err_up,dens_err_low])
                print('rho_sigma is {}'.format(rho_sig))
            else:
                print('Fit for {} will not include stellar density as a prior'.format(planet))
                rho = 0
                rho_sig = 0

            # If data is not float (e.g., empty values), reject system:
            if (type(tdur) is float) and (type(tdepth) is float) and (type(period) is float) and (type(period_err) is float) and \
               (type(t0) is float) and (type(t0_err) is float):
                has_data = True
            else:
                print('Something is wrong with ',planet,' data. Skipping.')
                file2 = open("problem_planets.txt", "a")  # append mode
                string="{} \n".format(planet)
                file2.write(string)
                file2.close()
                has_data = False

            # Now check eccentricity and omega. If no data, set to 0 and 90:
            ecc, omega = planet_data['eccentricity'], planet_data['omega']
            if (type(ecc) is not float or type(omega) is not float):
                ecc = 0.
                omega = 90.
        else:
            has_data = False

    except:
        print('No planetary data for ',planet)
        file2 = open("problem_planets.txt", "a")  # append mode
        string="{} \n".format(planet)
        file2.write(string)
        file2.close()
        has_data = False

    # If it has data, we move ahead with the analysis:
    if has_data:# and not os.path.exists(planet):

        # First, load data for each sector:
        t, f, ferr = juliet.utils.get_all_TESS_data('TIC ' + target_list[planet]['ticid'])

        # If it has planetary data, go through each sector. First, estimate the transit depth precision 
        # (assuming a box-shaped transit) we would obtain with this 2-min TESS data. If this gives rise 
        # to a 5-sigma "detection" of the depth using all the (phased) transits, we analyze it:
        nsectors = len(list(t.keys()))
        good_sectors = []
        for sector in t.keys():

            # Estimate number of transits we should expect in this dataset:
            total_time = np.max(t[sector]) - np.min(t[sector])
            Ntransits = int(total_time/period) * nsectors

            # Estimate number of datapoints in-transit in the phased lightcurve:
            Nin = int(tdur/exp_time) * Ntransits

            # Estimate transit depth precision:
            I = 2. * Nin
            sigma = np.median(ferr[sector])
            sigma_depth = (2. * sigma)/np.sqrt(I)

            # If initial SNR estimate is larger than 5-sigma, perform the fit:
            if tdepth/sigma_depth > 5:

                print('\t >> Performing fit for '+sector+'; expected depth precision FOR ALL SECTORS: ',sigma_depth*1e6,' giving SNR:', tdepth/sigma_depth)

                if not os.path.exists(planet):
                    os.mkdir(planet)

                full_path = planet+'/'+sector

                utils.fit(t[sector], f[sector], ferr[sector], sector, period, period_err, t0, t0_err, ecc, omega, GPmodel = 'ExpMatern', outpath = full_path, \
                          method = method, in_transit_length = factor*tdur)
                utils.fit(t[sector], f[sector], ferr[sector], sector, period, period_err, t0, t0_err, ecc, omega, GPmodel = 'QP', outpath = full_path, \
                          method = method, in_transit_length = factor*tdur)

                if fit_catwoman:

                    utils.fit(t[sector], f[sector], ferr[sector], sector, period, period_err, t0, t0_err, ecc, omega, GPmodel = 'ExpMatern', outpath = full_path, \
                              method = method, in_transit_length = factor*tdur, fit_catwoman = fit_catwoman)
                    utils.fit(t[sector], f[sector], ferr[sector], sector, period, period_err, t0, t0_err, ecc, omega, GPmodel = 'QP', outpath = full_path, \
                              method = method, in_transit_length = factor*tdur, fit_catwoman = fit_catwoman)

                good_sectors.append(sector)

            else:
                print('\t WARNING: ',sector, ' DOES NOT look good! Not doing the fit. Expected depth precision: ',sigma_depth*1e6,' giving SNR:', tdepth/sigma_depth)
                with open('problem_planets.txt') as f:
                    if str(planet) in f.read():
                        print("Already complete:",planet)
                        continue
                    else:
                        file2 = open("problem_planets.txt", "a")  # append mode
                        string="{} \n".format(planet)
                        file2.write(string)
                        file2.close()
            
        print('run_multisector:',run_multisector,' | Good sectors:', len(good_sectors))

        if run_multisector and len(good_sectors)>1:

            print('Running multisector fit for ExpMatern:')
            utils.multisector_fit(t, f, ferr, period, period_err, t0, t0_err, ecc, omega, rho, rho_sig,  GPmodel = 'ExpMatern', outpath = planet, good_sectors = good_sectors, method = method, \
                                 in_transit_length = factor*tdur, nthreads = nthreads)

            print('Running multisector fit for QP:')
            utils.multisector_fit(t, f, ferr, period, period_err, t0, t0_err, ecc, omega, rho, rho_sig, GPmodel = 'QP', outpath = planet, good_sectors = good_sectors, method = method, \
                                 in_transit_length = factor*tdur, nthreads = nthreads)
            if fit_catwoman:

                print('Repeating fits for catwoman model:')

                print('Running multisector fit for ExpMatern for catwoman:')
                utils.multisector_fit(t, f, ferr, period, period_err, t0, t0_err, ecc, omega, rho, rho_sig, GPmodel = 'ExpMatern', outpath = planet, good_sectors = good_sectors, method = method, \
                                     in_transit_length = factor*tdur, fit_catwoman = fit_catwoman, nthreads = nthreads)

                print('Running multisector fit for QP for catwoman:')
                utils.multisector_fit(t, f, ferr, period, period_err, t0, t0_err, ecc, omega, rho, rho_sig, GPmodel = 'QP', outpath = planet, good_sectors = good_sectors, method = method, \
                                     in_transit_length = factor*tdur, fit_catwoman = fit_catwoman, nthreads = nthreads)
        #complete_list.append(planet)
        # Append-adds at last
        file1 = open("completed_planets.txt", "a")  # append mode
        string="{} \n".format(planet)
        file1.write(string)
        file1.close()
        command="mv /home/lmiller/TSRC/final_run/tess-limb-darkening/{}/ /user/lmiller/TSRC_results/".format(str(planet))
        os.system(command)