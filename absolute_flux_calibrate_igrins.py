#Script to generate absolute flux calibrated IGRINS data and output
#it as a .flux_a0v.fits

import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table
from muler import igrins
from muler.utilities import photometry
import glob
import copy
from matplotlib.backends.backend_pdf import PdfPages  #For outputting a pdf with multiple pages (or one page)
import os
from astroquery.simbad import Simbad
Simbad.add_votable_fields('ids')
import astropy.units as u
from astropy.coordinates import SkyCoord

date = sys.argv[1]

#USER MODIFY THESE PATHS
plp_path = '/Users/kk25239/Desktop/plp' #Path to the to IGRINS-2 Pipeline
outdata_path = plp_path + '/outdata/' + date 
plots_path = outdata_path
recipe_log_path = plp_path +'/recipe_logs/'+date+'.recipes'


def standard_star_name_check(std_star_name, coords):
	std_star_name_0 = std_star_name
	#standard name/coordinates preflight check:
	#check to see if the standard name we grabbed from the header is simbad searchable
	query_name = Simbad.query_object(std_star_name)
	#the name is not simbad searchable
	if len(query_name) == 0 and coords != '':
		#now we query the coordinates
		#query_region doesnt work well when the string has coords in HMS format, but making the string into a SkyCoord object seems to fix things
		if ':' in coords:
			sky_coord = SkyCoord(coords, unit = (u.hourangle, u.deg), frame = 'icrs')
		else:
			sky_coord = SkyCoord(coords, unit = (u.deg, u.deg), frame = 'icrs')

		#query the object by coordintes
		query_coords = Simbad.query_region(sky_coord, radius = '20 arcsec')

		#if any object is returned from querying the coordinates we can try to verify the names
		if len(query_coords) != 0:
			#print(f'this many objects to search through: {len(query_coords)}')
			#sometimes the object we want is not the first object simbad returns
			for i in range(0, len(query_coords)):
				#grab the simbad main id for the object
				main_id = query_coords['main_id'][i]

				#remove weird characters from the main id including white spaces 
				main_id_tg = ''.join(main_id.replace('*', '').split())

				#check to see if the main_id for SIMBAD is a substring in the standard star name
				#if it is in the standard star name
				#print(main_id_tg, ''.join(std_star_name.split()))
				if main_id_tg in ''.join(std_star_name.split()):
					std_star_name = main_id
					#end the loop because we found our match
					print(f'\033[38;5;{39}mTHE INPUT NAME FOR THE STANDARD STAR IS {std_star_name_0}\033[0m')
					print(f'\033[38;5;{39}mTHE SIMBAD MAIN ID FOR THE STANDARD STAR IS {std_star_name}\033[0m')
					
					return std_star_name

				#if there is only one object found within 20 arcsec it is usually the standard we want
				elif len(query_coords) == 1:
					#print the standard we are assuming
					print(f"\033[38;5;{196}mONLY ONE OBJECT FOUND THROUGH COORDINATE SEARCH FOR STANDARD, ASSUMING STANDARD IS {query_coords['main_id'][0]}\033[0m")
					#grab the main ID for the object returned
					main_id = query_coords['main_id'][0]
					#set the standard searchable name
					std_star_name = main_id

					#print some additional info about the standard coords and simbad coords for reference
					print(f'\033[38;5;{39}mTHE INPUT NAME FOR THE STANDARD STAR IS {std_star_name_0}; COORDINATES {sky_coord.ra.deg} {sky_coord.dec.deg}\033[0m')
					print(f"\033[38;5;{39}mTHE SIMBAD MAIN ID FOR THE STANDARD STAR IS {std_star_name}; COORDINATES {query_coords['ra'][0]}; {query_coords['dec'][0]}\033[0m")

					return std_star_name
				#if the main id is not in the standard star name
				else:
					#grab all the object IDs from simbad
					ids = query_coords['ids'][i].split('|')

					#check to see if any of the object IDs from simbad are in the standard star name
					for ID in ids:
						#remove weird characters and spaces
						id_tg = ''.join(ID.replace('*', '').split())
						#if the id for the object is in the standard star name
						#print(id_tg, ''.join(std_star_name.split()))
						if id_tg in ''.join(std_star_name.split()):
							#set the standard star_name to the main id for the object
							std_star_name = main_id
							#in this case we matched the standard name to something in the list of object IDs in simbad, so we dont have to check other objects returned from the query
							print(f'\033[38;5;{39}mTHE SIMBAD SEARCHABLE MAIN ID FOR THE STANDARD STAR IS {std_star_name}\033[0m')
							
							return std_star_name

				#check to see if we matched the object to this main ID with using the IDs
				if std_star_name == main_id:
					#in this case we matched the standard name to something in the list of object IDs in simbad, so we dont have to check other objects returned from the query
					print(f'\033[38;5;{39}mTHE INPUT NAME FOR THE STANDARD STAR IS {std_star_name_0}\033[0m')
					print(f'\033[38;5;{39}mTHE SIMBAD MAIN ID FOR THE STANDARD STAR IS {std_star_name}\033[0m')
					
					return std_star_name
				else:
					#in this case we did not match the standard name to something in the list of obejct IDs in simbad, so we check the next object in the query
					continue

			#check to see if none object names queried matched the standard star name (e.g. a main_id did not overwrite the std_star_name var)
			if std_star_name == std_star_name_0:
				print(f'\033[38;5;{196}m{std_star_name} IS NOT SIMBAD SEARCHABLE AND WE CANNOT MATCH THE COORDINATES TO AN OBJECT IN SIMBAD RELIABLY\033[0m\n')
				print(f'\033[38;5;{39}mTHESE ARE THE OBJECTS RETURNED FROM SIMBAD WHEN SEARCHING {sky_coord.ra.deg} {sky_coord.dec.deg}\n')
				print(query_coords[['main_id', 'ra', 'dec']])
				print('\033[0m\n')

				#ask the user to input the standard star name
				name = input(f"\033[38;5;{39}mINPUT THE SIMBAD SEARCHABLE STANDARD STAR NAME: \033[0m")
				#try querying the input name
				query_name = Simbad.query_object(name)

				#if the name input does not generate a simbad result
				while len(query_name) == 0:
					#tell the user the name is not searchable
					print(f"\033[38;5;{196}mTHE INPUT NAME IS NOT SIMBAD SEARCHABLE\033[0m")
					#ask the user for the searchable name
					name = input(f"\033[38;5;{39}mINPUT THE SIMBAD SEARCHABLE STANDARD STAR NAME: \033[0m")
					#query the new input name
					query_name = Simbad.query_object(name)

				#if the input name generates a simbad result
				if len(query_name) != 0:
					#tell the user the name is accepted
					print(f"\033[38;5;{41}m{name} ACCEPTED AS THE STANDARD STAR NAME\033[0m")
					#set the standard star name
					std_star_name = query_name['main_id'][0]

					return std_star_name
	#if there is no name match and there are no coordinates in the header
	elif len(query_name) == 0 and coords == '':
		#print an error
		print(f'\033[38;5;{196}m{std_star_name} IS NOT SIMBAD SEARCHABLE AND HAS NO COORDINATES IN THE IGRINS HEADER\033[0m')

		#ask the user to input the standard star name
		name = input(f"\033[38;5;{39}mINPUT THE SIMBAD SEARCHABLE STANDARD STAR NAME: \033[0m")
		#try querying the input name
		query_name = Simbad.query_object(name)

		#if the name input does not generate a simbad result
		while len(query_name) == 0:
			#tell the user the name is not searchable
			print(f"\033[38;5;{196}mTHE INPUT NAME IS NOT SIMBAD SEARCHABLE\033[0m")
			#ask the user for the searchable name
			name = input(f"\033[38;5;{39}mINPUT THE SIMBAD SEARCHABLE STANDARD STAR NAME: \033[0m")
			#query the new input name
			query_name = Simbad.query_object(name)

		#if the input name generates a simbad result
		if len(query_name) != 0:
			#tell the user the name is accepted
			print(f"\033[38;5;{41}m{name} ACCEPTED AS THE STANDARD STAR NAME\033[0m")
			#set the standard star name
			std_star_name = query_name['main_id'][0]

			return std_star_name

	#if the query returned from the standard name has length greater than 0
	else:
		#take the first object as the standard star match
		std_star_name = query_name['main_id'][0]
		print(f'\033[38;5;{39}mTHE INPUT NAME FOR THE STANDARD STAR IS {std_star_name_0}\033[0m')
		print(f'\033[38;5;{39}mTHE SIMBAD MAIN ID FOR THE STANDARD STAR IS {std_star_name}\033[0m')

	return std_star_name


#make the directory for the plots if it is not made yet
if not os.path.isdir(plots_path+'/'):
	os.makedirs(os.path.dirname(plots_path+'/'))


#Find all .spec_a0v.fits files in directory
found_spec_a0v_file_paths = sorted(glob.glob(outdata_path+"/*spec_a0v.fits"))
spec_a0v_file_paths = []
for filepath in found_spec_a0v_file_paths: #Isolate only H band for now
	filename = filepath.split('/')[-1]
	if ('SDCH_' in filename) or ('_H.spec_a0v.fits' in filename):
		spec_a0v_file_paths.append(filepath)

#If no stanards were used for a night, say so and end the script
if len(spec_a0v_file_paths) == 0:
	print(f'\033[38;5;{177}mNO STANDARD STARS WERE OBSERVED ON THIS NIGHT. FLUX CALIBRATION IS NOT POSSIBLE. ENDING SCRIPT.\033[0m')
	quit()

#Load recipe log
recipe_log = Table.read(recipe_log_path, format='ascii', delimiter=',')
obsids = []
for recipe_log_entry in recipe_log: #Grab obsid for each entry
	obsid = int(recipe_log_entry['OBSIDS'].split(' ')[0])
	obsids.append(obsid)
obsids = np.array(obsids)

#Loop through each target
last_std_filename = ''
for spec_a0v_file_path in spec_a0v_file_paths:

	hdr = fits.Header() #Dictionary that stores all the results that will later be put into the final outputted FITS header

	#Read in target and standard star info and calcualte total exposure time on each
	pdfobj = PdfPages(plots_path + '/' + spec_a0v_file_path.split('/')[-1].replace('.fits', '.pdf'))
	spec_a0v_hdul = fits.open(spec_a0v_file_path)

	#grab the telescope of the observations so we can set the correct slit length and guiding error
	telescope = spec_a0v_hdul[4].header['TELESCOP']

	#setting the slit length and guiding based on the telescope
	if telescope == 'Discovery Channel':
		slit_length = 9.3 #slit length in arcsec
		guiding_error = 0 #kyle rec
	elif telescope == 'Gemini South':
		slit_length = 4.9 #slit length in arcseconds
		guiding_error = 0 #kyle rec
	else:
		slit_length = 14.8 #slit length in arcsecs at McDonald
		guiding_error = 1 #in arcseconds

	#grab the obsid for the target
	tgt_obsid = int(spec_a0v_hdul[4].header['OBSID'])
	#grab the obsid for the standard
	std_obsid = int(spec_a0v_hdul[6].header['OBSID'])

	#finding the index in the recipe log corresponding to the target obsid
	tgt_indx = np.where(tgt_obsid == obsids)[0][0]
	#finding the index in the recipe log corresponding to the standard obsid
	std_indx = np.where(std_obsid == obsids)[0][0]

	#construct the target .spec filename
	tgt_filename = 'SDCH_' + date + '_' + str(tgt_obsid).zfill(4) + '.spec.fits'
	#construct the standard .spec filename
	std_filename = 'SDCH_' + date + '_' + str(std_obsid).zfill(4) + '.spec.fits'
	#construct .spec_a0v.fits filename
	spec_a0v_filename = 'SDCH_' + date + '_' + str(tgt_obsid).zfill(4) + '.spec_a0v.fits'

	#grab the standard star name from the header
	std_star_name = spec_a0v_hdul[6].header['OBJECT']
	#grab the target name from the header
	tgt_name = spec_a0v_hdul[4].header['OBJECT']

	print(f'\033[38;5;{125}m\nFLUX CALIBRATION FOR TARGET: \033[0m'+f'\033[38;5;{196}m{tgt_name}\033[0m'+f'\033[38;5;{125}m WITH STANDARD: \033[0m'+f'\033[38;5;{196}m{std_star_name}\033[0m')
	#orange color text, see: https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797
	print(f'\033[38;5;{210}mOBSERVED AT {telescope}: ASSUMING {slit_length} ARCSEC SLIT LENGTH AND {guiding_error} ARCSEC GUIDING ERROR\033[0m')
	#yellow color text, see: https://gist.github.com/fnky/458719343aabd01cfb17a3a4f7296797
	print(f'\033[38;5;{222}mTARGET FILE PATH: ', spec_a0v_file_path+"\n\033[0m")

	try:
		#grabbing the standard star RA
		std_ra = spec_a0v_hdul[6].header['TELRA']
		#grabbing the standard star DEC
		std_dec = spec_a0v_hdul[6].header['TELDEC']
		#creating a string with the RA and DEC seperated by a space
		std_coords = str(std_ra) + ' ' + str(std_dec)
	except KeyError:
		print(f'\033[38;5;{196}mWARNING: NO STANDARD COORDINATES FROM HEADER\033[0m\n')
		std_coords = ''

	try:
		#grabbing the target position angle
		tgt_pa = spec_a0v_hdul[4].header['PASTART']
		#grabbing the standard position angle
		std_pa = spec_a0v_hdul[6].header['PASTART']
	except KeyError:
		#grabbing the target position angle
		tgt_pa = spec_a0v_hdul[4].header['TEL_PA']	
		#grabbing the standard position angle
		std_pa = spec_a0v_hdul[6].header['TEL_PA']	


	#grabbing the target exposure time
	tgt_exptime = float(spec_a0v_hdul[4].header['EXPTIME']) 
	#grabbing the standard exposure time
	std_exptime =  float(spec_a0v_hdul[6].header['EXPTIME'])

	#grabbing the number of target exposures
	tgt_nexp = len(recipe_log['FRAMETYPES'][tgt_indx].split(' '))
	#grabbing the number of standard exposures
	std_nexp = len(recipe_log['FRAMETYPES'][std_indx].split(' ')) 

	#calculating exposure time differently if the target/standard is nodded on/off slit instead of A/B pairs
	#see if the target is nodded on/off
	tgt_nod_off_slit = '_ONOFF' in recipe_log['RECIPE'][tgt_indx]
	#see if the standard is nodded on/off
	std_nod_off_slit = '_ONOFF' in recipe_log['RECIPE'][std_indx]
	if tgt_nod_off_slit: #If nod-off-slit
		tgt_nexp = tgt_nexp / 2 #Cut number of exposures in half since we are only on source half the time
	if std_nod_off_slit: #If nod-off-slit
		std_nexp = std_nexp / 2 #Cut number of exposures in half since we are only on source half the time

	#calculating the time spent on target	
	tgt_total_on_time = tgt_exptime * tgt_nexp
	#calculating the time spent on standard
	std_total_on_time = std_exptime * std_nexp

	#adding the time spent on target to the header
	hdr['FTTOTEXP'] = (tgt_total_on_time, "Total exp. time (s) on target.")
	#adding the time spent on standard to the header
	hdr['FSTOTEXP'] = (std_total_on_time, "Total exp. time (s) on std. star.")

	#close the spec_a0v.fits file
	spec_a0v_hdul.close()

	#Read target spectrum as a muler IGRINSSpectrumList object
	tgt_spec = igrins.readIGRINS(outdata_path+'/'+tgt_filename) 
	#check to see if the standard name has changed for this target (if it is the same, we dont have to read in a new standard file)
	if std_filename != last_std_filename:
		#Read standard spectrum as a muler IGRINSSpectrumList object
		std_spec = igrins.readIGRINS(outdata_path+'/'+std_filename)


	print(f'\033[38;5;{79}mESTIMATING TARGET SLIT THROUGHPUT\033[0m')
	#Apply throughput corrections if the target is not an extended object
	if 'STELLAR' in recipe_log['RECIPE'][tgt_indx]:
		#Calculate slit throughput for target
		tgt_throughput, tgt_m, tgt_b, target_flux_corr, target_flux_m, target_flux_b = tgt_spec.getSlitThroughput(name=tgt_name, name_prefix='Target:', slit_length=slit_length, PA=tgt_pa, guiding_error=guiding_error, plot=True, pdfobj=pdfobj, nod_off_slit=tgt_nod_off_slit)
		#divide target spec by the throughput (normalize by the throughput)
		tgt_spec = tgt_spec / tgt_throughput

		#Fit Target ThroughPut M (coeff from linear fit)
		hdr['FTTPM'] = (tgt_m, "Coeff m for target throughput m*(1/wave[um]+b")
		#Fit Target ThroughPut B (coeff from linear fit)
		hdr['FTTPB'] = (tgt_b, "Coeff b for target throughput m*(1/wave[um]+b")

		print(f'\033[38;5;{79}mTHE ESTIMATED TARGET SLIT THROUGHPUT IS {np.nanmedian(tgt_throughput)*100:.2f} PERCENT\033[0m')
		print(f'\033[38;5;{79}mTHE MEDIAN DIFFERENCE IN FLUX BETWEEN THE NODS FOR THE TARGET IS {(np.nanmedian(target_flux_corr)-1)*100:.2f} PERCENT\033[0m')

	#for printing aesthetics
	print(' ')
	if std_filename != last_std_filename: #If standard has changed
		#pre flight check--verify that the standard star name is SIMBAD searchable
		std_star_name = standard_star_name_check(std_star_name, std_coords)

		print(f'\033[38;5;{33}m\nESTIMATING STANDARD SLIT THROUGHPUT\033[0m')
		#Calculate slit throughput for standard
		std_throughput, std_m, std_b, std_flux_corr, std_flux_m, std_flux_b = std_spec.getSlitThroughput(name=std_star_name, name_prefix='Standard Star:', slit_length=slit_length, PA=std_pa, guiding_error=guiding_error, plot=True, pdfobj=pdfobj, nod_off_slit=std_nod_off_slit)		#divide standard spec by the throughput (normalize by the throughput)
		#divide standard spec by the throughput (normalize by the throughput)
		std_spec = std_spec / std_throughput

	#Fit Standard ThroughPut M (coeff from linear fit)
	hdr['FSTPM'] = (std_m, "Coeff m for stdstar throughput m*(1/wave[um])+b")
	#Fit Standard ThroughPut B (coeff from linear fit)
	hdr['FSTPB'] =  (std_b, "Coeff b for stdstar throughput m*(1/wave[um])+b")

	print(f'\033[38;5;{33}mTHE ESTIMATED STANDARD SLIT THROUGHPUT IS {np.nanmedian(std_throughput)*100:.2f} PERCENT\033[0m')
	print(f'\033[38;5;{33}mTHE MEDIAN DIFFERENCE IN FLUX BETWEEN THE NODS FOR THE STANDARD IS {(np.nanmedian(std_flux_corr)-1)*100:.2f} PERCENT\033[0m')

	print(f'\033[38;5;{63}m\nMODELING STANDARD STAR\033[0m')
	#Run standard star crude telluric correction and stellar atmosphere model fitting to H I line profiles
	if std_filename != last_std_filename: #If standard has changed
		#model the star using PHOENIX template

		#Get tellurics estimate from the PLP's standard star continuum estimate
		spec_a0v_ext1 = igrins.readIGRINS(outdata_path+'/'+spec_a0v_filename)
		spec_a0v_ext9 = igrins.readIGRINS(outdata_path+'/'+spec_a0v_filename, extension=9)
		tellurics_estimate = spec_a0v_ext9 / spec_a0v_ext1
		total_trans = tellurics_estimate.fullFitTellurics(plot=True, pdfobj=pdfobj)

		std_model, resampled_std_model, std_simbad_phot, std_fit_results = std_spec.fitStandardStar(name=std_star_name, coords=std_coords, name_prefix='Standard Star:', plot=True, pdfobj=pdfobj, total_trans=total_trans)
		std_model_phot = photometry() #Estimate photometry from best fit model
		std_model_phot.set_photometry(std_model)

	#Fit Model Temperature (Standard)
	hdr['FMTEFF'] = (std_fit_results['TEFF'], "Std star model fit Teff (K)")
	#Fit Model logg (Standard)
	hdr['FMLOGG'] = (std_fit_results['LOGG'], "Std star model fit log(g)")
	#Fit Model metallicity (Standard)
	hdr['FMZ'] = (std_fit_results['Z'], "Std star model fit metallicity")
	#Fit Model rotational velocity (Standard)
	hdr['FMROTV'] = (std_fit_results['ROTV'], "Std star model fit rotational velocity (km/s)")
	#Fit Model radial velocity(Standard)
	hdr['FMRADV'] = (std_fit_results['RADV'], "Std star model fit radial velocity (km/s)")
	#Fit Model alpha (Standard)
	hdr['FMALPHA'] = (std_fit_results['ALPHA'], "Std star model fit alpha, H I scale power")

	#this the file for a bunch of diagnostic plots for the fitting (and other things)
	pdfobj.close()
	
	#Do the absolute flux calibration of the target
	tgt_spec_flux_calibrated = tgt_spec * (std_total_on_time/tgt_total_on_time) * (resampled_std_model/std_spec) 
	#this stitched spectrum is just for the synth photometry check
	#these functions must be done in this specific order to work for IGRINS/IGRINS-2
	flattened_tgt_spec_flux_calibrated = tgt_spec_flux_calibrated.trim_overlap().remove_nans().remove_outliers(threshold = 20).stitch() #trim overlap here default is 0.5

	#SNR sanity check
	#read in the spec_a0v.fits file
	tgt_old_spec_a0v = igrins.readIGRINS(spec_a0v_file_path)
	print(f'\033[38;5;{57}m\nCHECK SNR FOR TARGET\033[0m')
	print(f'\033[38;5;{99}m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m')
	for order in [10, 20, 30 ,40]:
		total_flux_uncalibrated = np.nansum(tgt_old_spec_a0v[order].flux.value)
		total_var_uncalibrated = np.nansum(tgt_old_spec_a0v[order].uncertainty.array**2)
		total_flux_calibrated = np.nansum(tgt_spec_flux_calibrated[order].flux.value)
		total_var_calibrated = np.nansum(tgt_spec_flux_calibrated[order].uncertainty.array**2)
		print('Sanity SNR check for order '+str(order)+'.  SNR should be similar.')
		print('.spec_a0v.fits SNR = ', total_flux_uncalibrated / total_var_uncalibrated**0.5)
		print('New calibrated SNR = ', total_flux_calibrated / total_var_calibrated**0.5)
	print(f'\033[38;5;{99}m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m')

	#Do synthetic photomtery on the absolute flux calibrated target spectrum to calculate H and K magnitudes of target
	#and print results to command line as a sanity check and for comparison to 2MASS photometry
	tgt_data_phot = photometry() 
	tgt_H_mag = tgt_data_phot.get(flattened_tgt_spec_flux_calibrated, band='H')
	tgt_K_mag = tgt_data_phot.get(flattened_tgt_spec_flux_calibrated, band='K')
	print(f'\033[38;5;{201}m\nCOMPARE SYNTH PHOTMETRY OF ABS FLUX CALIBRATION TO 2MASS\033[0m')
	print(f'\033[38;5;{207}m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\033[0m')
	#print('Synthetic photometry of absolute flux calibrated target spectrum, compare to 2MASS')
	print(tgt_name +' H mag: '+str(tgt_H_mag))
	print(tgt_name +' K mag: '+str(tgt_K_mag))
	print(f'\033[38;5;{207}m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\033[0m')

	#Fit Target H-band Mag
	hdr['FTHMAG'] = (tgt_H_mag, 'Target synthetic H band mag.')
	#Fit Target K-band Mag
	hdr['FTKMAG'] = (tgt_K_mag, 'Target synthetic K band mag.')


	#Construct new .flux_spec.fits files to store the calibrated target spectrum and best fit model spectrum in almost the same format as the .spec_a0v.fits files
	for band in ['H', 'K']:
		#the front string for .flux_spec.fits
		filename = spec_a0v_file_path.split('/')[-1]
		#place to write the new file to
		filepath = spec_a0v_file_path.split(filename)[0]
		#full filename to write to
		filename = filename.replace('SDCH_', 'SDC'+band+'_').replace('_H.spec_a0v', '_'+band+'.spec_a0v')

		#Generate 2D arrays for calibrated target flux, variance, and throughput (and carry the wavelength array)
		tgt_flux_arr = tgt_spec_flux_calibrated.get_plp_array(band=band, kind='flux')
		tgt_var_arr = tgt_spec_flux_calibrated.get_plp_array(band=band, kind='var')
		tgt_wave_arr = tgt_spec_flux_calibrated.get_plp_array(band=band, kind='wave')
		if 'STELLAR' in recipe_log['RECIPE'][tgt_indx]: #Only apply throughput correction for stellar targets
			tgt_throughput_arr = tgt_m*(1/tgt_wave_arr) + tgt_b
		else:
			tgt_throughput_arr = np.ones(np.shape(tgt_wave_arr))

		#Generate 2D array for standard throughput (and carry the wavelength)
		std_wave_arr = std_spec.get_plp_array(band=band, kind='wave')
		std_throughput_arr = std_m*(1/std_wave_arr) + std_b

		#Generate 2D array for standard star best fit model flux
		std_model_arr = resampled_std_model.get_plp_array(band=band, kind='flux') 

		#Open hdul object for spec_a0v.fits file
		hdul = fits.open(filepath + filename)

		#Append the new keywords from the fit to the primary header
		hdul[0].header += hdr

		#adding the data to the fits file
		hdul[1] = fits.ImageHDU(data=tgt_flux_arr, header=hdul[1].header)
		hdul[1].header['BUNIT'] = ('erg s-1 cm-2 micron-1', 'Flux Units')
		hdul[2] = fits.ImageHDU(data=tgt_var_arr, header=hdul[2].header)
		hdul[8] = fits.ImageHDU(data=std_model_arr, header=hdul[8].header)
		hdul[8].header["EXTNAME"] = "A0V_MODEL_SPEC"
		hdul[9] = fits.ImageHDU(data=tgt_throughput_arr, header=hdul[9].header)
		#hdul[9].header['EXTNAME'] = "SCI" #this is a GN IGRINS-2 thing
		hdul[9].header['EXTNAME'] = "TGT_THROUGHPUT"
		#hdul[9].header["EXTVER"] = 6 #this is a GN IGRINS-2 thing
		hdul[10] = fits.ImageHDU(data=std_throughput_arr, header=hdul[10].header)
		#hdul[10].header['EXTNAME'] = "SCI" #this is a GN IGRINS-2 thing
		hdul[10].header['EXTNAME'] = "A0V_THROUGHPUT"
		#hdul[10].header["EXTVER"] = 7 #this is a GN IGRINS-2 thing

		#Delete last extension since it is not used here
		hdul.pop(11)
		#Run astropy fits verification
		hdul.verify()
		#Save final file with the extension flux_a0v.fits
		hdul.writeto(filepath + filename.replace('spec_a0v.fits', 'flux_a0v.fits'), overwrite=True)
		#close the fits file
		hdul.close()

	print(f'\033[38;5;{219}mCOMPLETED FLUX CALIBRATION FOR \033[0m'+f'\033[38;5;{196}m{tgt_name}\033[0m'+f"\033[38;5;{219}m WITH STANDARD \033[0m"+f'\033[38;5;{196}m{std_star_name}\033[0m')

	#Keep track of the last standard filename so we don't need to re-run the fitting if we are using the same standard (same standard file number--standards observed multiple times have each observation sequence fit)
	last_std_filename = std_filename



