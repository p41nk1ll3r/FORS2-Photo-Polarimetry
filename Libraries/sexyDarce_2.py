import math
import sep
import random
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse as mplEllipse
from astroquery.ipac.ned import Ned
from astropy.coordinates import Angle
from astropy.io import fits
from photutils.isophote import EllipseGeometry
from photutils.isophote import Ellipse


def list_sources_info(data_m, sources, acp_file, rej_file, max_ecc, min_Area, max_Area, data_err):
	src_a = open(acp_file, "w")
	src_r = open(rej_file, "w")
	
	src_a.write("Obj_n;Seg_n;Eccentricity;a;b;theta;x_barycenter;y_barycenter;x_peak;y_peak;x_cpeak;y_cpeak;Flux;cFlux;elFlux;elFlux_err;Area\n")
	src_r.write("Obj_n;Seg_n;Eccentricity;a;b;theta;x_barycenter;y_barycenter;x_peak;y_peak;x_cpeak;y_cpeak;Flux;cFlux;elFlux;elFlux_err;Area\n")
	
	j = 0
	k = 0
	rej_obj = []
	
	for i in range(len(sources)):
		ecc = math.sqrt(1 - (sources['b'][i]/sources['a'][i])**2)
		area = sources['npix'][i]
		flux = sources['flux'][i]
		
		#x != x is a library independent simple test for NaNs
		flux_flag = flux == flux
		
		if (ecc <= max_ecc) and (area <= max_Area) and (area >= min_Area) and flux_flag:
			
			j = j + 1
			
			el_flux, el_flux_err, el_flags = sep.sum_ellipse(data = data_m, x = np.array(sources['x'][i], ndmin = 1), y = np.array(sources['y'][i], ndmin = 1), a = np.array(sources['a'][i], ndmin = 1), b = np.array(sources['b'][i], ndmin = 1), theta = np.array(sources['theta'][i], ndmin = 1), err = data_err, subpix = 0)
			el_flux_flag = (el_flux == el_flux) and (el_flux_err == el_flux_err)
			
			if(el_flux_flag):
				#python treats arrays differently than R, this file will be read by R, as such x and y (rows and columns) must be switched around 
				src_a.write(str(j)+";"+str(i+1)+";"+str(ecc)+";"+str(sources['a'][i])+";"+str(sources['b'][i])+";"+str(sources['theta'][i])+";"+str(sources['y'][i])+";"+str(sources['x'][i])+";"+str(sources['ypeak'][i])+";"+str(sources['xpeak'][i])+";"+str(sources['ycpeak'][i])+";"+str(sources['xcpeak'][i])+";"+str(sources['flux'][i])+";"+str(sources['cflux'][i])+";"+str(el_flux[0])+";"+str(el_flux_err[0])+";"+str(sources['npix'][i])+"\n")
		
		if (ecc > max_ecc) or (area > max_Area) or (area < min_Area) or (not flux_flag) or (not el_flux_flag):
			k = k + 1
			
			src_r.write(str(k)+";"+str(i+1)+";"+str(ecc)+";"+str(sources['a'][i])+";"+str(sources['b'][i])+";"+str(sources['theta'][i])+";"+str(sources['y'][i])+";"+str(sources['x'][i])+";"+str(sources['ypeak'][i])+";"+str(sources['xpeak'][i])+";"+str(sources['ycpeak'][i])+";"+str(sources['xcpeak'][i])+";"+str(sources['flux'][i])+";"+str(sources['cflux'][i])+";"+"NaN"+";"+"NaN"+";"+str(sources['npix'][i])+"\n")
			
			rej_obj.append(i+1)
	
	src_a.close()
	src_r.close()
	
	return(rej_obj)



def ens_max_snr_circ_aperture(data_m, sources, data_err, min_r=3, max_r=15, eps=0.01):
	
	N = len(sources)
	ens_r = 0
	ens_SNR = 0
	SNR_eps = 10 * eps
	i = 0
	my_step = 1
	min_ri = min_r
	max_ri = max_r
	SNR_rec = []#tagged elements are temporary and should be deleted
	flux_rec = []#
	r_rec = []#
	
	while eps < SNR_eps and i < 4:
		hold_SNR = ens_SNR
		radii = np.arange(min_ri, max_ri+my_step, my_step)
		
		for ri in radii:
			r_rec.append(ri)#
			
			temp_SNR = 0
			temp_sum_flux = 0#
			temp_N = 0
			
			for n in range(N):
				temp_flux, temp_err, temp_fl = sep.sum_circle(data = data_m, x = np.array(sources['x'][n], ndmin = 1), y = np.array(sources['y'][n], ndmin = 1), r = ri, err = data_err, subpix = 0)
				if temp_flux == temp_flux:
					temp_N += 1
					temp_sum_flux += temp_flux
					temp_SNR += (temp_flux / temp_err)
					#print('With aperture of r = '+str(ri)+', the flux of source #'+str(n)+' of '+str(N)+' is '+str(temp_flux)+'...')
			
			temp_SNR /= temp_N			
			SNR_rec.append(float(temp_SNR))#
			flux_rec.append(float(temp_sum_flux))#
			
			if temp_SNR > ens_SNR:
				ens_SNR = temp_SNR
				ens_r = ri
		
		SNR_eps = ens_SNR - hold_SNR
		min_ri = ens_r - my_step
		max_ri = ens_r + my_step
		
		print('Present best aperture r is '+str(ens_r)+' with an average SNR of '+str(ens_SNR)+'...')
		
		if min_ri < min_r:
			min_ri = min_r
		
		if max_ri > max_r:
			max_ri = max_r
		
		my_step /= 10
		i += 1 
		print('New step is '+str(my_step)+'...')
	
	return (ens_r, ens_SNR, flux_rec, r_rec, SNR_rec)



#RFR stands for Radius to Fwhm Ratio
def indiv_max_snr_circ_aperture(data_m, sources, file, data_err, max_r=10, eps=0.001, max_it = 5, gain = 1):
	src_f = open(file, "w")
	src_f.write("******* OPTIMAL RADII COMPUTATION HISTORY *******\n")
	
	N = len(sources['x'])
	indiv_r = []
	indiv_flux = []
	indiv_err = []
	indiv_SNR = []
	min_r = 1 
	
	for n in range(N):
		src_f.write("##### Source "+str(n+1)+"/"+str(N)+" #####\n")
		
		# This function was thought of has being called after source selection has been performed
		# Whitin that selection sources with a NA pixel close to their center have been rejected
		# As such setting opt_r as 1 is the safest option in an agnostic setting 
		# regarding the actual apperture radius used before to discard sources
		opt_r = min_r
		temp_flux, temp_err, temp_fl = sep.sum_circle(data = data_m, x = np.array(sources['x'][n], ndmin = 1), y = np.array(sources['y'][n], ndmin = 1), r = opt_r, err = data_err, subpix = 0, gain = gain)
		opt_flux = float(temp_flux)
		opt_err = float(temp_err)
		opt_r_SNR = opt_flux / opt_err
		
		src_f.write("+++ Initial values: r = "+str(opt_r)+"; flux = "+str(opt_flux)+"; SNR = "+str(opt_r_SNR)+" +++\n")
		
		SNR_eps = eps + 1
		i = 0
		my_step = 1
		min_ri = min_r
		max_ri = max_r
		del_ri = max_ri - min_ri
		
		while SNR_eps > eps or i < max_it:
			t_SNR = opt_r_SNR
			radii = np.arange(min_ri, max_ri+my_step, my_step)
			
			src_f.write("««« Iteration "+str(i+1)+"/"+str(max_it)+", testing r = "+str(radii)+" »»»\n")
			
			for ri in radii:
				temp_flux, temp_err, temp_fl = sep.sum_circle(data = data_m, x = np.array(sources['x'][n], ndmin = 1), y = np.array(sources['y'][n], ndmin = 1), r = ri, err = data_err, subpix = 0, gain = gain)
				
				if(temp_err != temp_err):
					temp_err = np.sqrt(np.nanmedian(data_err)**2 * np.pi * ri**2 + temp_flux / gain)
				
				temp_SNR = temp_flux / temp_err
				
				if (temp_SNR == temp_SNR) and (temp_SNR > opt_r_SNR):
					opt_r_SNR = temp_SNR
					opt_flux = float(temp_flux)
					opt_err = float(temp_err)
					opt_r = ri
					src_f.write("!! A radius that yields a higher SNR has been found !!\n-")
					
				src_f.write("--> r = "+str(ri)+"; flux = "+str(float(temp_flux))+"; SNR = "+str(float(temp_SNR))+"\n")
			
			#simple nan test
			if opt_flux == opt_flux:
				SNR_eps = opt_r_SNR - t_SNR
				del_ri = del_ri / 2
				min_ri = opt_r - del_ri / 2
				max_ri = opt_r + del_ri / 2
				my_step = my_step / 2
			
			if min_ri < min_r:
				min_ri = min_r
			
			if max_ri > max_r:
				max_ri = max_r
			
			i += 1
		
		indiv_r.append(opt_r)
		indiv_SNR.append(opt_r_SNR)
		indiv_flux.append(opt_flux)
		indiv_err.append(opt_err)
		
		src_f.write("\n")
	
	src_f.close()
	
	return (indiv_err, indiv_flux, indiv_r, indiv_SNR)#



#Expects avg_el to have parameters 'theta' in position 0, and 'ax_ratio' in position 3
def ens_max_snr_ellp_aperture(data_m, sources, avg_el, data_err, min_a=3, max_a=15, eps=0.01):
	
	N = len(sources)
	ens_a = 0
	ens_SNR = 0
	SNR_eps = 10 * eps
	i = 0
	my_step = 1
	min_ai = min_a
	max_ai = max_a
	SNR_rec = []#tagged elements are temporary and should be deleted
	flux_rec = []#
	a_rec = []#
	
	while eps < SNR_eps and i < 4:
		hold_SNR = ens_SNR
		ax_ai = np.arange(min_ai, max_ai+my_step, my_step)
		
		for ai in ax_ai:
			#print('testing a = '+str(ai))
			a_rec.append(ai)#
			temp_SNR = 0
			temp_sum_flux = 0#
			
			for n in range(N):
				bi = ai / avg_el[3]
				#print('and b = '+str(bi)+' for star '+str(n+1))
				temp_flux, temp_err, temp_fl = sep.sum_ellipse(data = data_m, x = np.array(sources['x'][n], ndmin = 1), y = np.array(sources['y'][n], ndmin = 1), a = np.array(ai, ndmin = 1), b = np.array(bi, ndmin = 1), theta = np.array(avg_el[0], ndmin = 1), err = data_err, subpix = 0)
				#print('flux = '+str(temp_flux))
				#print('err_flux = '+str(temp_err)+'\n')
				temp_sum_flux += temp_flux
				temp_SNR += 1 / N * temp_flux / temp_err
			
			SNR_rec.append(float(temp_SNR))#
			flux_rec.append(float(temp_sum_flux))#
			
			if temp_SNR > ens_SNR:
				ens_SNR = temp_SNR
				ens_a = ai
		
		SNR_eps = ens_SNR - hold_SNR
		max_ai = ens_a + my_step
		min_ai = ens_a - my_step
		
		if min_ai <= 3:
			min_ai = 3
		
		my_step /= 10
		
		i += 1 
	
	return (ens_a, ens_SNR, flux_rec, a_rec, SNR_rec)



#RFR stands for Radius to Fwhm Ratio
def indiv_max_snr_ellp_aperture(data_m, sources, data_err, min_a=3, max_a=15, eps=0.01):
	
	N = len(sources)
	indiv_a = []
	indiv_flux = []
	indiv_err = []
	indiv_SNR = []
	a_rec = []#tagged elements are temporary and should be deleted
	flux_rec = []#
	SNR_rec = []#
	
	for n in range(N):
		a_rec.append([])#
		flux_rec.append([])#
		SNR_rec.append([])#
		opt_a = float("nan")
		opt_r_SNR = 0
		opt_flux = float("nan")
		opt_err = float("nan")
		opt_SNR = 0
		SNR_eps = 10 * eps
		i = 0
		my_step = 1
		min_ai = min_a
		max_ai = max_a
		
		while eps < SNR_eps and i < 4:
			t_SNR = opt_SNR
			ax_ai = np.arange(min_ai, max_ai+my_step, my_step)
			
			for ai in ax_ai:
				bi = ai / sources['ax_ratio'][n]
				temp_flux, temp_err, temp_fl = sep.sum_ellipse(data = data_m, x = np.array(sources['x'][n], ndmin = 1), y = np.array(sources['y'][n], ndmin = 1), a = np.array(ai, ndmin = 1), b = np.array(bi, ndmin = 1), theta = np.array(sources['theta'][n], ndmin = 1), err = data_err, subpix = 0)
				temp_SNR = temp_flux / temp_err
				a_rec[n].append(ai)#
				flux_rec[n].append(temp_flux)#
				SNR_rec[n].append(temp_SNR)#
				
				if (temp_SNR > opt_SNR)  and (temp_flux == temp_flux):
					opt_SNR = temp_SNR
					opt_flux = float(temp_flux)
					opt_a = ai
					opt_err = float(temp_err)
			
			if opt_flux == opt_flux:
				SNR_eps = opt_SNR - t_SNR
				min_ai = opt_a - my_step
				max_ai = opt_a + my_step
			
			if min_ai <= 3:
				min_ai = 3
			
			my_step /= 10
			i += 1
		
		indiv_a.append(opt_a)
		indiv_SNR.append(opt_SNR)
		indiv_flux.append(opt_flux)
		indiv_err.append(opt_err)
	
	return (indiv_err, indiv_flux, indiv_a, indiv_SNR, flux_rec, a_rec, SNR_rec)#



def plot_sources_map(data_m, sources, file, max_ecc):
	fig, ax = plt.subplots()
	m, s = np.mean(data_m), np.std(data_m)
	datas = np.transpose(data_m)
	im = ax.imshow(datas, interpolation='nearest', cmap='gray', aspect='equal', vmin=m-s, vmax=m+s, origin='lower')
	
	for i in range(len(sources)):
		ecc = math.sqrt(1 - (sources['b'][i]/sources['a'][i])**2)
		
		if ecc <= max_ecc:
			e = mplEllipse(xy=(sources['y'][i], sources['x'][i]), width=sources['a'][i], height=sources['b'][i], angle=(sources['theta'][i]*180./np.pi))
			e.set_facecolor('none')
			e.set_edgecolor('red')
			ax.add_artist(e)
	
	plt.savefig(file)



#Converts a coordinate value to decimal degrees	
def sexag_to_deg(value, unit, N):
	
	if unit == 'Decimal Degree':
		return round(float(value), N)
	elif unit == 'Sexagesimal':
		angle = Angle(value, unit='hourangle' if ':' in value else 'deg')
		return round(angle.degree, N)
	else:
		print(f"ERROR: Unsupported coordinate unit: {unit}. Returning NaN.")
		return float('nan')



#Retrieves the most recent Equatorial J2000 coordinates of the given object from NED in decimal degrees
def get_NED_coord(object_name, N = 5):
	
	# Query NED for object coordinates
	table = Ned.get_table(object_name, table='positions')

	# Find the first (most recent) row with Equatorial J2000 coordinates
	j2000_row = None
	for row in table:
		if row['Published System Coordinate'] == 'Equatorial' and row['Published Equinox'] in ['J2000', 'J2000.0']:
			j2000_row = row
			break

	if j2000_row is None:
		print("ERROR: Could not retrieve Equatorial J2000 coordinates for the object. Returning NaN.")
		return float('nan')

	# Extract RA and Dec values and their units
	ra_value = j2000_row['RA']
	dec_value = j2000_row['DEC']
	ra_unit = j2000_row['Published Unit']
	dec_unit = j2000_row['Published Unit']
	
	# Convert to decimal degrees if necessary
	ra_decimal = sexag_to_deg(ra_value, ra_unit, N)
	dec_decimal = sexag_to_deg(dec_value, dec_unit, N)
	
	return round(ra_decimal, N), round(dec_decimal, N)



#Retrieves the nth most recent major axis of the given object from NED and returns that value in arcsecs.
def get_NED_major_axis(object_name, n=1):
	
	# Query NED for object diameters
	table = Ned.get_table(object_name, table='diameters')
	
	if len(table) == 0:
		print("ERROR: There are no measurements, could not retrieve major axis for the object. Returning NaN.")
		return float('nan')
	
	# Sort the table by Refcode in descending order
	sorted_table = table[np.argsort(table['Refcode'])[::-1]]
	
	if n == 1:
		card = 'st'
	elif n == 2:
		card = 'nd'
	elif n == 3:
		card = 'rd'
	else:
		card = 'th'

	if len(sorted_table) < n:
		print(f"ERROR: There is no {n}{card} measurement, could not retrieve major axis for the object. Returning NaN.")
		return float('nan')

	# Find the row with the major axis
	the_row = None
	count = 0
	for row in sorted_table:
		count += 1
		if count == n:
			the_row = row
			break

	if the_row is None:
		print("ERROR: Could not retrieve major axis for the object. Returning NaN.")
		return float('nan')
	
	# Extract major axis value and its unit
	major_axis_value = the_row['NED Major Axis']
	
	return float(major_axis_value)



#Retrieves the nth most recent axis ratio of the given object from NED and returns that value in arcsecs.
def get_NED_axis_ratio(object_name, n=1):
	
	# Query NED for object diameters
	table = Ned.get_table(object_name, table='diameters')
	
	if len(table) == 0:
		print("ERROR: There are no measurements, could not retrieve minor axis for the object. Returning NaN.")
		return float('nan')
	
	# Sort the table by Refcode in descending order
	sorted_table = table[np.argsort(table['Refcode'])[::-1]]
	
	if n == 1:
		card = 'st'
	elif n == 2:
		card = 'nd'
	elif n == 3:
		card = 'rd'
	else:
		card = 'th'
	
	if len(sorted_table) < n:
		print(f"ERROR: There is no {n}{card} measurement, could not retrieve axis ratio for the object. Returning NaN.")
		return float('nan')

	# Find the row with the axis ratio
	the_row = None
	count = 0
	for row in sorted_table:
		count += 1
		if count == n:
			the_row = row
			break

	if the_row is None:
		print("ERROR: Could not retrieve axis ratio for the object. Returning NaN.")
		return float('nan')
	
	# Extract minor axis value and its unit
	axis_ratio_value = the_row['NED Axis Ratio']
	
	return float(axis_ratio_value)



#Retrieves the nth most recent position angle of the given object from NED and returns that value in arcsecs.
def get_NED_position_angle(object_name, n=1):
	
	# Query NED for object diameters
	table = Ned.get_table(object_name, table='diameters')
	
	if len(table) == 0:
		print("ERROR: There are no measurements, could not retrieve minor axis for the object. Returning NaN.")
		return float('nan')
	
	# Sort the table by Refcode in descending order
	sorted_table = table[np.argsort(table['Refcode'])[::-1]]
	
	if n == 1:
		card = 'st'
	elif n == 2:
		card = 'nd'
	elif n == 3:
		card = 'rd'
	else:
		card = 'th'
	
	if len(sorted_table) < n:
		print(f"ERROR: There is no {n}{card} measurement, could not retrieve position angle for the object. Returning NaN.")
		return float('nan')

	# Find the row with the axis ratio
	the_row = None
	count = 0
	for row in sorted_table:
		count += 1
		if count == n:
			the_row = row
			break

	if the_row is None:
		print("ERROR: Could not retrieve position angle for the object. Returning NaN.")
		return float('nan')
	
	# Extract minor axis value and its unit
	PA_value = the_row['NED Position Angle']
	
	return float(PA_value)



def get_closest_iso_for_target(DATA, X, Y, SMA, EPS, PA):
	TSMA = SMA
	
	while True:		
		try:
			with warnings.catch_warnings(record=True) as w:
				warnings.simplefilter("always")
				
				test_circ = EllipseGeometry(x0 = X, y0 = Y, sma = TSMA, eps = EPS, pa = PA * np.pi / 180.0)
				ellipse = Ellipse(DATA, test_circ)
				isolist = ellipse.fit_image()
				
				if w:
					TSMA = abs(random.gauss(SMA, SMA / 3))
					
					for warning in w:
						print(warning.message)
						print("There was a warning, attempting a new ellipse fit with a smaller input circle.")
						continue
				else:
					break
		except Exception as e:
			TSMA = abs(random.gauss(SMA, SMA / 2))
			print("Error:", str(e))
			print("Attempting a new ellipse fit with a smaller input circle.")
			continue
	
	iso = isolist.get_closest(SMA)
	
	return iso
	


def load_ctr_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Holder variables for contour data
    data = []
    current_contour = []
    reading_data = False
    
    # Iterate through each line in the file
    for line in lines:
        line = line.strip()
        # Check if the line is part of the header
        if line.startswith('#') or line == 'image':
            if current_contour:  # If we were reading data, reset for the next block
                data.append(current_contour)
                current_contour = []
            reading_data = False
        elif line.startswith('level='):
            if current_contour:  # Start of a new contour block, save the previous
                data.append(current_contour)
                current_contour = []
            reading_data = True  # Start reading contour data
        elif reading_data:
            if line.startswith('(') or line == ')':
                continue  # Ignore these lines, as they don't contain coordinates
            else:
                # Parse the x, y coordinates and append to the current contour
                x, y = map(float, line.split())
                current_contour.append((x, y))
    
    # Add the last contour data block if not empty
    if current_contour:
        data.append(current_contour)
    
    # Return list of contours
    return data
    

    
def create_2d_array_from_ctr(ctr_data, max_x, max_y):
    my_arr = np.zeros((max_y, max_x), dtype=int)
    
    for contour in ctr_data:
        for coord in contour:
            x, y = map(int, coord)
            my_arr[y - 1, x - 1] = 1
    
    return my_arr
    


def save_arr_as_fits(arr, output_file):
    arr_int16 = arr.astype(np.int16)
    hdu = fits.PrimaryHDU(arr_int16)
    hdu.writeto(output_file, overwrite=True)
    


