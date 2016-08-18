#Writes dcr corrected galaxy lists which will be rearranged and then input to IDL program

#Only r-band

#Run separately for each chip

#Input correct paths in 2 places initially. Then, just change chip number every time.

import os
import math
import numpy as np
import operator
from bisect import bisect_left
import matplotlib.pyplot as plt

#INPUT DETECTOR NUMBER***
chip = '1'


obs_info = np.loadtxt("r-band_day2_obs_info.txt")
rot_ang = obs_info[:,4] #rotator angles
z_ang = obs_info[:,5] #zenith angles
all_lst = obs_info[:,6]

conv_info = np.loadtxt("r%s_deg-per-pixel.txt" %chip)
ra_conv = conv_info[:,1]
dec_conv = conv_info[:,2]

star_tabl = np.loadtxt("star_shift.txt")
r_min_i = star_tabl[:,0] #is sorted already
star_b_off = star_tabl[:,1]


x_diffs = []

y_diffs = []

colors = []

offsets = []

x_rats = []

y_rats = []

starcolors = []

#Function to be used later
def takeClosest(myList, myNumber):
				    """
				    Assumes myList is sorted. Returns closest value to myNumber.

				    If two numbers are equally close, return the smallest number.
				    """
				    pos = bisect_left(myList, myNumber)
				    if pos == 0:
				        return myList[0]
				    if pos == len(myList):
				        return myList[-1]
				    before = myList[pos - 1]
				    after = myList[pos]
				    if after - myNumber < myNumber - before:
				       return after
				    else:
				       return before


# Empty list to be filled by select galaxies

gal_list = []


path = '/home/martine/Documents/Astr_research/Extended_catalogs/'  #***CHANGE PATH to wherever extended catalogs are****
for datafile in os.listdir(path):
	if datafile.startswith("bluer") and datafile[8] == chip:
		
		Catalog = open(path+datafile, 'r')
		for source in Catalog:
			data = source.split() #splits each property up by the spaces in between
			SDSS_objid = long(float(data[73]))
			if SDSS_objid != 0.0: #if the LBT source has a closest matching SDSS source
				LBT_number = float(data[0])
				Spread_model = float(data[13])
				probPSF_g = float(data[90])
				probPSF_r = float(data[91])
				probPSF_i = float(data[92])
				SDSS_dec = float(data[76])
				SDSS_decErr_given = float(data[77]) #in arcseconds
				SDSS_decErr_final = math.sqrt((SDSS_decErr_given*1.186)**2+0.01721**2) #in arcseconds
				SDSS_decErr_degrees = SDSS_decErr_final / 3600 #converting " into degrees
				SDSS_ra = float(data[74])
				SDSS_raErr_given = float(data[75]) # in "
				SDSS_raErr_final = math.sqrt((SDSS_raErr_given*1.186)**2+0.01721**2) #final SDSS ra error in ", as done in Pal 5 paper
				SDSS_raErr_degrees = SDSS_raErr_final / 3600 / math.cos(math.radians(SDSS_dec))#converting " into degrees
				Distance = float(data[71]) #The distance to the nearest source
				XWIN_IMAGE = float(data[30]) #In pixels
				YWIN_IMAGE = float(data[31]) #In pixels
				ERRX2WIN_IMAGE = float(data[34]) #Pixels^2
				ERRY2WIN_IMAGE = float(data[35])
				transf_ra = float(data[69]) #Tobias's original transformation from LBT pixel position to RA/Dec
				transf_dec = float(data[70])
				LBT_adj_mag = float(data[95])
				LBT_flux = float(data[11])
				LBT_fluxerr = float(data[12])
				if LBT_flux > 0:
					LBT_orig_mag = -2.5*math.log10(LBT_flux)
				#the following lines calculate the +/- LBT magnitude errors
				if LBT_flux > 0 and (LBT_flux-LBT_fluxerr) > 0:
					zero_point = (LBT_adj_mag - LBT_orig_mag)
					lower_mag = (-2.5*math.log10(LBT_flux + LBT_fluxerr) + zero_point)
					higher_mag = (-2.5*math.log10(LBT_flux - LBT_fluxerr) + zero_point)
					LBT_magerr_minus = LBT_adj_mag - lower_mag
					LBT_magerr_plus = higher_mag - LBT_adj_mag
				else:
					LBT_magerr_plus = 9999
					LBT_magerr_minus = 9999

				SDSS_mag_g = float(data[80])
				SDSS_magerr_g = float(data[81])
				SDSS_mag_r = float(data[82])
				SDSS_magerr_r = float(data[83])
				SDSS_mag_i = float(data[84])
				SDSS_magerr_i = float(data[85])
				image_number = datafile[5:7]


				row_num = int(image_number) - 13
				

				#convert LBT transformed ra/dec positions to Alt/Az coordinates
				lst = all_lst[row_num]
				HA = lst - transf_ra #hour angle
				alt = math.degrees(math.asin(math.sin(math.radians(transf_dec)) * math.sin(math.radians(32.7014))
				 	+ math.cos(math.radians(transf_dec)) * math.cos(math.radians(32.7014)) * math.cos(math.radians(HA)))) # in degrees
				az = math.degrees(math.acos((math.sin(math.radians(transf_dec))-math.sin(math.radians(32.7014))*math.sin(math.radians(alt)))/(math.cos(math.radians(32.7014))*math.cos(math.radians(alt)))))  # in degrees

				color = SDSS_mag_r - SDSS_mag_i

				if abs(Spread_model) > 0.003: #if source is a galaxy
					if color > 1.1:
						color_g = 1.1
					elif color < 0.1:
						color_g = 0.1
					else:
						color_g = color
					beta_offset = 11.17-35.76*(color_g)*math.tan(math.radians(90-alt)) #mas

				else: #if source is a star, find closest color to actual color in starshift table and corresponding beta offset.
					clr = takeClosest(r_min_i, color)	
					row = np.where(r_min_i==clr)
					beta_offset = star_b_off[row] * math.tan(math.radians(90-alt)) #mas
				
				beta_off_deg = beta_offset / 3600 / 1000
				alt_prime = alt - beta_off_deg #applying dcr correction, all in degrees
				
				#Convert back into ra/dec
				corr_dec = math.degrees(math.asin(math.sin(math.radians(alt_prime))*math.sin(math.radians(32.7014))
						+ math.cos(math.radians(alt_prime))*math.cos(math.radians(32.7014))*math.cos(math.radians(az))))

				HA_prime = math.degrees(math.acos((math.sin(math.radians(alt_prime))-math.sin(math.radians(corr_dec))*math.sin(math.radians(32.7014)))
						 / (math.cos(math.radians(corr_dec))*math.cos(math.radians(32.7014)))))
				corr_ra = lst - HA_prime

				ra_diff = corr_ra - transf_ra #in degrees
				dec_diff = corr_dec - transf_dec
				
				x_diff = ra_diff / ra_conv[row_num] #in pixels
				y_diff = dec_diff / dec_conv[row_num] #in pixels

				if datafile[8] == '4':
					XWIN_IMAGE_prime = XWIN_IMAGE + y_diff
					YWIN_IMAGE_prime = YWIN_IMAGE + x_diff

				else:
					XWIN_IMAGE_prime = XWIN_IMAGE + x_diff
					YWIN_IMAGE_prime = YWIN_IMAGE + y_diff
				
				#Checking Tobias's calculations
				# if abs(Spread_model) < 0.003:
				# 	x_angle_calc = beta_off_deg * math.sin(math.radians(rot_ang[row_num])) / ra_conv[row_num] #Tobias's calculation
				# 	y_angle_calc = beta_off_deg * math.cos(math.radians(rot_ang[row_num])) / dec_conv[row_num]
				# 	x_rats.append(x_diff / x_angle_calc)
				# 	y_rats.append(y_diff / y_angle_calc)
				# 	starcolors.append(color)

				#Plotting
				# offsets.append(math.sqrt(x_diff**2 + y_diff**2))
				x_diffs.append(x_diff)
				y_diffs.append(y_diff)
				colors.append(color)




				# Now add sources to list if they have a matching SDSS object within 3".
				
 	 			if (SDSS_mag_r < 21.5  and SDSS_magerr_r < 0.2 and np.mean([LBT_magerr_minus, LBT_magerr_plus]) < .02 and 
 	 			(SDSS_mag_r - .3) < LBT_adj_mag < (SDSS_mag_r + .3) and abs(Spread_model) > 0.003 and probPSF_g == 0 and 
 	 			probPSF_r == 0 and probPSF_i == 0 and Distance < 3 and SDSS_decErr_given < .15 and SDSS_raErr_given < .15 
				and ERRX2WIN_IMAGE < math.exp(1.2*LBT_adj_mag-30)+.006 and ERRY2WIN_IMAGE < math.exp(1.2*LBT_adj_mag-30)+.006 
				and 55 < XWIN_IMAGE < 2094 and 5 < YWIN_IMAGE < 4604):
 					gal_list.append([XWIN_IMAGE_prime, YWIN_IMAGE_prime, SDSS_ra, SDSS_dec, SDSS_raErr_degrees, SDSS_decErr_degrees, float(image_number)])



 		Catalog.close()



plt.plot(colors,x_diffs, 'r.')
plt.plot(colors, y_diffs, 'g.')
plt.show()

# plt.plot(colors, offsets, 'b.')
# plt.show()

#following plots were used to compare rotator angle method with horizontal coordinates method.

# plt.plot(starcolors, x_rats, 'b.')
# plt.ylim(-1.03, -1.05)
# plt.show()

# plt.plot(starcolors, y_rats, 'r.')
# plt.show()



final_list = sorted(gal_list, key = operator.itemgetter(6)) #Puts the sources in order by image number, lowest to highest.



#Creating a file for each list and writing the list to it in nice format.

#****CHANGE PATH****

Gal_list_for_program = open("/home/martine/Documents/Astr_research/DCR_corrections/Galaxies_DCR_corrected/r_filter_%s.txt" %chip, 'w')
for galaxy in final_list:
	for aspect in galaxy:
		Gal_list_for_program.write(str(aspect) + " ") #writes each of the chosen properties of the galaxy separated by spaces
	Gal_list_for_program.write("\n") #moves to a new line before the next galaxy's properties are written
Gal_list_for_program.close()
