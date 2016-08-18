# Arranges lists of galaxies or all sources into the form necessary for IDL distortion correction code, then rearranges
# results into a readable form.
# Needs to be run independently for each band and chip combination, both for "galaxies" and "all objects" data sets (separately)
# For example, run it for g band chip 1 galaxies


import numpy as np

path = "All_sources_DCR_corrected/" #Change path (all sources vs just galaxies)

#Change filename for every chip+filter combo
dat = open(path + 'r_filter_4.txt', 'r')

data_lol = [] # is a list of lists (hence "lol"): a list of every line in dat, with every line being a list of the numbers inside

for line in dat:
  line = line.split() #Turns line into list by separating at spaces
  thisline = []
  for item in line:
    if item == '\n':
      pass
    else:
      thisline.append(float(item))
  data_lol.append(thisline)


#stat2 is a list in which the first elements (or first column if you think of it as an array) are
#consecutive numbers starting with the lowest image number for that filter and going until the highest. The
#second elements (second column) is how many objects are in that image. For example, for
#g filter image 11 there are always 0 because this image has no data.

stat2 = []

image_number = [galaxy[6] for galaxy in data_lol] #list of all image numbers

for i in range(int(max(image_number))-int(min(image_number))+1):
  stat2.append([min(image_number) + i, 0])

for element in stat2:
  q = 0
  for number in image_number:
    if number == element[0]:
      q+=1
  element[1] = q

# print stat2

x = 0


#the following makes a list called "stat," which is the same as stat2 except it has images with 0 objects removed
#(ex. g filter image 11 no longer appears)

for element in stat2:
  if element[1] == 0:
    x += 1

if x == 0:
  stat = stat2 #If there are no 0 object images, then stat is the same as stat2

else:
  stat = []
  for element in stat2:
    if element[1] != 0:
      stat.append(element)

# print stat



#makes "ar2", an array of zeros to be filled later

rows = 2*len(data_lol)
columns = 4
ar2 = np.zeros((rows, columns))


#The following rearranges the data to be used as input for IDL, puts it into ar2
#column 1: LBT x or y
#column 2: Ra or Dec (with constants subtracted)
#column 3: errRA or errDec (SDSS modified)
#column 4: image number
#rows first get filled in with x/ra/errRa/image number of image 1
#then after these, all objects' y/dec/errDec/image number of image 1
#continue for image 2, etc.

obj_amounts = [item[1] for item in stat] #list of all "number of objects in image" values

i = 1

for element in stat:
  c = 0
  for dat in data_lol:
    if dat[6] == element[0]:
      c += 1
      if i == 1:
        y = 0
      if i > 1:
        y = int(2*(sum(obj_amounts[0:(i-1)])))

 #subtracts 51 from Ra and 16 from dec to avoid a numeric problem in idl     

      ar2[c+y-1][0] = dat[0]
      ar2[c+y-1][1] = dat[2] - 151
      ar2[c+y-1][2] = dat[4]
      ar2[c+y-1][3] = dat[6]

      ar2[c+element[1]+y-1][0] = dat[1]
      ar2[c+element[1]+y-1][1] = dat[3] - 16
      ar2[c+element[1]+y-1][2] = dat[5]
      ar2[c+element[1]+y-1][3] = dat[6]
  i += 1


#Writes the galaxy or allob file to be used in IDL code
#Change path and filename
newfile = open('Galaxies_DCR_corrected/gal_r_c4_v3.txt', 'w')
for line in ar2:
	for number in line:
		newfile.write(str(number)+' ')
	newfile.write('\n')
newfile.close()

#Writes the length file to be used in IDL code
#Change path and filename
len_file = open("Lengths/Gal_lengths/len_gal_r_c4_v3.txt", 'w')
m = 0
for sublist in stat:
	len_file.write("len[%s]=%s" %(m,sublist[1]) + "\n" )
	m += 1
len_file.close()


#IDL runs here


#Arranges the results back into a readable form
#Change path and filename
res = open("Distortion Solution Outputs/res_allob_r_c4_v3dcr.txt", 'r')

nrows = len(ar2)

ncolumns = 5
ndata = np.zeros((nrows, ncolumns))

ndata[:, 0:4] = ar2



z = 0
for line in res:
	line = line.split()
	for item in line:
		ndata[z, 4] = item
		z += 1


ndatab = np.zeros((nrows/2 ,11))

w = 0

for element in stat:
	x = sum(obj_amounts[0:w])
	w += 1
	ndatab[x:x+element[1], 0] = ndata[2*x:2*x+element[1],0]
	ndatab[x:x+element[1], 2] = ndata[2*x:2*x+element[1],1] + 151
	ndatab[x:x+element[1], 4] = ndata[2*x:2*x+element[1],2]
	ndatab[x:x+element[1], 6] = ndata[2*x:2*x+element[1],3]
	ndatab[x:x+element[1], 8] = ndata[2*x:2*x+element[1],4] + 151

	ndatab[x:x+element[1],1] = ndata[2*x+element[1]:2*x+2*element[1],0]
	ndatab[x:x+element[1],3] = ndata[2*x+element[1]:2*x+2*element[1],1] + 16
	ndatab[x:x+element[1],5] = ndata[2*x+element[1]:2*x+2*element[1],2]
	ndatab[x:x+element[1],7] = ndata[2*x+element[1]:2*x+2*element[1],3]
	ndatab[x:x+element[1],9] = ndata[2*x+element[1]:2*x+2*element[1],4] + 16

ndatab = np.ndarray.tolist(ndatab) #Turns it into list

#Makes rearranged file of results
#Change path and filename
rearranged = open("Distortion Solution Outputs/fullresults_all_r_c4_v3dcr.txt", 'w')
for line in ndatab:
	for number in line:
		rearranged.write(str(number) + " ")
	rearranged.write("\n")
rearranged.close() 

