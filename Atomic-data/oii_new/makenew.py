#!/usr/bin/python

# script to get updated recombination coefficients for the 208 OII lines we have in Roii.dat
# the new data file is at http://cdsarc.u-strasbg.fr/ftp/VI/150/DataFiles/OIIlines_ABC
# LSBC94 provided coefficients for a function fitted to the data for each line. Storey et al 2017 just give the calculated coefficients
# NEAT will do bilinear interpolation in log(Te),log(ne) to get the data

# read the old data

with open("../Roii.dat") as f:
    content = f.readlines()

del content[0] # remove comment row

wlens_old=[]
data_old=[]

for line in content:
    wlens_old.append(float(line[6:15]))
    data_old.append(line)


# read the new data

with open("OIIlines_ABC") as f:
    content=f.readlines()

wlens_new=[]
data_new=[]

newlinelist=content[18:8907] # this is the index
newlinedata=content[8908:] # this is the data

for line in newlinelist:
    wlens_new.append(float(line[92:102]))
    data_new.append(line)


# for the new data

newdata=[]

# cross correlate. find nearest wavelength and check if it's the right transition

for wlen in wlens_old:
    wlen_new=min(wlens_new, key=lambda x:abs(x-wlen))

    iold=wlens_old.index(wlen)
    inew=wlens_new.index(wlen_new)

    term_old=data_old[iold][77:79]+data_old[iold][108:110]
    term_new=data_new[inew][81:83]+data_new[inew][62:64]

# if term_old = term_new then it's the right line.
# inew=lineid-1
# block of data starts at inew*84+28 (case B)
# copy the multiplet name, terms and stat weights

    if term_old==term_new:
        newlinedata[inew*84+28]=newlinedata[inew*84+28].strip()+" "+data_old[iold].strip()+"\n"
        newdata.extend(newlinedata[inew*84+28:(inew*84)+55])

# save the new data

with open('Roii_storey2017.dat', 'w') as f:
    for item in newdata:
# fix some wavelengths for backward compatibility when writing
        f.write("%s" % item.replace("4283.72","4283.73").replace("4676.23","4676.24").replace("4285.69","4285.68").replace("4491.22","4491.23"))
