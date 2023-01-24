###################################################################################################
## MAKE MISCELLANEOUS PLOTS #######################################################################
## For Pleistocene paper 2020-02-03 ###############################################################
###################################################################################################

## CW's version 3 of original extra plot code, i.e. to make the multi-forcing figures, Which does the same as v2 but now for first pathway in linear factorisation i.e.
## 00000, 10000, 11000, 11100, 11110, 11111 (so 6 in total). Working (as of /23)

.libPaths(c("C:/Users/cw18831/OneDrive - University of Bristol/Documents/R/win-library/4.1", .libPaths()))

library(maps)
library(ggplot2)
library(viridis)
library(plyr)
library(plotrix)
library(gstat)
library(fields)
library(ncdf4)
library(mapproj)

require(colorspace)
require(graphics)

par_orig = par()

# Read in emulated data ----------------------------------------------------------------------------

# Set up data

print("Setting up")

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_insol_Laskar_jul_65N_800kyr_BP_cjrw.R", sep=""))

forcings_included = c("00000", "10000", "11000", "11100", "11110", "11111")

num_exps = length(forcings_included)


site_list = c("DomeC", "ODP982", "ODP722" ,"ODP1143" ,"ODP846", "ODP1090")

num_sites = length(site_list)


times = matrix(seq(-799, 0, 1))

num_times = length(times)


# load data

print("Loading input data (from emulator)")

data = array(0, c(num_times, num_exps, num_sites))
n = 1
m = 1

for (n in 1:num_sites){
  
  for (m in 1:num_exps){

    file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_E",forcings_included[m],"_",site_list[n],"_natural_0.8_to_0MyrBP.txt", sep="")
    data_all = read.table(file_nam, sep=" ", header=TRUE, numerals = "no.loss")

    data[,m,n] = data_all[,2] 
      
  }
}


single_forcing_paleodata_DomeC_800kyr_BP = matrix(data[,,1], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_DomeC_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_DomeC_800kyr_BP) <- c(times)

single_forcing_paleodata_ODP982_800kyr_BP = matrix(data[,,2], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_ODP982_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_ODP982_800kyr_BP) <- c(times)

single_forcing_paleodata_ODP722_800kyr_BP = matrix(data[,,3], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_ODP722_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_ODP722_800kyr_BP) <- c(times)

single_forcing_paleodata_ODP1143_800kyr_BP = matrix(data[,,4], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_ODP1143_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_ODP1143_800kyr_BP) <- c(times)

single_forcing_paleodata_ODP846_800kyr_BP = matrix(data[,,5], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_ODP846_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_ODP846_800kyr_BP) <- c(times)

single_forcing_paleodata_ODP1090_800kyr_BP = matrix(data[,,6], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_ODP1090_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_ODP1090_800kyr_BP) <- c(times)


print("Loading SST data")

# Read in proxy data and interpolate to 1 kyr intervals ----------------------------------------------------------------------------

# Read in SST data at ODP-982 for 800-0 kyr BP [Lawrence et al 2009]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_SST_ODP982_Lawrenceetal09_800kyr_BP_cjrw.R") # North Atlantic

data_interp = approx(data_SST_ODP982_Lawrenceetal09_800kyr_BP[,2],data_SST_ODP982_Lawrenceetal09_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_SST_ODP982_Lawrenceetal09_800kyr_BP[,2],data_SST_ODP982_Lawrenceetal09_800kyr_BP[,3], type="l")
# plot(data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp[,1],data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp[,2], type="l")


# Read in SST data at ODP-722 for 800-0 kyr BP [Herbert et al 2010]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_SST_ODP722_Herbertetal10_800kyr_BP_cjrw.R") # Arabian Sea

data_interp = approx(data_SST_ODP722_Herbertetal10_800kyr_BP[,2],data_SST_ODP722_Herbertetal10_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_SST_ODP722_Herbertetal10_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_SST_ODP722_Herbertetal10_800kyr_BP[,2],data_SST_ODP722_Herbertetal10_800kyr_BP[,3], type="l")
# plot(data_SST_ODP722_Herbertetal10_800kyr_BP_interp[,1],data_SST_ODP722_Herbertetal10_800kyr_BP_interp[,2], type="l")


# Read in SST data at ODP-1143 for 800-0 kyr BP [Li et al 2011]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_SST_ODP1143_Lietal11_800kyr_BP_cjrw.R") # South China Sea

data_interp = approx(data_SST_ODP1143_Lietal11_800kyr_BP[,2],data_SST_ODP1143_Lietal11_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_SST_ODP1143_Lietal11_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_SST_ODP1143_Lietal11_800kyr_BP[,2],data_SST_ODP1143_Lietal11_800kyr_BP[,3], type="l")
# plot(data_SST_ODP1143_Lietal11_800kyr_BP_interp[,1],data_SST_ODP1143_Lietal11_800kyr_BP_interp[,2], type="l")


# Read in SST data at ODP-846 for 800-0 kyr BP [Herbert et al 2010]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_SST_ODP846_Herbertetal10_800kyr_BP_cjrw.R") # Equatorial East Pacific

data_interp = approx(data_SST_ODP846_Herbertetal10_800kyr_BP[,2],data_SST_ODP846_Herbertetal10_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_SST_ODP846_Herbertetal10_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_SST_ODP846_Herbertetal10_800kyr_BP[,2],data_SST_ODP846_Herbertetal10_800kyr_BP[,3], type="l")
# plot(data_SST_ODP846_Herbertetal10_800kyr_BP_interp[,1],data_SST_ODP846_Herbertetal10_800kyr_BP_interp[,2], type="l")


# Read in SST data at ODP-1090 for 800-0 kyr BP [Martinez-Garcia et al 2010]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_cjrw.R") # Central Subarctic

data_interp = approx(data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,2],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,2],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,3], type="l")
# plot(data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,1],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,2], type="l")


# Read in Deuterium temperature data at Dome C for 800-0 kyr BP [Jouzel et al 2007]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_D_Jouzeletal07_800kyr_BP_cjrw.R") # Dome C + EPICA

data_interp = approx(data_D_Jouzeletal07_800kyr_BP[,2],data_D_Jouzeletal07_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_D_Jouzeletal07_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_D_Jouzeletal07_800kyr_BP[,2],data_D_Jouzeletal07_800kyr_BP[,3], type="l")
# plot(data_D_Jouzeletal07_800kyr_BP_interp[,1],data_D_Jouzeletal07_800kyr_BP_interp[,2], type="l")


# Read in SST stack data for 800-0 kyr BP [Shakun et al 2015]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_SST_stack_Shakunetal15_800kyr_BP_cjrw.R") # Global stack

data_interp = approx(data_SST_stack_Shakunetal15_800kyr_BP[,2],data_SST_stack_Shakunetal15_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_SST_stack_Shakunetal15_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_SST_stack_Shakunetal15_800kyr_BP[,2],data_SST_stack_Shakunetal15_800kyr_BP[,3], type="l")
# plot(data_SST_stack_Shakunetal15_800kyr_BP_interp[,1],data_SST_stack_Shakunetal15_800kyr_BP_interp[,2], type="l")



# Pre-industrial SST observations (HadISST) -----------------------------------------------------

# Load data

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_HadISST_SST_mon_1870_1900_BP_cjrw.R")

data_HadISST_SST_mon_1870_1900_BP_all = data_HadISST_SST_mon_1870_1900_BP_all[,2:361]

lats_HadISST = seq(-89.5,89.5,1)
lons_HadISST = seq(-179.5,179.5,1)


# Split data into separate months

mx = 180
my = 360
mmon = 12
myr = 31
me = mmon*myr

nrow_start = 2
nrow_end = mx+1

data_HadISST_SST_mon_1870_1900_BP = array(0,c(mx,my,me))
for (month in c(1:me)){
  data_HadISST_SST_mon_1870_1900_BP[,,month] = data_HadISST_SST_mon_1870_1900_BP_all[nrow_start:nrow_end,]
  data_HadISST_SST_mon_1870_1900_BP[,,month] = data_HadISST_SST_mon_1870_1900_BP[,,month]/100 # Convert to degC as currently degC*100
  nrow_start = nrow_start+mx+1
  nrow_end = nrow_end+mx+1
}


for (month in 1:me) { # Change land grid boxes (-327.68) to NA
  for (row in 1:mx) { 
    for (col in 1:my) {
      if (data_HadISST_SST_mon_1870_1900_BP[row,col,month] == -327.68) {
        data_HadISST_SST_mon_1870_1900_BP[row,col,month] = NA
      }
    }
  }
}


# Calculate annual averages for period

nyr_start = 1
nyr_end = 12

data_HadISST_SST_ann_1870_1900_BP = array(0,c(mx,my,myr))
for (year in 1:myr) {
  data_HadISST_SST_ann_1870_1900_BP[,,year] = apply(data_HadISST_SST_mon_1870_1900_BP[,,nyr_start:nyr_end], c(1,2), mean)
  nyr_start = nyr_start+mmon
  nyr_end = nyr_end+mmon
}


# Calculate average for entire period

data_HadISST_SST_av_1870_1900_BP = array(0,c(mx,my))
data_HadISST_SST_av_1870_1900_BP[,] = apply(data_HadISST_SST_ann_1870_1900_BP[,,], c(1,2), mean)


data_HadISST_SST_av_1870_1900_BP_flipped = t(data_HadISST_SST_av_1870_1900_BP)



# Extract proxy data climate for each site ----------------------------------------------------------------------------

# Select site data

print("Extracting locations")

lats_HadISST_ODP982 = 33 # ODP-982, North Atlantic
lons_HadISST_ODP982 = 165
lats_HadISST_ODP722 = 74 # ODP-722, Arabian Sea
lons_HadISST_ODP722 = 240
lats_HadISST_ODP1143 = 81 # ODP-1143, South China Sea
lons_HadISST_ODP1143 = 294
lats_HadISST_ODP846 = 94 # ODP-846, Equatorial East Pacific
lons_HadISST_ODP846 = 90
lats_HadISST_ODP1090 = 133 # ODP-1090, Central Subantarctic
lons_HadISST_ODP1090 = 189

data_site_982_HadISST_SST_av_1870_1900_BP = data_HadISST_SST_av_1870_1900_BP_flipped[lons_HadISST_ODP982,lats_HadISST_ODP982] # ODP-982, North Atlantic
data_site_722_HadISST_SST_av_1870_1900_BP = data_HadISST_SST_av_1870_1900_BP_flipped[lons_HadISST_ODP722,lats_HadISST_ODP722] # ODP-722, Arabian Sea
data_site_1143_HadISST_SST_av_1870_1900_BP = data_HadISST_SST_av_1870_1900_BP_flipped[lons_HadISST_ODP1143,lats_HadISST_ODP1143] # ODP-1143, South China Sea
data_site_846_HadISST_SST_av_1870_1900_BP = data_HadISST_SST_av_1870_1900_BP_flipped[lons_HadISST_ODP846,lats_HadISST_ODP846] # ODP-846, Equatorial East Pacific
data_site_1090_HadISST_SST_av_1870_1900_BP = data_HadISST_SST_av_1870_1900_BP_flipped[lons_HadISST_ODP1090,lats_HadISST_ODP1090] # ODP-1090, Central Subarctic


# image.plot(lats_HadISST, lons_HadISST, data_HadISST_SST_av_1870_1900_BP[,], zlim=c(-30,30), xlab="Longitude", ylab="Latitude", main="HadISST observations 1870-1900")
# 
# image.plot(lons_HadISST, lats_HadISST, data_HadISST_SST_av_1870_1900_BP_flipped[,c(180:1)], zlim=c(-30,30), xlab="Longitude", ylab="Latitude", main="HadISST observations 1870-1900")
# map(add=T, interior=F, col="black")


# Plot map showing all sites

data_sites_HadISST_SST_av_1870_1900_BP = matrix(0,nrow=360, ncol=180)

data_sites_HadISST_SST_av_1870_1900_BP[lons_HadISST_ODP982,lats_HadISST_ODP982] = 5 # ODP 982, North Atlantic
data_sites_HadISST_SST_av_1870_1900_BP[lons_HadISST_ODP722,lats_HadISST_ODP722] = 5 # ODP-722, Arabian Sea
data_sites_HadISST_SST_av_1870_1900_BP[lons_HadISST_ODP1143,lats_HadISST_ODP1143] = 5 # ODP-1143, South China Sea
data_sites_HadISST_SST_av_1870_1900_BP[lons_HadISST_ODP846,lats_HadISST_ODP846] = 5 # ODP-846, Equatorial East Pacific
data_sites_HadISST_SST_av_1870_1900_BP[lons_HadISST_ODP1090,lats_HadISST_ODP1090] = 5 # ODP-1090, Central Subarctic


col <- colorRampPalette(c("white","red"))(5)

image.plot(lons_HadISST, lats_HadISST, data_sites_HadISST_SST_av_1870_1900_BP[,c(180:1)], zlim=c(0,5), xlab="Longitude", ylab="Latitude", main="HadISST observations 1870-1900",col=col)
map(add=T, interior=F, col="black")



# Calculate anomaly between proxy data and PI in proxy data -----------------------------------------------------
# and proxy data and PI SST observations (HadISST) ---

print("Calculating anomalies")

data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP = array(c(data_SST_ODP982_Lawrenceetal09_800kyr_BP[,2],data_SST_ODP982_Lawrenceetal09_800kyr_BP[,3],data_SST_ODP982_Lawrenceetal09_800kyr_BP[,3],data_SST_ODP982_Lawrenceetal09_800kyr_BP[,3]),c(dim(data_SST_ODP982_Lawrenceetal09_800kyr_BP)[1],4))
data_SST_anom_ODP722_Herbertetal10_800kyr_BP = array(c(data_SST_ODP722_Herbertetal10_800kyr_BP[,2],data_SST_ODP722_Herbertetal10_800kyr_BP[,3],data_SST_ODP722_Herbertetal10_800kyr_BP[,3],data_SST_ODP722_Herbertetal10_800kyr_BP[,3]),c(dim(data_SST_ODP722_Herbertetal10_800kyr_BP)[1],4))
data_SST_anom_ODP1143_Lietal11_800kyr_BP = array(c(data_SST_ODP1143_Lietal11_800kyr_BP[,2],data_SST_ODP1143_Lietal11_800kyr_BP[,3],data_SST_ODP1143_Lietal11_800kyr_BP[,3],data_SST_ODP1143_Lietal11_800kyr_BP[,3]),c(dim(data_SST_ODP1143_Lietal11_800kyr_BP)[1],4))
data_SST_anom_ODP846_Herbertetal10_800kyr_BP = array(c(data_SST_ODP846_Herbertetal10_800kyr_BP[,2],data_SST_ODP846_Herbertetal10_800kyr_BP[,3],data_SST_ODP846_Herbertetal10_800kyr_BP[,3],data_SST_ODP846_Herbertetal10_800kyr_BP[,3]),c(dim(data_SST_ODP846_Herbertetal10_800kyr_BP)[1],4))
data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP = array(c(data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,2],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,3],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,3],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP[,3]),c(dim(data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP)[1],4))

data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp = array(c(data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp[,1],data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp[,2],data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp[,2],data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp[,2]),c(dim(data_SST_ODP982_Lawrenceetal09_800kyr_BP_interp)[1],4))
data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp = array(c(data_SST_ODP722_Herbertetal10_800kyr_BP_interp[,1],data_SST_ODP722_Herbertetal10_800kyr_BP_interp[,2],data_SST_ODP722_Herbertetal10_800kyr_BP_interp[,2],data_SST_ODP722_Herbertetal10_800kyr_BP_interp[,2]),c(dim(data_SST_ODP722_Herbertetal10_800kyr_BP_interp)[1],4))
data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp = array(c(data_SST_ODP1143_Lietal11_800kyr_BP_interp[,1],data_SST_ODP1143_Lietal11_800kyr_BP_interp[,2],data_SST_ODP1143_Lietal11_800kyr_BP_interp[,2],data_SST_ODP1143_Lietal11_800kyr_BP_interp[,2]),c(dim(data_SST_ODP1143_Lietal11_800kyr_BP_interp)[1],4))
data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp = array(c(data_SST_ODP846_Herbertetal10_800kyr_BP_interp[,1],data_SST_ODP846_Herbertetal10_800kyr_BP_interp[,2],data_SST_ODP846_Herbertetal10_800kyr_BP_interp[,2],data_SST_ODP846_Herbertetal10_800kyr_BP_interp[,2]),c(dim(data_SST_ODP846_Herbertetal10_800kyr_BP_interp)[1],4))
data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp = array(c(data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,1],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,2],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,2],data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,2]),c(dim(data_SST_ODP1090_MartinezGarciaetal10_800kyr_BP_interp)[1],4))


data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP[,3] = data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP[,2] - data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP[1,2]
data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP[,4] = data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP[,2] - data_site_982_HadISST_SST_av_1870_1900_BP
data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp[,3] = data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp[,2] - data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp[1,2]
data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp[,4] = data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp[,2] - data_site_982_HadISST_SST_av_1870_1900_BP

data_SST_anom_ODP722_Herbertetal10_800kyr_BP[,3] = data_SST_anom_ODP722_Herbertetal10_800kyr_BP[,2] - data_SST_anom_ODP722_Herbertetal10_800kyr_BP[1,2]
data_SST_anom_ODP722_Herbertetal10_800kyr_BP[,4] = data_SST_anom_ODP722_Herbertetal10_800kyr_BP[,2] - data_site_722_HadISST_SST_av_1870_1900_BP
data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp[,3] = data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp[,2] - data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp[9,2]
data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp[,4] = data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp[,2] - data_site_722_HadISST_SST_av_1870_1900_BP

data_SST_anom_ODP1143_Lietal11_800kyr_BP[,3] = data_SST_anom_ODP1143_Lietal11_800kyr_BP[,2] - data_SST_anom_ODP1143_Lietal11_800kyr_BP[1,2]
data_SST_anom_ODP1143_Lietal11_800kyr_BP[,4] = data_SST_anom_ODP1143_Lietal11_800kyr_BP[,2] - data_site_1143_HadISST_SST_av_1870_1900_BP
data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp[,3] = data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp[,2] - data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp[7,2]
data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp[,4] = data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp[,2] - data_site_1143_HadISST_SST_av_1870_1900_BP

data_SST_anom_ODP846_Herbertetal10_800kyr_BP[,3] = data_SST_anom_ODP846_Herbertetal10_800kyr_BP[,2] - data_SST_anom_ODP846_Herbertetal10_800kyr_BP[1,2]
data_SST_anom_ODP846_Herbertetal10_800kyr_BP[,4] = data_SST_anom_ODP846_Herbertetal10_800kyr_BP[,2] - data_site_846_HadISST_SST_av_1870_1900_BP
data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp[,3] = data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp[,2] - data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp[7,2]
data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp[,4] = data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp[,2] - data_site_846_HadISST_SST_av_1870_1900_BP

data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP[,3] = data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP[,2] - data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP[1,2]
data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP[,4] = data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP[,2] - data_site_1090_HadISST_SST_av_1870_1900_BP
data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,3] = data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,2] - data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[1,2]
data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,4] = data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,2] - data_site_1090_HadISST_SST_av_1870_1900_BP

########################################################################################################################################################
## FIGURES #############################################################################################################################################
########################################################################################################################################################

print("Creating figures")

SST_anom_PI_opt = 3 # Anomaly for proxy data from PI proxy data SST

SST_anom_PI_opt_nam = "PIprox"

# Set some constants, not including y-axes because these will change

fs = 1 # Size of font
fs2 = 0.5 # Size of font
fs3 = 0.75 # Size of font
fs4 = 1.75 # Size of font
lns = 1.5 # Width of lines
lns2 = 5 # Width of lines
pts = 0.75 # Size of points

lb=0.01 # Left border
rb=0.985 # Right border
bb=0.04 # Bottom border
ht=(1-bb)/1 # Height of panel

col_proxy = "black"
col_sim = c("green", "blue", "purple", "pink", "red", "black")

# Dome C, Antarctica

xlimits = c(0, 800)
ylimits = c(-12,8)

xticks = seq(0, 800, by = 50)
yticks=seq(-12, 8, by=1)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels = c("-12", "", "-10", "", "-8", "", "-6", "", "-4", "", "-2", "", "0", "", "2", "", "4", "", "6", "", "8")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/SAT_linfac_DomeC_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(data_D_Jouzeletal07_800kyr_BP_interp[,1]),data_D_Jouzeletal07_800kyr_BP_interp[,2], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_DomeC_800kyr_BP[,m],type="l", cex=lns, col=col_sim[m])
}
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(2, at = yticks, labels=ylabels, tck=0.01, lwd=lns, cex.axis=fs)
#text((xlimits[1]+10), (((ylimits[2]-ylimits[1])*0.96)+ylimits[1]), "a)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Temperature anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("("~ degree*C*")", side=2, line=2.25, font=1, cex=fs)

legend((xlimits[2]-267), (((ylimits[2]-ylimits[1])*1.08)+ylimits[1]), c("Proxy record", "E00000", "E10000", "E11000", "E11100", "E11110", "E11111"), pch=NA, lty=c(2,1,1,1,1,1,1), lwd=lns, cex=fs3, pt.cex=pts, adj = 0.05, col=c(col_proxy, col_sim), ncol=2, bty="n")

par(las=1)

dev.off()

# ODP-982, North Atlantic

xlimits = c(0, 800)
ylimits=c(-6,4)

xticks = seq(0, 800, by = 50)
yticks=seq(-6, 4, by = 1)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels=c("-6", "", "-4", "", "-2", "", "0", "", "2", "", "4")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/SST_linfac_ODP982_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp[,1]),data_SST_anom_ODP982_Lawrenceetal09_800kyr_BP_interp[,SST_anom_PI_opt], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_ODP982_800kyr_BP[,m],type="l", cex=lns, col=col_sim[m])
}
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(2, at = yticks, labels=ylabels, tck=0.01, lwd=lns, cex.axis=fs)
#text((xlimits[1]+10), (((ylimits[2]-ylimits[1])*0.96)+ylimits[1]), "b)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Temperature anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("("~ degree*C*")", side=2, line=2.25, font=1, cex=fs)

par(las=1)

dev.off()

# ODP-722, Arabian Sea

xlimits = c(0, 800)
ylimits=c(-4,2)

xticks = seq(0, 800, by = 50)
yticks=seq(-4, 2, by = 1)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels=c("-4", "", "-2", "", "0", "", "2")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/SsT_linfac_ODP722_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp[,1]),data_SST_anom_ODP722_Herbertetal10_800kyr_BP_interp[,SST_anom_PI_opt], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_ODP722_800kyr_BP[,m],type="l", cex=lns, col=col_sim[m])
}
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(2, at = yticks, labels=ylabels, tck=0.01, lwd=lns, cex.axis=fs)
#text((xlimits[1]+10), (((ylimits[2]-ylimits[1])*0.96)+ylimits[1]), "c)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Temperature anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("("~ degree*C*")", side=2, line=2.25, font=1, cex=fs)

par(las=1)

dev.off()

# ODP-1143, South China Sea

xlimits = c(0, 800)
ylimits=c(-4,2)

xticks = seq(0, 800, by = 50)
yticks=seq(-4, 2, by = 1)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels=c("-4", "", "-2", "", "0", "", "2")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/SST_linfac_ODP1143_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp[,1]),data_SST_anom_ODP1143_Lietal11_800kyr_BP_interp[,SST_anom_PI_opt], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_ODP1143_800kyr_BP[,m],type="l", cex=lns, col=col_sim[m])
}
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(2, at = yticks, labels=ylabels, tck=0.01, lwd=lns, cex.axis=fs)
#text((xlimits[1]+10), (((ylimits[2]-ylimits[1])*0.96)+ylimits[1]), "d)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Temperature anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("("~ degree*C*")", side=2, line=2.25, font=1, cex=fs)

par(las=1)

dev.off()

# ODP-846, Equatorial Atlantic

xlimits = c(0, 800)
ylimits=c(-4,2)

xticks = seq(0, 800, by = 50)
yticks=seq(-4, 2, by = 1)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels=c("-4", "", "-2", "", "0", "", "2")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/SST_linfac_ODP846_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp[,1]),data_SST_anom_ODP846_Herbertetal10_800kyr_BP_interp[,SST_anom_PI_opt], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_ODP846_800kyr_BP[,m],type="l", cex=lns, col=col_sim[m])
}
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(2, at = yticks, labels=ylabels, tck=0.01, lwd=lns, cex.axis=fs)
#text((xlimits[1]+10), (((ylimits[2]-ylimits[1])*0.96)+ylimits[1]), "e)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Temperature anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("("~ degree*C*")", side=2, line=2.25, font=1, cex=fs)

par(las=1)

dev.off()

# ODP-1090, Central Subantarctic

xlimits = c(0, 800)
ylimits = c(-8, 6)

xticks = seq(0, 800, by = 50)
yticks=seq(-8, 6, by=1)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels = c("-8", "", "-6", "", "-4", "", "-2", "", "0", "", "2", "", "4", "", "6")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/SST_linfac_ODP1090_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,1]),data_SST_anom_ODP1090_MartinezGarciaetal10_800kyr_BP_interp[,SST_anom_PI_opt], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_ODP1090_800kyr_BP[,m],type="l", cex=lns, col=col_sim[m])
}
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")
.
axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(2, at = yticks, labels=ylabels, tck=0.01, lwd=lns, cex.axis=fs)
#text((xlimits[1]+10), (((ylimits[2]-ylimits[1])*0.96)+ylimits[1]), "f)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Temperature anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("("~ degree*C*")", side=2, line=2.25, font=1, cex=fs)

par(las=1)

dev.off()


print("Finished")

#########################################################################################