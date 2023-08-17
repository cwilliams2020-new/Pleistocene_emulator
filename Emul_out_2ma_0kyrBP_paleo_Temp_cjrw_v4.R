###################################################################################################
## USE TWO EMULATORS TO PREDICT PAST TEMP FOR 800-0 kyr BP AT VARIOUS SITES #######################
## modice+tdab emulator (gl) AND modlowice emulator (ig), calibrated on temp ######################
## Using palaeo-proxy CO2 data [Luthi et al 2008] #################################################
## For Pleistocene paper 2020-02-03 ###############################################################
###################################################################################################

## CW's version 4 of original palaeo-emulator, Which is the same as version 3 but now removes all the reading in of proxy data and plotting (because this is now
## done in cjrw3*) and includes writing out of global variances
## Working (as of 17/8/23)

.libPaths(c("C:/Users/cw18831/OneDrive - University of Bristol/Documents/R/win-library/4.1", .libPaths()))

# Need to run this once the first time using the emulator:
#install.packages("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/Packages/GP_0.1.4/GP",repos = NULL, type="source")
# Will not be able to load plotrix library, so separately (at command line), do install.packages("plotrix")
#install.packages("plotrix")

library(maps)
library(ggplot2)
library(viridis)
library(plyr)
library(plotrix)
library(gstat)
library(fields)

require(colorspace)
require(graphics)
require(corrplot)

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_parula.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_bwr.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_wr.R')

par_orig = par()

## Load and calibrate emulator #########################################################################################

print("Loading emulator")

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/Loading/Emul_in_modice_tdab_+_modlowice_LT2000ppm_optd_on_Temp_sim_Temp_cjrw.R')

## Set up #########################################################################################

# CW: Because there is only one forcing experiment for this, no need to have forcings_included in either input or output files

# Load data

print("Loading input data for emulator")

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_forcing_paleodata_2Ma_BP_cjrw.R", sep=""))
source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_Rnum_paleodata_Focs_800kyr_BP_cjrw.R", sep="")) # Leave here, but not actually used apart from plotting

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/Data - Climate NC/Tdst_temp_mm_1_5m_cjrw.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/Data - Climate NC/Tdst_temp_mm_uo_cjrw.R')

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/GSL_CGSLM_120kyr_BP_cjrw.R", sep="")) # Leave here, but not actually used apart from plotting

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_insol_Laskar_jul_65N_800kyr_BP_cjrw.R", sep="")) # Leave here, but not actually used apart from plotting

# Create 2 lots of time arrays, one for 2.58 Ma and one for 800 kya
  
times = matrix(seq(-2579, 0, 1))
times2 = matrix(seq(-799, 0, 1))

K = 273.15


# Run emulator with LGM conditions

print("Running emulator")
Sys.time()
check1 = Sys.time()

paleoclimate_LGM_time = times[which(times == -21)]
paleoclimate_LGM_pos = which(times == -21)

model_input_lgm = matrix(rep(1,5),c(1,5),dimnames = list(NULL, colnames(model_input_paleodata_800kyr_BP)))
model_input_lgm[1] = model_input_paleodata_800kyr_BP[paleoclimate_LGM_pos,1] # CO2 (190.9) at -20.96 kyr [Bereiter et al 2014]
model_input_lgm[2] = model_input_paleodata_800kyr_BP[paleoclimate_LGM_pos,2] # Orbital (22.96) at -21 kyr [Laskar et al 2004]
model_input_lgm[3] = model_input_paleodata_800kyr_BP[paleoclimate_LGM_pos,3] # Orbital (0.017037) at -21 kyr [Laskar et al 2004]
model_input_lgm[4] = model_input_paleodata_800kyr_BP[paleoclimate_LGM_pos,4] # Orbital (-0.008031) at -21 kyr [Laskar et al 2004]
model_input_lgm[5] = model_input_paleodata_800kyr_BP[paleoclimate_LGM_pos,5] #c(-129.481) # GSL at -21 kyr [Spratt + Lisiecki]


Scenario_paleodata <- data.frame(model_input_paleodata_800kyr_BP)
Rnum_rcp <- data.frame(data_input_Rnum_paleodata_800kyr_BP)

Scenario_lgm <- data.frame(model_input_lgm)



# set CO2 to be log (natural) of CO2
Scenario_paleodata[,1] = log(Scenario_paleodata[,1])

Scenario_lgm[1] = log(Scenario_lgm[1])


# rescaled to make it compatible with the 'X' used to calibrate the emulator

XScenario_paleodata <- sweep(Scenario_paleodata,  2, attr(X_all, "scaled:center"), '-')

XScenario_paleodata <- sweep(XScenario_paleodata, 2, attr(X_all, "scaled:scale"), '/')

XScenario_lgm <- sweep(Scenario_lgm,  2, attr(X_all, "scaled:center"), '-')

XScenario_lgm <- sweep(XScenario_lgm, 2, attr(X_all, "scaled:scale"), '/')


## Run emulator using input data for last 800 kyr #################################################

# run emulator (feed every row of XScenario into emulator, put this in a list, 
# and repack this into an array.

OUT_gl_paleodata = lapply(seq(nrow(XScenario_paleodata)), function(i) pe_p((XScenario_paleodata[i,, drop=FALSE]) , E_gl ) )

OUT_ig_paleodata = lapply(seq(nrow(XScenario_paleodata)), function(i) pe_p((XScenario_paleodata[i,, drop=FALSE]) , E_ig ) )

OUT_gl_lgm = lapply(seq(nrow(XScenario_lgm)), function(i) pe_p((XScenario_lgm[i,, drop=FALSE]) , E_gl ) )

#rm(list=ls()[(ls() %in% c('E_gl', 'E_ig'))]) # CW: Not entirely sure why these are being removed (other than to save space), but this needs to be commented out when running all sensitivity experiments together


# extract mean and variance

scenario_means_gl_paleodata_anom <- simplify2array( lapply(OUT_gl_paleodata, function(i) i$mean))

scenario_means_ig_paleodata_anom <- simplify2array( lapply(OUT_ig_paleodata, function(i) i$mean))

scenario_var_gl_paleodata_anom <- simplify2array( lapply(OUT_gl_paleodata, function(i) i$var))

scenario_var_ig_paleodata_anom <- simplify2array( lapply(OUT_ig_paleodata, function(i) i$var))

scenario_means_gl_lgm_anom <- simplify2array( lapply(OUT_gl_lgm, function(i) i$mean))

rm(list=ls()[(ls() %in% c('OUT_gl_paleodata'))])
rm(list=ls()[(ls() %in% c('OUT_ig_paleodata'))])


# attach times as an attribute

attr(scenario_means_gl_paleodata_anom, 'times') <- times

attr(scenario_means_ig_paleodata_anom, 'times') <- times

attr(scenario_var_gl_paleodata_anom, 'times') <- times

attr(scenario_var_ig_paleodata_anom, 'times') <- times


# Extract data depending on glacial state to create full timeseries

glacial_state_paleodata = array(0,c(length(Scenario_paleodata[,1]),1))

for (n in 1:length(Scenario_paleodata[,1])){
  if (Scenario_paleodata[n,5] < 0){ # GLACIAL
    glacial_state_paleodata[n,1] = 0
  } else { # INTERGLACIAL
    glacial_state_paleodata[n,1] = 1
  }
}


scenario_means_paleodata_anom = array(0,c(dim(scenario_means_gl_paleodata_anom)[1],dim(scenario_means_gl_paleodata_anom)[2],dim(scenario_means_gl_paleodata_anom)[3]))

scenario_var_paleodata_anom = array(0,c(dim(scenario_var_gl_paleodata_anom)[1],dim(scenario_var_gl_paleodata_anom)[2],dim(scenario_var_gl_paleodata_anom)[3]))

for (n in 1:dim(scenario_means_gl_paleodata_anom)[3]){
  if (glacial_state_paleodata[n,1] == 0){ # GLACIAL
    #  message(n, " Glacial")
    scenario_means_paleodata_anom[,,n] = scenario_means_gl_paleodata_anom[,,n]
    scenario_var_paleodata_anom[,,n] = scenario_var_gl_paleodata_anom[,,n]
  } else { # INTERGLACIAL
    #  message(n, " Interglacial")
    scenario_means_paleodata_anom[,,n] = scenario_means_ig_paleodata_anom[,,n]
    scenario_var_paleodata_anom[,,n] = scenario_var_ig_paleodata_anom[,,n]
  }
}

# CW: For unknown reasons, possibly because of the removal line (123) or possibly because of call to Emulator_remove_vars2.R, various colour bars need to be sourced again
# as otherwise it does not recognise bwr function

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_parula.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_bwr.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_wr.R')

col_bwr = rgb(bwr(20))


# Use pythag. to calculate eccentricity value (hypotenuse)

model_input_paleodata_800kyr_BP_ep = model_input_paleodata_800kyr_BP
colnames(model_input_paleodata_800kyr_BP_ep, do.NULL = TRUE, prefix = "col")
colnames(model_input_paleodata_800kyr_BP_ep) = c("co2","obliquity","eccentricity","l. of perihelion","ice")

for (i in 1:length(model_input_paleodata_800kyr_BP_ep[,1])) {
  model_input_paleodata_800kyr_BP_ep[i,3] = sqrt((model_input_paleodata_800kyr_BP[i,3]^2) + (model_input_paleodata_800kyr_BP[i,4]^2)) # Calculate Ecc
  model_input_paleodata_800kyr_BP_ep[i,4] = atan2(model_input_paleodata_800kyr_BP[i,3]*pi/180,model_input_paleodata_800kyr_BP[i,4]*pi/180) # Calculate Per
  model_input_paleodata_800kyr_BP_ep[i,4] = model_input_paleodata_800kyr_BP_ep[i,4]*180/pi
}

for (i in 1:length(model_input_paleodata_800kyr_BP_ep[,1])) {
  if (model_input_paleodata_800kyr_BP_ep[i,4] < 0) {
    model_input_paleodata_800kyr_BP_ep[i,4] = model_input_paleodata_800kyr_BP_ep[i,4] + 360; # Correct Per to be between 0 and 360 (from -180 to 180)
  }
}



## Format results from emulator #####################################################################

# Convert SAT to Deg C

model_output_tdst = model_output_tdst - K


# Calculate true global mean annual temperature evolution for next 200 kyr ------------------------
# i.e grid box area weighted mean

nx = 96 # lon
ny = 73 # lat
nx2 = nx/2
ne = length(times) # timeslice exps
nc = 5 # control exps

lons <- array(0:(nx-1))*3.75 # lons for atmospheric data
lons[1:(nx+1)] <- c(lons[(nx2+1):nx]-360,lons[1:nx2],lons[nx2+1]) # lons for atmospheric data

lats <- array(0:(ny-1))*2.5-90 # lats for atmospheric data


lons_o <- array(0:(288-1))*1.25 # lons for ocean data
lons_o[1:(288+1)] <- c(lons_o[(144+1):288]-360,lons_o[1:144],lons_o[144+1]) # lons for ocean data

lats_o <- array(0:(144-1))*1.25-89.375 # lats for ocean data (not wrong was round!)


gb_lat_length = 2.5 # degrees
gb_lat_length_pole = gb_lat_length/2 # degrees
gb_lon_length = 3.75 # degrees

lats_true = array(0,c((ny+1),1))
lats_true[1,] = lats[1];
lats_true[ny+1,] = lats[ny];
lats_true[2,] = lats_true[1,]+gb_lat_length_pole;
lats_true[ny,] = lats_true[ny+1,]-gb_lat_length_pole;


for (ncol in 3:(ny-1)){
  lats_true[ncol,] = lats_true[ncol-1,]+gb_lat_length;
}


lons_true = array(0,c((nx+2),1))
lons_true[1,] = lons[1];
lons_true[2,] = lons[1]+(gb_lon_length/2);
lons_true[nx+2,] = lons[nx+1];
lons_true[nx+1,] = lons[nx+1]-(gb_lon_length/2);

for (nrow in 3:nx){
  lons_true[nrow,] = lons_true[nrow-1,]+gb_lon_length;
}


# Convert latitude and longitude to radians

lats_true_rad = lats_true*(pi/180);

lons_true_rad = lons_true*(pi/180);


# Calculate grid box areas manually

gb_area_manual_split = array(0,c((nx+1),ny))
gb_area_manual = array(0,c(nx,ny))

for (ncol in 1:ny){
  for (nrow in 1:nx+1){
    gb_area_manual_split[nrow,ncol] = abs(sin(lats_true_rad[ncol,])-sin(lats_true_rad[ncol+1,]))*abs(lons_true_rad[nrow,]-lons_true_rad[nrow+1,])/(4*pi)
  }
}

gb_area_manual_split[1,] = gb_area_manual_split[(nx+1),]

gb_area_manual_total = sum(sum(gb_area_manual_split))

gb_area_manual[1,] = gb_area_manual_split[1,]+gb_area_manual_split[nx+1,]
gb_area_manual[2:nx,] = gb_area_manual_split[2:nx,]


# Calculate true global mean annual temperature (anom compared to tdstb)

my_tim_paleodata_e_truemean = array(0,c(ne,1))

my_tim_lgm_e_truemean = array(0,c(ne,1))

for (y in 1:ne){
  my_tim_paleodata_e_truemean[y,] = sum(sum((scenario_means_paleodata_anom[,,y]*gb_area_manual[,])))
}

my_tim_lgm_e_truemean = sum(sum((scenario_means_gl_lgm_anom[,,1]*gb_area_manual[,])))

my_tim_tdstb_e_truemean = sum(sum((model_output_tdst[,,1]*gb_area_manual[,])))


# Calculate non-weighted global mean annual temperature (anom compared to tdstb)

my_tim_paleodata_e_mean = apply(scenario_means_paleodata_anom, c(3), mean)


# Convert data to correct format for mapping
# CW: This changes array to be firstly Eurocentric (i.e. Europe/Africa in the middle) and secondly to be slightly larger (hence "big") in longitudes/rows (97 instead of 96)
# This is where the anomalies come from, so need to use "big" arrays

scenario_means_paleodata_big_anom = array(0,c((nx+1),ny,ne))
for (y in c(1:ne))
  for (x in c(1:ny))
    scenario_means_paleodata_big_anom[1:(nx+1),x,y] = c(scenario_means_paleodata_anom[(nx2+1):nx,x,y],scenario_means_paleodata_anom[1:nx2,x,y],scenario_means_paleodata_anom[(nx2+1),x,y])

scenario_means_lgm_big_anom = array(0,c((nx+1),ny,1))
for (y in c(1:1))
  for (x in c(1:ny))
    scenario_means_lgm_big_anom[1:(nx+1),x,y] = c(scenario_means_gl_lgm_anom[(nx2+1):nx,x,y],scenario_means_gl_lgm_anom[1:nx2,x,y],scenario_means_gl_lgm_anom[(nx2+1),x,y])


model_output_tdst_big = array(0,c((nx+1),ny,nc))
for (y in c(1:nc))
  for (x in c(1:ny))
    model_output_tdst_big[1:(nx+1),x,y] = c(model_output_tdst[(nx2+1):nx,x,y],model_output_tdst[1:nx2,x,y],model_output_tdst[(nx2+1),x,y])

model_output_tdst_SST_big = array(0,c((288+1),144,nc))
for (y in c(1:nc))
  for (x in c(1:144))
    model_output_tdst_SST_big[1:(288+1),x,y] = c(model_output_tdst_SST[(144+1):288,x,y],model_output_tdst_SST[1:144,x,y],model_output_tdst_SST[(144+1),x,y])

model_output_modice_tdab_big = array(0,c((nx+1),ny,dim(model_input_modice_tdab)[1]))
for (y in c(1:dim(model_input_modice_tdab)[1]))
  for (x in c(1:ny))
    model_output_modice_tdab_big[1:(nx+1),x,y] = c(model_output_modice_tdab[(nx2+1):nx,x,y],model_output_modice_tdab[1:nx2,x,y],model_output_modice_tdab[(nx2+1),x,y])


scenario_var_paleodata_big_anom = array(0,c((nx+1),ny,ne))
for (y in c(1:ne))
  for (x in c(1:ny))
    scenario_var_paleodata_big_anom[1:(nx+1),x,y] = c(scenario_var_paleodata_anom[(nx2+1):nx,x,y],scenario_var_paleodata_anom[1:nx2,x,y],scenario_var_paleodata_anom[(nx2+1),x,y])


# Calculate actual global annual temperature (add pre-industrial back on) ------------------------
# CW: Emulator works with anomalies, not absolutes i.e. GCM data going into emulator = anomalies, so direct output = anomalies.  If absolutes are wanted from emulator, need to take
# first value from GCM (of which there are 5 values, where 1st = PI) and add this back on.  As it says below, need to use "big" arrays for both emulatted anomalies and GCM

scenario_means_paleodata_orig = array(0,c(nx,ny,ne))
scenario_means_lgm_orig = array(0,c(nx,ny,1))

scenario_means_paleodata_orig_big = array(0,c((nx+1),ny,ne))
scenario_means_lgm_orig_big = array(0,c((nx+1),ny,1))

# CW: If wanting to extract proxy locations from emulated data, need to use "big" arrays

for (e in c(1:ne)){
  scenario_means_paleodata_orig[,,e] = scenario_means_paleodata_anom[,,e]+model_output_tdst[,,1]
}

for (e in c(1:ne)){
  scenario_means_paleodata_orig_big[,,e] = scenario_means_paleodata_big_anom[,,e]+model_output_tdst_big[,,1]
}

scenario_means_lgm_orig_big[,,1] = scenario_means_lgm_big_anom[,,1]+model_output_tdst_big[,,1]

# Select GCM reproduction of LGM temperature anomaly compared to tdstb ----------

tdab_exps_count = (seq(1, Exp_num_tdab))-1
tdab_LGM_time = tdab_exps_count[which(tdab_exps_count == 21)]
tdab_LGM_pos = tdab_exps[which(tdab_exps_count == 21)]

model_output_tdab_lgm_big_anom = model_output_modice_tdab_big[,,tdab_LGM_pos]



# Calculate individual timeslice maps for emulated experiments anomaly compared to tdstb ----------
# (800-0 kyr BP)

timeslices = array(c(800,782,700,300,1))

plot_timeslices_paleodata = array(0,c((nx+1),ny,length(timeslices)))


for (i in 1:length(timeslices)){
  plot_timeslices_paleodata[,,i] = scenario_means_paleodata_big_anom[,,timeslices[i]]
}



# Timeslices at 0, 18, 100, 500, 799kyr BP
# SLR (m) is 0, -130, 42, 44, 105 m

nam_timeslice = c(0, 18, 100, 500, 799)
nam_SLR = c(0, -130, 42, 44, 105)

plot_timeslices_paleodata_mod = plot_timeslices_paleodata

for (n in 1:length(plot_timeslices_paleodata)){
  if (plot_timeslices_paleodata[n] < -20){
    plot_timeslices_paleodata_mod[n] = -20
  } else {
  }
  if (plot_timeslices_paleodata[n] > 20){
    plot_timeslices_paleodata_mod[n] = 20
  } else {
  }
}


col_bwr = rgb(bwr(20))

i_order = c(1, 2, 3, 4)
labels_plot = c("<-20","-16","-12","-8","-4","0","4","8","12","16",">20")

timeslices_tdab = array(c(61, 65, 68, 115, 101, 91, 79, 117))

nam_timeslice_tdab = c(0, 4, 7, 92, 58, 38, 18, 100)
nam_SLR_tdab = c(0, -0.2, -9.5, -54.2, -88.9, -103.7, -123.8, -41)

model_output_modice_tdab_big_mod = model_output_modice_tdab_big

for (n in 1:length(model_output_modice_tdab_big)){
  if (model_output_modice_tdab_big[n] < -20){
    model_output_modice_tdab_big_mod[n] = -20
  } else {
  }
  if (model_output_modice_tdab_big[n] > 20){
    model_output_modice_tdab_big_mod[n] = 20
  } else {
  }
}


i_order_tdab = c(8, 9, 10, 12, 13, 14, 15, 11)


# Timeslices every kyr

col_bwr = rgb(bwr(20))

num = sprintf("%04d",1:length(times))
  

print("Finished running emulator")

check2 = Sys.time()
print(check2-check1)

################## CW: Write out netcdf of global absolutes, anomalies and variances #################################################

print("Writing out global absolutes")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "2m air temperature"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_temp_absolutes_Global_2.58_to_0MyrBP.nc", sep="")

# define dimensions
lon_dim <- ncdim_def("longitude","degrees_east",as.double(lons))
lat_dim <- ncdim_def("latitude","degrees_north",as.double(rev(lats))) # CW: Need to reverse the latitudes, as otherwise data are upside down
time_dim <- ncdim_def("time","",as.double(times))

# define variables
miss_value <- 2e20
nc_var_name <- "SAT"
nc_var_unit <- "degC"
nc_var_lname <- "2m air temperature"

# create variable and ncdf file with first variable so dimensions are set
var_def <- ncvar_def(nc_var_name,nc_var_unit,list(lon_dim,lat_dim,time_dim),miss_value,nc_var_lname,prec="double")
nc_file <- nc_create(nc_file_name,list(var_def))

# define data to save
data_to_save = scenario_means_paleodata_orig_big[,,]

# add data to variable in ncdf file
ncvar_put(nc_file,var_def,data_to_save)

# put additional attributes into dimension and data variables
ncatt_put(nc_file,"longitude","axis","X")
ncatt_put(nc_file,"latitude","axis","Y")
ncatt_put(nc_file,"time","axis","T")

# add global attributes
ncatt_put(nc_file,0,"title",nc_title)
ncatt_put(nc_file,0,"institution",nc_institution)
ncatt_put(nc_file,0,"source",nc_datasource)

# close the file, writing data to disk
nc_close(nc_file)

############

print("Writing out global anomalies")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "2m air temperature"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_temp_anomalies_Global_2.58_to_0MyrBP.nc", sep="")

# define dimensions
lon_dim <- ncdim_def("longitude","degrees_east",as.double(lons))
lat_dim <- ncdim_def("latitude","degrees_north",as.double(lats))
time_dim <- ncdim_def("time","",as.double(times))

# define variables
miss_value <- 2e20
nc_var_name <- "SAT"
nc_var_unit <- "degC"
nc_var_lname <- "2m air temperature"

# create variable and ncdf file with first variable so dimensions are set
var_def <- ncvar_def(nc_var_name,nc_var_unit,list(lon_dim,lat_dim,time_dim),miss_value,nc_var_lname,prec="double")
nc_file <- nc_create(nc_file_name,list(var_def))

# define data to save
data_to_save = scenario_means_paleodata_big_anom[,,]

# add data to variable in ncdf file
ncvar_put(nc_file,var_def,data_to_save)

# put additional attributes into dimension and data variables
ncatt_put(nc_file,"longitude","axis","X")
ncatt_put(nc_file,"latitude","axis","Y")
ncatt_put(nc_file,"time","axis","T")

# add global attributes
ncatt_put(nc_file,0,"title",nc_title)
ncatt_put(nc_file,0,"institution",nc_institution)
ncatt_put(nc_file,0,"source",nc_datasource)

# close the file, writing data to disk
nc_close(nc_file)

#############

print("Writing out global variances")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "2m air temperature"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_temp_variances_Global_2.58_to_0MyrBP.nc", sep="")

# define dimensions
lon_dim <- ncdim_def("longitude","degrees_east",as.double(lons))
lat_dim <- ncdim_def("latitude","degrees_north",as.double(rev(lats))) # CW: Need to reverse the latitudes, as otherwise data are upside down
time_dim <- ncdim_def("time","",as.double(times))

# define variables
miss_value <- 2e20
nc_var_name <- "SAT"
nc_var_unit <- "degC"
nc_var_lname <- "2m air temperature"

# create variable and ncdf file with first variable so dimensions are set
var_def <- ncvar_def(nc_var_name,nc_var_unit,list(lon_dim,lat_dim,time_dim),miss_value,nc_var_lname,prec="double")
nc_file <- nc_create(nc_file_name,list(var_def))

# define data to save
data_to_save = scenario_var_paleodata_big_anom[,,]

# add data to variable in ncdf file
ncvar_put(nc_file,var_def,data_to_save)

# put additional attributes into dimension and data variables
ncatt_put(nc_file,"longitude","axis","X")
ncatt_put(nc_file,"latitude","axis","Y")
ncatt_put(nc_file,"time","axis","T")

# add global attributes
ncatt_put(nc_file,0,"title",nc_title)
ncatt_put(nc_file,0,"institution",nc_institution)
ncatt_put(nc_file,0,"source",nc_datasource)

# close the file, writing data to disk
nc_close(nc_file)

############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################

print("Finished")
