###################################################################################################
## USE TWO EMULATORS TO PREDICT PAST TEMP FOR 800-0 kyr BP AT VARIOUS SITES #######################
## modice+tdab emulator (gl) AND modlowice emulator (ig), calibrated on temp ######################
## Using palaeo-proxy CO2 data [Luthi et al 2008] #################################################
## For Pleistocene paper 2020-02-03 ###############################################################
###################################################################################################

## CW's version 8 of original palaeo-emulator, Which is the same as version 7 but instead of looping over selected forcing combinations, does ALL possible combinations
## ie 32 simulations.  Also removes unnecessary stuff e.g. creating figures.  So only saves emulated anomalies/absolutes, either globally or at each proxy location.
## Working (as of /22)

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

# Setup 5 loops, which can either be 0 or 1.  Originally tried to use i,j,k,l,m, but these are already use lower down (where i becomes 800 and j 
# becomes 6) so to use something original

counter=1

for (forcing1 in 0:1){             # CO2 loop

  for (forcing2 in 0:1){           # Obliquity loop
    
    for (forcing3 in 0:1){         # Eccentricity (esinw) loop
      
      for (forcing4 in 0:1){       # Precession (ecosw) loop
        
        for (forcing5 in 0:1){     # Ice loop
          
          print(forcing1)
          print(forcing2)
          print(forcing3)
          print(forcing4)
          print(forcing5)
          
# Load data.  This needs to be done each and every time (rather than just once, before the loops, because otherwise constant values from previous)
# iteration will still exist.  Therefore need to start every time with fresh all forcings file
          
print("Loading all forcings for emulator")

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_forcing_paleodata_Focs_800kyr_BP_cjrw.R", sep=""))

# Setup a switch so that when any of the above loops = 0, the corresponding column is set to a constant value, as the last value (i.e. PI).  Otherwise, 
# i.e. when they = 1, don't do anything i.e. it will stay as varying
          
          if (forcing1 == 0){
            model_input_paleodata_800kyr_BP[,1] = model_input_paleodata_800kyr_BP[800,1]
          }
            
          if (forcing2 == 0){
            model_input_paleodata_800kyr_BP[,2] = model_input_paleodata_800kyr_BP[800,2]
          }
          
          if (forcing3 == 0){
            model_input_paleodata_800kyr_BP[,3] = model_input_paleodata_800kyr_BP[800,3]
          }
          
          if (forcing4 == 0){
            model_input_paleodata_800kyr_BP[,4] = model_input_paleodata_800kyr_BP[800,4]
          }
          
          if (forcing5 == 0){
            model_input_paleodata_800kyr_BP[,5] = model_input_paleodata_800kyr_BP[800,5]
          }
          
#stop()

print("Loading rest of input data")
          
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/Data - Climate NC/Tdst_temp_mm_1_5m_cjrw.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/Data - Climate NC/Tdst_temp_mm_uo_cjrw.R')

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/GSL_CGSLM_120kyr_BP_cjrw.R", sep=""))

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_insol_Laskar_jul_65N_800kyr_BP_cjrw.R", sep=""))


times = matrix(seq(-799, 0, 1))

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
#Rnum_rcp <- data.frame(data_input_Rnum_paleodata_800kyr_BP)

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

#rm(list=ls()[(ls() %in% c('scenario_means_gl_paleodata_anom', 'scenario_means_ig_paleodata_anom', 'scenario_var_gl_paleodata_anom', 'scenario_var_ig_paleodata_anom'))])

# CW: For unknown reasons, possibly because of the removal line (123) or possibly because of call to Emulator_remove_vars2.R, various colour bars need to be sourced again
# as otherwise it does not recognise bwr function

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_parula.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_bwr.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_wr.R')

col_bwr = rgb(bwr(20))

#image.plot(scenario_means_paleodata_anom[,,1]-scenario_var_gl_paleodata_anom[,,1],zlim=c(-0.1,0.1))
#image.plot(scenario_means_paleodata_anom[,,1]-scenario_means_ig_paleodata_anom[,,1],zlim=c(-0.1,0.1))

#image.plot(scenario_var_paleodata_anom[,,1]-scenario_var_gl_paleodata_anom[,,1],zlim=c(-0.1,0.1))
#image.plot(scenario_var_paleodata_anom[,,1]-scenario_means_ig_paleodata_anom[,,1],zlim=c(-0.1,0.1))


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
  
for (i in 1:length(times)){
  #png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/Other\\SAT\\gif\\",num[i],".png", sep=""),width=500,height=300)
  #image.plot(lons, lats, scenario_means_rcp45_big_anom[,c(73:1),i], zlim=c(-20,20), col=col_bwr, xlab="", ylab="", main=paste("Upd RCP4.5 SAT anomaly at ",(i-1)," kyr AP (SLR ",model_input_rcp45_800kyr_BP[i,5]," m)", sep=""),
   #          axis.args=list(at=seq(-20, 20, 4),labels=seq(-20, 20, 4),cex.axis=1), 
    #         legend.args=list(text='', side=4, font=4, line=2.5, cex=0.8))
  #map(add=T, interior=F, col="black")
  #dev.off()
}

# CW: Write out netcdf of global anomalies and absolutes #################################################################################################

print("Writing out global absolutes")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "2m air temperature"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Global/Em_output_data_temp_absolutes_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_Global_0.8_to_0MyrBP.nc", sep="")

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

print("Writing out global anomalies")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "2m air temperature"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Global/Em_output_data_temp_anomalies_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_Global_0.8_to_0MyrBP.nc", sep="")

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

##########################################################################################################################################################

# Extract emulated climate for each site ----------------------------------------------------------------------------

# ODP-982, North Atlantic (lat=57.8, lon=-15.9) -----------------------------------------------------------
# Equivalent lat GB=57.5, #14, lon GB=-15, #45 (#37 = 0 lat, #49 = 0 lon)

lats_ODP982 = 14
lons_ODP982 = 45

# Extract site data

my_tim_anom_ODP982_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_ODP982,lats_ODP982,])  # CW: These are the anomalies coming out of the emulator, using "big" array

my_var_anom_ODP982_paleodata_e = scenario_var_paleodata_big_anom[lons_ODP982,lats_ODP982,]            # CW: These are the variability coming out of the emulator, using "big" array

my_tim_means_ODP982_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_ODP982,lats_ODP982,]) # CW: These are the absolutes, calculated above as anomalies + PI, using "big" array

my_tim_pi_ODP982_paleodata_e = matrix(model_output_tdst[lons_ODP982,lats_ODP982,])                    # CW: These are the GCM values (where ,1 = PI), using original array

my_tim_pi_ODP982_paleodata_e_big = matrix(model_output_tdst_big[lons_ODP982,lats_ODP982,])            # CW: These are the GCM values (where ,1 = PI), using "big" array

# Using "big" array, my_tim_means* = my_time_anom* + my_tim_pi_*_big i.e. absolutes = anomalies from big array + PI from big array

# ODP-722, Arabian Sea (lat=16.6, lon=59.8) -----------------------------------------------------------
# Equivalent lat GB=17.5, #30, lon GB=60, #65 (#37 = 0 lat, #49 = 0 lon)

lats_ODP722 = 30
lons_ODP722 = 65

# Extract site data

my_tim_anom_ODP722_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_ODP722,lats_ODP722,])

my_var_anom_ODP722_paleodata_e = scenario_var_paleodata_big_anom[lons_ODP722,lats_ODP722,]

my_tim_means_ODP722_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_ODP722,lats_ODP722,]) # CW: These are the absolutes, calculated above as anomalies + PI, using "big" array

my_tim_pi_ODP722_paleodata_e_big = matrix(model_output_tdst_big[lons_ODP722,lats_ODP722,])            # CW: These are the GCM values (where ,1 = PI), using "big" array

# ODP-1143, South China Sea (lat=9.4, lon=113.3) -----------------------------------------------------------
# Equivalent lat GB=10, #33, lon GB=112.5, #79 (#37 = 0 lat, #49 = 0 lon)

lats_ODP1143 = 33
lons_ODP1143 = 79

# Extract site data

my_tim_anom_ODP1143_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_ODP1143,lats_ODP1143,])

my_var_anom_ODP1143_paleodata_e = scenario_var_paleodata_big_anom[lons_ODP1143,lats_ODP1143,]

my_tim_means_ODP1143_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_ODP1143,lats_ODP1143,]) # CW: These are the absolutes, calculated above as anomalies + PI, using "big" array

my_tim_pi_ODP1143_paleodata_e_big = matrix(model_output_tdst_big[lons_ODP1143,lats_ODP1143,])            # CW: These are the GCM values (where ,1 = PI), using "big" array

# ODP-846, Equatorial East Pacific (lat=-3.1, lon=-90.07) -----------------------------------------------------------
# Equivalent lat GB=-2.5, #38, lon GB=-90, #25 (#37 = 0 lat, #49 = 0 lon)

lats_ODP846 = 38
lons_ODP846 = 25
  
# Extract site data

my_tim_anom_ODP846_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_ODP846,lats_ODP846,])

my_var_anom_ODP846_paleodata_e = scenario_var_paleodata_big_anom[lons_ODP846,lats_ODP846,]

my_tim_means_ODP846_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_ODP846,lats_ODP846,]) # CW: These are the absolutes, calculated above as anomalies + PI, using "big" array

my_tim_pi_ODP846_paleodata_e_big = matrix(model_output_tdst_big[lons_ODP846,lats_ODP846,])            # CW: These are the GCM values (where ,1 = PI), using "big" array

# ODP-1090, Central Subantarctic (lat=-42.9, lon=8.9) -----------------------------------------------------------
# Equivalent lat GB=-42.5, #54, lon GB=7.5, #51 (#37 = 0 lat, #49 = 0 lon)

lats_ODP1090 = 54
lons_ODP1090 = 51

# Extract site data

my_tim_anom_ODP1090_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_ODP1090,lats_ODP1090,])

my_var_anom_ODP1090_paleodata_e = scenario_var_paleodata_big_anom[lons_ODP1090,lats_ODP1090,]

my_tim_means_ODP1090_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_ODP1090,lats_ODP1090,]) # CW: These are the absolutes, calculated above as anomalies + PI, using "big" array

my_tim_pi_ODP1090_paleodata_e_big = matrix(model_output_tdst_big[lons_ODP1090,lats_ODP1090,])            # CW: These are the GCM values (where ,1 = PI), using "big" array


# Dome C, Antarctica (lat=-75.1, lon=123.3) -----------------------------------------------------------
# Equivalent lat GB=-75, #67, lon GB=123.75, #82 (#37 = 0 lat, #49 = 0 lon)

lats_DomeC = 67
lons_DomeC = 82

# Extract site data

my_tim_anom_DomeC_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_DomeC,lats_DomeC,])

my_var_anom_DomeC_paleodata_e = scenario_var_paleodata_big_anom[lons_DomeC,lats_DomeC,]

my_tim_means_DomeC_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_DomeC,lats_DomeC,]) # CW: These are the absolutes, calculated above as anomalies + PI, using "big" array

my_tim_pi_DomeC_paleodata_e_big = matrix(model_output_tdst_big[lons_DomeC,lats_DomeC,])            # CW: These are the GCM values (where ,1 = PI), using "big" array


# Hulu Cave, China (lat=32.5, lon=119.2) -----------------------------------------------------------
# Equivalent lat GB=32.5, #24, lon GB=120, #81 (#37 = 0 lat, #49 = 0 lon)

lats_Hulu = 24
lons_Hulu = 81

# Extract site data

my_tim_anom_Sanbao_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_Hulu,lats_Hulu,])

my_var_anom_Sanbao_paleodata_e = scenario_var_paleodata_big_anom[lons_Hulu,lats_Hulu,]



# Sanbao Cave, China (lat=31.7, lon=110.4) -----------------------------------------------------------
# Equivalent lat GB=32.5, #24, lon GB=108.75, #78 (#37 = 0 lat, #49 = 0 lon)

lats_Sanbao = 24
lons_Sanbao = 78

# Extract site data

my_tim_anom_Sanbao_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_Sanbao,lats_Sanbao,])

my_var_anom_Sanbao_paleodata_e = scenario_var_paleodata_big_anom[lons_Sanbao,lats_Sanbao,]



# Dongge Cave, China (lat=25.3, lon=108.1) -----------------------------------------------------------
# Equivalent lat GB=25, #27, lon GB=108.75, #78 (#37 = 0 lat, #49 = 0 lon)

lats_Dongge = 25
lons_Dongge = 78

# Extract site data

my_tim_anom_Dongge_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_Dongge,lats_Dongge,])

my_var_anom_Dongge_paleodata_e = scenario_var_paleodata_big_anom[lons_Dongge,lats_Dongge,]


# Plot map showing all sites (emulated experiments anomaly compared to tdstb)

scenario_means_sites_paleodata_big = array(0,c((nx+1),ny,1))

scenario_means_sites_paleodata_big[lons_ODP982,lats_ODP982,1] = 20 # ODP-982, North Atlantic
scenario_means_sites_paleodata_big[lons_ODP722,lats_ODP722,1] = 20 # ODP-722, Arabian Sea
scenario_means_sites_paleodata_big[lons_ODP1143,lats_ODP1143,1] = 20 # ODP-1143, South China Sea
scenario_means_sites_paleodata_big[lons_ODP846,lats_ODP846,1] = 20 # ODP-846, Equatorial East Pacific
scenario_means_sites_paleodata_big[lons_ODP1090,lats_ODP1090,1] = 20 # ODP-1090, Central Subarctic
scenario_means_sites_paleodata_big[lons_DomeC,lats_DomeC,1] = 20 # Dome C
scenario_means_sites_paleodata_big[lons_Hulu,lats_Hulu,1] = 20 # Hulu, China
scenario_means_sites_paleodata_big[lons_Sanbao,lats_Sanbao,1] = 20 # Sanbao, China
scenario_means_sites_paleodata_big[lons_Dongge,lats_Dongge,1] = 20 # Dongge, China


image.plot(lons, lats, scenario_means_sites_paleodata_big[,c(73:1),1], zlim=c(-25,25), xlab="Longitude", ylab="Latitude", main="Sites")
map(add=T, interior=F, col="black")

print("Reading SST data")

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
# This is already in anomaly space, so if absolute values are wanted they need to be calculated by adding first valid year (i.e. year 1) back onto anomalies

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_D_Jouzeletal07_800kyr_BP_cjrw.R") # Dome C + EPICA

data_interp = approx(data_D_Jouzeletal07_800kyr_BP[,2],data_D_Jouzeletal07_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_D_Jouzeletal07_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

data_D_Jouzeletal07_800kyr_BP_interp_ab = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2)) # Create another array identical to above
data_D_Jouzeletal07_800kyr_BP_interp_ab[,2] = data_D_Jouzeletal07_800kyr_BP_interp_ab[,2] + data_D_Jouzeletal07_800kyr_BP_interp_ab[2,2] # Add first valid year (2nd row) onto all anomalies

# plot(data_D_Jouzeletal07_800kyr_BP[,2],data_D_Jouzeletal07_800kyr_BP[,3], type="l")
# plot(data_D_Jouzeletal07_800kyr_BP_interp[,1],data_D_Jouzeletal07_800kyr_BP_interp[,2], type="l")


# Read in SST stack data for 800-0 kyr BP [Shakun et al 2015] #  CW: This doesn't appear to be used anywhere else?

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_SST_stack_Shakunetal15_800kyr_BP_cjrw.R") # Global stack

data_interp = approx(data_SST_stack_Shakunetal15_800kyr_BP[,2],data_SST_stack_Shakunetal15_800kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_SST_stack_Shakunetal15_800kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_SST_stack_Shakunetal15_800kyr_BP[,2],data_SST_stack_Shakunetal15_800kyr_BP[,3], type="l")
# plot(data_SST_stack_Shakunetal15_800kyr_BP_interp[,1],data_SST_stack_Shakunetal15_800kyr_BP_interp[,2], type="l")

###############################################################################################################################################################

# Pre-industrial SST observations (HadISST) -----------------------------------------------------

# Load data

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_HadISST_SST_mon_1870_1900_BP_cjrw.R')

data_HadISST_SST_mon_1870_1900_BP_all = data_HadISST_SST_mon_1870_1900_BP_all[,2:361]

lats_HadISST = seq(-89.5,89.5,1)
lons_HadISST = seq(-179.5,179.5,1)

print("Finished reading SST data")

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

print("Extracting proxy data")

# Select site data

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

check2 = Sys.time()
print(check2-check1)

# Save data ---------------------------------------------------------

# SAT ANOMALY AT PROXY LOCATIONS (TS) ---------------------------------------------------------------

print("Saving emulated anomalies")

ncol = 3

data_to_save = array(c(times, round(my_tim_anom_ODP982_paleodata_e, digits = 5), round(my_var_anom_ODP982_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC", "Uncertainty_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP982_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_anom_ODP722_paleodata_e, digits = 5), round(my_var_anom_ODP722_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC", "Uncertainty_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP722_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_anom_ODP1143_paleodata_e, digits = 5), round(my_var_anom_ODP1143_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC", "Uncertainty_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP1143_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_anom_ODP846_paleodata_e, digits = 5), round(my_var_anom_ODP846_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC", "Uncertainty_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP846_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_anom_ODP1090_paleodata_e, digits = 5), round(my_var_anom_ODP1090_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC", "Uncertainty_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP1090_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_anom_DomeC_paleodata_e, digits = 5), round(my_var_anom_DomeC_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC", "Uncertainty_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_DomeC_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

# SAT ABSOLUTE VALUES AT PROXY LOCATIONS (TS) ---------------------------------------------------------------

print("Saving emulated absolutes")

ncol = 2

data_to_save = array(c(times, round(my_tim_means_ODP982_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_absolutes_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP982_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_means_ODP722_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_absolutes_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP722_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_means_ODP1143_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_absolutes_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP1143_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_means_ODP846_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_absolutes_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP846_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_means_ODP1090_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_absolutes_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_ODP1090_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

data_to_save = array(c(times, round(my_tim_means_DomeC_paleodata_e, digits = 5)),c(length(times),ncol))
colnames(data_to_save) = c("Time_kyrAP", "Temp_degC")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/All 32 simulations/Data/Em_output_data_temp_absolutes_E",forcing1,forcing2,forcing3,forcing4,forcing5,"_DomeC_natural_0.8_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

print(counter)
counter=counter+1

print("Onto next simulation")

        } # Close forcing5 loop
        
      } # Close forcing4 loop
      
    } # Close forcing3 loop
    
  } # Close forcing2 loop
  
} # Close forcing1 loop
        
print("Finished everything")

############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################

