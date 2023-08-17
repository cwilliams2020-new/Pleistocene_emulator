###################################################################################################
## USE TWO EMULATORS TO PREDICT PAST TEMP FOR 800-0 kyr BP AT VARIOUS SITES #######################
## modice+tdab emulator (gl) AND modlowice emulator (ig), calibrated on temp ######################
## Using palaeo-proxy CO2 data [Luthi et al 2008] #################################################
## For Pleistocene paper 2020-02-03 ###############################################################
###################################################################################################

## CW's version 2 of original palaeo-emulator, Which is the same as version 1 but now includes writing out of global variances
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

# CW: Just to speed up plotting, create a switch here so don't need to run emulator each time.  So only runs emulator if runemul = 1, and only saves figures if savefigs = 1

runemul = 1
savefigs = 0

par_orig = par()

if (runemul == 1){

## Load and calibrate emulator #########################################################################################

print("Loading emulator")

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/Loading/Emul_in_modice_tdab_+_modlowice_LT2000ppm_optd_on_Temp_sim_Precip_cjrw.R')
  
## Set up #########################################################################################

# CW: Because there is only one forcing experiment for this, no need to have forcings_included in either input or output files

# Load data

print("Loading input data for emulator")

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_forcing_paleodata_2Ma_BP_cjrw.R", sep=""))
source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_Rnum_paleodata_Focs_800kyr_BP_cjrw.R", sep="")) # Leave here, but not actually used apart from plotting

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/Data - Climate NC/Tdst_precip_mm_srf_cjrw.R')
  
source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/GSL_CGSLM_120kyr_BP_cjrw.R", sep="")) # Leave here, but not actually used apart from plotting

source(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_insol_Laskar_jul_65N_800kyr_BP_cjrw.R", sep="")) # Leave here, but not actually used apart from plotting

# Create 2 lots of time arrays, one for 2.58 Ma and one for 800 kya
  
times = matrix(seq(-2579, 0, 1))
times2 = matrix(seq(-799, 0, 1))

Scenario_paleodata <- data.frame(model_input_paleodata_800kyr_BP)
Rnum_rcp <- data.frame(data_input_Rnum_paleodata_800kyr_BP)

# set CO2 to be log (natural) of CO2

Scenario_paleodata[,1] = log(Scenario_paleodata[,1])

# rescaled to make it compatible with the 'X' used to calibrate the emulator

XScenario_paleodata <- sweep(Scenario_paleodata,  2, attr(X_all, "scaled:center"), '-')

XScenario_paleodata <- sweep(XScenario_paleodata, 2, attr(X_all, "scaled:scale"), '/')

## Run emulator using input data for last 800 kyr #################################################

# run emulator (feed every row of XScenario into emulator, put this in a list, 
# and repack this into an array.

print("Running emulator")
Sys.time()
check1 = Sys.time()

OUT_gl_paleodata = lapply(seq(nrow(XScenario_paleodata)), function(i) pe_p((XScenario_paleodata[i,, drop=FALSE]) , E_gl ) )

OUT_ig_paleodata = lapply(seq(nrow(XScenario_paleodata)), function(i) pe_p((XScenario_paleodata[i,, drop=FALSE]) , E_ig ) )

#rm(list=ls()[(ls() %in% c('E_gl', 'E_ig'))]) # CW: Not entirely sure why these are being removed (other than to save space), but this needs to be commented out when running all sensitivity experiments together


# extract mean and variance

scenario_means_gl_paleodata_anom <- simplify2array( lapply(OUT_gl_paleodata, function(i) i$mean))

scenario_means_ig_paleodata_anom <- simplify2array( lapply(OUT_ig_paleodata, function(i) i$mean))

scenario_var_gl_paleodata_anom <- simplify2array( lapply(OUT_gl_paleodata, function(i) i$var))

scenario_var_ig_paleodata_anom <- simplify2array( lapply(OUT_ig_paleodata, function(i) i$var))

rm(list=ls()[(ls() %in% c("OUT_gl_paleodata"))])
rm(list=ls()[(ls() %in% c("OUT_ig_paleodata"))])


# attach times as an attribute

attr(scenario_means_gl_paleodata_anom, "times") <- times

attr(scenario_means_ig_paleodata_anom, "times") <- times

attr(scenario_var_gl_paleodata_anom, "times") <- times

attr(scenario_var_ig_paleodata_anom, "times") <- times


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

#rm(list=ls()[(ls() %in% c("scenario_means_gl_paleodata_anom", "scenario_means_ig_paleodata_anom", "scenario_var_gl_paleodata_anom", "scenario_var_ig_paleodata_anom"))])

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

# Convert from mm per sec to mm per month (mm same as kg m-2)

model_output_tdst[,,] = model_output_tdst[,,]*60*60*24*30

scenario_means_paleodata_anom[,,] = scenario_means_paleodata_anom[,,]*30

scenario_var_paleodata_anom[,,] = scenario_var_paleodata_anom[,,]*30

# Calculate true global mean annual temperature evolution for next 200 kyr ------------------------
# i.e grid box area weighted mean

nx = 96 # lon
ny = 73 # lat
nx2 = nx/2
ne = length(times) # timeslice exps
nc = 5 # control exps

lons <- array(0:(nx-1))*3.75
lons[1:(nx+1)] <- c(lons[(nx2+1):nx]-360,lons[1:nx2],lons[nx2+1])

lats <- array(0:(ny-1))*2.5-90

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

for (y in 1:ne){
  my_tim_paleodata_e_truemean[y,] = sum(sum((scenario_means_paleodata_anom[,,y]*gb_area_manual[,])))
}

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


model_output_tdst_big = array(0,c((nx+1),ny,nc))
for (y in c(1:nc))
  for (x in c(1:ny))
    model_output_tdst_big[1:(nx+1),x,y] = c(model_output_tdst[(nx2+1):nx,x,y],model_output_tdst[1:nx2,x,y],model_output_tdst[(nx2+1),x,y])

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

scenario_means_paleodata_orig_big = array(0,c((nx+1),ny,ne))

# CW: If wanting to extract proxy locations from emulated data, need to use "big" arrays

for (e in c(1:ne)){
  scenario_means_paleodata_orig[,,e] = scenario_means_paleodata_anom[,,e]+model_output_tdst[,,1]
}

for (e in c(1:ne)){
  scenario_means_paleodata_orig_big[,,e] = scenario_means_paleodata_big_anom[,,e]+model_output_tdst_big[,,1]
}


scenario_means_paleodata_orig[which(scenario_means_paleodata_orig<0, arr.ind = T)] = 0

scenario_means_paleodata_orig_big[which(scenario_means_paleodata_orig_big<0, arr.ind = T)] = 0



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

# for (i in 1:length(timeslices)){
#   #png(file=paste("C:\\Users\\nl6806\\OneDrive - University of Bristol\\PostDoc\\2017-02-15 Posiva + SKB\\5. Output\\Plots\\2020-02-03 Pleistocene paper\\Other\\SAT\\Mp_em_SAT_rcp8.5_upd_",i_order[i],"_SLR",nam_SLR[i],"m_",nam_timeslice[i],"kyrAP.png", sep=""),width=1000,height=600)
#   image.plot(lons, lats, plot_timeslices_paleodata_mod[,c(73:1),i], zlim=c(-20,20), col=col_bwr, xlab="", ylab="", main=paste("Upd RCP8.5 SAT anomaly at ",nam_timeslice[i]," kyr AP (SLR ",nam_SLR[i]," m)", sep=""),
#              axis.args=list(at=seq(-20, 20, 4),labels=labels_plot,cex.axis=1), 
#              legend.args=list(text="", side=4, font=4, line=2.5, cex=0.8))
#   map(add=T, interior=F, col="black")
#   #dev.off()
# }


timeslices_tdab = array(c(61, 65, 68, 115, 101, 91, 79, 117))

nam_timeslice_tdab = c(0, 4, 7, 92, 58, 38, 18, 100)
nam_SLR_tdab = c(0, -0.2, -9.5, -54.2, -88.9, -103.7, -123.8, -41)

model_output_modice_tdab_big_mod = model_output_modice_tdab_big

for (n in 1:length(model_output_modice_tdab_big)){
  if (model_output_modice_tdab_big[n] < -3){
    model_output_modice_tdab_big_mod[n] = -3
  } else {
  }
  if (model_output_modice_tdab_big[n] > 3){
    model_output_modice_tdab_big_mod[n] = 3
  } else {
  }
}


i_order_tdab = c(8, 9, 10, 12, 13, 14, 15, 11)

#for (i in 1:length(timeslices_tdab)){
#  png(file=paste("C:\\Users\\nl6806\\OneDrive - University of Bristol\\PostDoc\\2017-02-15 Posiva + SKB\\5. Output\\Plots\\2020-02-03 Pleistocene paper\\Other\\Ppt\\Mp_tdab_Ppt_",i_order_tdab[i],"_SLR",nam_SLR_tdab[i],"m_",nam_timeslice_tdab[i],"kyrBP.png", sep=""),width=1000,height=600)
#  image.plot(lons, lats, model_output_modice_tdab_big_mod[,c(73:1),timeslices_tdab[i]], zlim=c(-3,3), col=col_bwr, xlab="", ylab="", main=paste("tdab SAT anomaly at ",nam_timeslice_tdab[i]," kyr BP (SLR ",nam_SLR_tdab[i]," m)", sep=""),
#             axis.args=list(at=seq(-3, 3, 1),labels=labels_plot,cex.axis=1), 
#             legend.args=list(text="", side=4, font=4, line=2.5, cex=0.8))
#  map(add=T, interior=F, col="black")
#  dev.off()
#}


# Timeslices every kyr

col_bwr = rgb(bwr(20))

num = sprintf("%04d",1:length(times))

for (i in 1:length(times)){
  #png(file=paste("C:\\Users\\nl6806\\OneDrive - University of Bristol\\PostDoc\\2017-02-15 Posiva + SKB\\5. Output\\Plots\\2020-02-03 Pleistocene paper\\Other\\SAT\\gif\\",num[i],".png", sep=""),width=500,height=300)
  #image.plot(lons, lats, scenario_means_rcp45_big_anom[,c(73:1),i], zlim=c(-20,20), col=col_bwr, xlab="", ylab="", main=paste("Upd RCP4.5 SAT anomaly at ",(i-1)," kyr AP (SLR ",model_input_rcp45_800kyr_BP[i,5]," m)", sep=""),
  #          axis.args=list(at=seq(-20, 20, 4),labels=seq(-20, 20, 4),cex.axis=1), 
  #         legend.args=list(text="", side=4, font=4, line=2.5, cex=0.8))
  #map(add=T, interior=F, col="black")
  #dev.off()
}

print("Finished running emulator")

check2 = Sys.time()
print(check2-check1)

################## CW: Write out netcdf of global absolutes, anomalies and variances #################################################

print("Writing out global absolutes")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "Surface precipitation"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_precip_absolutes_Global_2.58_to_0MyrBP.nc", sep="")

# define dimensions
lon_dim <- ncdim_def("longitude","degrees_east",as.double(lons))
lat_dim <- ncdim_def("latitude","degrees_north",as.double(rev(lats))) # CW: Need to reverse the latitudes, as otherwise data are upside down
time_dim <- ncdim_def("time","",as.double(times))

# define variables
miss_value <- 2e20
nc_var_name <- "Precip"
nc_var_unit <- "mm/mo"
nc_var_lname <- "Surface precipitation"

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

#############

print("Writing out global anomalies")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "Surface precipitation"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_precip_anomalies_Global_2.58_to_0MyrBP.nc", sep="")

# define dimensions
lon_dim <- ncdim_def("longitude","degrees_east",as.double(lons))
lat_dim <- ncdim_def("latitude","degrees_north",as.double(lats))
time_dim <- ncdim_def("time","",as.double(times))

# define variables
miss_value <- 2e20
nc_var_name <- "Precip"
nc_var_unit <- "mm/mo"
nc_var_lname <- "Surface precipitation"

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

print("Writing out global variances")

# define global attributes
nc_institution = "University of Bristol, 2022"
nc_datasource = "Climate emulator"

nc_title = "Surface precipitation"
nc_file_name <- paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_precip_variances_Global_2.58_to_0MyrBP.nc", sep="")

# define dimensions
lon_dim <- ncdim_def("longitude","degrees_east",as.double(lons))
lat_dim <- ncdim_def("latitude","degrees_north",as.double(lats))
time_dim <- ncdim_def("time","",as.double(times))

# define variables
miss_value <- 2e20
nc_var_name <- "Precip"
nc_var_unit <- "mm/mo"
nc_var_lname <- "Surface precipitation"

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

stop()

####################

##########################################################################################################################################################

# Extract emulated climate for each site ----------------------------------------------------------------------------

# Hulu Cave, China (lat=32.5, lon=119.2) -----------------------------------------------------------
# Equivalent lat GB=32.5, #24, lon GB=120, #81 (#37 = 0 lat, #49 = 0 lon)

lats_Hulu = 24
lons_Hulu = 81

# Extract site data

my_tim_anom_Hulu_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_Hulu,lats_Hulu,])  # CW: These are the anomalies coming out of the emulator, using "big" array

my_var_anom_Hulu_paleodata_e = scenario_var_paleodata_big_anom[lons_Hulu,lats_Hulu,]            # CW: These are the variability coming out of the emulator, using "big" array

my_tim_means_Hulu_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_Hulu,lats_Hulu,]) # CW: These are the absolutes, calculated above as anomalies + PI, using "big" array

my_tim_pi_Hulu_paleodata_e = matrix(model_output_tdst[lons_Hulu,lats_Hulu,])                    # CW: These are the GCM values (where ,1 = PI), using original array

my_tim_pi_Hulu_paleodata_e_big = matrix(model_output_tdst_big[lons_Hulu,lats_Hulu,])            # CW: These are the GCM values (where ,1 = PI), using "big" array

# Using "big" array, my_tim_means* = my_time_anom* + my_tim_pi_*_big i.e. absolutes = anomalies from big array + PI from big array

# Sanbao Cave, China (lat=31.7, lon=110.4) -----------------------------------------------------------
# Equivalent lat GB=32.5, #24, lon GB=108.75, #78 (#37 = 0 lat, #49 = 0 lon)

lats_Sanbao = 24
lons_Sanbao = 78

# Extract site data

my_tim_anom_Sanbao_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_Sanbao,lats_Sanbao,])

my_var_anom_Sanbao_paleodata_e = scenario_var_paleodata_big_anom[lons_Sanbao,lats_Sanbao,]

my_tim_means_Sanbao_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_Sanbao,lats_Sanbao,]) 

my_tim_pi_Sanbao_paleodata_e = matrix(model_output_tdst[lons_Sanbao,lats_Sanbao,])

my_tim_pi_Sanbao_paleodata_e_big = matrix(model_output_tdst_big[lons_Sanbao,lats_Sanbao,])

# Dongge Cave, China (lat=25.3, lon=108.1) -----------------------------------------------------------
# Equivalent lat GB=25, #27, lon GB=108.75, #78 (#37 = 0 lat, #49 = 0 lon)

lats_Dongge = 25
lons_Dongge = 78

# Extract site data

my_tim_anom_Dongge_paleodata_e = matrix(scenario_means_paleodata_big_anom[lons_Dongge,lats_Dongge,])

my_var_anom_Dongge_paleodata_e = scenario_var_paleodata_big_anom[lons_Dongge,lats_Dongge,]

my_tim_means_Dongge_paleodata_e = matrix(scenario_means_paleodata_orig_big[lons_Dongge,lats_Dongge,]) 

my_tim_pi_Dongge_paleodata_e = matrix(model_output_tdst[lons_Dongge,lats_Dongge,])

my_tim_pi_Dongge_paleodata_e_big = matrix(model_output_tdst_big[lons_Dongge,lats_Dongge,])

# Read in proxy data and interpolate to 1 kyr intervals ----------------------------------------------------------------------------

# Read in d18O data for 641-0 kyr BP [Wang et al 2001, 2008]
# Original D + H data has been made 1.6 %o more negative to account for higher values at D + H than SB [Wang et al 2008, Cheng et al. 2016]

print("Reading proxy data")

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_d18O_Chengetal16_641kyr_BP_cjrw.R") # Composite record, China (Sanbao, Hulu, Dongge)

data_interp = approx(data_d18O_Chengetal16_641kyr_BP[,2],data_d18O_Chengetal16_641kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_d18O_Chengetal16_641kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_d18O_Chengetal16_641kyr_BP[,2],data_d18O_Chengetal16_641kyr_BP[,3], type="l", xlim=c(min(times), max(times)), ylim=c(-12,-4))
# plot(data_d18O_Chengetal16_641kyr_BP_interp[,1],data_d18O_Chengetal16_641kyr_BP_interp[,2], type="l", ylim=c(-12,-4))

###############################################################################################################################################################

## Print out 2.58 Ma emulated data, absolutes and anomalies

print("Saving emulated precipitation anomalies")

ncol = 3

data_to_save = array(c(times, round(my_tim_anom_Hulu_paleodata_e, digits = 5), round(my_var_anom_Hulu_paleodata_e, digits = 5)),c(length(rev(times)),ncol))
colnames(data_to_save) = c("Time_kyrAP", "dPrecip_mm/day", "Uncertainty_mm/day")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Em_output_data_precip_Hul_natural_2.58_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)


data_to_save = array(c(times, round(my_tim_anom_Sanbao_paleodata_e, digits = 5), round(my_var_anom_Sanbao_paleodata_e, digits = 5)),c(length(rev(times)),ncol))
colnames(data_to_save) = c("Time_kyrAP", "dPrecip_mm/day", "Uncertainty_mm/day")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Em_output_data_precip_San_natural_2.58_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)


data_to_save = array(c(times, round(my_tim_anom_Dongge_paleodata_e, digits = 5), round(my_var_anom_Dongge_paleodata_e, digits = 5)),c(length(rev(times)),ncol))
colnames(data_to_save) = c("Time_kyrAP", "dPrecip_mm/day", "Uncertainty_mm/day")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Em_output_data_precip_Don_natural_2.58_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

print("Saving emulated precipitation absolutes")

ncol = 3

data_to_save = array(c(times, round(my_tim_means_Hulu_paleodata_e, digits = 5), round(my_var_anom_Hulu_paleodata_e, digits = 5)),c(length(rev(times)),ncol))
colnames(data_to_save) = c("Time_kyrAP", "dPrecip_mm/day", "Uncertainty_mm/day")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Em_output_data_precip_absolutes_Hul_natural_2.58_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)


data_to_save = array(c(times, round(my_tim_means_Sanbao_paleodata_e, digits = 5), round(my_var_anom_Sanbao_paleodata_e, digits = 5)),c(length(rev(times)),ncol))
colnames(data_to_save) = c("Time_kyrAP", "dPrecip_mm/day", "Uncertainty_mm/day")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Em_output_data_precip_absolutes_San_natural_2.58_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)


data_to_save = array(c(times, round(my_tim_means_Dongge_paleodata_e, digits = 5), round(my_var_anom_Dongge_paleodata_e, digits = 5)),c(length(rev(times)),ncol))
colnames(data_to_save) = c("Time_kyrAP", "dPrecip_mm/day", "Uncertainty_mm/day")
file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Em_output_data_precip_absolutes_Don_natural_2.58_to_0MyrBP.txt", sep="")
write.table(data_to_save, file = file_nam, append = FALSE, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)

} # Close emulator and data saving

############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################

# Save figures ------------------------------------------------------------

if (savefigs == 1){
  
print("Saving figures")
 
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_parula.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_bwr.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_wr.R')

# Calculate error bars for mean annual precipitation anomaly at sites for -2.58-0 Ma AP

my_tim_anom_Hulu_paleodata_e_lower_1sd = my_tim_anom_Hulu_paleodata_e-sqrt(my_var_anom_Hulu_paleodata_e)
my_tim_anom_Hulu_paleodata_e_upper_1sd = my_tim_anom_Hulu_paleodata_e+sqrt(my_var_anom_Hulu_paleodata_e)
  
my_tim_anom_Sanbao_paleodata_e_lower_1sd = my_tim_anom_Sanbao_paleodata_e-sqrt(my_var_anom_Sanbao_paleodata_e)
my_tim_anom_Sanbao_paleodata_e_upper_1sd = my_tim_anom_Sanbao_paleodata_e+sqrt(my_var_anom_Sanbao_paleodata_e)
  
my_tim_anom_Dongge_paleodata_e_lower_1sd = my_tim_anom_Dongge_paleodata_e-sqrt(my_var_anom_Dongge_paleodata_e)
my_tim_anom_Dongge_paleodata_e_upper_1sd = my_tim_anom_Dongge_paleodata_e+sqrt(my_var_anom_Dongge_paleodata_e)

# Plot evolution of mean annual precipitation anomaly and d18O (precipitation) from speloethem data [Cheng et al. 2016] at China sites for -641-0 kyr AP

xlimits = c(0,2580)
ylimits_a = c(-30,30)
ylimits_b = c(-4.5,-14.8)

xticks = seq(0, 2580, by = 250)
yticks_a = seq(-30, 30, by = 5)
yticks_b = seq(-4.5, -14.8, by = -0.5)

xlabels = c("0", "", "500", "", "1000", "", "1500", "", "2000", "", "2500")
ylabels_a = c("-30", "", "-20", "", "-10", "", "0", "", "10", "", "20", "", "30")
ylabels_b = c("", "", "", "-6", "", "", "", "-8", "", "", "", "-10", "", "", "", "-12", "", "", "", "-14", "")

# Setup constants
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

opacity_shading = 0.5

col_proxy = "black"
col_Hulu = "lightskyblue4"
col_Hulu_800 = "orange"
col_Sanbao = "blue"
col_Sanbao_800 = "red"
col_Dongge = "darkcyan"
col_Dongge_800 = "pink"
col_shade_lightteal = c(col2rgb("lightskyblue")/255, opacity_shading)
col_shade_lightblue = c(col2rgb("royalblue1")/255, opacity_shading)
col_shade_lightcyan = c(col2rgb("cyan")/255, opacity_shading)

# Extract first 800 kya (i.e. last 800 values from timeseries)

my_tim_anom_Hulu_paleodata_e_800=my_tim_anom_Hulu_paleodata_e[1781:2580,]
my_tim_anom_Sanbao_paleodata_e_800=my_tim_anom_Sanbao_paleodata_e[1781:2580,]
my_tim_anom_Dongge_paleodata_e_800=my_tim_anom_Dongge_paleodata_e[1781:2580,]

### Plot evolution of mean annual precipitation anomaly at the 3 sites for 2.58-0 Ma AP.  Need to do this twice in order to get the right legends (and then manipulate image in PowerPoint)

# Plot, with first 800 kya shown in a different colour, and with legend showing first 800 kya as well as proxy data

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Plots/Precip_d18OC_SEAsia_2.5MaBP_v1.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(times), my_tim_anom_Hulu_paleodata_e,type="l", cex=lns, col=col_Hulu, xlab="", ylab="", xlim=xlimits, ylim=ylimits_a, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

lines(-(times), my_tim_anom_Sanbao_paleodata_e, lwd=lns, col=col_Sanbao)
lines(-(times), my_tim_anom_Dongge_paleodata_e, lwd=lns, col=col_Dongge)
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")
polygon(c(-(times), rev(-(times))), c(my_tim_anom_Hulu_paleodata_e_lower_1sd,rev(my_tim_anom_Hulu_paleodata_e_upper_1sd)), col = rgb(col_shade_lightteal[1],col_shade_lightteal[2],col_shade_lightteal[3],col_shade_lightteal[4]), border = NA)
polygon(c(-(times), rev(-(times))), c(my_tim_anom_Sanbao_paleodata_e_lower_1sd,rev(my_tim_anom_Sanbao_paleodata_e_upper_1sd)), col = rgb(col_shade_lightblue[1],col_shade_lightblue[2],col_shade_lightblue[3],col_shade_lightblue[4]), border = NA)
polygon(c(-(times), rev(-(times))), c(my_tim_anom_Dongge_paleodata_e_lower_1sd,rev(my_tim_anom_Dongge_paleodata_e_upper_1sd)), col = rgb(col_shade_lightcyan[1],col_shade_lightcyan[2],col_shade_lightcyan[3],col_shade_lightcyan[4]), border = NA)
lines(-(times), my_tim_anom_Hulu_paleodata_e, lwd=lns, col=col_Hulu)
lines(-(times2), my_tim_anom_Hulu_paleodata_e_800, lwd=lns, col=col_Hulu_800)
lines(-(times), my_tim_anom_Sanbao_paleodata_e, lwd=lns, col=col_Sanbao)
lines(-(times2), my_tim_anom_Sanbao_paleodata_e_800, lwd=lns, col=col_Sanbao_800)
lines(-(times), my_tim_anom_Dongge_paleodata_e, lwd=lns, col=col_Dongge)
lines(-(times2), my_tim_anom_Dongge_paleodata_e_800, lwd=lns, col=col_Dongge_800)

axis(2, at = yticks_a, labels=ylabels_a, tck=0.01, lwd=lns, cex.axis=fs)
par(new=TRUE)

plot(-(data_d18O_Chengetal16_641kyr_BP_interp[,1]), data_d18O_Chengetal16_641kyr_BP_interp[,2],type="l", lty=2, cex=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits_b, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)
lines(xlimits,c(data_d18O_Chengetal16_641kyr_BP_interp[1,2],data_d18O_Chengetal16_641kyr_BP_interp[1,2]), lty=3, lwd=lns, col="grey49")

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(4, at = yticks_b, tck=0.01, lwd=lns, cex.axis=fs, padj=0.4)
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Precipitation anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("(mm/month)", side=2, line=2.25, font=1, cex=fs)
mtext(expression(paste(delta^{18},"O (\u2030, VPDB)")), side=4, line=2.5, font=1, cex=fs3, col=col_proxy)

legend((xlimits[2]-900), (((ylimits_b[2]-ylimits_b[1])*0.32)+ylimits_b[1]), c("Proxy record","Emulated data, Hulu cave (800 kya)","Emulated data, Sanbao cave (800 kya)","Emulated data, Dongge cave (800 kya)"), pch=NA, lty=c(2,1,1,1), lwd=lns, cex=fs3, pt.cex=pts, col=c(col_proxy, col_Hulu_800, col_Sanbao_800, col_Dongge_800), bty="n")

par(las=1)

dev.off()

# Plot again, still with first 800 kya shown in a different colour but with legend showing 2.58 Ma

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Plots/Precip_d18OC_SEAsia_2.5MaBP_v2.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(times), my_tim_anom_Hulu_paleodata_e,type="l", cex=lns, col=col_Hulu, xlab="", ylab="", xlim=xlimits, ylim=ylimits_a, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

lines(-(times), my_tim_anom_Sanbao_paleodata_e, lwd=lns, col=col_Sanbao)
lines(-(times), my_tim_anom_Dongge_paleodata_e, lwd=lns, col=col_Dongge)
lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")
polygon(c(-(times), rev(-(times))), c(my_tim_anom_Hulu_paleodata_e_lower_1sd,rev(my_tim_anom_Hulu_paleodata_e_upper_1sd)), col = rgb(col_shade_lightteal[1],col_shade_lightteal[2],col_shade_lightteal[3],col_shade_lightteal[4]), border = NA)
polygon(c(-(times), rev(-(times))), c(my_tim_anom_Sanbao_paleodata_e_lower_1sd,rev(my_tim_anom_Sanbao_paleodata_e_upper_1sd)), col = rgb(col_shade_lightblue[1],col_shade_lightblue[2],col_shade_lightblue[3],col_shade_lightblue[4]), border = NA)
polygon(c(-(times), rev(-(times))), c(my_tim_anom_Dongge_paleodata_e_lower_1sd,rev(my_tim_anom_Dongge_paleodata_e_upper_1sd)), col = rgb(col_shade_lightcyan[1],col_shade_lightcyan[2],col_shade_lightcyan[3],col_shade_lightcyan[4]), border = NA)
lines(-(times), my_tim_anom_Hulu_paleodata_e, lwd=lns, col=col_Hulu)
lines(-(times2), my_tim_anom_Hulu_paleodata_e_800, lwd=lns, col=col_Hulu_800)
lines(-(times), my_tim_anom_Sanbao_paleodata_e, lwd=lns, col=col_Sanbao)
lines(-(times2), my_tim_anom_Sanbao_paleodata_e_800, lwd=lns, col=col_Sanbao_800)
lines(-(times), my_tim_anom_Dongge_paleodata_e, lwd=lns, col=col_Dongge)
lines(-(times2), my_tim_anom_Dongge_paleodata_e_800, lwd=lns, col=col_Dongge_800)

axis(2, at = yticks_a, labels=ylabels_a, tck=0.01, lwd=lns, cex.axis=fs)
par(new=TRUE)

plot(-(data_d18O_Chengetal16_641kyr_BP_interp[,1]), data_d18O_Chengetal16_641kyr_BP_interp[,2],type="l", lty=2, cex=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits_b, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)
lines(xlimits,c(data_d18O_Chengetal16_641kyr_BP_interp[1,2],data_d18O_Chengetal16_641kyr_BP_interp[1,2]), lty=3, lwd=lns, col="grey49")

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(4, at = yticks_b, tck=0.01, lwd=lns, cex.axis=fs, padj=0.4)
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Precipitation anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("(mm/month)", side=2, line=2.25, font=1, cex=fs)
mtext(expression(paste(delta^{18},"O (\u2030, VPDB)")), side=4, line=2.5, font=1, cex=fs3, col=col_proxy)

legend((xlimits[2]-900), (((ylimits_b[2]-ylimits_b[1])*0.32)+ylimits_b[1]), c("Emulated data, Hulu cave (1780 kya)","Emulated data, Sanbao cave (1780 kya)","Emulated data, Dongge cave (1780 kya)"), pch=NA, lty=c(1,1,1), lwd=lns, cex=fs3, pt.cex=pts, col=c(col_Hulu, col_Sanbao,col_Dongge), bty="n")

par(las=1)

dev.off()

#########################################################################################

} # Close figure saving

print("Finished")
