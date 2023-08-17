## Program to read in new proxy data (from Rachel) and interpolate to 1000 kyr resolution, read in emulated temperature, find relevant locations, and plot (anomalies only)
## Working (as of /8/23)

.libPaths(c("C:/Users/cw18831/OneDrive - University of Bristol/Documents/R/win-library/4.1", .libPaths()))

# Need to run this once the first time:
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
library(ncdf4)

require(colorspace)
require(graphics)
require(corrplot)

source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_parula.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_bwr.R')
source('C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/R/col_wr.R')

#################### Read in proxy data and interpolate to 1 kyr intervals ##############################

print("Reading data")

nloc = 23

# Files contain: Site [,1], age [,2], proxy type [,3], latitude [,4], longitude [,5], temporal resolution [,6], SST [,7], SST 1 SD [,8], SST anomaly [,9]

times = matrix(seq(-2579, 0, 1)) # Create times for the last 2.58 Myr, starting at -2.58
times2 = matrix(seq(-799, 0, 1)) # Create times for the last 800 kyr, starting at -800
proxy_anom_temp = array(0,c(2580,nloc)) # Create array to hold interpolated (to 1 kyr) SST anomalies for all timeslices for the last 2.58 Myr (rows) by nloc locations (columns)

all_filenames <- c("DSDP594_MAT", "geoB1008_UK37", "geoB1016_UK37", "geoB1028_UK37", "geob3327_UK37", "lpaz21p_UK37", "MD01-2414_tex86", "MD01-2414_UK37", "MD03-2611_UK37", "ODP1012_UK37", "ODP1082_UK37", "ODP1094_tex86", "ODP1123_MAT", "ODP1143_UK37", "ODP1239_UK37a", "ODP1239_UK37b", "ODP806_MgCa", "ODP871_MgCa", "PC1_UK37", "PS75_UK37", "S0164_MgCa", "U1313_UK37", "U1385_UK37")

for (loc in 1:nloc){         # Begin filename loop
  
  print(loc)
  print(all_filenames[loc])
  
  filename=file.path(paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/Papers/My papers/Charlie/Work/New proxy from Rachel/Individual locations/",all_filenames[loc],".txt", sep=""))
  cont_paramdat <- read.table(filename, sep="")
  print(filename)
  
  data = data.matrix(cont_paramdat)
  data[,2]=data[,2]*-1.0 # Multiply years by -1 so they are negative
  
  data_interp = approx(data[,2],data[,9], xout=rev(times)) # Interpolate proxy data anomalies to 1 kyr resolution
  data_anom_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))
  data_anom_interp = data_anom_interp[,2]
  
  proxy_anom_temp[,loc] = data_anom_interp
  
  rm(list=ls()[(ls() %in% c('cont_paramdat','filename'))])

} # Close filename loop

#################### Read in emulated data and find relevant locations ##############################

print("Reading emulated data")

# Open emulated SAT anomalies.  Here, initially, data are upside down.  Also, time slices are reversed i.e. time 1 = -2579 kyr

nc_file = nc_open("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_temp_anomalies_Global_2.58_to_0MyrBP.nc")
longs = ncvar_get(nc_file,"longitude")
lats = ncvar_get(nc_file,"latitude")
time_emul = ncvar_get(nc_file,"time")
emul_anom_temp = ncvar_get(nc_file,"SAT")

nx = length(longs)
ny = length(lats)
nt = length(time_emul)

# Open emulated SAT variances.  Here, initially, data are upside down.  Also, time slices are reversed i.e. time 1 = -2579 kyr

nc_file = nc_open("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Data/Global/Em_output_data_temp_variances_Global_2.58_to_0MyrBP.nc")
emul_var_temp = ncvar_get(nc_file,"SAT")

# Flip data to be right way up, for all timesteps

emul_anom_temp_new = array(0,c(nx,ny,nt))
emul_var_temp_new = array(0,c(nx,ny,nt))

n = ny

for (y in c(1:ny)){
  emul_anom_temp_new[1:nx,y,] = emul_anom_temp[(1:nx),n,]
  emul_var_temp_new[1:nx,y,] = emul_var_temp[(1:nx),n,]
  n = n-1
}

# Read locations (GB) at emulator resolution, i.e. GB corresponding to location closest to actual proxy location (manually calculated and listed in
# ~/New proxy from Rachel/New proxy locations.txt, contained in proxy_gb.txt)

filename=file.path("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/Papers/My papers/Charlie/Work/New proxy from Rachel/proxy_gb.txt")
cont_paramdat <- read.table(filename, sep="")
data = data.matrix(cont_paramdat)

emul_loc = array(0,c(nloc,2)) # Create array to hold all GB i.e. nloc locations (rows) by 2 GB (lon, lat)

emul_loc[,1] = data[,1]
emul_loc[,2] = data[,2]

#for (loc in c(1:nloc)){ # Check that GB match lon/lat at emulator resolution (manually calculated in proxy_lonlat_HadCM3res.txt).  And indeed they do!
# print(longs[emul_loc[loc,1]])
# print(lats[emul_loc[loc,2]])
#}

#emul_loc[1,1] = 45 # Compare with existing figure using original data e.g. ODP 982, where lon,lat = 45,14 (because latitudes are flipped here).  Should be identical, and indeed it is!
#emul_loc[1,2] = 60

# Extract data

emul_anom_temp_loc = array(0,c(nt,nloc)) # Create array to hold emulated SAT anomalies/variances for all timeslices for the last 2.58 Myr (rows) by nloc locations (columns)
emul_var_temp_loc = array(0,c(nt,nloc)) 

for (loc in c(1:nloc)){ 
  
  emul_anom_temp_loc[,loc] = emul_anom_temp_new[emul_loc[loc,1],emul_loc[loc,2],] # Save data at every location for all timeslices
  emul_var_temp_loc[,loc] = emul_var_temp_new[emul_loc[loc,1],emul_loc[loc,2],]
  
}

# Extract first 800 kya (i.e. last 800 values from timeseries, because timeslices are reversed), anomalies only

emul_anom_temp_loc_800 = array(0,c(800,nloc)) # Create array to hold emulated SAT for all timeslices for the last 800 kyr (rows) by nloc locations (columns)

for (loc in c(1:nloc)){ 
  
  emul_anom_temp_loc_800[,loc] = emul_anom_temp_loc[1781:2580,loc] # Save data at every location for all timeslices
  
}

####### Calculate error bars at every location (where lower = anomalies - square root of variances, upper = anomalies + square root of variances) ##########

lower_1sd_loc = array(0,c(nt,nloc)) # Create array to hold lower error bars for all timeslices for the last 2.58 Myr (rows) by nloc locations (columns)
upper_1sd_loc = array(0,c(nt,nloc)) # Create array to hold upper error bars for all timeslices for the last 2.58 Myr (rows) by nloc locations (columns)

for (loc in c(1:nloc)){ 
  
  lower_1sd_loc[,loc] = emul_anom_temp_loc[,loc]-sqrt(emul_var_temp_loc[,loc])
  upper_1sd_loc[,loc] = emul_anom_temp_loc[,loc]+sqrt(emul_var_temp_loc[,loc])
  
}

#################### Visualise ##############################

print("Creating plots")

# Setup constants for all figures

xlimits = c(0, 2580)
ylimits = c(-12, 20)

xticks = seq(0, 2580, by = 250)
yticks=seq(-12, 20, by=2)

xlabels = c("0", "", "500", "", "1000", "", "1500", "", "2000", "", "2500")
ylabels = c("-12", "", "-8", "", "-4", "", "0", "", "4", "", "8", "", "12", "", "16", "", "20")

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
col_emul = "blue"
col_emul_800 = "orange"
col_shade_lightblue = c(col2rgb("royalblue1")/255, opacity_shading)

# Plot, for every location individually, with first 800 kya shown in a different colour

for (loc in 1:nloc){         # Begin filename loop
  
png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output 2Ma/Plots/SST_anom_",all_filenames[loc],".png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(times), emul_anom_temp_loc[,loc], type="l", cex=lns, col=col_emul, xlab="", ylab="", xlim=xlimits, ylim=ylimits, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")
polygon(c(-(times), rev(-(times))), c(lower_1sd_loc[,loc],rev(upper_1sd_loc[,loc])), col = rgb(col_shade_lightblue[1],col_shade_lightblue[2],col_shade_lightblue[3],col_shade_lightblue[4]), border = NA)
lines(-(times), emul_anom_temp_loc[,loc], lwd=lns, col=col_emul)
lines(-(times2), emul_anom_temp_loc_800[,loc], lwd=lns, col=col_emul_800)

lines(-(rev(times)),proxy_anom_temp[,loc], lty=2, lwd=lns, col=col_proxy)

axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(2, at = yticks, labels=ylabels, tck=0.01, lwd=lns, cex.axis=fs)
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Temperature anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("("~ degree*C*")", side=2, line=2.25, font=1, cex=fs)
par(las=1)

dev.off()

print(loc)

}

print("Finished")

