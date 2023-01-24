###################################################################################################
## MAKE MISCELLANEOUS PLOTS #######################################################################
## For Pleistocene paper 2020-02-03 ###############################################################
###################################################################################################

## CW's version 2 of original extra plot code, i.e. to make the multi-forcing figures, now with unnecessary figures cut and including 3 more
## sensitivity experiments (so 10 in total). Working (as of 4/11/22)

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

forcings_included = c("F", "Foo", "Foew", "Fo", "Foepc", "Foeps", "Focs", "Fc", "Fs", "Fcs")

exp_name = c(expression(paste(bold("F"))), expression(paste(bold("F"[Oo]))), expression(paste(bold("F"[Oew]))), expression(paste(bold("F"[O]))), expression(paste(bold("F"[Oepc]))), expression(paste(bold("F"[Oeps]))), expression(paste(bold("F"[Ocs]))), expression(paste(bold("F"[C]))), expression(paste(bold("F"[S]))), expression(paste(bold("F"[Cs]))))

num_exps = length(exp_name)


site_list = c("Hul", "San", "Don")

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

    file_nam = paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Data/Em_output_data_precip_",forcings_included[m],"_",site_list[n],"_natural_0.8_to_0MyrBP.txt", sep="")
    data_all = read.table(file_nam, sep=" ", header=TRUE, numerals = "no.loss")

    data[,m,n] = data_all[,2] 
      
  }
}


single_forcing_paleodata_Hulu_800kyr_BP = matrix(data[,,1], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_Hulu_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_Hulu_800kyr_BP) <- c(times)

single_forcing_paleodata_Sanbao_800kyr_BP = matrix(data[,,2], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_Sanbao_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_Sanbao_800kyr_BP) <- c(times)

single_forcing_paleodata_Dongge_800kyr_BP = matrix(data[,,3], nrow = num_times, ncol = num_exps)
colnames(single_forcing_paleodata_Dongge_800kyr_BP) <- c(forcings_included)
rownames(single_forcing_paleodata_Dongge_800kyr_BP) <- c(times)

print("Loading proxy data")

# Read in proxy data and interpolate to 1 kyr intervals ----------------------------------------------------------------------------

# Read in d18O data for 641-0 kyr BP [Wang et al 2001, 2008]
# Original D + H data has been made 1.6 %o more negative to account for higher values at D + H than SB [Wang et al 2008, Cheng et al. 2016]

source("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Emulator/2015_Bristol_5D_v001/Data/2018-08-01 Final report/Import_d18O_Chengetal16_641kyr_BP_cjrw.R") # Composite record, China (Sanbao, Hulu, Dongge)

data_interp = approx(data_d18O_Chengetal16_641kyr_BP[,2],data_d18O_Chengetal16_641kyr_BP[,3], xout=rev(times)) # Interpolate proxy data to 1 kyr resolution for calculating correlation coefficient
data_d18O_Chengetal16_641kyr_BP_interp = array(c(data_interp[1]$x,data_interp[2]$y),c(length(rev(times)),2))

# plot(data_d18O_Chengetal16_641kyr_BP[,2],data_d18O_Chengetal16_641kyr_BP[,3], type="l", xlim=c(min(times), max(times)), ylim=c(-12,-4))
# plot(data_d18O_Chengetal16_641kyr_BP_interp[,1],data_d18O_Chengetal16_641kyr_BP_interp[,2], type="l", ylim=c(-12,-4))

########################################################################################################################################################
## FIGURES #############################################################################################################################################
########################################################################################################################################################

print("Creating figures")

# Hulu cave, China

xlimits = c(0, 800)
ylimits_a = c(-30,30)
#ylimits_b = c(-4.5,-12)
ylimits_b = c(-4.5,-14.8)

xticks = seq(0, 800, by = 50)
yticks_a = seq(-30, 30, by = 5)
#yticks_b = seq(-4.5, -12, by = -0.5)
yticks_b = seq(-4.5, -14.8, by = -0.5)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels_a = c("-30", "", "-20", "", "-10", "", "0", "", "10", "", "20", "", "30")
#ylabels_b = c("", "", "", "-6", "", "", "", "-8", "", "", "", "-10", "", "", "", "-12")
ylabels_b = c("", "", "", "-6", "", "", "", "-8", "", "", "", "-10", "", "", "", "-12", "", "", "", "-14", "")


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
col_Focs = "blue"
col_Fooewcs_list = c("green", "blue", "purple", "pink", "red", "cyan", "black", "orange", "yellow", "chocolate4")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/Fig7_Ts_paleo_Ppt_Fsep6_Hulu_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(times), single_forcing_paleodata_Hulu_800kyr_BP[,1],type="l", cex=lns, col=col_Focs, xlab="", ylab="", xlim=xlimits, ylim=ylimits_a, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

## CW: Instead of plotting Focs first (in blue) and then plotting all the others, plots all of them together with new colouring system

#lines(-(times), single_forcing_paleodata_Hulu_800kyr_BP[,1],type="l", cex=lns, col=col_Focs)
#for (m in 2:num_exps){
#  lines(-(times), single_forcing_paleodata_Hulu_800kyr_BP[,m],type="l", cex=lns, col=col_Fooewcs_list[m-1])
#}
#lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")
#lines(-(times), single_forcing_paleodata_Hulu_800kyr_BP[,1],type="l", cex=lns, col=col_Focs)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_Hulu_800kyr_BP[,m],type="l", cex=lns, col=col_Fooewcs_list[m])
}

lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(2, at = yticks_a, labels=ylabels_a, tck=0.01, lwd=lns, cex.axis=fs)
par(new=TRUE)
plot(-(data_d18O_Chengetal16_641kyr_BP_interp[,1]),data_d18O_Chengetal16_641kyr_BP_interp[,2], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits_b, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)
lines(xlimits,c(data_d18O_Chengetal16_641kyr_BP_interp[1,2],data_d18O_Chengetal16_641kyr_BP_interp[1,2]), lty=3, lwd=lns, col="grey49")
axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(4, at = yticks_b, labels=ylabels_b, tck=0.01, lwd=lns, cex.axis=fs, padj=0.4)
text((xlimits[1]+10), (((ylimits_b[2]-ylimits_b[1])*0.96)+ylimits_b[1]), "a)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Precipitation anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("(mm/month)", side=2, line=2.25, font=1, cex=fs)
mtext(expression(paste(delta^{18},"O (\u2030, VPDB)")), side=4, line=2.5, font=1, cex=fs, col=col_proxy)

legend((xlimits[2]-387), (((ylimits_b[2]-ylimits_b[1])*1.08)+ylimits_b[1]), c("Proxy", "E", "Eob", "Eep", "Eob,ep", "Eob,ep,co", "Eob,ep,ic", "Eob,ep,co,ic", "Eco", "Eic", "Eco,ic"), pch=NA, lty=c(2,1,1,1,1,1,1,1), lwd=lns, cex=fs3, pt.cex=pts, adj = 0.05, col=c(col_proxy, col_Fooewcs_list), ncol=2, bty="n")

par(las=1)

dev.off()

# Sanbao cave, China

xlimits = c(0, 800)
ylimits_a = c(-30,30)
#ylimits_b = c(-4.5,-12)
ylimits_b = c(-4.5,-14.8)

xticks = seq(0, 800, by = 50)
yticks_a = seq(-30, 30, by = 5)
#yticks_b = seq(-4.5, -12, by = -0.5)
yticks_b = seq(-4.5, -14.8, by = -0.5)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels_a = c("-30", "", "-20", "", "-10", "", "0", "", "10", "", "20", "", "30")
#ylabels_b = c("", "", "", "-6", "", "", "", "-8", "", "", "", "-10", "", "", "", "-12")
ylabels_b = c("", "", "", "-6", "", "", "", "-8", "", "", "", "-10", "", "", "", "-12", "", "", "", "-14", "")


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
col_Focs = "blue"
col_Fooewcs_list = c("green", "blue", "purple", "pink", "red", "cyan", "black", "orange", "yellow", "chocolate4")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/Fig7_Ts_paleo_Ppt_Fsep6_Sanbao_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(times), single_forcing_paleodata_Sanbao_800kyr_BP[,1],type="l", cex=lns, col=col_Focs, xlab="", ylab="", xlim=xlimits, ylim=ylimits_a, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_Sanbao_800kyr_BP[,m],type="l", cex=lns, col=col_Fooewcs_list[m])
}

lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(2, at = yticks_a, labels=ylabels_a, tck=0.01, lwd=lns, cex.axis=fs)
par(new=TRUE)
plot(-(data_d18O_Chengetal16_641kyr_BP_interp[,1]),data_d18O_Chengetal16_641kyr_BP_interp[,2], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits_b, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)
lines(xlimits,c(data_d18O_Chengetal16_641kyr_BP_interp[1,2],data_d18O_Chengetal16_641kyr_BP_interp[1,2]), lty=3, lwd=lns, col="grey49")
axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(4, at = yticks_b, labels=ylabels_b, tck=0.01, lwd=lns, cex.axis=fs, padj=0.4)
text((xlimits[1]+10), (((ylimits_b[2]-ylimits_b[1])*0.96)+ylimits_b[1]), "b)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Precipitation anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("(mm/month)", side=2, line=2.25, font=1, cex=fs)
mtext(expression(paste(delta^{18},"O (\u2030, VPDB)")), side=4, line=2.5, font=1, cex=fs, col=col_proxy)

# CW: Don't show legend here
#legend((xlimits[2]-267), (((ylimits_b[2]-ylimits_b[1])*1.08)+ylimits_b[1]), c("Proxy", "E", "Eob", "Eep", "Eob,ep", "Eob,ep,co", "Eob,ep,ic", "Eob,ep,co,ic", "Eco", "Eic", "Eco,ic"), pch=NA, lty=c(2,1,1,1,1,1,1,1), lwd=lns, cex=fs3, pt.cex=pts, adj = 0.05, col=c(col_proxy, col_Fooewcs_list), ncol=2, bty="n")

par(las=1)

dev.off()

# Dongge cave, China

xlimits = c(0, 800)
ylimits_a = c(-30,30)
#ylimits_b = c(-4.5,-12)
ylimits_b = c(-4.5,-14.8)

xticks = seq(0, 800, by = 50)
yticks_a = seq(-30, 30, by = 5)
#yticks_b = seq(-4.5, -12, by = -0.5)
yticks_b = seq(-4.5, -14.8, by = -0.5)

xlabels = c("0", "", "100", "", "200", "", "300", "", "400", "", "500", "", "600", "", "700", "", "800")
ylabels_a = c("-30", "", "-20", "", "-10", "", "0", "", "10", "", "20", "", "30")
#ylabels_b = c("", "", "", "-6", "", "", "", "-8", "", "", "", "-10", "", "", "", "-12")
ylabels_b = c("", "", "", "-6", "", "", "", "-8", "", "", "", "-10", "", "", "", "-12", "", "", "", "-14", "")


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
col_Focs = "blue"
col_Fooewcs_list = c("green", "blue", "purple", "pink", "red", "cyan", "black", "orange", "yellow", "chocolate4")

png(file=paste("C:/Users/cw18831/OneDrive - University of Bristol/Documents/Research/SKB/SKB_Alan_Code_DJL/PosivaSKB/PosivaSKB/Palaeo output/Plots/Fig7_Ts_paleo_Ppt_Fsep6_Dongge_800kyrBP.png", sep=""),width=3000,height=850,res=300)
par(fig=c(lb,rb,(bb),(bb+ht)), mar=c(2,4.5,0.5,3), pty="m", mgp=c(2,0.5,0), las=1, lwd=lns)
plot(-(times), single_forcing_paleodata_Dongge_800kyr_BP[,1],type="l", cex=lns, col=col_Focs, xlab="", ylab="", xlim=xlimits, ylim=ylimits_a, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)

for (m in 1:num_exps){
  lines(-(times), single_forcing_paleodata_Dongge_800kyr_BP[,m],type="l", cex=lns, col=col_Fooewcs_list[m])
}

lines(xlimits,c(0,0), lty=3, lwd=lns, col="grey49")

axis(2, at = yticks_a, labels=ylabels_a, tck=0.01, lwd=lns, cex.axis=fs)
par(new=TRUE)
plot(-(data_d18O_Chengetal16_641kyr_BP_interp[,1]),data_d18O_Chengetal16_641kyr_BP_interp[,2], type="l", lty=2, lwd=lns, col=col_proxy, xlab="", ylab="", xlim=xlimits, ylim=ylimits_b, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE)
lines(xlimits,c(data_d18O_Chengetal16_641kyr_BP_interp[1,2],data_d18O_Chengetal16_641kyr_BP_interp[1,2]), lty=3, lwd=lns, col="grey49")
axis(1, at = xticks, labels=xlabels, tck=0.01, lwd=lns, cex.axis=fs)
axis(4, at = yticks_b, labels=ylabels_b, tck=0.01, lwd=lns, cex.axis=fs, padj=0.4)
text((xlimits[1]+10), (((ylimits_b[2]-ylimits_b[1])*0.96)+ylimits_b[1]), "c)", font=2, cex=fs, adj = c(0,NA))
par(las=0, xpd=TRUE,lwd=lns)
mtext("Time BP (kyr)", side=1, line=1.5, font=1, cex=fs)
mtext("Precipitation anomaly", side=2, line=3.5, font=1, cex=fs)
mtext("(mm/month)", side=2, line=2.25, font=1, cex=fs)
mtext(expression(paste(delta^{18},"O (\u2030, VPDB)")), side=4, line=2.5, font=1, cex=fs, col=col_proxy)

# CW: Don't show legend here
#legend((xlimits[2]-267), (((ylimits_b[2]-ylimits_b[1])*1.08)+ylimits_b[1]), c("Proxy", "E", "Eob", "Eep", "Eob,ep", "Eob,ep,co", "Eob,ep,ic", "Eob,ep,co,ic", "Eco", "Eic", "Eco,ic"), pch=NA, lty=c(2,1,1,1,1,1,1,1), lwd=lns, cex=fs3, pt.cex=pts, adj = 0.05, col=c(col_proxy, col_Fooewcs_list), ncol=2, bty="n")

par(las=1)

dev.off()

print("Finished")

#########################################################################################