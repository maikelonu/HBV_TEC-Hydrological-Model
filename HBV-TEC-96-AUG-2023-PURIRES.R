#///////////////////////////////////////////////////////////////////////////////////////////
# HBV-TEC-96 HYDROLOGICAL MODEL
#///////////////////////////////////////////////////////////////////////////////////////////
# VERSION: 1.0.2016
# Instituto Tecnologico de Costa Rica (www.tec.ac.cr)
# Maikel Mendez-M (mamendez@itcr.ac.cr);(maikel.mendez@gmail.com)
# Luis Alexander Calvo-V (lcalvo@itcr.ac.cr);(lualcava.sa@gmail.com)
# This script is structured in R (www.r-project.org)
# General purpose: Implementation of the SMHI-HBV-96 Hydrological Model on R
# Input files: "hbvtecptq.txt", "hbvtecpar.txt", "hbvtecevap.txt", "hbvtecattri.txt" 
# Output files: "hbvtecptq_desc_hbv.csv", "hbvtecout_hbv.csv", "hbvtecout_desc_hbv.csv",
# "hbvqtecsim.csv", "hbvqteceff.csv"
#
#///////////////////////////////////////////////////////////////////////////////////////////

# REFERENCES:
#
# Aghakouchak, A., Habib, E. 2010. Application of a conceptual hydrologic model in 
# teaching hydrologic processes. Int J Engng Ed 26: 963-973.
# 
# Bergström, S. 1995. The HBV Model. In: V.P. Singh (editor) Computer Models of Watershed
# Hydrology. Water Resources Publications, Highlands Ranch, Colorado pp. 443- 470.
#
# Herman, J.D., Reed, P.M., Wagener, T. 2013. Time-varying sensitivity analysis 
# clarifies the effects of watershed model formulation on model behavior.
# Water Resour Res 49: 1400-1414.
#
# Lindström, G., Johansson, B., Persson, M., Gardelin, M., Bergström, S. 1997. Development
# and test of the distributed HBV-96 hydrological model. J Hydrol 201: 272-288.
# 
# Seibert, J. Vis, MJP. 2012. Teaching hydrological modeling with a user-friendly 
# catchment-runoff-model software package. Hydrol Earth Syst Sci 16: 3315-3325.
#
#///////////////////////////////////////////////////////////////////////////////////////////

# Workspace is cleared
# rm(list = ls())

# CRAN libraries are loaded
require(DescTools)
require(ggplot2)
require(grid)
require(gridExtra)
require(lubridate)
require(pastecs)
require(RColorBrewer)
require(reshape)
require(visreg)

# Working directory is defined
setwd("C:/DATOS/MPE_2018/R_Workshop_UCR_2018/TOPIC11/HBV_TEC_CARTAGO/PURIRES")

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Creating and organizing input data.frames
# ///////////////////////////////////////////////////////////////////////////////

# Precipitation, temperature and Q-observed data.frame is loaded
df.ptq <- read.delim ("hbvtecptq.txt", header = TRUE, sep = "\t")

# Model input-parameters data.frame is loaded
df.param <- read.delim ("hbvtecpar.txt", header = TRUE, sep = "\t")

# Monthly potential evapotranspiration data.frame is loaded
df.evap <- read.delim ("hbvtecevap.txt", header = TRUE, sep = "\t")

# Watershed attributes data.frame is loaded
df.attri <- read.delim ("hbvtecattri.txt", header = TRUE, sep = "\t")

# Desc {DescTools} function is requested
Desc(df.ptq, plotit = TRUE) # Optional !!

# Descriptive statistics are requested and rounded to 4 decimals
df.ptq.desc <- round((as.data.frame(stat.desc(df.ptq[, 2 : 4]))),4)

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Resetting and defining and data containers 
# ///////////////////////////////////////////////////////////////////////////////

# Output data.frame containers are reset to NULL values
df.out01 <- NULL
df.out02 <- NULL
df.out03 <- NULL

# Model state-variables containers are defined
Rstore <- vector()    # mass storage [mm/T]
SM <- vector()        # soil moisture [mm/T]
AET <- vector()       # actual evapotranspiration [mm/T]
PET <- vector()       # potential evapotranspiration [mm/T]
R <- vector()         # recharge [mm/T]
SUZ <- vector()       # storage in the upper zone [mm/T]
SLZ <- vector()       # storage in the lower zone [mm/T]
Q0 <- vector()        # mass coming from SUZ [mm/T] (if SUZ > uzl)
Q1 <- vector()        # mass coming from SUZ [mm/T]
Q2 <- vector()        # mass coming from SLZ [mm/T]
QSR <- vector()       # mass coming from surface runoff [mm/T] if (SM >= fc)
MASSIN <- vector()    # mass entering the loop as precipitation [mm/T]
MASSOUT <- vector()   # mass exiting the loop [mm/T]
DELTAMASS <- vector() # delta of mass within the loop (MASSIN - MASSOUT) [mm/T]
QOT <- vector()       # total mass summation (Q0 + Q1 + Q2) [mm/T]
QOUL <- vector()      # mass summation of upper and lower zones (Q1 + Q2) [mm/T]
QOSR <- vector()      # mass summation of surface runoff and uzl (Q0 + QSR) [mm/T]

# Model transformation-variables containers are defined
mxrelweisum <- 0     # maxbas relative weigths summation [unitless]
mxabsweisum <- 0     # maxbas absolute weigths summation [unitless]
mxrelwei <- vector() # maxbas relative weigth [unitless]
mxabswei <- vector() # maxbas absolute weigth [unitless]
QSIM <- vector()     # model simulated mass [mm/T]
QSIM3 <- vector()    # model simulated mass [m3/sec]

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Loading input variables and model parameters
# ///////////////////////////////////////////////////////////////////////////////

# Model input parameters are loaded
fc <- df.param [1, 2]      # soil field capacity [mm]
lp <- df.param [2, 2]      # threshold at which AET reaches PET [unitless]
beta <- df.param [3, 2]    # soil shape calibration parameter [unitless]
perc <- df.param [4, 2]    # percolation rate [mm/T]
uzl <- df.param [5, 2]     # threshold parameter for quick flow [mm]
k0 <- df.param [6, 2]      # recession coefficient for upper zone [1/T]
k1 <- df.param [7, 2]      # recession coefficient for upper zone [1/T]
k2 <- df.param [8, 2]      # recession coefficient for lower zone [1/T]
maxbas <- df.param [9, 2]  # length of weighting  transformation function [T]

# Monthly potential evapotranspiration containers are loaded [mm/month]
evap.jan <- df.evap [1, 2]
evap.feb <- df.evap [2, 2]
evap.mar <- df.evap [3, 2]
evap.apr <- df.evap [4, 2]
evap.may <- df.evap [5, 2]
evap.jun <- df.evap [6, 2]
evap.jul <- df.evap [7, 2]
evap.aug <- df.evap [8, 2]
evap.set <- df.evap [9, 2]
evap.oct <- df.evap [10, 2]
evap.nov <- df.evap [11, 2]
evap.dic <- df.evap [12, 2]

# Watershed attributes are loaded
warea <- as.numeric(df.attri [1, 2])      # catchment area (km2)
fract_SMi <- as.numeric(df.attri [2, 2])  # initial SM fraction [unitless]
fract_SLZi <- as.numeric(df.attri [3, 2]) # initial SLZ fraction [unitless]
resolution <- as.numeric(df.attri [4, 2]) # model temporal resolution (1 = daily; 2 = hourly)
routing <- as.numeric(df.attri [5, 2])    # maxbas selection-control parameter 

# Numeric counters are defined
counterlength <- seq(1, (length(df.ptq$QOBS)), by = 1)
countermain <- length (counterlength)

# A SEQ variable is created at df.ptq data.frame
df.ptq$SEQ <- seq(1, (length(df.ptq$QOBS)), by = 1)

# df.ptq$DATE factor class is converted to date class 
DATEtemp01 <- df.ptq$DATE
df.ptq$DATE <- as.Date(DATEtemp01, format = "%d/%m/%Y")

# lubridate Library functions are loaded
df.ptq$YEAR <- year(df.ptq$DATE)                    # year component of DATE
df.ptq$YEAR_CH <- as.character(year(df.ptq$DATE))   # year component of DATE as character
df.ptq$MONTH <- month(df.ptq$DATE, label = FALSE)   # month component of DATE
df.ptq$MONTH_CH <- month(df.ptq$DATE, label = TRUE) # month component of DATE as character
df.ptq$WEEK <- week(df.ptq$DATE)                    # week component of DATE
df.ptq$DAY <- yday(df.ptq$DATE)                     # day component of a DATE
df.ptq$DAY_MONTH <- days_in_month(df.ptq$DATE)      # number of days in the month of DATE

# A selection data.frame is created based on df.ptq data.frame
df.selection <- df.ptq[c("SEQ", "DATE", "YEAR", "YEAR_CH", "MONTH", "MONTH_CH", "DAY")]

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Model state variables are initialized
# ///////////////////////////////////////////////////////////////////////////////
SMtemp  <- (fract_SMi * fc)           # temporal soil moisture [mm/T]
SLZtemp <- (fract_SLZi * (perc / k2)) # temporal storage in the lower zone [mm/T]
SUZtemp <- 0                          # temporal storage in the upper zone [mm/T]

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Main mass-balance loop
# ///////////////////////////////////////////////////////////////////////////////

# Main mass-balance loop is initialized
for (i in 1 : countermain) {
  
  # ----------------------------------------------------------
  # Evapotranspiration sorting
  # ----------------------------------------------------------
  
  # PET values are selected [mm/month] (according to month of the year and date)
  if (df.ptq$MONTH[i] == 1) { PET[i] <- evap.jan } else if (df.ptq$MONTH[i] == 2) {
    PET[i] <- evap.feb } else if (df.ptq$MONTH[i] == 3) {
      PET[i] <- evap.mar } else if (df.ptq$MONTH[i] == 4) {
        PET[i] <- evap.apr } else if (df.ptq$MONTH[i] == 5) {
          PET[i] <- evap.may } else if (df.ptq$MONTH[i] == 6) {
            PET[i] <- evap.jun } else if (df.ptq$MONTH[i] == 7) {
              PET[i] <- evap.jul } else if (df.ptq$MONTH[i] == 8) {
                PET[i] <- evap.aug } else if (df.ptq$MONTH[i] == 9) {
                  PET[i] <- evap.set } else if (df.ptq$MONTH[i] == 10) {
                    PET[i] <- evap.oct } else if (df.ptq$MONTH[i] == 11) {
                      PET[i] <- evap.nov } else if (df.ptq$MONTH[i] == 12) {
                        PET[i] <- evap.dic }
  
  
  # ----------------------------------------------------------
  # Precipitation routine
  # ----------------------------------------------------------
  
  # MASSIN is defined [mm/T] (entering the loop as precipitation/snow) 
  PRECtemp  <- df.ptq$PREC[i]    
  MASSIN[i] <- df.ptq$PREC[i]  
  
  # ----------------------------------------------------------
  # Soil moisture routine
  # ----------------------------------------------------------
  
  # This is for i = 1; meaning day = 1; meaning initialization
  if (i == 1) {
    SM[i] <- SMtemp
  } else { # This is for i > 1; meaning day > 1
    SMtemp2 <- SM[i - 1]
    SM[i] <- SMtemp2 
  }
  
  # The process continues...
  if (SM[i] >= fc) {
    # R[i] <- PRECtemp + (SM[i] - fc) # Original HBV statement (deprecated!!)
    QSR[i] <- SM[i] - fc # New HBV-TEC statement
    R[i] <- PRECtemp + QSR[i] # New HBV-TEC statement
    SM[i] <- fc
  } else {
    QSR[i] <- 0 # New HBV-TEC statement
    Rstore[i] <- PRECtemp * (1 - ((SM[i] / fc) ^ beta))
    SM[i] <- SM[i] + Rstore[i]
    R[i] <- PRECtemp - Rstore[i]
    if (SM[i] > fc) {
      R[i] <- SM[i] - fc
      SM[i] <- fc
    }
  }
  
  # AET compensation within soil moisture routine
  if ((SM[i] / fc) > lp) {
    AET[i] <- PET[i]
  } else {
    AET[i] <- (SM[i] / (fc * lp)) * PET[i]
  }
  if (AET[i] < 0) {
    AET[i] <- 0
  }
  if (SM[i] > AET[i]) {
    SM[i] <- SM[i] - AET[i]
  } else {
    AET[i] <- SM[i]
    SM[i] <- 0
  } 
  
  # ----------------------------------------------------------
  # Responce function routine
  # ----------------------------------------------------------
  
  # Storage in the Upper Zone (SUZ)
  # This is for i = 1; meaning day = 1; meaning initialization
  if(i == 1) {
    SUZ[i] <- SUZtemp + R[i]
  } else { # This is for i > 1; meaning day > 1
    SUZ[i] <- SUZ[i - 1] + R[i]
  }
  
  # The process continues...
  if (SUZ[i] > uzl) {
    Q0[i] <- k0 * (SUZ[i] - uzl)
    SUZ[i] <- SUZ[i] - Q0[i]
  } else {
    Q0[i] <- 0
  }
  
  # Percolation threshold
  if (SUZ[i] > perc) {
    SUZ[i] <- SUZ[i] - perc
    Q1[i] <- k1 * (SUZ[i])
    SUZ[i] <- SUZ[i] - Q1[i]
    
    # Storage in the Lower Zone (SLZ)
    # This is for i = 1; meaning day = 1; meaning initialization
    if (i == 1) {
      SLZ[i] <- SLZtemp + perc
    } else { 
      SLZ[i] <- SLZ[i - 1] + perc
    }
  } else {
    Q1[i] <- 0 # New HBV-TEC statement
    # This is for i = 1; meaning day = 1; meaning initialization
    if (i == 1) {
      SLZ[i] <- SLZtemp + SUZ[i]
    } else { 
      SLZ[i] <- SLZ[i - 1] + SUZ[i]
    }
    SUZ[i] <- 0
  }
  
  if (SLZ[i] > 0) {
    Q2[i] <- k2 * SLZ[i]
    SLZ[i] <- SLZ[i] - Q2[i]  
  } else {
    Q2[i] <- 0
    SLZ[i] <- SLZ[i - 1] # a warning should be issued here!!!!! 
  }
  
  # ----------------------------------------------------------
  # Check mass-balance loop 
  # ----------------------------------------------------------
  
  # This is for i = 1; meaning day = 1; meaning initialization
  if (i == 1) {
    
    # Mass exiting the loop is calculated [mm/T]
    MASSOUT[i] <- ((SM[i] - 0) +
                     (SUZ[i] - SUZ[i]) +
                     (SLZ[i] - SLZ[i]) +
                     (Q0[i]) +
                     (Q1[i]) +
                     (Q2[i]) +
                     (AET[i]))
  } else {
    # Mass exiting the loop is calculated [mm/T]
    MASSOUT[i] <- ((SM[i] - SM[i - 1]) +
                     (SUZ[i] - SUZ[i - 1]) +
                     (SLZ[i] - SLZ[i - 1]) +
                     (Q0[i]) +
                     (Q1[i]) +
                     (Q2[i]) +
                     (AET[i]))
  }
  
  # Delta of mass within the loop (MASSIN - MASSOUT) is calculated [mm/T]  
  DELTAMASS[i] <- MASSIN[i] - MASSOUT[i]
  
  # Total mass summation (Q0 + Q1 + Q2) is calculated [mm/T]      
  QOT[i] <- Q0[i] + Q1[i] + Q2[i] # original HBV statement
  
  # Mass summation of upper and lower zones (Q1 + Q2) is calculated [mm/T]
  QOUL[i] <- Q1[i] + Q2[i] # new HBV-TEC statement
  
  # Mass summation of surface runoff and uzl (Q0 + QSR) is calculated [mm/T]
  QOSR[i] <- Q0[i] + QSR[i] # new HBV-TEC statement
}
# ----------------------------------------------------------
# Main mass-balance loop is closed
# ----------------------------------------------------------

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Transformation function routine (maxbas-t) Q0+Q1+Q2
# ///////////////////////////////////////////////////////////////////////////////

# maxbas counters are defined
mx2 <- (maxbas / 2) # half maxbas
maxbascount <- ceiling(maxbas) # maxbas counter based on ceiling

# ----------------------------------------------------------
# maxbas relative weights are calculated
# ----------------------------------------------------------

# If maxbas < 2
if (maxbas < 2) {
  for (i in 1 : (maxbascount)) {
    if (i <= maxbas) {
      mxrelwei[i] <- (maxbas / 2) 
    } else {
      mxrelwei[i] <- (((maxbas - i) + 1.0)) 
    }
    mxrelwei2 <- (mxrelwei[i])
    mxrelweisum <- (mxrelwei[i]) + mxrelweisum
  }
}

# If maxbas >= 2
if (maxbas >= 2) {
  for (i in 1 : (maxbascount)) {    
    if (i <= mx2) {
      mxrelwei[i] <- (i) 
    } else {
      mxrelwei[i] <- (((maxbas - i) + 1.0))
    }
    mxrelwei2 <- (mxrelwei[i])
    mxrelweisum <- (mxrelwei[i]) + mxrelweisum             
  }
}

# ----------------------------------------------------------
# maxbas relative weights normalization loop
# ----------------------------------------------------------

# maxbas relative weights normalization loop is initialized
for (i in 1 : (maxbascount)) {
  mxabswei[i] <- (mxrelwei[i]) / mxrelweisum
  mxabsweisum <- (mxabswei[i]) + mxabsweisum
}

# Total duration of hydraulic routing + maxbas is defined
routingdur <- ((length(QOT)) + maxbascount - 1 + 1)

# A "qresult[i]" vector inizialized to cero for the duration of "routingdur" is created
qresult <- c(rep(0.0, routingdur))

# Duration of hydraulic routing ONLY is defined
QOTdur <- (length(QOT))

# ----------------------------------------------------------
# Flow integration over-time loop
# ----------------------------------------------------------

# Flow integration over-time loop is initialized
for (n in 1 : (QOTdur)) { # day of the year external Loop
  for (i in 1 : (maxbascount)) { # internal maxbas Loop                 
    qrouting <- ((QOT[n]) * (mxabswei[i]))
    qresult[n + i - 1] <- qresult[n + i - 1] + qrouting
  }
} 
# Flow integration over-time loop is closed

# Define QSIM based on "qresult[i]" for the duration of countermain ONLY 
QSIM <- (qresult[1 : countermain])

# Convert QSIM from mm/day to m3/sec
QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400

# If df.attri data.frame routing = 0; maxbas is IGNORED and NO tranformation is calculated
if (routing == 0) {
  QSIM <- Q0 + Q1 + Q2 + QSR
  QSIM3 <- (QSIM / 1000) * (warea * 1000000) / 86400
}

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Setup model outputs
# ///////////////////////////////////////////////////////////////////////////////

# Output variables are defined and rounded to 4 decimals
outQOBS <- round(df.ptq$QOBS, 4)
outPREC <- round(df.ptq$PREC, 4)
outTEMP <- round(df.ptq$TEMP, 4)
outPET <- round(PET, 4)
outAET <- round(AET, 4)
outMASSIN <- round(MASSIN, 4)
outRstore <- round(Rstore, 4)
outSM <- round(SM, 4)
outR <- round(R, 4)
outSUZ <- round(SUZ, 4)
outSLZ <- round(SLZ, 4)
outDELTAMASS <- round(DELTAMASS, 4) 
outQ0 <- round(Q0, 4)
outQ1 <- round(Q1, 4)
outQ2 <- round(Q2, 4)
outQSR <- round(QSR, 4)
outQOT <- round(QOT, 4)
outQSIM <- round(QSIM, 4)
outQRES <- round((df.ptq$QOBS - QSIM), 4)
outQSIM3 <- round(QSIM3, 4)

# Generalizaed mass-balance outputs are defined and rounded to 4 decimals
sumQOBS <- round((sum(outQOBS)), 4)
sumQSIM <- round((sum(outQSIM)), 4)
deltaQ <- round((sumQOBS - sumQSIM ), 4)
sumPREC <- round((sum(outPREC)), 4)
sumPET <- round((sum(outPET)), 4)
sumAET <- round((sum(outAET)), 4)
fractSUZ <- round(((sum(outQ0) + sum(outQ1)) / sum(outQOT)), 4)
fractSLZ <- round((sum(outQ2) / sum(outQOT)), 4)
fractUZL <- round((sum(outQ0) / sum(outQOT)), 4)
MAXBASabswei <- round(mxabswei, 4)
MAXBASabsweisum <- round(mxabsweisum, 4)

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Model goodness of fit (objective functions)
# ///////////////////////////////////////////////////////////////////////////////

# Vector containers are reset to NULL values
sumA <- NULL
sumB <- NULL
sumC <- NULL
sumD <- NULL
sumE <- NULL
sumF <- NULL
sumG <- NULL

# Vector containers are defined
sumA <- vector()
sumB <- vector()
sumC <- vector()
sumD <- vector()
sumE <- vector()
sumF <- vector()
sumG <- vector()

# ----------------------------------------------------------
# Generic objective-functions loop is initialized
# ----------------------------------------------------------

for (i in 1 : countermain) {
  tempsum <- ((outQOBS[i] - outQSIM[i]) ^ 2)
  sumA[i] <- tempsum 
  
  tempsum2 <- ((outQOBS[i] - (mean(outQOBS))) ^ 2)
  sumB[i] <- tempsum2
  
  tempsum3 <- (outQOBS[i] - (mean(outQOBS)))
  sumC[i] <- tempsum3
  
  tempsum4 <- (outQSIM[i] - (mean(outQSIM)))
  sumD[i] <- tempsum4
  
  tempsum5 <- ((outQSIM[i] - (mean(outQSIM))) ^ 2)
  sumE[i] <- tempsum5
  
  tempsum6 <- ((log(outQOBS[i]) - log(outQSIM[i])) ^ 2)
  sumF[i] <- tempsum6 
  
  tempsum7 <- ((log(outQOBS[i]) - (log(mean(outQOBS)))) ^ 2)
  sumG[i] <- tempsum7
}

# -----------------------------------------------------------------------------
# The Nash and Sutcliffe efficiency criterion - NSeff [fraction]
# -----------------------------------------------------------------------------
NSeff <- ( 1- (sum(sumA) / sum(sumB)))
NSeff <- round (NSeff, 4)

# -----------------------------------------------------------------------------
# The Nash and Sutcliffe efficiency logarithmic values - LNNSeff [fraction]
# -----------------------------------------------------------------------------
LNNSeff <- (1 - (sum(sumF) / sum(sumG)))
LNNSeff <- round (NSeff, 4)

# -----------------------------------------------------------------------------
# Correlation coefficient - R2 [fraction]
# -----------------------------------------------------------------------------
R2 <- ((sum(sumC * sumD)) ^ 2) / (sum(sumB) * sum(sumE))
R2 <- round (R2, 4)

# -----------------------------------------------------------------------------
# Percent Bias - PBIAS [%]
# -----------------------------------------------------------------------------
PBIAS <- ((sum(outQOBS - outQSIM)) / sum(outQOBS)) * 100
PBIAS <- round (PBIAS, 4)

# -----------------------------------------------------------------------------
# Absolute Percent Bias APB - [%]
# -----------------------------------------------------------------------------
APB <- ((sum(abs(outQOBS - outQSIM))) / sum(outQOBS)) * 100
APB <- round (APB, 4)

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Filling up output data containers 
# ///////////////////////////////////////////////////////////////////////////////

# df.out01 data.frame is created
df.out01 <- rbind(df.out01, data.frame(outQOBS,
                                       outPREC,
                                       outTEMP,
                                       outPET,
                                       outAET,
                                       outMASSIN,
                                       outRstore,
                                       outSM,
                                       outR,
                                       outSUZ,
                                       outSLZ,
                                       outDELTAMASS,
                                       outQSR,
                                       outQ0,
                                       outQ1,
                                       outQ2,
                                       outQOT,
                                       outQSIM,
                                       outQRES,
                                       outQSIM3))

# cbind function is applied to df.selection and df.out01 data.frames
df.out01 <- cbind(df.selection, df.out01)

# Descriptive statistics are applied on df.out01 data.frame
df.out01.desc <- round((as.data.frame(stat.desc(df.out01[, 8 : 27]))),4)

# A summary data.frame is created based on df.out01.desc data.frame
df.summ <- df.out01.desc[7 , ]

# df.summ data.frame is sorted
df.summ <- df.summ [c("outPREC",
                      "outQOBS",
                      "outMASSIN",
                      "outPET",
                      "outAET",
                      "outQSIM",
                      "outQSR",
                      "outQ0",
                      "outQ1",
                      "outQ2",
                      "outQRES")]

# df.summ data.frame is melted
df.summ <- melt(df.summ)

# df.summ data.frame is normalized
df.summ$norm <- round((df.summ$value) / (df.summ[3 , 2]), 3)

# df.out02 data.frame is created for calibration purposes
df.out02 <- rbind(df.out02, data.frame(outQSIM))

# df.out03 data.frame is created
df.out03 <- data.frame(NSeff, LNNSeff, R2, PBIAS, APB)

# df.out03 data.frame is melted
df.out03 <- melt(df.out03)

# /////////////////////////////////////////////////////////
# BLOCK: ggplot2 time series analysis 
# /////////////////////////////////////////////////////////

# A input precipitation boxplot (by month) is created
gghbvtec01 <- ggplot(aes(y = outPREC,x = MONTH_CH),data = df.out01) +
  geom_boxplot(aes(colour = MONTH_CH),size = 0.75, outlier.size = 4.5) +
  stat_boxplot(aes(colour = MONTH_CH),geom = "errorbar", size = 0.75) +
  scale_fill_brewer(guide = guide_legend(),palette = 'Spectral') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0, min.n = 20.0)) +
  ggtitle("HBV-TEC - Input precipitation boxplot (by month)") +
  xlab("Month") +
  ylab("Precipitation (mm/T)") +
  theme_grey(base_size = 22.0)

# A input precipitation boxplot (by month) is requested
gghbvtec01

# A input precipitation boxplot (by year) is created
gghbvtec02 <- ggplot(aes(y = outPREC,x = YEAR_CH), data = df.out01) +
  geom_boxplot(aes(colour = YEAR_CH), size = 0.75, outlier.size = 4.5) +
  stat_boxplot(aes(colour = YEAR_CH), geom ="errorbar", size = 0.75) +
  scale_fill_brewer(guide = guide_legend(), palette = 'Spectral') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
  ggtitle("HBV-TEC - Input precipitation boxplot (by year)") +
  xlab("Year") +
  ylab("Precipitation (mm/T)") +
  theme_grey(base_size = 22.0)

# A input precipitation boxplot (by year) is requested
gghbvtec02

# A Qsimulated vs Qobserved plot is created
gghbvtec03 <- ggplot(aes(x = outQOBS,y = outQSIM),data = df.out01) +
  geom_point(aes(colour = outQRES,size = abs(outQRES)),alpha = 0.95) +
  geom_abline(data=df.attri,colour = '#666666',size = 0.95,alpha = 0.95) +
  scale_colour_gradient(guide = guide_colourbar(direction = 'vertical'),
                        breaks = scales::pretty_breaks(min.n = 5.0),low = '#ff0000',high = '#00ff00') +
  #scale_size(guide = c(pos = 'none'), range = c(2,10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("HBV-TEC - Qsimulated vs Qobserved") +
  xlab("Qobserved (mm/T)") +
  ylab("Qsimulated (mm/T)") +
  theme_grey(base_size = 22.0)

# A Qsimulated vs Qobserved plot is requested
gghbvtec03

# A residuals histogram is created
gghbvtec04 <- ggplot() +
  geom_histogram(aes(y = ..density..,x = outQRES),
                 data = df.out01, colour = '#000000',fill = '#cc0000') +
  geom_density(aes(x = outQRES,y = ..density..),
               data = df.out01,colour = '#000099',size = 2.0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  ggtitle("HBV-TEC - Residuals Histogram") +
  xlab("Qresidual (mm/T)") +
  ylab("Density") +
  theme_grey(base_size = 22.0)

# A residuals histogram is requested
gghbvtec04

# A Model residuals boxplot (by month) is created
gghbvtec05 <- ggplot(aes(y = outQRES,x = MONTH_CH), data=df.out01) +
  geom_boxplot(aes(colour = MONTH_CH), size = 0.75,outlier.size = 4.50) +
  stat_boxplot(aes(colour = MONTH_CH), geom ="errorbar", size = 0.75) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0 ,min.n = 10.0)) +
  ggtitle("HBV-TEC - Model residuals boxplot (by month)") +
  xlab("Month") +
  ylab("Qresidual (mm/T)") +
  theme_grey(base_size = 22.0)

# A Model residuals boxplot (by month) is requested
gghbvtec05

# A Model residuals boxplot (by year) is created
gghbvtec06 <- ggplot(aes(y = outQRES,x = YEAR_CH), data=df.out01) +
  geom_boxplot(aes(colour = YEAR_CH), size = 0.75,outlier.size = 4.50) +
  stat_boxplot(aes(colour = YEAR_CH), geom ="errorbar", size = 0.75) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0 ,min.n = 10.0)) +
  ggtitle("HBV-TEC - Model residuals boxplot (by year)") +
  xlab("Year") +
  ylab("Qresidual (mm/T)") +
  theme_grey(base_size = 22.0)

# A Model residuals boxplot (by year) is requested
gghbvtec06

# A mass balance plot is created 
gghbvtec07 <- ggplot() +
  stat_identity(aes(x = variable,y = value,fill = variable), data = df.summ,
                geom = 'bar',position = position_identity()) +
  scale_fill_brewer(guide = guide_legend(),palette = 'Spectral') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0, min.n = 10.0)) +
  geom_text(aes(x = variable,y = value,label = value), data = df.summ,
            size = 5.0, vjust = 0.5, parse = FALSE) +
  ggtitle("HBV-TEC - Mass balance summary (by entire period)") +
  geom_hline(data = df.attri,size = 0.75, yintercept = 0.0) +
  xlab("Mass Balance Component") +
  ylab("Mass (mm)") +
  theme_grey(base_size = 22.0)

# A mass balance plot is requested
gghbvtec07

# A mass balance plot is created (normalized)
gghbvtec08 <- ggplot() +
  stat_identity(aes(x = variable,y = norm,fill = variable), data = df.summ,
                geom = 'bar',position = position_identity()) +
  scale_fill_brewer(guide = guide_legend(), palette = 'Spectral') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0, min.n = 10.0)) +
  geom_text(aes(x = variable,y = norm,label = norm), data = df.summ,
            size = 5.0, vjust = 0.5, parse = FALSE) +
  ggtitle("HBV-TEC - Normalized Mass Balance Summary (by entire period)") +
  geom_hline(data = df.attri,size = 0.75, yintercept = 0.0) +
  xlab("Mass Balance Component") +
  ylab("Fraction (unitless)") +
  theme_grey(base_size = 22.0)

# A mass balance plot is requested (normalized)
gghbvtec08

# A objective functions summary plot is created
gghbvtec09 <- ggplot(aes(x = variable,y = value),data=df.out03) +
  geom_point(aes(shape = variable,colour = variable, fill = variable), size = 7.5) +
  facet_wrap(facets = ~variable, scales = 'free_y') +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0, min.n = 10.0)) +
  ggtitle("HBV-TEC - Objective functions summary") +
  xlab("Objective Function") +
  ylab("Value (-/%)") +
  theme_grey(base_size = 22.0)

# A objective functions summary plot is requested
gghbvtec09

# A parameter summary plot is created
gghbvtec10 <- ggplot(aes(x = PARAMETER,y = VALUE), data=df.param) +
  geom_point(aes(colour = PARAMETER, shape = PARAMETER),size = 8.0) +
  facet_wrap(facets = ~ PARAMETER, scales = 'free_y') +
  scale_shape_manual(values = 8: (length(unique(df.param$PARAMETER)) + 8) ) +
  ggtitle("HBV-TEC - Model parameter summary") +
  xlab("Model Parameter") +
  ylab("Value (-)") +
  theme_grey(base_size = 22.0)

# A parameter summary plot is created
gghbvtec10

# -----------------------------------------------------------------------------
# Main Hydrograph
# -----------------------------------------------------------------------------

# Top main hydrograph is created
gghbvtec10.top <- ggplot(aes(x = DATE,y = outPREC),data = df.out01) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0),trans = 'reverse') +
  stat_identity(size = 0.01,geom = 'bar',position = position_identity()) +
  scale_x_date(breaks = scales::pretty_breaks(n = 8.0, min.n = 8.0), expand = c(0.05,0.5)) +
  theme_gray(base_size = 16.0) +     
  ggtitle(label = 'HBV-TEC. Yearly Hydrograph') +
  ylab(label = 'Rainfall (mm/T)') #+
#theme(plot.margin = unit(c(10,10,-((max(df.out01$outPREC))*.55),10),units="points"))

# Top main hydrograph is requested
gghbvtec10.top

# Bottom main hydrograph is created
gghbvtec10.bottom <- ggplot(aes(x = DATE,y = outQSIM),data=df.out01) +
  geom_line(colour = '#333333',size = 1.25) +
  geom_point(aes(x = DATE,y = outQOBS, colour = abs(outQRES),size = abs(outQRES))) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10.0,min.n = 10.0)) +
  scale_x_date(breaks = scales::pretty_breaks(n = 8.0,min.n = 8.0),
               expand = c(0.05,0.5)) +
  scale_colour_gradient(guide = guide_colourbar(),low = '#0000ff',high = '#ff0000') +
  scale_size(guide = guide_legend(),range = c(2,7)) +
  theme_grey(base_size = 16.0) +
  ylab(label = 'Flow (mm/d)') +
  xlab(label = 'Date') +
  theme(legend.position="bottom") #+
#theme(plot.margin = unit(c(10,10,10,10),units="points")) 

# Bottom main hydrograph is requested
# gghbvtec10.bottom

# Main hydrograph is generated
grid.arrange(gghbvtec10.top, gghbvtec10.bottom, heights = c(0.25, 0.75)) 

# ///////////////////////////////////////////////////////////////////////////////
# BLOCK: Export and display of ouput data containers
# ///////////////////////////////////////////////////////////////////////////////

# ----------------------------------------------------------
# Summary LogFile:
# ----------------------------------------------------------
LogFile <- file("hbvtec_logfile.txt", "w")
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("HBV-TEC-96. V.1.0. Summary of Activities."), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("Instituto Tecnologico de Costa Rica. (www.tec.ac.cr)"), LogFile)
writeLines(c("Maikel Mendez-M (mamendez@itcr.ac.cr);(maikel.mendez@gmail.com)"), LogFile)
writeLines(c("Luis Alexander Calvo (lcalvo@itcr.ac.cr);(lualcava.sa@gmail.com)"), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("Cacthment Attributes"), LogFile)
writeLines(c("CathmentArea [km2] =",warea), LogFile)
writeLines(c("fract_SMi [unitless] =",fract_SMi), LogFile)
writeLines(c("fract_SLZi [unitless] =",fract_SLZi), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("Mass-Balance [(mm/TimeStep)/TimePeriod]"), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("sumQOBS =",sumQOBS), LogFile)
writeLines(c("sumQSIM =",sumQSIM), LogFile)
writeLines(c("deltaQ = (if negative, QSIM overestimates QOBS)",deltaQ), LogFile)
writeLines(c("sumPREC =",sumPREC), LogFile)
writeLines(c("sumPET =",sumPET), LogFile)
writeLines(c("sumAET =",sumAET), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("Relative contribution fractions [unitless]"), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("fractSUZ =",fractSUZ), LogFile)
writeLines(c("fractSLZ =",fractSLZ), LogFile)
writeLines(c("fractUZL = (UZL threshold ONLY)",fractUZL), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("MAXBAS relative weights/TimeStep for HBV [unitless]"), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("MAXBASabswei =",MAXBASabswei), LogFile)
writeLines(c("MAXBASabsweisum = (must be equal to 1)",MAXBASabsweisum), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("Goodness of fit (objective functions)/TimePeriod for HBV"), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("Nash and Sutcliff efficiency criterion; NSeff [unitless] =",NSeff), LogFile)
writeLines(c("Nash and Sutcliffe efficiency with logarithmic values; LNNSeff [unitless]=",LNNSeff), LogFile)
writeLines(c("Percent Bias; expressed as percentage; PBIAS [%] =",PBIAS), LogFile)
writeLines(c("Absolute Percent Bias; expressed as percentage; APB [%] =",APB), LogFile)
writeLines(c("Coefficient of determination squared; R2 [unitless]=",R2), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
writeLines(c("End of program"), LogFile)
writeLines(c("-----------------------------------------------------------------"), LogFile)
close(LogFile)

# df.ptq.desc is displayed
#View(df.ptq.desc)

# df.out01 is displayed
#View(df.out01)

# df.out01.desc is displayed
#View(df.out01.desc)

# df.out02 is displayed
#View(df.out02)

# df.out03 is displayed
#View(df.out03)

# df.ptq.desc data.frame is exported
write.csv(df.ptq.desc, file = "hbvtecptq_desc_hbv.csv")

# df.out01 data.frame is exported
write.csv(df.out01, file = "hbvtecout_hbv.csv")

# df.out01.desc data.frame is exported
write.csv(df.out01.desc, file = "hbvtecout_desc_hbv.csv")

# df.out02 data.frame is exported
write.csv(df.out02, file = "hbvqtecsim.csv")

# df.out03 data.frame is exported
write.csv(df.out03, file = "hbvqteceff.csv")

# ----------------------------------------------------------
# END OF PROGRAM
# ----------------------------------------------------------

#////////////////////////////////////////////////////////////////////////////
