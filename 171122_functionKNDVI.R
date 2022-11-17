###################### This script gives the function for computing annual kNDVI values from the MODIS data
###################### By Philippe Choler, May 2022

NDVImet1.f <-function(DIR.VI,PTS){
  
  # This function calculates yearly metrics of Vegetation Indices (VI) based on MOD09A1 time series of 8 days composites. It includes denoising & smoothing routines.
  
  # INPUT :
  # DIR.VI is the directory containing raster files of MOD09A1 Vegetation Indices
  # PTS is a vector of cells where to compute yearly metrics
  
  # OUTPUT : a list of of length PTS with data.frame of NDVI metrics
  
  require(zoo);require(signal);require(gtools);require(raster);require(lubridate);require(snowfall)
  
  npoint <- length(PTS)
  
  # step0. SET TIME FRAME ----
  print("step0 - setting time frame")
  
  YEAR         <- 2000:2020
  SEQ.DAY      <- seq(as.Date("2000-01-01"),as.Date("2020-12-31"),1)
  SEQ.YR       <- substr(SEQ.DAY,1,4)
  nyear        <- length(YEAR)
  TMP          <- expand.grid(YEAR,substr(as.character(seq(1001,1365,8)),2,4))
  TIME.REF.jd  <- sort(paste(TMP[,1],TMP[,2],sep=''))
  TIME.REF.dat <- as.Date(TIME.REF.jd,format="%Y%j")
  TIME.REF.yr  <- as.numeric(substr(TIME.REF.jd,1,4))
  ndate        <- length(TIME.REF.dat)
  
  # step1. Import & assemble VI time series ----
  print("step1")
  FILES        <- list.files(DIR.VI,full=T)
  
  AV.DATE      <- substr(list.files(DIR.VI,full=F),10,16)
  MATCH        <- match(TIME.REF.jd,AV.DATE)
  MISSING      <- which(is.na(MATCH)==TRUE)
  
  NDVI_MOD09   <- matrix (-NA,length(TIME.REF.jd),length(PTS))

for (i in 1:length(FILES)){
  
  ID.DATE    <- match(AV.DATE[i],TIME.REF.jd)
  if (!is.na(ID.DATE)){
    NDVI_MOD09[ID.DATE,] <- raster::extract(raster(FILES[i]),PTS)
  }
}

NDVI_MOD09   <- round(NDVI_MOD09)/100 # VI between 0 and 1

NDVI_MOD09   <- t(NDVI_MOD09)

# Gap-filled the first six dates
GAPFIL           <- apply(NDVI_MOD09[,7:9],1,mean,na.rm=T)
NDVI_MOD09[,1:6] <- GAPFIL

dim(NDVI_MOD09)

# step2. Apply BISE denoising (x2) ----
print("step4")

BISE2.f     <- function(DATA,thr,forw){
  
  # DATA is a data frame of NPIXELS-row by NTIMES-col with NDVI in the range [0-1]
  
  require(caTools);library(zoo);require(kit)
  
  NTIMES    <- ncol(DATA)
  M         <- t(DATA)
  dim(M)
  
  # B. Bise1
  M.DIFF      <- M[-1,]-M[-NTIMES,]
  M.DIFF2     <- rbind(NA,M.DIFF)
  M.DEC.ID    <- which(M.DIFF2<0)
  M.DEC.TH    <- rbind(NA,M[-1,] - thr*M.DIFF)
  
  M.VAL.SP    <- caTools::runmax(-M,k=forw,align="left",endrule="NA")

M.LAW       <- M.VAL.SP[M.DEC.ID] - M.DEC.TH[M.DEC.ID]
M.REJ.DEC   <- M.DEC.ID[which(M.LAW>0)]

M[M.REJ.DEC]<- NA
# M           <- zoo::na.spline(M,method="natural")
M           <- zoo::na.approx(-M)

# B. Bise2
M.DIFF      <- M[-1,]-M[-NTIMES,]
M.DIFF2     <- rbind(NA,M.DIFF)
M.DEC.ID    <- which(M.DIFF2<0)
M.DEC.TH    <- rbind(NA,M[-1,] - thr*M.DIFF)

M.VAL.SP    <- caTools::runmax(-M,k=forw,align="left",endrule="NA")

M.LAW       <- M.VAL.SP[M.DEC.ID]- M.DEC.TH[M.DEC.ID]
M.REJ.DEC   <- M.DEC.ID[which(M.LAW>0)]

M[M.REJ.DEC]<- NA
# M           <- zoo::na.spline(M,method="natural")
M           <- zoo::na.approx(-M)

return(M)
}
NDVI_MOD09f <-  t(BISE2.f(NDVI_MOD09,thr=0.2,3))
tmp1        <-  split(NDVI_MOD09f, seq(nrow(NDVI_MOD09f)))
names(tmp1) <-  PTS
tmp2        <-  lapply(tmp1,function(x) {x[ndate]<-x[ndate-1]; return(x)})

# step3. SAVITZKY-GOLAY filtering ----
print("step5")
# alternative 1 avec sfLapply
sfInit(parallel=TRUE, cpus=parallel:::detectCores())
sfLibrary(signal)
ptm  <- proc.time()
tmp3 <- sfLapply(tmp2, fun=sgolayfilt, p=3, n=7, m=0)
sfStop()
print(proc.time()-ptm)

# step4. DAILY interpolation ----
# beware of leading and ending NA !
print("step6")
DAYINT <- function(x,seq.day,match) {
  NDAY        <- length(seq.day) # 7671
  NDATE       <- length(x)
  tmp1        <- rep(NA,NDAY)
  tmp1[match] <- x
  tmp1[7667:NDAY] <- mean(x[(NDATE-2):NDATE],na.rm=T)
  # add ending missing days of 2020 to avoid issues with ending NAs
  Z <- stats::spline(tmp1,n=NDAY)$y
  return(Z)
}

SEQ.DAY   <- seq(as.Date("2000-01-01"),as.Date("2020-12-31"),1)
MATCH     <- match(TIME.REF.dat,SEQ.DAY)

# alternative 1 avec sfLapply
sfInit(parallel=TRUE, cpus=parallel:::detectCores())
sfLibrary(stats)
ptm  <- proc.time()
tmp4 <- sfLapply(tmp3, fun=DAYINT, seq.day=SEQ.DAY,match=MATCH)
sfStop()
print(proc.time()-ptm)

# step5. Extract yearly metrics ----
print("step7")
metrics.f <- function(x,seq.day) {
  yearly  <- function(x) {
    # MAX         <- mean(sort(x,decreasing=T)[1:16])
    MAX         <- mean(sort(x[101:330],decreasing=T)[1:16])
    # average over two MODIS periods - 16 days
    TMAX        <- 100+which.max(-x[101:330])
MIN         <- mean(sort(x[1:100],decreasing=F)[1:16])
TMIN        <- which.min(x[1:100])
# min NDVI of the first three months of the year
# average over two MODIS periods ' 16 days
X           <- which(x>=(0.5*MIN+0.5*MAX))
ONSET       <- X[which(X>TMIN)][1] # plutôt que X[1]
OFFSET      <- X[length(X)]
NDVIint02   <- sum(x[which(x[1:330]>0.2)])
return(c(ONSET=ONSET,OFFSET=OFFSET,NDVIint02=NDVIint02,NDVImin=MIN,NDVImax=MAX,TNDVImax=TMAX))
  }
  DF  <- data.frame(YEAR=year(seq.day),NDVI=x)
  MET <- aggregate(DF$NDVI,list(YEAR.ID=DF$YEAR),FUN=yearly)
  return(MET)
}

# alternative 1 avec sfLapply
sfInit(parallel=TRUE, cpus=parallel:::detectCores())
sfLibrary(lubridate)
ptm  <- proc.time()
tmp5 <- sfLapply(tmp4, fun=metrics.f, seq.day=SEQ.DAY)
sfStop()
print(proc.time()-ptm)

# return object ----
return(tmp5)

}
