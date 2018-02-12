
#########################################################################################

# Program in R - Recommended times for SUNCAL calibrations 
# Project
# P. Altube Vazquez - Feb 2018

## USER SETTINGS ########################################################################

work_path <- "/home/pav/Desktop/PEDESTAL/SC_procedures"
year <- "2018"

radars <- c("CDV", "LMI", "PBE", "PDA")

# Nominal measurement angles [deg]
el0 <- 25
az0 <- 120
tolerance <- 5

# Precision for the time
precision <- 1 # minutes

## SCRIPT SETTINGS ######################################################################

coords <- list("CDV"=c(41.60192, 1.40283), "LMI"=c(41.09175, 0.86348), 
               "PBE"=c(41.37335, 1.88197), "PDA"=c(41.88882, 2.99717))

heights <- list("CDV"=825, "LMI"=910, "PBE"=631, "PDA"=542) #[m]

date_i <- as.Date(paste(year, "-01-21", sep=""))
date_f <- as.Date(paste(year, "-12-20", sep=""))

date_seq <- seq(date_i, date_f, by="day")

out_path <- paste(work_path, year, sep="/")
if (!file.exists(out_path)){dir.create(out_path, recursive=TRUE)}

## FUNCTIONS ############################################################################

# Gives decimal time in format "HH:MM:SS:mS"
time_format <- function(hour_dec, precision="sec"){
  
  t_dec <- hour_dec

  prec <- list("hour"=1, "min"=2, "sec"=3, "msec"=4)
  loop <- seq(1:prec[[precision]])
  t_str <- ""
  sep <- ""
  
  for (i in loop){
    
    tt <- floor(t_dec)
    
    XX <- as.character(tt)
    XX[tt<10] <- paste("0", XX[tt<10], sep="")
    
    t_str <- paste(t_str, XX, sep=sep)
    dg <- 4-i
    sep <-":"
    
    if (i==length(loop)-1){
      dg <- 0
    }
    t_dec <- round((t_dec-tt)*60, digits=dg)
    
  }
  
  return(t_str)
  
}

# Calculates reference Julian Day
prime_JD <- function(date, utc){
  
  # Prime Julian Day function: days since 1st Jan 2000, 12UTC
  # (JD is days since 1st Jan 1949, 00UTC)
  year <- as.integer(format(date, format="%Y"))
  dof <- as.integer(format(date, format="%j"))
  # n is the number of days since 12UT of 1st January, 2000
  n <- dof + (year-2000)*365 + (year - 2000)%/%4 + utc/24 - 0.5
  # -0.5 is the offset in days between 0UT and 12UT
  return(n)
}

# Solar elevation and azimuth (local), atm. refraction not considered:
solar_position <- function(lat, lon, date, utc){
  
  rad <- pi/180
  
  # Prime Julian Day with reference at noon 1st Jan 2000 
  n <- prime_JD(date, utc)
  
  # Ecliptic coordinates (in degrees)
  L <- (280.460 + 0.9856474*n)%%360 # mean longitude
  g <- (357.528 + 0.9856003*n)%%360 # mean anomaly
  ll <- (L + 1.915*sin(g*rad) + 0.02*sin(2*g*rad))%%360 # ecliptic longitude
  ep <- 23.439 - 0.0000004*n # obliquity of ecliptic
  
  # Right Ascension (in degrees!)
  ra <- atan2(cos(ep*rad)*sin(ll*rad), cos(ll*rad))*(1/rad)
  # Declination (in degrees)
  dec <- asin(sin(ep*rad)*sin(ll*rad))*(1/rad)
  
  # Greenwich Mean Sidereal Time (in degrees: 15deg per hour) 
  GMSTdeg <- ((6.697375 + 0.0657098242*n + utc)%%24)*15
  # Local Mean Sidereal Time (in degrees) 
  LMSTdeg <- GMSTdeg + lon
  
  # Hour angle (in degrees) 
  ha <- LMSTdeg - ra
  
  # Solar elevation
  sinEL <- sin(dec*rad)*sin(lat*rad)  + cos(dec*rad)*cos(lat*rad)*cos(ha*rad)
  EL <- asin(sinEL)*(1/rad)
  
  # Solar azimuth
  sinAZ <- -cos(dec*rad)*sin(ha*rad)/cos(EL*rad)
  cosAZ <- (sin(dec*rad) - (sin(EL*rad)*sin(lat*rad)))/(cos(EL*rad)*cos(lat*rad))
  AZ <- (atan2(sinAZ, cosAZ)*(1/rad))%%360

  return(data.frame("el"=EL, "az"=AZ))
  
}

# Atmospheric refraction model (Holleman & Huuskonen, 2013):
ref_kmodel <- function(el, k=1.25){
  
  N0 <- 313
  H0 <- 8.4 #[km]
  Re <- 6378 #[km]
  k <- 1.25
  
  a <- N0*1e-6
  b <- ((k-1)/(2*k-1))*cos(el*pi/180)
  c <- (sin(el*pi/180)^2) + 2*((2*k-1)/(k-1))*a + 2*((2*k-1)/(k*Re))*H0
  ref_ang <- b*(sqrt(c) - sin(el*pi/180))*180/pi # refraction angle in degrees
  return(ref_ang)
}

# Gives the time interval that minimizes a difference within an specified tolerance
find_meas_time <- function(diff_df, tolerance){
  
  # diff_df data-frame needs at least two columns: time & diff
  
    sel <- diff_df[diff_df$diff<=tolerance,]
    
    m <- which(sel$diff==min(sel$diff))
    
    tt <- time_format(sort(unique(c(min(sel$time), sel$t[m], max(sel$time)))), precision="min")
    return(paste(tt, collapse = "-"))

}

# Main function
SC_meas_times <- function(date, az0, el0, lat=0, lon=0, masl=0, tolerance=5, precision=1){

  precision <- precision/60
  # Time sequence to calculate positions
  times <- seq(0, 24-precision, precision)
  # Target angles
  angles <- list("az"=az0, "el"=el0)
  
  # Solar position
  sun_pos <- solar_position(lat, lon, date, times)
  # Include refraction effect in elevation
  sun_pos$el <- sun_pos$el + ref_kmodel(sun_pos$el)
  sun_pos$time <- times
  
  # Separate East/West
  data <- list("E"=sun_pos[sun_pos$az<180,], "W"=sun_pos[sun_pos$az>=180,])
  # Translation of az values
  data[["W"]]$az <- data[["W"]]$az - 180
  
  tt <- list()
  tt["date"] <- format(date, format="%Y-%m-%d")

  for (d in names(data)){
    
    for (a in names(angles)){
    
      n <- paste(a, d, sep="_")
      data[[d]]$diff <- abs(data[[d]][[a]] - angles[[a]])
      tt[n] <- find_meas_time(data[[d]], tolerance=tolerance)
    }
    
  }
  
  return(tt)
  
}

#########################################################################################


# Common Headers for output files:
# Headers of output files
H1 <- paste("# Recommended SUNCAL measurement times -", year, sep=" ")
H4 <- paste("# Target azimuth:", az0, "deg", sep=" ")
H5 <- paste("# Target elevation:", el0, "deg", sep=" ")
H6 <- paste("# Tolerance:", tolerance, "deg", sep=" ")
H7 <- paste("# Time precision:", precision, "min", sep=" ")
H8 <- rep("#", 8)

for (r in radars){
  
  meas_times <- list()
  for (i in seq_along(date_seq)){
    
    date <- date_seq[i]
    tmp <- SC_meas_times(date, az0, el0, lat=coords[[r]][1], lon=coords[[r]][2], 
                         masl=heights[[r]], tolerance=tolerance, precision=precision)
    meas_times[[i]] <- as.data.frame(tmp)

  }
  
  SC_times <- do.call(rbind.data.frame, meas_times)
  
  H2 <- paste("# File created:", Sys.time(), sep=" ")
  H3 <- paste("# Radar abbrv.:", r, sep=" ")
  hdr <- as.character(unlist(mget(paste("H", seq(1:8), sep=""))))
  
  file_out <- paste(out_path, "/SC_times_", year, "_", r, ".txt", sep="")
  write(hdr, file=file_out, ncolumns=length(hdr), sep="\n", append=FALSE)
  write.table(SC_times, file=file_out, quote=FALSE, row.names = FALSE, col.names = TRUE, append=TRUE) 
  
}
    




