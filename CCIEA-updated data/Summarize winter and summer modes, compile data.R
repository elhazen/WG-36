rm(list=ls())

library(reshape2)
library(dplyr)
library(tidyr)
library(lubridate) #https://rstudio-pubs-static.s3.amazonaws.com/28038_1bcb9aa80ca84f27ace07d612872861a.html
library(ggplot2)
library(gridExtra)
library(ggrepel)

base.dir <-"~/Dropbox/CCIEA/2018/Council report"
setwd(base.dir)

### NPGO
npgo <- read.csv("cciea_OC_NPGO_97bf_314a_7e4e.csv", header=TRUE)[-1,]
head(npgo)
class(npgo$time)

npgo$date <- as_date(ymd_hms(npgo$time))
npgo$year <- year(npgo$date)
npgo$month <- month(npgo$date)

npgo_out <- npgo %>%
  filter(year > 1995 & year < 2017) %>%
  group_by(year)  %>%
  summarise(npgo.winter = mean(NPGO[month==1 | month==2 | month==3]))
#View(npgo_out)
npgo_out2 <- melt(npgo_out, id.vars="year", variable.name = "timeseries")
names(npgo_out2)[1] <- "Year"

### NOI
noi <- read.csv("cciea_OC_NOI_40fa_4687_c4da.csv", header=TRUE)[-1,]
head(noi)
class(noi$time)

noi$date <- as_date(ymd_hms(noi$time))
noi$year <- year(noi$date)
noi$month <- month(noi$date)
noi$NOI <- as.numeric(as.character(noi$NOI))

noi_out <- noi %>%
  filter(year > 1995 & year < 2017) %>%
  group_by(year)  %>%
  summarise(noi.summer = mean(NOI[month==6 | month==7 | month==8]))
#View(noi_out)
noi_out2 <- melt(noi_out, id.vars="year", variable.name = "timeseries")
names(noi_out2)[1] <- "Year"

### COPEPODS
summer.copepods <- read.csv("summer copepods.csv", header=TRUE)
summer.copepods2 <- melt(summer.copepods, id.vars="Year", variable.name = "timeseries")

### SEA LION PUP COUNTS
pup.counts <- read.csv("pup counts.csv", header=TRUE)
pup.counts2 <- melt(pup.counts, id.vars="Year", variable.name = "timeseries")

### COMPILE INTO SINGLE DATA FRAME
# columns are year, value, time series

out.df <- rbind(
  npgo_out2,
  noi_out2,
  summer.copepods2,
  pup.counts2
)

write.csv(out.df, "Compiled data for 2018 report.csv", row.names=FALSE)


### LOOK AT THE DATA
plot.df <- dcast(out.df, Year~timeseries)

### COMPARE CURRENT COMPILED DATA TO DATA USED IN 2015 WORKSHOP

df.2015 <- read.csv("compiled data used in 2015 workshop.csv", header = TRUE)
plot.df.2015 <- dcast(df.2015, year~timeseries)

plot(plot.df.2015$`CA sea lion pup production`,plot.df$live_pup_count[-c(20:21)])
plot(plot.df.2015$`Copepod anomaly summer`,plot.df$AvgOfNorthernBiomassAnomaly[-c(20:21)])
plot(plot.df.2015$`NOI summer`,plot.df$noi.summer[-c(20:21)])
plot(plot.df.2015$`NPGO winter`,plot.df$npgo.winter[-c(20:21)])

### MAKE SOME PLOTS

p.2014.copepods <- ggplot(data=plot.df[-c(20:21),], aes(x=npgo.winter, y=AvgOfNorthernBiomassAnomaly))+
  geom_point()+
  geom_text_repel(aes(label=Year))+
  geom_smooth(method="loess")+
  ggtitle("1996-2014: nonlinear model is best")+
  xlab("NPGO Winter\n(Jan-Mar)")+
  ylab("Summer Northern Copepod Anomaly\n(May-September)")+
  theme_bw()
p.2014.copepods

p.2016.copepods <- ggplot(data=plot.df, aes(x=npgo.winter, y=AvgOfNorthernBiomassAnomaly))+
  geom_point()+
  geom_text_repel(aes(label=Year))+
  geom_smooth(method="lm")+
  ggtitle("1996-2016: linear model is best")+
  xlab("NPGO Winter\n(Jan-Mar)")+
  ylab("Summer Northern Copepod Anomaly\n(May-September)")+
  theme_bw()
p.2016.copepods


grid.arrange(p.2014.copepods, p.2016.copepods, ncol=2)

quartz(file = "Copepods plots.pdf",type="pdf",dpi=300, height=9,width=12)
print(grid.arrange(p.2014.copepods, p.2016.copepods, ncol=2))
dev.off()
  

p.2014.pups <- ggplot(data=plot.df[-c(20:21),], aes(x=noi.summer, y=live_pup_count))+
  geom_point()+
  geom_text_repel(aes(label=Year))+
  geom_smooth()+
  ggtitle("1997-2014: nonlinear model is best")+
  xlab("NOI Summer\n(Jun-Aug)")+
  ylab("California Sea Lion Pup Counts")+
  theme_bw()
p.2014.pups

p.2016.pups <- ggplot(data=plot.df, aes(x=noi.summer, y=live_pup_count))+
  geom_point()+
  geom_text_repel(aes(label=Year))+
  geom_smooth()+
  ggtitle("1997-2016: nonlinear model is best")+
  xlab("NOI Summer\n(Jun-Aug)")+
  ylab("California Sea Lion Pup Counts")+
  theme_bw()
p.2016.pups

grid.arrange(p.2014.pups, p.2016.pups, ncol=2)

quartz(file = "Sea lion pups plots.pdf",type="pdf",dpi=300, height=9,width=12)
print(grid.arrange(p.2014.pups, p.2016.pups, ncol=2))
dev.off()


### DO in ml/L
dO2a <- read.csv("cciea_OC_DO_8084_32b4_96ff.csv", header=TRUE)[-1,]
head(dO2a)
class(dO2a$time)

dO2a$date <- as_date(ymd_hms(dO2a$time))
dO2a$year <- year(dO2a$date)
dO2a$month <- month(dO2a$date)
dO2a$dissolved_oxygen <- as.numeric(as.character(dO2a$dissolved_oxygen))

# dO2a_out <- dO2a %>%
#   filter(year > 1995 & year < 2017) %>%
#   group_by(year)  %>%
#   summarise(npgo.winter = mean(NPGO[month==1 | month==2 | month==3]))
#View(npgo_out)
#dO2a_out2 <- melt(dO2a[,c('year','dissolved_oxygen')], id.vars="year", variable.name = "timeseries")
#names(dO2a_out2)[1] <- "Year"
ggplot(data=dO2a, aes(x=year, y=dissolved_oxygen))+
  geom_point()+
  geom_line()+
  theme_bw()
