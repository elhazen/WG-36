library(MARSS)
library(xtable)
library(reshape2)
tabledir="~/Dropbox/Documents/R/IEA/IEA-thresholds/data/output"
setwd(tabledir)

# prepare the data
ccALL.old<-read.csv("~/Dropbox/Documents/R/IEA/IEA-thresholds/data/Indicator data for workshop 2015- coastwide.csv")
ccALL<-read.csv("~/Dropbox/Documents/R/IEA/IEA-thresholds/data/coastwide data for reference points.csv")
#UNIX DIR
#ccALL<-read.csv("/home/ehazen/R/IEA/data/Indicator data for workshop 2015- coastwide.csv")
#ccALL<-read.csv("/home/ehazen/R/IEA/data/coastwide data for reference points.csv")

ccALL$year <- as.numeric(ccALL$year)
ccALL$timeseries <- gsub("\\(", "", ccALL$timeseries)
ccALL$timeseries <- gsub("\\)", "", ccALL$timeseries)
ccALL$timeseries <- gsub(" ", "_", ccALL$timeseries)
ccALL$timeseries <- gsub("_-_", "_", ccALL$timeseries)
#
# Wide df with columns as variables
dat.full <- dcast(ccALL, year ~ timeseries, value.var = "value")
#dat.types <- dcast(ccALL, Ecological.Oceanographic.or.Human.activities ~ timeseries, value.var = "value")
#dat.red1 <- dat.full[,c("year","Coastal_engineering","Nutrients","Inorganic_pollution","Commercial_landings_of_all_fisheries","NOI_summer","NOI_winter","PDO_summer","PDO_winter","NPGO_summer","NPGO_winter","GF_spp_richness_coastwide","GF-Simp_coastwide","Scav_ratio_coastwide","GF_MTL_coastwide","Cop_Spp_Rich_Anom_summer","Cop_Spp_Rich_Anom-_winter","California_sea_lion_pup_production")]
dat.red1 <- dat.full
for (i in 1:dim(dat.red1)[2]) dat.red1[,i]<-as.numeric(dat.red1[,i])

save(dat.red1,file="CC.CW.redv2.0.RData")

# load the data
load("CC.CW.redv2.0.RData")

indat=dat.red1
#subset data to only drivers and pressures
#indat=indat[,1:11]

envind<-c(12:15,17,18)
humind<-c(3,4,5,9,16)
ecosysind<-c(2,6,7,8,10,11,19)
cantcontrolind<-c(3,13,17,18)

envind<-c(grep("NOI",names(dat.red1)),grep("NPGO",names(dat.red1)),grep("PDO",names(dat.red1)))
ecosysind<-c(grep("lion",names(dat.red1)),grep("GF",names(dat.red1)),grep("Scav",names(dat.red1)),grep("Cop",names(dat.red1)))
humind<-c(grep("Commercial",names(dat.red1)),grep("Coastal",names(dat.red1)),grep("pollution",names(dat.red1)),grep("Nutrient",names(dat.red1)),grep("landings",names(dat.red1)),grep("Dredging",names(dat.red1)),grep("Habitat",names(dat.red1)))

#A - Human 
indat<-indat[,c(1,humind)]; datagroup="A_Human"
#B - Env
#indat<-indat[,c(1,envind)]; datagroup="B_Env"
#C - Human and Env
#indat<-indat[,c(1,envind,humind)]; datagroup="C_Human_Env"
#D - Eco States
#indat<-indat[,c(1,ecosysind)]; datagroup="D_Ecosystem_State"
#E - Can't Control
#indat<-indat[,c(1,cantcontrolind)]; datagroup="E_forcings"


names<-names(indat)[-1]
##############################################################################################
####USING THIS WILL NORMALIZE DATA ACROSS THE PERIOD and SPATIAL SCALE OF INTEREST (E.G. 1994-2008 OR 1985-2011)
##############################################################################################

#transpose data so time goes across columns
start=1985
end=2014
yrdiff=end-start

yrs<-indat$year>=start & indat$year<=end

dat = t(indat[yrs,])
dat = dat[-1,]
#dat=dat[-1,] #removes teh year row
# get number of time series
N.ts = dim(dat)[1]
# get length of time series
TT = dim(dat)[2] 

Sigma = sqrt(apply(dat, 1, var, na.rm=TRUE))
y.bar = apply(dat, 1, mean, na.rm=TRUE)
dat.z = (dat - y.bar) * (1/Sigma)
rownames(dat.z) = rownames(dat)
#################END OF NORMALIZING ROUTINES################


#######Plotting time series#################################


spp = rownames(dat)
par(mfrow=c(5,5), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in spp){
    plot(dat.z[i,],xlab="",ylab="Index", bty="L", xaxt="n", pch=16, col="blue", type="b")
    axis(1,(0:dim(dat)[2])+1,start+0:dim(dat)[2])
    title(i)
}
##############################################################
#####DFA models with no covariates##############
#############################################
  n<-dim(dat.z)[1]
  Rnew<-matrix(0.01,n,n);diag(Rnew)<-paste("r",1:n,sep="")

cntl.list = list(minit=200, maxit=100000, allow.degen=FALSE, safe=TRUE, trace=1)
 # set up forms of R matrices
 levels.R = list("identity",
              "diagonal and equal",
              "equalvarcov",
              "diagonal and unequal",
              "unconstrained",Rnew)
 
 
  model.data=data.frame()
  for (j in 1:5){
  for (i in 1:5){
  

	model.list = list(m=i, R=levels.R[[j]])

  kemz = MARSS(dat.z, model=model.list, control=cntl.list,
  		z.score=TRUE, form="dfa")
	
  model.data = rbind(model.data, 
	 data.frame(R=levels.R[[j]],
           m=i,
           covariates="none",
           logLik=kemz$logLik,
           K=kemz$num.params,
           AICc=kemz$AICc,
           Convergence=kemz$convergence,
           NumIter=kemz$numIter,
           stringsAsFactors=FALSE))
    assign(paste("kemz", i, levels.R[[j]], "none", sep="."), kemz)
    write.csv(model.data, file=paste("model.resultsv2_",datagroup,".csv",sep=''))
    save.image(paste("modelOutputv2",datagroup,".Rdata",sep='')) 
}
}

###################################################
### Kir's plotting code
###################################################
year<-indat$year[yrs]

get_loadings<-function(tmp_model=H1.3){
    # get loadings
    # get the inverse of the rotation matrix
    if(dim(coef(tmp_model, type="matrix")$Z)[2]==1){
      H.inv =1
    }else{
      H.inv = varimax(coef(tmp_model, type="matrix")$Z)$rotmat
    }
    # Finally, we use H􀀀1 to rotate the factor loadings and H to rotate the trends
    # as in Equation 9.10.
    # rotate factor loadings
    Z.rot = coef(tmp_model, type="matrix")$Z %*% H.inv
    # rotate trends
    trends.rot = solve(H.inv) %*% tmp_model$states
    return(list(H.inv=H.inv,Z.rot=Z.rot,trends.rot=trends.rot))
}   

col1<-colorRampPalette(colors()[c(71,73)])
col1<-colorRampPalette(c("blue","light blue","light green","dark green"))
  # plot results function

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

plot_MARSS_load<-function(tmp_MARSS=kemzMod,lab=""){
    dev.new(height=6,width=9)
    par(mfrow=c(1,2))
    par(mar=c(10,2,0,0)) # margins of graph: (bottom,left, top, right)
    par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
    par(oma=c(1,1,1,1))# outer margins of graph: (bottom,left, top, right)
    trends<-get_loadings(tmp_MARSS)$trends.rot
    legendText=paste("Trend",1:dim(trends)[1])
    plot(year,trends[1,],type="l",ylim=c(min(trends),max(trends)))
    for(i in 1:dim(trends)[1]){
      lines(year,trends[i,],col=col1(dim(trends)[1])[i],lwd=2)
    }
    legend("topright",legendText,col=col1(dim(trends)[1]),box.lty=0,lwd=2)
    mtext(side=3,adj=0.02,paste("conv =",kemzMod$convergence))
    mtext(side=1,line=5,lab)
    abline(h=0)
    tmpZ<-get_loadings(tmp_MARSS)$Z.rot
    barplot(t(tmpZ),beside=T,names.arg=rownames(kemzMod$model$data),las=2,col=col1(dim(trends)[1]))
    abline(h=0)
}
plot_MARSS_load2<-function(tmp_MARSS=kemzMod,lab="",cutoff=0.2){
    dev.new(height=6,width=9)
    par(mfrow=c(1,2))
    par(mar=c(10,2,0,0)) # margins of graph: (bottom,left, top, right)
    par(mgp=c(2,.5,0)) # axis margins: (distance of lab, distance of numbers, axis line)
    par(oma=c(1,1,1,1))# outer margins of graph: (bottom,left, top, right)
    trends<-get_loadings(tmp_MARSS)$trends.rot
    legendText=paste("Trend",1:dim(trends)[1])
    plot(year,trends[1,],type="l",ylim=c(min(trends),max(trends)))
    for(i in 1:dim(trends)[1]){
      lines(year,trends[i,],col=col1(dim(trends)[1])[i],lwd=2)
    }
    legend("topright",legendText,col=col1(dim(trends)[1]),box.lty=0,lwd=2)
    mtext(side=3,adj=0.02,paste("conv =",kemzMod$convergence))
    mtext(side=1,line=5,lab)
    abline(h=0)
    tmpZ<-t(get_loadings(tmp_MARSS)$Z.rot)
    tmpZlow<-tmpZ
    tmpZlow[which(abs(tmpZ)>cutoff)]<-0
    tmpZ[which(abs(tmpZ)<cutoff)]<-0
    barplot(tmpZ,beside=T,names.arg=rownames(kemzMod$model$data),las=2,ylim=c(min(tmpZ)-.05,max(tmpZ)+.05),col=col1(dim(trends)[1]))
    par(new=T); barplot(tmpZlow,beside=T,names.arg=rownames(kemzMod$model$data),las=2,ylim=c(min(tmpZ)-.05,max(tmpZ)+.05),angle=45,density=20,border=TRUE,col=add.alpha(col1(dim(trends)[1]),0.3),axes=FALSE)
    abline(h=0)
}

###################################################
### makemodeltable
###################################################
# calculate delta-AICc
model.data$delta.AICc = model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt = exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt = wt/sum(wt)
# sort results
model.tbl = model.data[order(model.data$AICc),-4]
# drop AICc from table
# calculate cumulative wts
model.tbl$Ak.wt.cum = cumsum(model.tbl$Ak.wt)
model.tbl = model.tbl[,-4]
tmpaln = "c" #figure out the number of cols automatically
for(i in 1:ncol(model.tbl)) {tmpaln = paste(tmpaln,"c",sep="")}
thetable = xtable(model.tbl,caption='Model selection results.', label='tab:tablefits', align=tmpaln, digits=c(1,1,1,1,1,1,1,1,2,2))
align(thetable) = "ccccp{3.5cm}p{0.7cm}p{1.5cm}p{1.75cm}cc"
print(thetable, type = "latex", file = paste(tabledir,"tablefit.tex",sep=""), include.rownames=FALSE,include.colnames=TRUE, caption.placement="top",table.placement="htp", sanitize.text.function = function(x){x},hline.after = c(-1,0,nrow(model.data)))

###### WEIGHTED AICs
#aics<-data.frame(paste("m",1:3,sep=""),c(m1$aic,m2$aic,m3$aic),row.names=NULL) 
#colnames(aics)<-c("model","AIC") 
#aics<-aics[order(-aics$AIC),] 
#for(i in 1:dim(aics)[1]){ 
#aics$diff[i]<-aics$AIC[1]-aics$AIC[i]} 
#aics$wi<-2.71828182845904523536^(-0.5*aics$diff) 
#aics$aic.weights<-aics$wi/sum(aics$wi) 

###################################################
### getbestmodel
###################################################
# get the "best" model
###use when doing this after running all models together
best.model = model.tbl[1,]
fitname = paste("kemz",best.model$m,best.model$R,best.model$covariates,sep=".")
best.fit = get(fitname)

###################################################
### varimax rotations
###################################################
# get the inverse of the rotation matrix
H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat

###################################################
### rotations
###################################################
# rotate factor loadings
Z.rot = coef(best.fit, type="matrix")$Z %*% H.inv   

# rotate trends
trends.rot = solve(H.inv) %*% best.fit$states

##only one trend so just Z
Z.rot=coef(best.fit, type='matrix')$Z

###################################################
### plotfacloadings
###################################################
###VARIMAX###
spp = rownames(dat.z)
minZ = 0.2
ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
pdf(paste("loadings.",start,".",end,".pdf"))
par(mfrow=c(3,2), mar=c(2,4,1.5,0.5), oma=c(0.4,1,1,1))
 for(i in 1:best.model$m) {
    plot(c(1:N.ts)[abs(Z.rot[,i])>minZ], as.vector(Z.rot[abs(Z.rot[,i])>minZ,i]),
         type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1))
    for(j in 1:N.ts) {
        if(Z.rot[j,i] > minZ) {text(j, -0.05, spp[j], srt=90, adj=1, cex=0.9)}
        if(Z.rot[j,i] < -minZ) {text(j, 0.05, spp[j], srt=90, adj=0, cex=0.9)}
        abline(h=0, lwd=1, col="gray")
    } # end j loop
    mtext(paste("Factor loadings on trend",i,sep=" "),side=3,line=.5)
} # end i loop
dev.off()


###################################################
### plottrends
###################################################

#VARIMAX rotation: Use the following if there was only one trend in the best model###
pdf(paste("trends.",start,".",end,".pdf"))
plot(seq(start,end),best.fit$states, type='l', lwd=2, xlab="year", ylab="index")
dev.off()


#VARIMAX rotation: Use the following if there was more than one trend in the best model###
# get ts of trends
ts.trends = t(trends.rot)
#ts.trends = t(get_loadings(kemz.1.unconstrained.none)$trends.rot)
pdf(paste("trends.",start,".",end,".pdf"))
par(mfrow=c(ceiling(dim(ts.trends)[2]/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
# loop over each trend
for(i in 1:dim(ts.trends)[2]) {
    # set up plot area
    plot(seq(start,end),ts.trends[,i],
         ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
         type="n", lwd=2, bty="L", 
         xlab="", ylab="",  yaxt="n")
    # draw zero-line
    abline(h=0, col="gray")
    # plot trend line
    par(new=TRUE)
    plot(seq(start,end),ts.trends[,i],
         ylim=c(-1.1,1.1)*max(abs(ts.trends)), 
         type="l", lwd=2, bty="L", 
         xlab="", ylab="", xaxt="n")
    # add panel labels
    mtext(paste("Trend",i,sep=" "), side=3, line=-2,cex=1.5,font=2)
    #axis(1,12*(0:dim(dat.1980)[2])+1,start+0:dim(dat.1980)[2])
} # end i loop (trends)

dev.off()





###################################################
### plotbestfits
###################################################
par.mat=coef(best.fit, type="matrix")

###USE WHEN NO COVARIATES ARE INVOLVED IN THE BEST MODEL###
fit.b = par.mat$Z %*% best.fit$states + matrix(par.mat$A, nrow=N.ts, ncol=TT)

###USE WHEN COVARIATES ARE INVOLVED IN THE BEST MODEL###
#fit.b = par.mat$Z %*% best.fit$states + matrix(best.fit$par$A, 20,1) %*% pop.z


#########Plot all pressures on one pdf page#######
#names = c("Coastal_engineering","Nutrients","Inorganic_pollution","Commercial_landings_of_all_fisheries","NOI_summer","NOI_winter","PDO_summer","PDO_winter","NPGO_summer","NPGO_winter","GF_spp_richness_coastwide","GF-Simp_coastwide","Scav_ratio_coastwide","GF_MTL_coastwide","Cop_Spp_Rich_Anom_summer","Cop_Spp_Rich_Anom-_winter","California_sea_lion_pup_production")
#names = c("Coastal_engineering","Nutrients","Inorganic_pollution","Commercial_landings_of_all_fisheries","NOI_summer","NOI_winter","PDO_summer","PDO_winter","NPGO_summer","NPGO_winter","GF_spp_richness_coastwide","GF-Simp_coastwide","Scav_ratio_coastwide","Cop_Spp_Rich_Anom_summer","Cop_Spp_Rich_Anom-_winter","California_sea_lion_pup_production")

pdf(paste("best.fits.",start,".",end,".pdf"))
par(mfrow=c(3,5), mar=c(3,4,1.5,1), oma=c(0.4,1,1,1))
for(i in 1:length(names)){
    plot(seq(start,end), dat.z[i,],xlab="",ylab="Normalized index",bty="L", ylim=c(-3,3), pch=16, col="blue")
    #this plots the lines only across years where there is real data
    lines(seq(start,end)[is.na(dat.z[i,])==F],fit.b[i,][is.na(dat.z[i,])==F], lwd=2)
    title(names[i],cex.main=0.8)
    }
dev.off()


####################################################
### KIR PLOTS
####################################################


kemzMod<-kemz.1.equalvarcov.none
kemzMod<-get(fitname)
m<-2

# get_loadings(kemzMod)$Z.rot  # get loadings 

# plot_MARSS_load(
#       tmp_MARSS=kemzMod,
#       lab="kemz.2.equalvarcov.none")  

# plot_MARSS_load2(
#       tmp_MARSS=get("kemz.1.equalvarcov.none"),
#       lab="kemz.1.equalvarcov.none")  

# plot_MARSS_load2(
#       tmp_MARSS=get("kemz.2.equalvarcov.none"),
#       lab="kemz.2.equalvarcov.none")  

plot_MARSS_load2(
      tmp_MARSS=get(fitname),
      lab=fitname)  

quartz.save(file=file.path(tabledir,datagroup,fitname,".pdf",sep=""),type="pdf",dpi=500)



