
setwd("C:/Users/thoma/Rprojects/BSMdata")
we.data=read.csv("hg_age_lat_yr_WE.csv")
we.data$sizecat=cut(we.data$Total_Length, c(min(we.data$Total_Length-10), 372, 435, 508, max(we.data$Total_Length)),
                    + labels=c("minnow", "small","medium", "large"))
###Gompertz model fitted
gomp.we.hg=(nls(VALUE~SSgompertz(Assessed_Fish_Age, Asym, b2, b3)))
summary(gomp.we.hg)

logis.we.hg=(nls(VALUE~SSlogis(Assessed_Fish_Age, Asym, xmid, sccal)))
summary(logis.we.hg)

png("Hg_age_WE.png", width=5, height=4, units="in", res=1200)
plot(Assessed_Fish_Age, VALUE, xlab = "Age",ylab = "Mercury")
curve(SSlogis(x, 1.66, 9.91, 4.91), add=TRUE, col="green", lwd=2)
curve(SSgompertz(x, 1.99, 2.54, 0.89), add=TRUE, col="red", lwd=2)
abline(lm(VALUE~Assessed_Fish_Age), col="grey", lwd=2)
legend("topleft", c("Logistic","Gompertz", "Linear"), col=c("green", "red", "grey"), lty=c(1,1,1), cex=1, bty = "n")
dev.off()

#########################################
AIC(lm.we.hg, gomp.we.hg, logis.we.hg)  
#             df    AIC
#lm.we.hg     3 2379.161
#gomp.we.hg   4 2367.458
#logis.we.hg  4 2371.548
#########################################

##########################################################################################################################################################
#####################################################TESTING GROWTH DILUTION IN WALLEYE AND NORTHERN PIKE#################################################
########################################################using linear mixed effects regression models######################################################


WE.bsm=read.csv("hg_age_lat_yr_WE.csv")
head(WE.bsm)

WE.age.len_Hg=subset(WE.bsm[,c(8,9,13,15,16,18,28,29,30)])
head(WE.age.len_Hg)
WE.age.len_Hg$Age=WE.bsm$Assessed_Fish_Age
write.csv(WE.age.len_Hg, file="WE.age.len_HG")

we.growth.lmer=lmer(log(WE.age.len_Hg$Total_Length)~ log(WE.age.len_Hg$Age+1)+(1+log(WE.age.len_Hg$Age+1)|WE.age.len_Hg$LatZones))
summary(we.growth.lmer)
coef(we.growth.lmer)

we.biomag.lmer=lmer(log(WE.age.len_Hg$VALUE)~log(WE.age.len_Hg$Total_Length)+(1+log(WE.age.len_Hg$Total_Length)|WE.age.len_Hg$LatZones))
summary(we.biomag.lmer)
coef(we.biomag.lmer)


library(RColorBrewer)
library(rgdal)
myColors = brewer.pal(3,"Set1")
colScale = scale_colour_manual(name = "LatZones",values = myColors)
WE.growth.lats=ggplot(WE.age.len_Hg,aes(x=log(WE.age.len_Hg$Age), y=log(WE.age.len_Hg$Total_Length),colour = LatZones)) + geom_point(shape="o")+
geom_smooth(method='lm', col= "black", lwd=1.25, lty=2)
WE.growth.lats+geom_line(aes(y = predict(we.growth.lmer)),size=1.5) + colScale+ylab("Length")+xlab("Age")+theme_bw()+theme(legend.title=element_blank())+
theme(legend.justification=c(1,0), legend.position=c(1,0))

png("WE.grow.lats.png", width = 5, height = 4, units = 'in', res = 900)



we.growth.lmer.lat=lmer(log(WE.age.len_Hg$Total_Length)~ log(WE.age.len_Hg$Age+1)+(1+log(WE.age.len_Hg$Age+1)|WE.age.len_Hg$LAT))
summary(we.growth.lmer.lat)
ranef(we.growth.lmer.lat)
fixef(we.growth.lmer.lat)
coef(we.growth.lmer.lat)
growth.coefs=coef(we.growth.lmer.lat)$`WE.age.len_Hg$LAT`[,2]

we.biomag.lmer.lat=lmer(log(WE.age.len_Hg$VALUE)~log(WE.age.len_Hg$Total_Length)+(1+log(WE.age.len_Hg$Total_Length)|WE.age.len_Hg$LAT))
summary(we.biomag.lmer.lat)
ranef(we.biomag.lmer.lat)
fixef(we.biomag.lmer.lat)
coef(we.biomag.lmer.lat)
biomag.coefs=coef(we.biomag.lmer.lat)$`WE.age.len_Hg$LAT`[,2]

lats.we=sort(unique(WE.age.len_Hg$LAT))
WE.growth_biomag.coefs=as.data.frame(cbind(lats.we, growth.coefs, biomag.coefs))

ggplot(we.growth.biomag,aes(x=growth.coefs, y=biomag.coef ,colour = latzones)) + geom_point(shape=16)+
geom_smooth(method='lm', col= "black", lwd=1.25, lty=1)

#Randomize latitudes to test significance of correlations
rand.biomag.mods=rlply(10, lmer(log(WE.age.len_Hg$VALUE)~log(WE.age.len_Hg$Total_Length)+(1+log(WE.age.len_Hg$Total_Length)|sample(WE.age.len_Hg$LAT))))
coef.biomag.rands=lapply(rand.mods, function(f) coef(f)$`sample(WE.age.len_Hg$LAT)`[,2])
list.means=sapply(coef.rands, mean)
list.sd=sapply(coef.rands, sd)
median.biomag=sapply(coef.biomag.rands, quantile, probs =0.500)
up.biomag=sapply(coef.bioag.rands, quantile, probs =0.999)
low.biomag=sapply(coef.biomag.rands, quantile, probs =0.001)

P=ggplot(we.growth.biomag,aes(x=growth.coefs, y=biomag.coef)) + geom_point(shape="o")+
geom_smooth(method='lm', col= "black", lwd=1.25, lty=1)+
geom_abline(intercept = 1.71, slope = -0.003, lty=2, lwd=1.25)
P+theme_bw()
