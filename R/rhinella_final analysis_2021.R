setwd("C:/Users/USER/Dropbox/Producoes/Rhinella/rhinella amph-reptilia 2020-2021/rhinella amph-reptilia 2020-2021 Revisao 1/analises para publicacao 2021")

# install.packages("bbmle")
library(bbmle)

# install.packages("MuMIn")
library(MuMIn)

# install.packages("emmeans")
library(emmeans)


rdata=read.csv(file="rhinella_final data.csv",header=TRUE, sep=";")

#================================================
#==============Body condition index==============

regres_rhinella=lm(rdata$totalmass ~ rdata$svl)
plot(rdata$totalmass ~ rdata$svl) 
abline(regres_rhinella) 

bci = (residuals(regres_rhinella))*2 #Body condition index (BCI)
rdata$bci = bci 

hist(bci)
head(rdata)

rdata$bcipad <- as.numeric(scale(rdata$bci, center= T, scale = T)) # scaling BCI values to standard deviation units, to use as a predictor variable below.
hist(rdata$bcipad)

length(rdata$bci[rdata$locality=="A"])
mean(rdata$bci[rdata$locality=="A"])
sd(rdata$bci[rdata$locality=="A"])

length(rdata$bci[rdata$locality=="B"])
mean(rdata$bci[rdata$locality=="B"])
sd(rdata$bci[rdata$locality=="B"])

length(rdata$bci[rdata$locality=="C"])
mean(rdata$bci[rdata$locality=="C"])
sd(rdata$bci[rdata$locality=="C"])

bcim0 <- lm(bci ~ 1, rdata)
bcim1 <- lm(bci ~ locality, rdata)

ICtab(bcim0, bcim1, type="AIC", delta=TRUE, base=TRUE, weights=TRUE)

#tiff("Fig2.tif", width = 6.5, height = 8, units = 'in', res = 300)
par(bty="l", cex=1.4, font.main=3, tcl=-0.5)
boxplot(bci~locality, rdata, ylab="Body condition index (BCI)", 
        xaxt="n", xlab="Site", col="white")
axis(1, labels=c("Polluted", "Intermediate", "Reference"), at=c(1:3))
#dev.off()



#=================================================
#==============Organ Somatic Indeces==============

#======Liver=======


rdata$figado = rdata$liver/rdata$totalmass
head(rdata)
hist(rdata$figado)

length(rdata$figado[rdata$locality=="A"])
mean(rdata$figado[rdata$locality=="A"])
sd(rdata$figado[rdata$locality=="A"])

length(rdata$figado[rdata$locality=="B"])
mean(rdata$figado[rdata$locality=="B"])
sd(rdata$figado[rdata$locality=="B"])

length(rdata$figado[rdata$locality=="C"])
mean(rdata$figado[rdata$locality=="C"])
sd(rdata$figado[rdata$locality=="C"])

fm0 = lm(figado ~ 1, data=rdata)
fm1 = lm(figado ~locality, data=rdata)
fm2 = lm(figado ~ bcipad, data=rdata)
fm3 = lm(figado ~ bcipad + locality, data=rdata)
fm4 = lm(figado ~ bcipad * locality, data=rdata)

ICtab(fm0, fm1, fm2, fm3, fm4, type="AIC", delta=TRUE, base=TRUE, weights=TRUE)


models <- list(fm3,fm4)
model.avg(models)
summary(model.avg(models)) # get averaged coefficients
confint(model.avg(models))

summary(fm3)
confint(fm3)

?confint

fig_emmeans <- emmeans(fm3, "locality")
pairs(fig_emmeans)
confint(pairs(fig_emmeans))


#tiff("Fig3_figado.tif", width = 6, height = 8, units = 'in', res = 300)
par(bty="l", cex=1.3, font.main=3, tcl=-0.5)
boxplot(figado~locality, rdata, ylab="Hepatosomatic index (%)", 
        xaxt="n", ylim = c(0.013, 0.03),
        xlab="Site", col="white")
text(x= 1, y= 0.0285, labels= "a")
text(x= 2, y= 0.021, labels= "b")
text(x= 3, y= 0.0182, labels= "b")
axis(1, labels=c("Polluted", "Intermediate", "Reference"), at=c(1:3))
#dev.off()


#======Kidneys=======

rdata$rins = rdata$kidneys/rdata$totalmass
colnames(rdata)
hist(rdata$rins)

length(rdata$rins[rdata$locality=="A"])
mean(rdata$rins[rdata$locality=="A"])
sd(rdata$rins[rdata$locality=="A"])

length(rdata$rins[rdata$locality=="B"])
mean(rdata$rins[rdata$locality=="B"])
sd(rdata$rins[rdata$locality=="B"])

length(rdata$rins[rdata$locality=="C"])
mean(rdata$rins[rdata$locality=="C"])
sd(rdata$rins[rdata$locality=="C"])

rm0 = lm(rins~1, data=rdata)
rm1 = lm(rins~locality, data=rdata)
rm2 = lm(rins~bcipad, data=rdata)
rm3 = lm(rins~bcipad + locality, data=rdata)
rm4 = lm(rins~bcipad * locality, data=rdata)

ICtab(rm0, rm1, rm2, rm3, rm4, type="AIC", delta=TRUE, base=TRUE, weights=TRUE)

models <- list(rm1, rm3,rm4)
model.avg(models)
summary(model.avg(models)) # get averaged coefficients
confint(model.avg(models))

summary(rm1)
confint(rm1)

rins_emmeans <- emmeans(rm1,  "locality")
pairs(rins_emmeans)
confint(pairs(rins_emmeans))


#tiff("Fig3_rins.tif", width = 6, height = 8, units = 'in', res = 300)
par(bty="l", cex=1.3, font.main=3, tcl=-0.5)
boxplot(rins~locality, rdata, ylab="Kidneys-somatic index (%)", 
        xaxt="n",ylim = c(0.0025, 0.0060),
        xlab="Site", col="white")
text(x= 1, y= 0.0056, labels= "a")
text(x= 2, y= 0.0052, labels= "b")
text(x= 3, y= 0.0040, labels= "b")
axis(1, labels=c("Polluted", "Intermediate", "Reference"), at=c(1:3))
#dev.off()

#======Spleen=======

rdata$baco = rdata$spleen/rdata$totalmass
colnames(rdata)
hist(rdata$baco)

length(rdata$baco[rdata$locality=="A"])
mean(rdata$baco[rdata$locality=="A"])
sd(rdata$baco[rdata$locality=="A"])

length(rdata$baco[rdata$locality=="B"])
mean(rdata$baco[rdata$locality=="B"])
sd(rdata$baco[rdata$locality=="B"])

length(rdata$baco[rdata$locality=="C"])
mean(rdata$baco[rdata$locality=="C"])
sd(rdata$baco[rdata$locality=="C"])

bm0 = lm(baco~1, data=rdata)
bm1 = lm(baco~locality, data=rdata)
bm2 = lm(baco~bcipad, data=rdata)
bm3 = lm(baco~bcipad + locality, data=rdata)
bm4 = lm(baco~bcipad * locality, data=rdata)

ICtab(bm0, bm1, bm2, bm3, bm4, type="AIC", delta=TRUE, base=TRUE, weights=TRUE)

summary(bm4)
confint(bm4)

baco_emmeans <- emmeans(bm4,  "locality")
pairs(baco_emmeans)
confint(pairs(baco_emmeans))


#tiff("Fig3_baco.tif", width = 6, height = 8, units = 'in', res = 300)
par(bty="l", cex=1.3, font.main=3, tcl=-0.5)
boxplot(baco~locality, rdata, ylab="Spleen-somatic index (%)",
        xaxt="n", xlab="Site", col="white")
text(x= 1, y= 0.00161, labels= "a")
text(x= 2, y= 0.00117, labels= "b")
text(x= 3, y= 0.000778, labels= "b")
axis(1, labels=c("Polluted", "Intermediate", "Reference"), at=c(1:3))
#dev.off()

#======Gonads (testes)=======

rdata$gonadas = rdata$gonads/rdata$totalmass 
colnames(rdata)
hist(rdata$gonadas)

length(rdata$gonadas[rdata$locality=="A"])
mean(rdata$gonadas[rdata$locality=="A"])
sd(rdata$gonadas[rdata$locality=="A"])

length(rdata$gonadas[rdata$locality=="B"])
mean(rdata$gonadas[rdata$locality=="B"])
sd(rdata$gonadas[rdata$locality=="B"])

length(rdata$gonadas[rdata$locality=="C"])
mean(rdata$gonadas[rdata$locality=="C"])
sd(rdata$gonadas[rdata$locality=="C"])

gm0 = lm(gonadas~1, data=rdata)
gm1 = lm(gonadas~locality, data=rdata)
gm2 = lm(gonadas~bcipad, data=rdata)
gm3 = lm(gonadas~bcipad + locality, data=rdata)
gm4 = lm(gonadas~bcipad * locality, data=rdata)

ICtab(gm0, gm1, gm2, gm3, gm4, type="AIC", delta=TRUE, base=TRUE, weights=TRUE)

summary (gm0)

#tiff("Fig3_gonadas.tif", width = 6, height = 8, units = 'in', res = 300)
par(bty="l", cex=1.3, font.main=3, tcl=-0.5)
boxplot(gonadas~locality, rdata, ylab="Gonadosomatic index (%)",
        xaxt="n", xlab="Site", col="white")
axis(1, labels=c("Polluted", "Intermediate", "Reference"), at=c(1:3))
#dev.off()



#=====================================================
#==============Haematological indicators==============

head(rdata)
length(rdata$locality[rdata$locality=="A"])
length(rdata$locality[rdata$locality=="B"])
length(rdata$locality[rdata$locality=="C"])

# Searching for NAs in data.frame
any(is.na(rdata))

nrow(rdata)
?complete.cases #search in all data.frame and remove lines with any NA, in any column
rdata2 = rdata[complete.cases(rdata),] #saving the data.frame without NAs in dependent variables
any(is.na(rdata2))
nrow(rdata2)

length(rdata2$locality[rdata2$locality=="A"])
length(rdata2$locality[rdata2$locality=="B"])
length(rdata2$locality[rdata2$locality=="C"])


#======Neutrophil by linfocite ratio (neutrophilia)======
(razaoNL <- cbind(rdata2$neutrophil, rdata2$lymphocyte))

nl0=glm(razaoNL~1, data=rdata2, family=binomial) 
nl1=glm(razaoNL~locality, data=rdata2, family=binomial)
nl2=glm(razaoNL~bcipad, data=rdata2, family=binomial)
nl3=glm(razaoNL~bcipad + locality, data=rdata2, family=binomial)
nl4=glm(razaoNL~bcipad * locality, data=rdata2, family=binomial)

ICtab(nl0, nl1, nl2, nl3, nl4,
      type="AIC", delta=TRUE, base=TRUE, weights=TRUE)

summary(nl3)
confint(nl3)

nlmodels <- list(nl1, nl3,nl4)
model.avg(nlmodels)
summary(model.avg(nlmodels)) # get averaged coefficients
confint(model.avg(nlmodels))

nl_emmeans <- emmeans(nl3, ~locality)
pairs(nl_emmeans)
confint(pairs(nl_emmeans))

# to understand emmeans results with the script above:
# https://cran.r-project.org/web/packages/emmeans/vignettes/sophisticated.html

nl = rdata2$neutrophil/rdata2$lymphocyte

#tiff("Fig4_NL.tif", width = 6, height = 8, units = 'in', res = 300)
par(bty="l", cex=1.2, font.main=3, tcl=-0.5)
boxplot(nl~rdata2$locality, ylab="Individual stress index (N/L)", 
        xaxt="n", cex=0.8, bty="l", las = 1 , ylim=c(0.0,0.75),
        xlab ="Site", col="white")
text(x= 1, y= 0.14, labels= "a")
text(x= 2, y= 0.55, labels= "b")
text(x= 3, y= 0.714, labels= "c")
axis(1, labels=c("Polluted", "Intermediate", "Reference"), at=c(1:3))
#dev.off()


#======Proportion of eosinophils======

colnames(rdata2)
head(rdata2)

plot(eosinophil~bcipad, data= rdata2, ylab="Total eosinophil count (E)", 
     xlab="Body conditon index", cex=1, bty="l", las = 1)

#

rhinella_eosinofilos0=glm(eosinophil~1, data=rdata2, family=poisson)
rhinella_eosinofilos1=glm(eosinophil~locality,data=rdata2, family=poisson)
rhinella_eosinofilos2=glm(eosinophil~totalparasites, data=rdata2, family=poisson)
rhinella_eosinofilos3=glm(eosinophil~bcipad,data=rdata2, family=poisson)

rhinella_eosinofilos4=glm(eosinophil~ totalparasites + locality,data=rdata2, family=poisson)
rhinella_eosinofilos5=glm(eosinophil~ totalparasites * locality, data=rdata2, family=poisson)

rhinella_eosinofilos6=glm(eosinophil~ bcipad + locality, data=rdata2, family=poisson)
rhinella_eosinofilos7=glm(eosinophil~ bcipad * locality, data=rdata2, family=poisson)


ICtab(rhinella_eosinofilos0, rhinella_eosinofilos1, 
      rhinella_eosinofilos2, rhinella_eosinofilos3, 
      rhinella_eosinofilos4, rhinella_eosinofilos5, 
      rhinella_eosinofilos6, rhinella_eosinofilos7, 
      type="AICc", delta=TRUE, base=TRUE, weights=TRUE) #-> we used the AICc due to small sample size. 

summary (rhinella_eosinofilos6)
confint(rhinella_eosinofilos6)

eosi_emmeans <- emmeans(rhinella_eosinofilos6, ~locality)
pairs(eosi_emmeans)
confint(pairs(eosi_emmeans))

length(rdata2$eosinophil[rdata2$locality=="A"])
mean(rdata2$eosinophil[rdata2$locality=="A"])
sd(rdata2$eosinophil[rdata2$locality=="A"])

length(rdata2$eosinophil[rdata2$locality=="B"])
mean(rdata2$eosinophil[rdata2$locality=="B"])
sd(rdata2$eosinophil[rdata2$locality=="B"])

length(rdata2$eosinophil[rdata2$locality=="C"])
mean(rdata2$eosinophil[rdata2$locality=="C"])
sd(rdata2$eosinophil[rdata2$locality=="C"])



#tiff("Fig4_eosi.tif", width = 6, height = 8, units = 'in', res = 300)
par(bty="l", cex=1.2, font.main=3, tcl=-0.5)
boxplot(eosinophil~locality, data= rdata, ylab="Percentual number of eosinophils (E)", 
        xaxt="n", cex=0.8, bty="l", las = 1,xlab="Site", col="white")
axis(1, labels=c("Polluted", "Intermediate", "Reference"), at=c(1:3))
#dev.off()


#=======================
# Removing ID 39, an outlier

nrow(rdata2)
rdata2 <- rdata2[(rdata2$id != "39"),]
nrow(rdata2)
