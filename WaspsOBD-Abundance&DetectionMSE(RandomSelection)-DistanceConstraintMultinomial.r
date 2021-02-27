#### Bayesian optimal design - Wasps Quadratic Loss - minimise the 1/(msqe)
#### Script 2 - optimal design in 2 visits, how many sites? MH algorithm
#### Chosen purely at random - land use as a random effect
##### Varying distance surveyed in each cell + number of sites
##### Effort within function
##### MSE: abundance per cell & detection coefficients
### Multinomial distribution of values

sum(rmultinom(1, max.km, dist))

### Load the libraries
library(coda)
library(MASS)
library(MCMCpack)
library(MuMIn)
library(jagsUI)
library(ggmcmc)
library(corrplot)
library(nimble)
library(dclone)
library(beepr)
library(parallel)
library(doParallel)
library(foreach)
library(DHARMa)
library(compiler)
library(lhs)
library(modeest)
library(truncdist)

### Load data
data.wasp<-read.table("e:/CONTAIN/Experimental Design/Wasps/Nest-Count-OptimalExperimental2.csv", header=TRUE, sep=",")

summary(data.wasp)

#### Nimble model to fit
model.wt<-nimbleCode({

    ### Priors
    ### Intercepts of the abundance model - land uses
    for (i in 1:n.use){

        alpha.ab[i]~dnorm(0, var=var.int)

    }

    var.int~dunif(0, 100)

    ### Slopes of the abundance model
    for (j in 1:slope.ab){

        beta.ab[j]~dnorm(0, var=100)

    }

    #### Detection
    alpha.det~dnorm(0, var=100)

    beta.det~dnorm(0, var=100)

    ### 1st session before control

    for (i in 1:n.site){

        ### Abundance
        log(mean.ab[i])<-alpha.ab[land.cov[i]]+beta.ab[1]*pop.dens[i]+beta.ab[2]*mean.ndvi[i]+beta.ab[3]*var.ndvi[i]+beta.ab[4]*water.length[i]

        n0[i]~dpois(mean.ab[i])

        ### Quadratic error
        quad.loss[i]<-pow(param[i]-n0[i], 2)

        #### Detection
        logit(p.det1[i])<-alpha.det+beta.det*effort[i]

            ### Subsequent occasions
            for (j in 1:n.occ){

                p.det2[i, j]<-p.det1[i]

                y0[i, j]~dbinom(p.det2[i, j], n0[i])

            }

    }

    ### Mean loss
    quad.loss[n.site+1]<-pow(param[n.site+1]-alpha.det, 2)

    quad.loss[n.site+2]<-pow(param[n.site+2]-beta.det, 2)

    loss.fun<-mean(quad.loss[1:(n.site+2)])

})


#### Utility function to determine the inverse of the determinant
ut.fun1<-function(n.site, av.eff, sd.eff, n.occ, data.wasp){

    ### Define effort
    effort<-rnorm(n.site, av.eff, sd.eff)

    effort<-ifelse(effort<=0, av.eff, effort)

    ### Select cells to survey at random
    random.samp<-sample(c(1:nrow(data.wasp)), n.site)

    #### Data-generation function
    data.wasp<-data.wasp[random.samp, ]

    #### Land uses
    data.wasp$land.use<-factor(data.wasp$land.use)

    data.wasp$land.cov<-as.numeric(data.wasp$land.use)

    ### Number of land cover types
    n.use<-max(data.wasp$land.cov)

    ### Population density & standardise
    density<-(data.wasp$pop.dens-mean(data.wasp$pop.dens))/sd(data.wasp$pop.dens)

    ### NDVI
    mean.ndvi<-(data.wasp$median.ndvi-mean(data.wasp$median.ndvi))/sd(data.wasp$median.ndvi)

    var.ndvi<-(data.wasp$var.ndvi-mean(data.wasp$var.ndvi))/sd(data.wasp$var.ndvi)

    #### Water body length in each cell
    water.length<-(data.wasp$lenght.water-mean(data.wasp$lenght.water))/sd(data.wasp$lenght.water)

    cov.data<-data.frame(land.cov=data.wasp$land.cov, pop.dens=density, mean.ndvi=mean.ndvi, var.ndvi=var.ndvi, water.length=water.length)

    ######## Data-generation part
    #### Initial abundance in each cell
    alpha.ab<-rnorm(n.use, 0, 2)

    beta.ab<-c(runif(3, 0, 2), runif(1, -1, 1))

    mean.pop<-exp(alpha.ab[cov.data$land.cov]+beta.ab[1]*cov.data$pop.dens+beta.ab[2]*cov.data$mean.ndvi+beta.ab[3]*cov.data$water.length+beta.ab[4]*cov.data$var.ndvi)

    n0<-rpois(n.site, lambda=mean.pop)

    ### Intercept and slope of the effect of survey effort
    alpha.det<-runif(1, -5, 5)

    beta.det<-runif(1, 0, 5)

    ### 1st session - generate detection histories
    y0<-matrix(0, nrow=n.site, ncol=n.occ)

    for (i in 1:n.site){

            ### Probability of detection by site and occasion
            p.det1<-plogis(alpha.det+beta.det*log10(effort[i]+1))

            y0[i, ]<-rbinom(2, n0[i], p.det1)

            }

   ## Bayesian modelling
    ### Data for the model
    data.tot<-list(y0=y0, effort=log10(effort+1), param=c(n0, alpha.det, beta.det),
    mean.ndvi=cov.data$mean.ndvi, var.ndvi=cov.data$var.ndvi, pop.dens=cov.data$pop.dens,
     water.length=cov.data$water.length, land.cov=cov.data$land.cov)

    inits<-list(alpha.det=0, beta.det=0, alpha.ab=rep(0, n.use), beta.ab=rep(0, 4),
                n0=apply(y0, 1, max)*2)

    ### Constants
    constants<-list(n.site=n.site, n.use=n.use, n.occ=n.occ, slope.ab=4)

    ###### THE MODEL
    Rmodel<-nimbleModel(code=model.wt, constants=constants, data=data.tot, inits=inits)

    wt.conf<-configureMCMC(Rmodel,  monitors=list('loss.fun'), thin=1)

    Rmcmc<-buildMCMC(wt.conf)

    Cmodel<-compileNimble(Rmodel)

    Cmcmc<-compileNimble(Rmcmc, project = Rmodel)

    ### Three chains and check for burnin
    m2.wt<-runMCMC(Cmcmc, niter=50000, nchains=1, nburnin=5000, inits=inits, samplesAsCodaMCMC=TRUE, progress=TRUE)

    ### Inverse of the mean squared error
    return(1/mean(m2.wt, na.rm=TRUE))
}

ut.fun<-cmpfun(ut.fun1)

### Check number of cores for parallel computing
n.core<-detectCores()

n.core<-5

#### Simulations
### Number of simulation steps
n.sim<-200

n.site<-max.dist<-rep(NA, n.sim)

#### Initial design
lh<-improvedLHS(1, 2, dup=5)								### Generate values

n.site[1]<-round(qunif(lh[,1], 2, 10), digits=0)						### Number of survey sites

max.dist[1]<-round(qunif(lh[, 2], n.site[1]*1*1000, n.site[1]*5*1000))  #### Maximum total distance walked

n.occ<-2

### Multinomial distribution of values
dist<-round(rnorm(nrow(data.wasp), 10000, 1000), digits=-3)

dist

seq.tot<-seq(min(dist), max(dist), by=1000)

seq.tot

#### Number of cells for each band of distance to CEHUM (the closeer the better)
cell.per.band<-rmultinom(1, n.site[1], 1/seq.tot)

sel.band<-seq.tot[which(cell.per.band!=0)]

sel.band

number.per.band<-cell.per.band[which(cell.per.band!=0)]

number.per.band

### Sampling at random from the selected bands
id.sel<-rep(NA, sum(number.per.band))

for (i in 1:length(id.sel)){

    id.sel[i]<-sample(which(dist==sel.band[i]), number.per.band[i])

}

#### Simulating stuff!!
### Number of repeats per step
n.rep<-20

ut.sim<-rep(NA, n.sim)      ### To store utility value

#### Timing the function
start.time<-Sys.time()

### Set a progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

#### First design
registerDoParallel(n.core)

rep.sim<-foreach(i=1:n.rep, .combine = c, .packages=c("nimble")) %dopar% ut.fun(n.site[1], av.eff=av.eff[1], sd.eff=sd.eff[1], n.occ=n.occ, data.wasp=data.wasp)

ut.sim[1]<-median(rep.sim, na.rm=TRUE)

setTxtProgressBar(pb, 1)

#### Stop cluster
stopImplicitCluster()

#### Rest of simulations

for (j in 2:n.sim){

    #### Propose a new design
    new.site<-round(runif(1, 2, 10), digits=0)						### Number of survey sites

    new.aveffort<-runif(1, 1000, 1500)                      ### Mean effort per cell (in metres walked)

    new.sdeffort<-runif(1, 100, 500)                      ### Standard error of the effort per site (metres walked)

    ### Estimate the utility of the new design
    registerDoParallel(n.core)

    new.ut<-foreach(i=1:n.rep, .combine = c, .packages=c("nimble")) %dopar% ut.fun(new.site, new.aveffort, new.sdeffort, n.occ=n.occ, data.wasp=data.wasp)

    curr.ut<-median(new.ut, na.rm=TRUE)

    #### Stop cluster
    stopImplicitCluster()

    #### Metropolis
    acc.thre<-curr.ut/ut.sim[j-1]      ### Acceptance criteria

     rd.num<-runif(1, 0, 1)           ### Random number generation

     if(rd.num<=acc.thre){                ### Accept & update

         ut.sim[j]<-curr.ut

         n.site[j]<-new.site

         av.eff[j]<-new.aveffort

         sd.eff[j]<-new.sdeffort

     }else{                              ### Reject proposal and keep previous iteration

         ut.sim[j]<-ut.sim[j-1]

         n.site[j]<-n.site[j-1]

         av.eff[j]<-av.eff[j-1]

         sd.eff[j]<-sd.eff[j-1]

     }

    setTxtProgressBar(pb, j)

}

### Close progress bar
close(pb)

## End time
end.time<-Sys.time()

end.time-start.time     ### Time to run

beep(sound=8)

#### Outputs of the designs
ut.sim

hist(ut.sim)

n.site

sd.eff

## Plotting
par(mfrow=c(2, 2))

plot(ut.sim~c(1:n.sim), type="l", lwd=2, xlab="Iteration", ylab="Utility", main="MH Search - Inverse of MSE")

#### Number of cells to survey
hist(n.site, xlab="Number of cells to survey", main="Number of cells to survey randomly")

#### Mean effort
hist(av.eff, xlab="Mean survey effort (metres)",
            main="Mean survey effort (metres walked)")

#### SD effort
hist(sd.eff, xlab="Standard deviation of the survey effort (metres)",
            main="Standard deviation of the survey effort (metres walked)")


### Optimal design based on the mean, median, and mode of each trap (n.trap = 12)
library(modeest)

### Mode
#### Sites
mode.site<-meanshift(n.site)

mode.site

## Mean effort
mode.meaneff<-meanshift(av.eff)

mode.meaneff

## Variance effort
mode.vareff<-meanshift(sd.eff)

mode.vareff

### Mean
mean.site<-mean(n.site)

mean.site

mean.meaneff<-mean(av.eff)

mean.meaneff

mean.sdeff<-sd(sd.eff)

mean.sdeff

### Median
median.site<-median(n.site)

median.site

median.meaneff<-median(av.eff)

median.meaneff

median.sdeff<-sd(sd.eff)

median.sdeff

#### Export table
data.exp<-data.frame(mode.site=meanshift(n.site), mode.meaneff=meanshift(av.eff), mode.vareff=meanshift(sd.eff),
    mean.site=mean(n.site), mean.meaneff=mean(av.eff), mean.sdeff=sd(sd.eff),
    median.site=median(n.site), median.meaneff=median(av.eff), median.sdeff=sd(sd.eff))

data.exp

write.table(data.exp, "e:/CONTAIN/Experimental Design/Wasps/OptimalMH-Distance(RandomMSE).csv", sep=",")
