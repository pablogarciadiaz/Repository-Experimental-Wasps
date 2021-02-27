#### Bayesian optimal design - Wasps Quadratic Loss - minimise the 1/(msqe)
#### Script 2 - optimal design in 2 visits, how many sites? MH algorithm
#### Trade-offs number of sites vs distance walked & maximum distance walked
##### Effort within function
##### MSE: abundance per cell & detection coefficients
### Multinomial distribution of values

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
data.wasp<-read.table("e:/CONTAIN/Experimental Design/Wasps/Nest-Count-OptimalExperimental.csv", header=TRUE, sep=",")

id.exc<-which(data.wasp$pop.dens==0 | is.na(data.wasp$pop.dens))

data.wasp<-data.wasp[-id.exc, ]

summary(data.wasp)

nrow(data.wasp)

#### Nimble model to fit
model.wt<-nimbleCode({

    ### Priors
    ### Intercept of the abundance model
    alpha.ab~dnorm(0, var=100)

    ### Quadratic error
    quad.loss[1]<-pow(param[1]-alpha.ab, 2)

    ### Slopes of the abundance model
    for (j in 1:slope.ab){

        beta.ab[j]~dnorm(0, var=100)

        quad.loss[j+1]<-pow(param[j+1]-beta.ab[j], 2)
    }

    ### Probability of detection
    p.det~dbeta(1, 1)

    quad.loss[slope.ab+2]<-pow(param[slope.ab+2]-p.det, 2)

    ### 1st session before control

    for (i in 1:n.site){

        ### Abundance
        log(mean.ab[i])<-alpha.ab+beta.ab[1]*mean.ndvi[i]+beta.ab[2]*var.ndvi[i]+beta.ab[3]*water.length[i]

        n0[i]~dpois(mean.ab[i])

            ### Survey occasions
            for (j in 1:n.occ){

                y0[i, j]~dbinom(p.det, n0[i])

            }

    }

    ### Mean squared error
    loss.fun<-mean(quad.loss[1:(slope.ab+2)])

})


#### Utility function to determine the inverse of the determinant
ut.fun1<-function(n.site, transect.reps, tot.dist, n.occ, data.wasp){

    ### Define effort per cell (walking distance)
    effort<-rep(500*transect.reps, n.site)

    ### Choosing cells based on their distance to CEHUM
    seq.tot<-seq(min(data.wasp$dist.cehum), max(data.wasp$dist.cehum), by=5000)

    seq.band<-seq.tot[-length(seq.tot)]

    weight.fun<-function(x) length(which(data.wasp$dist.cehum>=x & data.wasp$dist.cehum<x+5000))

    count.cells.per.band<-sapply(seq.band, weight.fun)

    #weight.vals<-count.cells.per.band*(1/seq.band)

    #### Number of cells for each band of distance to CEHUM (the closeer the better)
    cell.per.band<-rmultinom(1, n.site[1], 1/seq.band)

    sel.band<-seq.band[which(cell.per.band!=0)]

    number.per.band<-cell.per.band[which(cell.per.band!=0)]

    ### Sampling at random from the selected bands
    sel.fun<-function(x) sample(which(data.wasp$dist.cehum>=x[1] & data.wasp$dist.cehum<x[1]+5000), x[2])

    id.sel<-c(unlist(apply(data.frame(sel.band, number.per.band), 1, sel.fun)))

    #### Chosen data
    data.wasp<-data.wasp[id.sel, ]

    ### NDVI
    mean.ndvi<-(data.wasp$median.ndvi-mean(data.wasp$median.ndvi))/sd(data.wasp$median.ndvi)

    var.ndvi<-(data.wasp$var.ndvi-mean(data.wasp$var.ndvi))/sd(data.wasp$var.ndvi)

    #### Water body length in each cell
    water.length<-(data.wasp$lenght.water-mean(data.wasp$lenght.water))/sd(data.wasp$lenght.water)

    cov.data<-data.frame(mean.ndvi=mean.ndvi, var.ndvi=var.ndvi, water.length=water.length)

    #cov.data

    ######## Data-generation part
    #### Initial abundance in each cell
    alpha.ab<-rnorm(1, 0, 2)

    beta.ab<-c(runif(2, 0, 2), runif(1, -1, 1))

    mean.pop<-exp(alpha.ab+beta.ab[1]*cov.data$mean.ndvi+beta.ab[2]*cov.data$water.length+beta.ab[3]*cov.data$var.ndvi)

    n0<-rpois(n.site, lambda=mean.pop)

    ### Intercept and slope of the effect of survey effort
    alpha.det<-runif(1, -5, 0)

    beta.det<-runif(1, 0, 5)

    ### 1st session - generate detection histories
    y0<-matrix(0, nrow=n.site, ncol=n.occ)

    for (i in 1:n.site){

            ### Probability of detection by site and occasion
            p.det1<-plogis(alpha.det+beta.det*log10(effort[i]))

            y0[i, ]<-rbinom(2, n0[i], p.det1)

            }

   ## Bayesian modelling
    ### Data for the model
    data.tot<-list(y0=y0, param=c(alpha.ab, beta.ab, p.det1),
    mean.ndvi=cov.data$mean.ndvi, var.ndvi=cov.data$var.ndvi, water.length=cov.data$water.length)

    inits<-list(p.det=runif(1, 0, 1), alpha.ab=0, beta.ab=rep(0, 3),
                n0=apply(y0, 1, max)*2)

    ### Constants
    constants<-list(n.site=n.site, n.occ=n.occ, slope.ab=3)

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
n.sim<-20

n.site<-transect.reps<-tot.dist<-rep(NA, n.sim)

#### Initial design
lh<-improvedLHS(1, 2, dup=5)								### Generate values

transect.reps[1]<-round(qunif(lh[, 1], 1, 6))                #### Number of 500-metre sections per site

tot.dist[1]<-round(qunif(lh[, 2], 10000, 50000), digits=0)                #### Total distance walked to distribute across sites

n.site[1]<-round(tot.dist[1]/(500*transect.reps[1]), digits=0)

n.occ<-2

#### Simulating stuff!!
### Number of repeats per step
n.rep<-10

ut.sim<-rep(NA, n.sim)      ### To store utility value

#### Timing the function
start.time<-Sys.time()

### Set a progress bar
pb<-txtProgressBar(min = 0, max = n.sim, style = 3)

#### First design
registerDoParallel(n.core)

rep.sim<-foreach(i=1:n.rep, .combine = c, .packages=c("nimble")) %dopar% ut.fun(n.site[1], transect.reps=transect.reps[1],
                    tot.dist=tot.dist[1], n.occ=n.occ, data.wasp=data.wasp)

ut.sim[1]<-median(rep.sim, na.rm=TRUE)

setTxtProgressBar(pb, 1)

#### Stop cluster
stopImplicitCluster()

#### Rest of simulations

for (j in 2:n.sim){

    new.transect.reps<-round(runif(1, 1, 6), digits=0)                #### Number of 500-metre sections per site

    new.tot.dist<-round(runif(1, 10000, 50000), digits=0)                #### Total distance walked to distribute across sites

    new.site<-round(new.tot.dist/(500*new.transect.reps), digits=0)

    ### Estimate the utility of the new design
    registerDoParallel(n.core)

    new.ut<-foreach(i=1:n.rep, .combine = c, .packages=c("nimble")) %dopar% ut.fun(new.site, new.transect.reps,
                new.tot.dist, n.occ=n.occ, data.wasp=data.wasp)

    curr.ut<-median(new.ut, na.rm=TRUE)

    #### Stop cluster
    stopImplicitCluster()

    #### Metropolis
    acc.thre<-curr.ut/ut.sim[j-1]      ### Acceptance criteria

     rd.num<-runif(1, 0, 1)           ### Random number generation

     if(rd.num<=acc.thre){                ### Accept & update

         ut.sim[j]<-curr.ut

         n.site[j]<-new.site

         transect.reps[j]<-new.transect.reps

         tot.dist[j]<-new.tot.dist

     }else{                              ### Reject proposal and keep previous iteration

         ut.sim[j]<-ut.sim[j-1]

         n.site[j]<-n.site[j-1]

         transect.reps[j]<-transect.reps[j-1]

         tot.dist[j]<-tot.dist[j-1]
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

transect.reps

tot.dist

## Plotting
par(mfrow=c(2, 2))

plot(ut.sim~c(1:n.sim), type="l", lwd=2, xlab="Iteration", ylab="Utility", main="MH Search - Inverse of MSE")

#### Number of cells to survey
hist(n.site, xlab="Number of cells to survey", main="Number of cells to survey")

#### Sections per cell
hist(transect.reps, xlab="Number of 500-m transects per cell",
            main="Number of 500-m transects per cell")

#### Total walking effort across all cells
hist(tot.dist, xlab="Total distance walked (metres)",
            main="Total distance across all sits (metres walked)")


### Optimal design based on the mean, median, and mode of each trap
library(modeest)

### Mode
#### Sites
mode.site<-meanshift(n.site)

mode.site

## 500-m sections per cell
mode.reps<-meanshift(transect.reps)

mode.reps

## Total distance walked
mode.totdist<-meanshift(tot.dist)

mode.totdist

### Mean
mean.site<-mean(n.site)

mean.site

mean.reps<-mean(transect.reps)

mean.reps

mean.totdist<-mean(tot.dist)

mean.totdist

### Median
median.site<-median(n.site)

median.site

median.reps<-median(transect.reps)

median.reps

median.totdist<-median(tot.dist)

median.totdist

#### Export table
data.exp<-data.frame(mode.site=meanshift(n.site), mode.reps=meanshift(transect.reps), mode.totdist=meanshift(tot.dist),
    mean.site=mean(n.site), mean.reps=mean(transect.reps), mean.totdist=mean(tot.dist),
    median.site=median(n.site), median.reps=median(transect.reps), median.tot.dist=median(tot.dist))

data.exp

write.table(data.exp, "e:/CONTAIN/Experimental Design/Wasps/OptimalMH-Distance(RandomMSE).csv", sep=",")
