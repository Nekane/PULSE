### variables ######################################################################
species <- 5
sp <- 1:5
age <-3
cat <-1:3
numInputOptions <- 3
option <- c("a","b","c")
numInputSeasons <- 4
season <- as.character(1:4)


### data ######################################################################
catchMean <- array(NA, dim=c(species, age, numInputOptions, numInputSeasons), 
                   dimnames=list(sp= sp,cat=cat,season=season,option=option))
catchMean[1,,,] <- array(c(rep(10,12),rep(10,12),rep(8,12)),dim=c(3,4,3), dimnames=list(cat=1:3,season=as.character(1:4),option=c("a","b","c")))
catchMean[2,,,] <- array(c(rep(14,12),rep(12,12),rep(10,12)),dim=c(3,4,3), dimnames=list(cat=1:3,season=as.character(1:4),option=c("a","b","c")))

catchSigma <- array(NA, dim=c(species, age, numInputOptions, numInputSeasons), 
                    dimnames=list(sp= sp,cat=cat,season=season,option=option))
catchSigma[1,,,] <- array(0.1,                                  dim=c(3,4,3), dimnames=list(cat=1:3,season=as.character(1:4),option=c("a","b","c")))
catchSigma[2,,,] <- array(0.2,                                  dim=c(3,4,3), dimnames=list(cat=1:3,season=as.character(1:4),option=c("a","b","c")))


### Calculate length of vector of steps for separate spec ######################
control@spp1Incs    <- max(max(inputSpp1@catchMean[1,,] + (3*inputSpp1@catchSigma[1,,])) /(control@increments - 1 ),
                           max(inputSpp1@catchMean[2,,] + (3*inputSpp1@catchSigma[2,,])) /(control@increments - 1 ),
                           max(inputSpp1@catchMean[3,,] + (3*inputSpp1@catchSigma[3,,])) /(control@increments - 1 ))     
control@spp2Incs	 <- max(max(inputSpp2@catchMean[1,,] + (3*inputSpp2@catchSigma[1,,])) /(control@increments - 1 ),
                         max(inputSpp2@catchMean[2,,] + (3*inputSpp2@catchSigma[2,,])) /(control@increments - 1 ),
                         max(inputSpp2@catchMean[3,,] + (3*inputSpp2@catchSigma[3,,])) /(control@increments - 1 ))
control@spp3Incs       <- max(max(inputSpp3@catchMean[1,,] + (3*inputSpp3@catchSigma[1,,])) /(control@increments - 1 ),
                              max(inputSpp3@catchMean[2,,] + (3*inputSpp3@catchSigma[2,,])) /(control@increments - 1 ),
                              max(inputSpp3@catchMean[3,,] + (3*inputSpp3@catchSigma[3,,])) /(control@increments - 1 ))
control@spp4Incs	 <- max(max(inputSpp4@catchMean[1,,] + (3*inputSpp4@catchSigma[1,,])) /(control@increments - 1 ),
                         max(inputSpp4@catchMean[2,,] + (3*inputSpp4@catchSigma[2,,])) /(control@increments - 1 ),
                         max(inputSpp4@catchMean[3,,] + (3*inputSpp4@catchSigma[3,,])) /(control@increments - 1 ))
control@spp5Incs       <- max(max(inputSpp5@catchMean[1,,] + (3*inputSpp5@catchSigma[1,,])) /(control@increments - 1 ),
                              max(inputSpp5@catchMean[2,,] + (3*inputSpp5@catchSigma[2,,])) /(control@increments - 1 ),
                              max(inputSpp5@catchMean[3,,] + (3*inputSpp5@catchSigma[3,,])) /(control@increments - 1 ))



sppincs <- array(spp1inc,spp2inc,spp3inc,spp4inc,spp5inc, dim=c(rep(spp1inc,5)))
                             
discSteps <- array(NA, dim=c(control@spp1DiscardSteps+1,
                             control@spp2DiscardSteps+1,
                             control@spp3DiscardSteps+1,
                             control@spp4DiscardSteps+1,
                             control@spp5DiscardSteps+1))

eachterm <- array(NA, dim=c(5,1), dimnames=list(sp= 1:5,term=1))
eachterm [1,] <-

catchparms <- array(NA,dim=c(numInputOptions,
                             discSteps,
                             numInputSeasons,species,age,k))

norm <- array(NA, dim=c(5,3,4,3), dimnames=list(sp= 1:5,cat=1:3,season=as.character(1:4),option=c("a","b","c")))


### data ######################################################################
for (sp in 1:species){ # the specie
  for (ag in 1:age){ # the age
    for (ii in 1:numInputOptions){ # the area
      for (jj in 1:numInputSeasons){ # the season
        norm[sp,ag,jj,ii] <- pnorm(sppincs[sp], catchMean[sp,ag,jj,ii], catchSigma[sp,ag,jj,ii])
      }
    }
  }
}

    
    eachterm1 <- ((control@spp1DiscardSteps+1)^2*(control@spp2DiscardSteps+1)^3*(control@spp3DiscardSteps+1)^3*(control@spp4DiscardSteps+1)^3*(control@spp5DiscardSteps+1)^3)
    catchparms[ii,1,,,,,,,,,,,,,,,jj,1,1, ] <- rep(a[2:(k+1)] - a[1:k], each=eachterm1)
    catchparms[ii,,1,,,,,,,,,,,,,,jj,1,2, ] <- rep(b[2:(k+1)] - b[1:k], each=eachterm1)
    catchparms[ii,,,1,,,,,,,,,,,,,jj,1,3, ] <- rep(c[2:(k+1)] - c[1:k], each=eachterm1)

    eachterm2 <- ((control@spp1DiscardSteps+1)^3*(control@spp2DiscardSteps+1)^2*(control@spp3DiscardSteps+1)^3*(control@spp4DiscardSteps+1)^3*(control@spp5DiscardSteps+1)^3)
    catchparms[ii,,,,1,,,,,,,,,,,,jj,2,1, ] <- rep(d[2:(k+1)] - d[1:k], each=eachterm2)
    catchparms[ii,,,,,1,,,,,,,,,,,jj,2,2, ] <- rep(e[2:(k+1)] - e[1:k], each=eachterm2)
    catchparms[ii,,,,,,1,,,,,,,,,,jj,2,3, ] <- rep(f[2:(k+1)] - f[1:k], each=eachterm2)
    eachterm3 <-((control@spp1DiscardSteps+1)^3*(control@spp2DiscardSteps+1)^3*(control@spp3DiscardSteps+1)^2*(control@spp4DiscardSteps+1)^3*(control@spp5DiscardSteps+1)^3)
    catchparms[ii,,,,,,,1,,,,,,,,,jj,3,1, ] <- rep(g[2:(k+1)] - g[1:k], each=eachterm3)
    catchparms[ii,,,,,,,,1,,,,,,,,jj,3,2, ] <- rep(h[2:(k+1)] - h[1:k], each=eachterm3)
    catchparms[ii,,,,,,,,,1,,,,,,,jj,3,3, ] <- rep(i[2:(k+1)] - i[1:k], each=eachterm3)
    eachterm4 <- ((control@spp1DiscardSteps+1)^3*(control@spp2DiscardSteps+1)^3*(control@spp3DiscardSteps+1)^3*(control@spp4DiscardSteps+1)^2*(control@spp5DiscardSteps+1)^3)
    catchparms[ii,,,,,,,,,,1,,,,,,jj,4,1, ] <- rep(j[2:(k+1)] - j[1:k], each=eachterm4)
    catchparms[ii,,,,,,,,,,,1,,,,,jj,4,2, ] <- rep(l[2:(k+1)] - l[1:k], each=eachterm4)
    catchparms[ii,,,,,,,,,,,,1,,,,jj,4,3, ] <- rep(m[2:(k+1)] - m[1:k], each=eachterm4)
    eachterm5 <- ((control@spp1DiscardSteps+1)^3*(control@spp2DiscardSteps+1)^3*(control@spp3DiscardSteps+1)^3*(control@spp4DiscardSteps+1)^3*(control@spp5DiscardSteps+1)^2) 
    catchparms[ii,,,,,,,,,,,,,1,,,jj,5,1, ] <- rep(n[2:(k+1)] - n[1:k], each=eachterm5)
    catchparms[ii,,,,,,,,,,,,,,1,,jj,5,2, ] <- rep(o[2:(k+1)] - o[1:k], each=eachterm5)
    catchparms[ii,,,,,,,,,,,,,,,1,jj,5,3, ] <- rep(p[2:(k+1)] - p[1:k], each=eachterm5)
    
  }
}