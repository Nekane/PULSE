### class ######################################################################
createFLAccesors <- function(object, exclude=character(1)) {

	slots <- getSlots(class(object))[!names(getSlots(class(object)))%in%exclude]

	defined <- list()

	for (x in names(slots)) {
		# check method is defined already and signatures match
		eval(
		substitute(if(isGeneric(x) && names(formals(x)) != "object") {warning(paste("Accesor
			method for", x, "conflicts with a differently defined generic. Type", x,
			"for more information")); break}, list(x=x))
			)
		# create new generic and accesor method
		eval(
		substitute(if(!isGeneric(x)) setGeneric(x, function(object, ...) standardGeneric(x)),
		list(x=x))
		)
		eval(
		substitute(setMethod(x, signature(y), function(object) return(slot(object, x))), list(x=x,
			y=class(object)))
		)
		# create replacement method
		xr <- paste(x, "<-", sep="")
		eval(
		substitute(if(!isGeneric(x)) setGeneric(x,
			function(object, ..., value) standardGeneric(x)), list(x=xr))
		)
		eval(
		substitute(setMethod(x, signature(object=y, value=v), function(object, value)
			{slot(object, s) <- value; object}), list(x=xr, y=class(object), s=x,
			v=unname(slots[x])))
		)
		defined[[x]] <- c(x, xr, paste('alias{',x,',',class(object),'-method}', sep=''),
			paste('\alias{',xr,',',class(object),',',unname(slots[x]), '-method}', sep=''),
			paste('\alias{',x,'-methods}', sep=''),
			paste('\alias{"',xr, '"-methods}', sep='')
		)
	}
	return(defined)
}	# }}}




validDynState.control <- function(object){
	if (object@spp1LndQuota <= 0)
		return ("value of spp1LndQuota must be > 0")
        if (object@spp2LndQuota <= 0)
		  return("value of spp2LndQuota must be > 0")
	if (object@simNumber <= 0)
		return("value of Simumber must be > 0")
	if (object@increments <= 0)
		return("value of Increments must be > 0")
  # Everything is fine
	return(TRUE)
}

validDynStateInput <- function(object){
  if (!all(names(dimnames(object@catchMean))) == c("cat", "season", "option")) 
		return ("dimensions of arrays should be cat, season, option")
  # Everything is fine
	return(TRUE)
}

setClass("DynStateInput",
	representation(
        catchMean  ="array",
        catchSigma ="array"),
	prototype=prototype(
	  catchMean  =array(),
        catchSigma =array()),
      validity=validDynStateInput
)

invisible(createFLAccesors(new("DynStateInput")))


setClass("DynState.control",
	representation(
        spp1LndQuota          ="numeric",
        spp2LndQuota          ="numeric",
        spp1LndQuotaFine      ="numeric",
        spp2LndQuotaFine      ="numeric",
        simNumber             ="integer",
        sigma                 ="numeric",
        fuelUse               ="numeric",
        fuelPrice             ="numeric",
        gearMaintenance       ="numeric",
        landingCosts          ="numeric", 
        increments            ="integer",
        spp1Incs              ="numeric",
        spp2Incs              ="numeric",
        spp3Incs              ="numeric",
        spp4Incs              ="numeric",
        spp5Incs              ="numeric",
        addNoFishing	   ="logical",
        spp1DiscardSteps      ="numeric",
        spp2DiscardSteps      ="numeric",
        spp3DiscardSteps      ="numeric",
        spp4DiscardSteps      ="numeric",
        spp5DiscardSteps      ="numeric",
        choiceDist            ="integer",
        verbose               ="integer",
        numThreads            ="integer"),
	prototype=prototype(
        spp1LndQuota          =as.double(8),
        spp2LndQuota          =as.double(8),
        spp1LndQuotaFine      =as.double(1e8),
        spp2LndQuotaFine      =as.double(1e8),
        simNumber             =as.integer(10),
        sigma                 =as.double(30),
        fuelUse               =as.double(1),
        fuelPrice             =as.double(30),
        gearMaintenance       =as.double(30), 
        landingCosts          =as.double(0), 
        increments            =as.integer(8),
        spp1Incs              =as.double(1),
        spp2Incs              =as.double(1),
        spp3Incs              =as.double(1),
        spp4Incs              =as.double(1),
        spp5Incs              =as.double(1),
        addNoFishing          =TRUE,
        spp1DiscardSteps      =as.double(0),
        spp2DiscardSteps      =as.double(0),
        spp3DiscardSteps      =as.double(0),
        spp4DiscardSteps      =as.double(0),
        spp5DiscardSteps      =as.double(0),
        choiceDist            =as.integer(0),
        verbose               =as.integer(0),
        numThreads            =as.integer(4)),
	validity=validDynState.control
)

#invisible(createFLAccesors(new("DynState.control")))

setClass("Sim",
	representation(
        choice          = "array",
        spp1Landings    = "array",
        spp2Landings    = "array",
        spp3Landings    = "array",
        spp4Landings    = "array",
        spp5Landings    = "array",
        spp1LndHold     = "array",
        spp2LndHold     = "array",
        spp3LndHold     = "array",
        spp4LndHold     = "array",
        spp5LndHold     = "array",
        spp1Discards    = "array",
        spp2Discards    = "array",
        spp3Discards    = "array",
        spp4Discards    = "array",
        spp5Discards    = "array",
        spp1DisHold     = "array",
        spp2DisHold     = "array",
        spp3DisHold     = "array",
        spp4DisHold     = "array",
        spp5DisHold     = "array",
        spp1Rand        = "array",
        spp2Rand        = "array",
        spp3Rand        = "array",
        spp4Rand        = "array",
        spp5Rand        = "array",
        effort          = "array"),
	prototype=prototype(
        choice          = array(),
        spp1Landings    = array(),
        spp2Landings    = array(),
        spp3Landings    = array(),
        spp4Landings    = array(),
        spp5Landings    = array(),
        spp1LndHold     = array(),
        spp2LndHold     = array(),
        spp3LndHold     = array(),
        spp4LndHold     = array(),
        spp5LndHold     = array(),
        spp1Discards    = array(),
        spp2Discards    = array(),
        spp3Discards    = array(),
        spp4Discards    = array(),
        spp5Discards    = array(),
        spp1DisHold     = array(),
        spp2DisHold     = array(),
        spp3DisHold     = array(),
        spp4DisHold     = array(),
        spp5DisHold     = array(),
        spp1Rand        = array(),
        spp2Rand        = array(),
        spp3Rand        = array(),
        spp4Rand        = array(),
        spp5Rand        = array(),
        effort          = array())
)

invisible(createFLAccesors(new("Sim")))


setClass("DynState",
	representation(
        sim         ="Sim",
        control     ="DynState.control",
        choiceDist  ="array",
        spp1Price   ="array",
        spp2Price   ="array",
        spp3Price   ="array",
        spp4Price   ="array",
        spp5Price   ="array"),
	prototype=prototype(
        sim = new("Sim"),
        control     = new("DynState.control"),
        choiceDist  = array(),
        spp1Price   = array(),
        spp2Price   = array(),
        spp3Price   = array(),
        spp4Price   = array(),
        spp5Price   = array())
)

invisible(createFLAccesors(new("DynState")))

DynState.control <- function(spp1LndQuota=4000, spp2LndQuota=4000, spp1LndQuotaFine = 1E8, spp2LndQuotaFine=1E8, simNumber=30, fuelUse=1, sigma=0.1, fuelPrice=40, gearMaintenance=2, landingCosts=0, increments=25, addNoFishing=TRUE, spp1DiscardSteps=0, spp2DiscardSteps=0, spp3DiscardSteps=0, spp4DiscardSteps=0, spp5DiscardSteps=0, choiceDist=0, verbose=0,numThreads=4){
    res <- new("DynState.control", spp1LndQuota= as.double(spp1LndQuota), spp2LndQuota=as.double(spp2LndQuota),  spp1LndQuotaFine= as.double(spp1LndQuotaFine), spp2LndQuotaFine= as.double(spp2LndQuotaFine),  simNumber=as.integer(simNumber), sigma=as.double(sigma), fuelUse=as.double(fuelUse), fuelPrice=as.double(fuelPrice), gearMaintenance=as.double(gearMaintenance), landingCosts=as.double(landingCosts), increments=as.integer(increments), addNoFishing=addNoFishing, spp1DiscardSteps=as.double(spp1DiscardSteps), spp2DiscardSteps=as.double(spp2DiscardSteps), spp3DiscardSteps=as.double(spp3DiscardSteps), 
spp4DiscardSteps=as.double(spp4DiscardSteps), spp5DiscardSteps=as.double(spp5DiscardSteps), choiceDist=as.integer(choiceDist),verbose=as.integer(verbose), numThreads=as.integer(numThreads))
   return(res)
}

DynState <- function(inputSpp1, inputSpp2, inputSpp3, inputSpp4, inputSpp5, spp1Price, spp2Price, spp3Price, spp4Price, spp5Price, inputEffort, control){

  if (!inherits(inputSpp1,"DynStateInput"))   stop("inputSpp1 should be of class DynStateInput")
  if (!inherits(inputSpp2,"DynStateInput"))   stop("inputSpp2 should be of class DynStateInput")
  if (!inherits(inputSpp3,"DynStateInput"))   stop("inputSpp3 should be of class DynStateInput")
  if (!inherits(inputSpp4,"DynStateInput"))   stop("inputSpp4 should be of class DynStateInput")
  if (!inherits(inputSpp5,"DynStateInput"))   stop("inputSpp5 should be of class DynStateInput")
  if (control@increments > 65)                stop("number of increments should be =< 65")

  
  numInputSizes   <- dim(inputSpp1@catchMean)[1] #sizes
  numInputOptions <- dim(inputSpp1@catchMean)[3] #areas
  numInputSeasons <- dim(inputSpp1@catchMean)[2] #seasons
  
  
  # ----------------------------- Calculate length of vector of steps for separate spec (the *3 is to get a broad enough statistical distribution to discretize)

    #  control@spp1Incs 	      <- max(inputSpp1@catchMean[1,,] + (4.0*inputSpp1@catchSigma[1,,])) /(control@increments - 1 ) 
    #  control@spp2Incs        <- max(inputSpp2@catchMean[1,,] + (4.0*inputSpp2@catchSigma[1,,])) /(control@increments - 1 ) 
    #  control@spp3Incs        <- max(inputSpp3@catchMean[1,,] + (4.0*inputSpp3@catchSigma[1,,])) /(control@increments - 1 )
    #  control@spp4Incs        <- max(inputSpp4@catchMean[1,,] + (4.0*inputSpp4@catchSigma[1,,])) /(control@increments - 1 )
    #  control@spp5Incs        <- max(inputSpp5@catchMean[1,,] + (4.0*inputSpp5@catchSigma[1,,])) /(control@increments - 1 )
    #  for(si in 2:numInputSizes){
    #       control@spp1Incs <-max(control@spp1Incs ,max(inputSpp1@catchMean[si,,] + (4.0*inputSpp1@catchSigma[si,,])) /(control@increments - 1 )) 
    #       control@spp2Incs <-max(control@spp2Incs ,max(inputSpp2@catchMean[si,,] + (4.0*inputSpp2@catchSigma[si,,])) /(control@increments - 1 ))
    #       control@spp3Incs <-max(control@spp3Incs ,max(inputSpp3@catchMean[si,,] + (4.0*inputSpp3@catchSigma[si,,])) /(control@increments - 1 ))
    #       control@spp4Incs <-max(control@spp4Incs ,max(inputSpp4@catchMean[si,,] + (4.0*inputSpp4@catchSigma[si,,])) /(control@increments - 1 ))
    #       control@spp5Incs <-max(control@spp5Incs ,max(inputSpp5@catchMean[si,,] + (4.0*inputSpp5@catchSigma[si,,])) /(control@increments - 1 ))
    #  }                                                                        
 
 # Make th increments follow the probability quantile function of the negative binomial
 # Because sigma (theta) is the same everywhere we just take the first value of CatchSimga[1,1,1]  
  #for(si in 1:numInputSizes){
  control@spp1Incs 	      <- qnbinom(0.90,mu= max(inputSpp1@catchMean),size=inputSpp1@catchSigma[1,1,1]) /(control@increments - 1 ) 
  control@spp2Incs        <- qnbinom(0.90,mu= max(inputSpp2@catchMean),size=inputSpp2@catchSigma[1,1,1]) /(control@increments - 1 ) 
  control@spp3Incs        <- qnbinom(0.90,mu= max(inputSpp3@catchMean),size=inputSpp3@catchSigma[1,1,1]) /(control@increments - 1 )
  control@spp4Incs        <- qnbinom(0.90,mu= max(inputSpp4@catchMean),size=inputSpp4@catchSigma[1,1,1]) /(control@increments - 1 )
  control@spp5Incs        <- qnbinom(0.90,mu= max(inputSpp5@catchMean),size=inputSpp5@catchSigma[1,1,1]) /(control@increments - 1 )
  #}
  
  # ----------------------------- Check which vector is longest and use that number of steps
  k <- control@increments

  maxOptions <- numInputOptions*((control@spp1DiscardSteps+1)^numInputSizes*(control@spp2DiscardSteps+1)^numInputSizes*(control@spp3DiscardSteps+1)^numInputSizes*
			(control@spp4DiscardSteps+1)^numInputSizes*(control@spp5DiscardSteps+1)^numInputSizes)

  # ---------------------------- Build vector of equal length for three species, starting at very large neg number so that sum always =1
  spp1inc <- seq(0,k*control@spp1Incs,control@spp1Incs)
  spp2inc <- seq(0,k*control@spp2Incs,control@spp2Incs)
  spp3inc <- seq(0,k*control@spp3Incs,control@spp3Incs)
  spp4inc <- seq(0,k*control@spp4Incs,control@spp4Incs)
  spp5inc <- seq(0,k*control@spp5Incs,control@spp5Incs)
   
  spp1inc[1] <- spp2inc[1] <- spp3inc[1] <- spp4inc[1] <- spp5inc[1]<- -8000
  spp1inc[length(spp1inc)] <- spp2inc[length(spp2inc)] <- spp3inc[length(spp3inc)]<- spp4inc[length(spp4inc)] <- spp5inc[length(spp5inc)] <- 200000000
 
  # ----------------- The structure of catchparms in the model is patch, discardsteps, season, species,size, increment
  
  eachterm1 <- ((control@spp1DiscardSteps+1)^(numInputSizes-1)*(control@spp2DiscardSteps+1)^numInputSizes*(control@spp3DiscardSteps+1)^numInputSizes*
		(control@spp4DiscardSteps+1)^numInputSizes*(control@spp5DiscardSteps+1)^numInputSizes)
  eachterm2 <- ((control@spp1DiscardSteps+1)^numInputSizes*(control@spp2DiscardSteps+1)^(numInputSizes-1)*(control@spp3DiscardSteps+1)^numInputSizes*
		(control@spp4DiscardSteps+1)^numInputSizes*(control@spp5DiscardSteps+1)^numInputSizes)
  eachterm3 <-((control@spp1DiscardSteps+1)^numInputSizes*(control@spp2DiscardSteps+1)^numInputSizes*(control@spp3DiscardSteps+1)^(numInputSizes-1)*
		(control@spp4DiscardSteps+1)^numInputSizes*(control@spp5DiscardSteps+1)^numInputSizes)
  eachterm4 <- ((control@spp1DiscardSteps+1)^numInputSizes*(control@spp2DiscardSteps+1)^numInputSizes*(control@spp3DiscardSteps+1)^numInputSizes*
		(control@spp4DiscardSteps+1)^(numInputSizes-1)*(control@spp5DiscardSteps+1)^numInputSizes)
  eachterm5 <- ((control@spp1DiscardSteps+1)^numInputSizes*(control@spp2DiscardSteps+1)^numInputSizes*(control@spp3DiscardSteps+1)^numInputSizes*
		(control@spp4DiscardSteps+1)^numInputSizes*(control@spp5DiscardSteps+1)^(numInputSizes-1)) 

  steps <- max(control@spp1DiscardSteps,control@spp2DiscardSteps,control@spp3DiscardSteps,control@spp4DiscardSteps,control@spp5DiscardSteps)+1

  catchparms <- array(NA,dim=c(numInputOptions,steps,numInputSeasons,5,numInputSizes,k))

  Lparms <- array(NA,dim=c(numInputOptions,
                           rep((control@spp1DiscardSteps+1),numInputSizes),
                           rep((control@spp2DiscardSteps+1),numInputSizes),
                           rep((control@spp3DiscardSteps+1),numInputSizes), 
                           rep((control@spp4DiscardSteps+1),numInputSizes),
                           rep((control@spp5DiscardSteps+1),numInputSizes),numInputSeasons,5,numInputSizes,k))

  Dparms <- array(NA,dim=dim(Lparms))


  # ----------------- The structure of discparms and landparms in the model is patch, season, species,size, increment
  #-------------------------the same structure as catchparms
  
  discparms <- array(0,dim=dim(catchparms))
  landparms <- catchparms

  #normdist <- array(NA,dim=c(numInputOptions,numInputSeasons,5,numInputSizes,k+1))
  binomdist <- array(NA,dim=c(numInputOptions,numInputSeasons,5,numInputSizes,k+1))

  # ----------------- Funtion to put landparms/discparms in their positions in Lparms/Dparms
  # Need to change the args by hand if we change the number of sizes, because the positions of the specie, 
  # nsize change depending on the number of sizes dim[numInputOptions,sp1size1, sp1size2,sp1size3, sp2size1....,species,sizes, k]

	foo <- function(x, ii, jj, sp, si, dp, value) {

  		# list w/ all dims
  		args <- lapply(as.list(dim(x)), seq)
        posit <- length(args)
  		pos <- ((sp - 1) * length(args[[posit-1]])) + si + 1

  		# modify ii, jj, sp, si
  		args[[1]]  <- ii # numInputOptions
  		args[[posit-3]] <- jj # numInputSeasons
  		args[[posit-2]] <- sp # specie
  		args[[posit-1]] <- si # nsize

  		# modify dp
  		args[[pos]] <- dp

  		y <- do.call('[<-', c(list(x=x), args, list(value=value)))

  		return(y)
		} 


  for (sp in 1:5){ # number species
     for (si in 1:numInputSizes){ # number sizes
	for (ii in 1:numInputOptions){ # the area
	    for (jj in 1:numInputSeasons){ # the season

		# ----------------- Fill landparms/discparms when discarding is not allowed (landparms= catchparms)
		# ---------------------------------------------------(discparms= the prob of discarding 0 is 1)
        #		normdist [ii,jj,sp,si,] <- pnorm(get(paste0("spp",sp,"inc")),
        #						slot(get(paste0("inputSpp",sp)),"catchMean")[si,jj,ii],
        #						slot(get(paste0("inputSpp",sp)),"catchSigma")[si,jj,ii])
						
		binomdist [ii,jj,sp,si,] <- pnbinom(get(paste0("spp",sp,"inc")),
						mu=   slot(get(paste0("inputSpp",sp)),"catchMean")[si,jj,ii],
						size= slot(get(paste0("inputSpp",sp)),"catchSigma")[si,jj,ii])

		#catchparms[ii,1,jj,sp,si, ] <- normdist [ii,jj,sp,si,2:(k+1)]- normdist [ii,jj,sp,si,1:k]  # Those are the total catches probabilities
		catchparms[ii,1,jj,sp,si, ] <- binomdist [ii,jj,sp,si,2:(k+1)]- binomdist [ii,jj,sp,si,1:k]  # Those are the total catches probabilities

		landparms [ii,1,jj,sp,si, ] <- catchparms[ii,1,jj,sp,si, ]

		discparms [ii,1,jj,sp,si, ] <- c(1,rep(0,k-1))

		Lparms <- foo(x=Lparms,sp=sp,si=si,dp=1,ii=ii,jj=jj, value=rep(landparms[ii,1,jj,sp,si, ], each= get(paste0("eachterm",sp))))
		Dparms <- foo(x=Dparms,sp=sp,si=si,dp=1,ii=ii,jj=jj, value=rep(discparms[ii,1,jj,sp,si, ], each= get(paste0("eachterm",sp))))

		# ----------------- Fill landparms/discparms with the different discard options

		for (dp in 0:slot(get("control"),paste0("spp",sp,"DiscardSteps"))){ #for the different sp discardsteps

			# ----------------- Fill landparms/discparms with the different discard options
      			if (dp >0 && dp != slot(get("control"),paste0("spp",sp,"DiscardSteps"))) {

      			spltvec <- floor(k*(1- dp/(slot(get("control"),paste0("spp",sp,"DiscardSteps")))))

      			landparms[ii,dp+1,jj,sp,si, ] <- c(landparms[ii,1,jj,sp,si,1:spltvec],sum(landparms[ii,1,jj,sp,si,(spltvec+1):k]), rep(0,k-spltvec-1)) 
            
      			discparms[ii,dp+1,jj,sp,si, ] <- c(sum(landparms[ii,1,jj,sp,si,1:(spltvec+1)]),landparms[ii,1,jj,sp,si,(spltvec+2):k], rep(0,spltvec)) 

		Lparms <- foo(x=Lparms,sp=sp,si=si,dp=dp+1,ii=ii,jj=jj, value=rep(landparms[ii,dp+1,jj,sp,si, ], each= get(paste0("eachterm",sp))))
		Dparms <- foo(x=Dparms,sp=sp,si=si,dp=dp+1,ii=ii,jj=jj, value=rep(discparms[ii,dp+1,jj,sp,si, ], each= get(paste0("eachterm",sp))))
      			}

			# ----------------- Everything is discarded
      			if (dp > 0 && dp == slot(get("control"),paste0("spp",sp,"DiscardSteps"))){ 

      			landparms[ii,dp+1,jj,sp,si, ] <- c(1,rep(0,k-1))

      			discparms[ii,dp+1,jj,sp,si, ] <- landparms [ii,1,jj,sp,si, ] 

		Lparms <- foo(x=Lparms,sp=sp,si=si,dp=dp+1,ii=ii,jj=jj, value=rep(landparms[ii,dp+1,jj,sp,si, ], each= get(paste0("eachterm",sp))))
		Dparms <- foo(x=Dparms,sp=sp,si=si,dp=dp+1,ii=ii,jj=jj, value=rep(discparms[ii,dp+1,jj,sp,si, ], each= get(paste0("eachterm",sp))))
      			}

    		} 

	     }
	 }
      }
  }
 

  # ---------------------------- Generate effort array from effort vector ----------------------------- 
  #   inputEffort <- array(inputEffort, dim=c(length(inputEffort),1) )
  #   inputEffort <- rep(inputEffort,((control@spp1DiscardSteps+1)^numInputSizes*
  # (control@spp2DiscardSteps+1)^numInputSizes*(control@spp3DiscardSteps+1)^numInputSizes*
  # (control@spp4DiscardSteps+1)^numInputSizes*(control@spp5DiscardSteps+1)^numInputSizes))

  li <- ((control@spp1DiscardSteps+1)^numInputSizes*(control@spp2DiscardSteps+1)^numInputSizes*
        (control@spp3DiscardSteps+1)^numInputSizes*(control@spp4DiscardSteps+1)^numInputSizes*
        (control@spp5DiscardSteps+1)^numInputSizes)

  inputEffort <- do.call(rbind, lapply(1:li, function (x) inputEffort))
  rownames(inputEffort)<- NULL

  # ---------------------------- Generate lndparms and dscparms ----------------------------------------------------
  lndparms<- array(0,dim=dim(Lparms))
  dim(lndparms) <- c(maxOptions,numInputSeasons,5,numInputSizes,k) # collapse first three dims so we get each derivative area in first dim
  lndparms <- array(NA,dim=c(maxOptions + control@addNoFishing,numInputSeasons,5,numInputSizes,k))
  lndparms[1:maxOptions,,,,] <- Lparms
  dscparms<- array(0,dim=dim(Dparms))  
  dim(dscparms) <- c(maxOptions,numInputSeasons,5,numInputSizes,k) # collapse first three dims so we get each derivative area in first dim
  dscparms <- array(NA,dim=c(maxOptions + control@addNoFishing,numInputSeasons,5,numInputSizes,k))
  dscparms[1:maxOptions,,,,] <- Dparms

  if (slot(control, "addNoFishing")){
      # add zero catchrates in highest patch which is not fishing
      for (j in 1:numInputSeasons)  { # the season
        for (s in (1:5))  {      # the number of species
          for (si in (1:numInputSizes))  {      # the number of sizeclasses
            lndparms[maxOptions + 1,j,s,si, 1] <- 1
            lndparms[maxOptions + 1,j,s,si,-1] <- 0
            dscparms[maxOptions + 1,j,s,si, 1] <- 1
            dscparms[maxOptions + 1,j,s,si,-1] <- 0
          }  
        }
      }  
      # inputEffort <- c(inputEffort,0)
      inputEffort <- rbind(inputEffort,0)
  }


  # ---------------------------- Generate prices-costs for individual species that are used in c++ part of model 
  spp1PminC <- spp1Price - control@landingCosts
  spp2PminC <- spp2Price - control@landingCosts
  spp3PminC <- spp3Price - control@landingCosts
  spp4PminC <- spp4Price - control@landingCosts
  spp5PminC <- spp5Price - control@landingCosts

  # ---------------------------- merge prices into a single object that goes into cpp 
  sppPminC <- array(cbind( spp1PminC, spp2PminC, spp3PminC, spp4PminC, spp5PminC),dim=c( dim(spp1Price),5))

  sppPminC
    
     
  # ---------------------------- Test if patchnumber in costarray is equal to patchnumber catchparms --
  # if (length(inputEffort) != dim(lndparms)[1]) stop ("Effortarray contains different number patches than landarray")
  # if (length(inputEffort) != dim(dscparms)[1]) stop ("Effortarray contains different number patches than discarray")
  if (dim(inputEffort)[1] != dim(lndparms)[1]) stop ("Effortarray contains different number patches than landarray")
  if (dim(inputEffort)[1] != dim(dscparms)[1]) stop ("Effortarray contains different number patches than discarray")

  # ----------------------------- Check if MaxOptionss+1 > 320000 because that is not possible if choice in cpp part is unsigned int THIS IS NOT TRUE ANYMORE, NOW IS UNISGNED INT, BUT NOT SURE ABOUT RANGE
  if (dim(lndparms)[1]>320000) stop ("More than 320000 choices in landarray, which is not allowed")    
  if (dim(dscparms)[1]>320000) stop ("More than 320000 choices in discarray, which is not allowed") 
  
                                                                         
  # ---------------------------- Call c++ code --------------------------------------------------------
  res <- .Call("DynStateF",lndparms, dscparms,inputEffort, sppPminC, control)

 
  # ---------------------------- Calculate derivative choice ------------------------------------------
  res@sim@choice <-
   ifelse(res@sim@choice == (maxOptions + 1), (maxOptions + 1),ifelse((res@sim@choice %% numInputOptions)==0,numInputOptions,(res@sim@choice %% numInputOptions)))

 
  res@sim@choice                 <- array(dimnames(inputSpp1@catchMean)$option[res@sim@choice],dim=c(1,control@simNumber,numInputSeasons))
  
  dimnames(res@sim@choice)       <- list(size="all",simulation=(1:control@simNumber),season=(dimnames(inputSpp1@catchMean)$season))
  
  dimnames(res@sim@spp1Landings) <- dimnames(res@sim@spp2Landings) <- dimnames(res@sim@spp3Landings) <- dimnames(res@sim@spp4Landings) <- dimnames(res@sim@spp5Landings)<- list(cat=1:numInputSizes,simulation=(1:control@simNumber),season=(dimnames(inputSpp1@catchMean)$season))
  
  dimnames(res@sim@spp1LndHold)  <- dimnames(res@sim@spp2LndHold)  <- dimnames(res@sim@spp3LndHold)  <- dimnames(res@sim@spp4LndHold)  <- dimnames(res@sim@spp5LndHold) <- list(cat=1:numInputSizes,simulation=(1:control@simNumber),season=(dimnames(inputSpp1@catchMean)$season))
  
  dimnames(res@sim@spp1Discards) <- dimnames(res@sim@spp2Discards) <- dimnames(res@sim@spp3Discards) <- dimnames(res@sim@spp4Discards) <- dimnames(res@sim@spp5Discards)<- list(cat=1:numInputSizes,simulation=(1:control@simNumber),season=(dimnames(inputSpp1@catchMean)$season))
  
  dimnames(res@sim@spp1DisHold)  <- dimnames(res@sim@spp2DisHold)  <- dimnames(res@sim@spp3DisHold)  <- dimnames(res@sim@spp4DisHold)  <- dimnames(res@sim@spp5DisHold) <- list(cat=1:numInputSizes,simulation=(1:control@simNumber),season=(dimnames(inputSpp1@catchMean)$season))

  res@control <- control
  res@spp1Price <-  sppPminC [,,1] #spp1Price 
  res@spp2Price <-  sppPminC [,,2] #spp2Price
  res@spp3Price <-  sppPminC [,,3] #spp3Price
  res@spp4Price <-  sppPminC [,,4] #spp4Price
  res@spp5Price <-  sppPminC [,,5] #spp5Price

  # --------------------------- Add catchparms to results if required --------------------------------
  if (control@choiceDist == 1) res@choiceDist <- lndparms
  if (control@choiceDist == 2) res@choiceDist <- dscparms
  return(res)
}
     
     
plotChoice <- function(inputSpp,dynState){

        choice(sim(dynState))[is.na(choice(sim(dynState)))] <- "no"
        options <-  c(unlist(unique(dimnames(catchMean(inputSpp))[3])),"no")
        seasons <- dimnames(catchMean(inputSpp))[[2]]
        res <- matrix(0,nrow=length(seasons),ncol=length(options), dimnames=list(seasons=seasons,options=options))
        
        for (ss in 1:dim(choice(sim(dynState)))[3]){
                restab <- table(choice(sim(dynState))[,,ss])
                for (ii in 1:length(restab)){
                res[ss,names(restab[ii])] <- restab[ii] 
                } 
        } 
        
        print(image(x=1:length(seasons), y=1:length(options),res, axes=F, col=gray((32:0)/32),xlab="Seasons",ylab="Options"))
        axis(1, at = seq(1, length(seasons), by = 1), labels=seasons)
        axis(2, at = seq(1, length(options), by = 1), labels=options)
 }

grossRev <- function(dynstate){
   res <- apply(sweep(spp1Landings(sim(dynstate)),c(1,3), (dynstate@spp1Price - control(dynstate)@landingCosts), "*") +
                sweep(spp2Landings(sim(dynstate)),c(1,3), (dynstate@spp2Price - control(dynstate)@landingCosts), "*") + 
                sweep(spp3Landings(sim(dynstate)),c(1,3), (dynstate@spp3Price - control(dynstate)@landingCosts), "*") +
                sweep(spp4Landings(sim(dynstate)),c(1,3), (dynstate@spp4Price - control(dynstate)@landingCosts), "*") + 
                sweep(spp5Landings(sim(dynstate)),c(1,3), (dynstate@spp5Price - control(dynstate)@landingCosts), "*"), 2, sum) 

    return(res)
 }
 
annualFine <- function(dynstate){
    res <- ifelse(apply(spp1Landings(sim(dynstate)),2,sum) > control(dynstate)@spp1LndQuota,
                 (apply(spp1Landings(sim(dynstate)),2,sum) - control(dynstate)@spp1LndQuota) *control(dynstate)@spp1LndQuotaFine ,0)+
           ifelse(apply(spp2Landings(sim(dynstate)),2,sum) > control(dynstate)@spp2LndQuota,
                 (apply(spp2Landings(sim(dynstate)),2,sum) - control(dynstate)@spp2LndQuota) *control(dynstate)@spp2LndQuotaFine ,0)
 
    return(res)
 }


netRev<- function (dynstate){
        res <- grossRev(dynstate) - 
                apply(effort(sim(dynstate)) * control(dynstate)@fuelUse * control(dynstate)@fuelPrice, 2, sum)- #Fuel cost
                (grossRev(dynstate) * 0)    -                                                                   #Crew cost
                apply(effort(sim(dynstate)) * control(dynstate)@gearMaintenance, 2, sum)-                       #Gear mantainance cost
                annualFine(dynstate)
      return(res)
}
 
