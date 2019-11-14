#' Run the fish density simulation
#'
#' Executes the simulation for a given region and species.
#'
#'
#' @examples
#' ## TBD
#'
#' @param year Either a single year or 2 element vector with maximum and minimum year to retrieve.
#' @param survey String for intended survey (EBS_SHELF, EBS_SLOPE, NBS_SHELF, GOA, AI).
#'
#' @return The output is a list of years and densities. The output must be re-joined to the Predict_data
#' object for spatial analysis.
#'   \item{YEARS_USE}{Vector of available years}
#'   \item{SDM_result}{Matrix of densities per year, rows for a given year correspond to Predict_data rows.}
#' @export

	### The function belows returns a list containing the simulated species distribution (the size of list depends on the number of mcmc iteration if using the MCMC option, OR a list of 1 element containing the mean predicted distribution based on the frequentist approach) within a specific survey region boundary (defined using shapefile). The output object name is: SDM_result and the unit is kg/km2

	runAnalysis <- function(
		Mcmc_prediction = TRUE,   			# TRUE if prediction is based on MCMC, FALSE if based on predicted mean. Default = TRUE.
		Abundance_or_catch = "Catch",   	# choose among 'Abundance' or 'Catch'. This determines the type of final output: either i) the predictive distribution (option 'Abundance'), or ii) the sampled catch (option 'Catch')
		Mcmc_sample = 1, 					# only accounted for if Mcmc_prediction = TRUE
		Mcmc_seed=200, 						# seed for producing mcmc sample. For reproducibility.
		Which_species,	# Which species to focus on? Using species code for now as I do not have species name
		Which_region, 	# Which regions to focus the analysis? It can be any combination of: "EBS", "SLP", "AI", "GOA", "NBS" (NBS is not available yet)
		Years = NULL,						# If Years=NULL, the code automatically defines it as annual from the min and max range from the data. Default = NULL.
		Env_covariate = c('SURFACE_TEMPERATURE','GEAR_TEMPERATURE','BOTTOM_DEPTH'), # choose among 'SURFACE_TEMPERATURE','GEAR_TEMPERATURE','BOTTOM_DEPTH', only (the three that are available in the dataset). P.S: When more than 2 regions are combined, 'SURVEY_NAME' effect is automatically added to the model to account for regional differences
		Transform_covariate = c("ID", "ID", "LOG"),		# variable transformation method for each environmental covariate. Choice between "ID" (identity) and "LOG" (log) only at the moment. NEEDS to be the same length as 'Env_covariate' as it applies the transformation to each of the above variable
		Include_year_effect = FALSE,  		# Adding the 'Year' effect as factor. P.S: Not recomended. It is VERY difficult to estimate because counfounded with background field. Default =FALSE
		Include_2_order = TRUE,				# Adding the polynomial effect for the different covariates (recommended). Default = TRUE
		Standardize_variable = FALSE,		# Standardize the variable or not. Not recommended as it "complicates" interpretation. Default = FALSE.
		Spacetime = "AR",					# The type of spatio-temporal model to run. Choice between "AR" (spatial field is 1st order AR), "IID" (spatial field is sampled independently every year from the same underlying spatial covariance structure), "NO" (no time varying spatial correlation structure), "NADA" (no spatial nor temporal correlation structure)
		YEARS_pred = 1982:2016,				# Prediction year (DO NOT CHANGE). These are the only extent possible at the moment
		Overwrite = FALSE  					# If the model was already run, do you want to overwrite the results. Default = FALSE to not erase past run.
	)
	{

		### Basic condition checking for the input variables
  	  if (identical(YEARS_pred, 1982:2016) == FALSE)
  	  {
  	    print("WARNING: YEARS_pred was changed back to 1982:2016. \n These are the only extent possible at the moment")
  	    YEARS_pred = 1982:2016
  	  }
  	  if (!Which_region %in% c("EBS", "SLP", "AI", "GOA", "NBS"))
  	  {
  	    print("WARNING: WRONG region name. Please modify")
  	    break()
  	  }
  	  if (!all(Env_covariate %in% c('SURFACE_TEMPERATURE','GEAR_TEMPERATURE','BOTTOM_DEPTH')))
  	  {
  	    print("WARNING: 'Env_covariate' can ONLY be either \n  'SURFACE_TEMPERATURE','GEAR_TEMPERATURE','BOTTOM_DEPTH'. \n Please check for spelling error")
  	    break()
  	  }
  	  if (length(Transform_covariate) != length(Env_covariate))
  	  {
  	    print("WARNING: The length of 'Transform_covariate' should be \n EQUAL to the 'Env_covariate'. Please modify")
  	    break()
  	  }
  	  if (!all(Transform_covariate %in% c("ID", "LOG")))
  	  {
  	    print("WARNING: 'Transform_covariate' can ONLY be either \n 'ID' or 'LOG'. Please check for spelling error")
  	    break()
  	  }
  	  if (!is.logical(Mcmc_prediction) | !is.logical(Include_year_effect) | !is.logical(Include_2_order) | !is.logical(Standardize_variable) | !is.logical(Overwrite))
  	  {
  	    print("WARNING: 'Mcmc_prediction', 'Include_year_effect', 'Include_2_order', \n 'Standardize_variable', and 'Overwrite' can ONLY be boolean. Please check for error")
  	    break()
  	  }
  	  if (Standardize_variable == TRUE)
  	  {
  	    print("WARNING: Standardize_variable is not yet fully implemented \n Only FALSE available for now")
  	    break()
  	  }
  	  if (!all(Spacetime %in% c("AR", "IID", "NO", "NADA")))
  	  {
  	    print("WARNING: 'Spacetime' can ONLY be either \n 'AR', 'IID', 'NO', or 'NADA'. Please check for spelling error")
  	    break()
  	  }



		Data_work <- subset(All_data, subset=c(SPECIES_CODE == Which_species) & SURVEY_NAME %in% Which_region)

		Env_covariate_base <- Env_covariate
		select_col <- match(c('YEAR','LONGITUDE','LATITUDE','CPUE','SURVEY_NAME',Env_covariate_base), colnames(Data_work))
		Data_work <- Data_work[,select_col]
		Data_work <- Data_work[complete.cases(Data_work),] # remove observations with missing data

		if (any(Transform_covariate == "LOG"))
		{
		  which_log <- which(Transform_covariate == "LOG")
		  colnames(Data_work)[match(Env_covariate_base[which_log],colnames(Data_work))] <- paste0("LOG_",Env_covariate_base[which_log])
		  Env_covariate[which_log] <- paste0("LOG_",Env_covariate_base[which_log])
		  Data_work[,match(Env_covariate[which_log],colnames(Data_work))] <- log(Data_work[,match(Env_covariate[which_log],colnames(Data_work))])
		}

		### If polynomial effects need to be included
		if (Include_2_order == TRUE)
		{
		  select_var <- match(Env_covariate, colnames(Data_work))
		  Polynomial_dat <- Data_work[,select_var]^2
		  colnames(Polynomial_dat) <- paste0(colnames(Polynomial_dat), 2)
		  Data_work <- cbind(Data_work, Polynomial_dat)
		  Env_covariate = c(Env_covariate, paste0(Env_covariate, "2"))
		}

		YEARS = seq(range(Data_work$YEAR)[1], range(Data_work$YEAR)[2])
		if(!is.null(Years))
		{
		  if (all(Years %in% YEARS)) YEARS = Years
		  if (!all(Years %in% YEARS))
		  {
		    print("WARNING: 'Years' contains values outside the range of the data. \n Please change to 'NULL' or check for error")
		    break()
		  }
		}

		### New variable for presence/absence
		Data_work$Pres <- with(Data_work, ifelse(CPUE>0, 1, 0))


		######################################################
		############### Run the binomial model ###############
		######################################################

		NAME <- paste("Bin",paste(Which_region, collapse="_"),paste(Which_species, collapse="_"),sep="_")

		if (all(c("AI") %in% Which_region)) Data_work$LONGITUDE[which(Data_work$LONGITUDE>0)] <- Data_work$LONGITUDE[which(Data_work$LONGITUDE>0)]-360
		Coord <- cbind(Data_work$LONGITUDE, Data_work$LATITUDE)					# the coordinates
		bnd = inla.nonconvex.hull(Coord, convex=200000)		# I need to get the boundary limit of the region on interest
		if (all(sort(c("EBS", "SLP", "AI", "GOA")) == sort(Which_region))) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(sort(c("AI", "GOA")) == sort(Which_region))) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(sort(c("EBS", "GOA"))== sort(Which_region))) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(sort(c("EBS", "SLP"))== sort(Which_region))) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(c("EBS") == Which_region)) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(c("SLP") == Which_region)) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(1000000,2000000), cutoff = 75000, offset=c(50000,100000))
		if (all(c("GOA") == Which_region)) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(c("AI") == Which_region)) mesh = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(50000,100000))


		spde=inla.spde2.matern(mesh, alpha=3/2)	# exponential decaying spatial correlation

		### The design matrix for the observation model
		# Categorical variables
		if (Include_year_effect == TRUE)
		{
		  N.bin 	<- 	dim(Data_work)[1]
		  temp		<- matrix(0,N.bin,length(YEARS)-1)
		  if(length(YEARS)>1){
		    for(i in 2:(length(YEARS))){
		      these = which(Data_work$YEAR==YEARS[i])
		      temp[these,(i-1)] <- 1
		    }
		    colnames(temp)	<-	paste("X.",YEARS[2:length(YEARS)],sep="")
		  }
		  TEMP <- cbind(temp)
		} else { TEMP <- NULL }

		if (length(Which_region)>1)
		{
		  N.bin 	<- 	length(Which_region)-1
		  temp		<- matrix(0,nrow(Data_work),N.bin)
		  for(i in 2:length(Which_region)){
		    these = which(Data_work$SURVEY_NAME==Which_region[i])
		    temp[these,(i-1)] <- 1
		  }
		  colnames(temp)	<-	paste("SURVEY_NAME.",Which_region[2:length(Which_region)],sep="")
		  TEMP <- cbind(TEMP, temp)
		}

		# Continuous variables
		Mean_Bin <- apply(Data_work[,match(Env_covariate, names(Data_work))], 2, mean, na.rm=T)
		Sd_Bin <- apply(Data_work[,match(Env_covariate, names(Data_work))], 2, sd, na.rm=T)

		A<-	NULL
		if(Standardize_variable==TRUE)
		{
		  for( i in 1:length(Env_covariate)){
		    A.new <- Data_work[,Env_covariate[i]]
		    A.new <- (A.new - Mean_Bin[i])/(2*Sd_Bin[i])
		    A	<-	cbind(A,A.new)
		  }
		}
		if(Standardize_variable==FALSE)
		{
		  for( i in 1:length(Env_covariate)){
		    A.new <- Data_work[,Env_covariate[i]]
		    A	<-	cbind(A,A.new)
		  }
		}
		A<-	data.frame(A)
		colnames(A)	<-	(Env_covariate)
		B <-	A

		# Combine both
		if (is.null(TEMP)) X.1 <- B
		if (!is.null(TEMP)) X.1 <- data.frame(cbind(TEMP,B))

		# Create the final output design matrix
		X.1.all <- as.data.frame(X.1)

		Covar.names <- colnames(X.1.all)										# the list of Env_covariate for the estimation model
		XX.list <- as.list(X.1.all)

		### Creating the design matrix for prediction (by year)
		DESIGN.pred	<-	list()

		if (any(Transform_covariate == "LOG"))
		{
		  which_log <- which(Transform_covariate == "LOG")
		  to_change <- match(Env_covariate_base[which_log],colnames(Predict_data))
		  colnames(Predict_data)[to_change] <- paste0("LOG_",Env_covariate_base[which_log])
		  # dealing with depth<=0 (i.e. land masses and erroneous bathymetry). Here we replace these entries with the minimum haul depth for the region.
		  Predict_data[!is.na(Predict_data[,to_change]) & Predict_data[,to_change]<0,to_change] <- min(Data_work[,to_change])
		  Predict_data[!is.na(Predict_data[,to_change]),to_change] <- log(Predict_data[!is.na(Predict_data[,to_change]),to_change])
		}

		for(j in 1:length(YEARS_pred))
		{
		  selected <- grep(YEARS_pred[j], colnames(Predict_data))
		  new_pred_dat <- Predict_data[,c(1,2,3,selected)]
		  if (length(Which_region)>1)
		  {
		    temp <- matrix(0,nrow(new_pred_dat),(length(Which_region)-1))
		    for (i in 2:length(Which_region))
		    {
		      selected <- grep(Which_region[i], colnames(Predict_data))
		      temp[,i-1] <- Predict_data[,selected]
		    }
		    colnames(temp) <- paste("SURVEY_NAME.",Which_region[2:length(Which_region)],sep="")
		    new_pred_dat <- data.frame(new_pred_dat, temp)
		  }

		  if (Include_2_order == TRUE)
		  {
		    select_var <- sort(unlist(sapply(1:length(Env_covariate), function(x) grep(Env_covariate[x], colnames(new_pred_dat)))))
		    Polynomial_dat <- new_pred_dat[,select_var]^2
		    colnames(Polynomial_dat) <- paste0(colnames(Polynomial_dat), 2)
		    new_pred_dat <- cbind(new_pred_dat, Polynomial_dat)
		  } else {
		    select_var <- sort(unlist(sapply(1:length(Env_covariate), function(x) grep(Env_covariate[x], colnames(new_pred_dat)))))
		  }

		  colnames(new_pred_dat) <- gsub(YEARS_pred[j], "", colnames(new_pred_dat))

		  # For continuous covariates
		  A<-	NULL
		  if(Standardize_variable==TRUE)
		  {
		    for( i in 1:length(Env_covariate)){
		      A.new <- new_pred_dat[,Env_covariate[i]]
		      A.new <- (A.new - Mean_Bin[i])/(2*Sd_Bin[i])
		      A	<-	cbind(A,A.new)
		    }
		  }
		  if(Standardize_variable==FALSE)
		  {
		    for( i in 1: length(Env_covariate)){
		      A.new <- new_pred_dat[,Env_covariate[i]]
		      A	<-	cbind(A,A.new)
		    }
		  }
		  A			<-	data.frame(A)
		  colnames(A)	<-	(Env_covariate)

		  B	<-	A

		  # For the year effect
		  if (j>1)
		  {
		    Intercept <- 1
		    if (Include_year_effect == TRUE)
		    {
		      X.1.pred		<-	data.frame(cbind(Intercept,1,B))
		      colnames(X.1.pred) <- c("Intercept", paste("X",YEARS_pred[j],sep="."), colnames(B))
		    } else {
		      X.1.pred		<-	data.frame(cbind(Intercept, B))
		      colnames(X.1.pred) <- c("Intercept", colnames(B))
		    }
		  }
		  if (j==1)
		  {
		    Intercept <- 1
		    X.1.pred		<-	data.frame(cbind(Intercept, B))
		    colnames(X.1.pred) <- c("Intercept", colnames(B))
		  }

		  DESIGN.pred[[j]] <- X.1.pred

		}


		### Setting the spatio-temporal model

		###### Specifies how to allocate the ID for each mesh
		if (Spacetime=="NO") field.indices = inla.spde.make.index("field", n.spde=mesh$n)
		if (Spacetime=="AR") field.indices = inla.spde.make.index("field", n.spde=mesh$n, n.group=length(YEARS))
		if (Spacetime=="IID") field.indices = inla.spde.make.index("field", n.spde=mesh$n, n.repl=length(YEARS))

		field.group <- rep(1:length(YEARS), as.numeric(table(factor(Data_work$YEAR, levels=YEARS))))
		field.repl <- rep(1:length(YEARS), as.numeric(table(factor(Data_work$YEAR, levels=YEARS))))
		if (Spacetime=="AR") A <- inla.spde.make.A(mesh, loc=Coord, group=field.group)
		if (Spacetime=="IID") A <- inla.spde.make.A(mesh, loc=Coord, repl=field.repl)
		if (Spacetime=="NO") A <- inla.spde.make.A(mesh, loc=Coord)
		if(Spacetime != "NADA")
		{
		  ## creating the list of object to be used to define covariates that I want to keep for the estimation model
		  effect.list <- list()
		  effect.list[[1]] <- c(field.indices, list(Intercept=1))
		  for (i in 1:ncol(X.1)) effect.list[[i+1]] <- XX.list[[i]]
		  names(effect.list) <- c("1", Covar.names)

		  ## creating the list of predictors (for the estimation model) that match the one specified above
		  A.list <- list()
		  A.list[[1]] <- A
		  for (i in 1:ncol(X.1)) A.list[[i+1]] <- 1

		  ## model frame specification
		  stack <- inla.stack(data=list(resp=(Data_work$Pres)), A=A.list, effects=effect.list, tag='est')
		}

		## model formula specification --> this is the main part to test different models. For example, assuming different spatial structure every year, or serial autocorrelation
		formula0 = as.formula(paste0("resp ~ ", paste(Covar.names, collapse="+")))
		formula = as.formula(paste0("resp ~ -1 + Intercept +",  paste(Covar.names, collapse="+"), "+ f(field, model=spde)"))
		formula2 = as.formula(paste0("resp ~ -1 + Intercept +",  paste(Covar.names, collapse="+"), "+ f(field, model=spde, group=field.group, control.group=list(model='ar1'))"))
		formula3 = as.formula(paste0("resp ~ -1 + Intercept +",  paste(Covar.names, collapse="+"), "+ f(field, model=spde, replicate=field.repl, control.group=list(model='exchangeable'))"))

		# Run model
		if (Overwrite == TRUE | file.exists(NAME)==FALSE)
		{
		  unlink(x=paste0(getwd(), "/", NAME), recursive = TRUE, force = TRUE)
		  if (Spacetime=="NADA") result = inla(formula0, family="binomial", data=cbind(resp=dat$y, X.1), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), keep=T, working.directory= NAME)
		  if (Spacetime=="NO") result = inla(formula, family="binomial", data=inla.stack.data(stack, spde=spde), verbose=TRUE, control.family = list(link = "logit"), control.predictor=list(A=inla.stack.A(stack), compute=TRUE), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), keep=T, working.directory= NAME)
		  if (Spacetime=="AR") result = inla(formula2, family="binomial", data=inla.stack.data(stack, spde=spde), verbose=TRUE, control.family = list(link = "logit"), control.predictor=list(A=inla.stack.A(stack), compute=TRUE), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), keep=T, working.directory= NAME)
		  if (Spacetime=="IID") result = inla(formula3, family="binomial", data=inla.stack.data(stack, spde=spde), verbose=TRUE, control.predictor=list(A=inla.stack.A(stack), compute=TRUE), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE, openmp.strategy="large"), keep=FALSE, num.threads=2, working.directory= NAME)
		}

		# if you have already run the model
		if (Overwrite == FALSE)
		{
		  result <- inla.collect.results(paste(NAME, "/results.files", sep=""))
		}

		summary(result)

		spde.result = inla.spde.result(result, "field", spde, do.transf=TRUE)

		# Calculate the estimated expected surface (for the expected value of parameters, replace by NA all values that are outside the range of the Predict_data)
		mesh_proj = inla.mesh.projector(mesh, xlim=range(Predict_data$LONG), ylim=range(Predict_data$LAT), dims=c(length(unique(Predict_data$LONG)),length(unique(Predict_data$LAT))))

		link_y_hat_bin <- list()
		LATENT_FIELD = matrix(result$summary.ran$field$mean, mesh$n, length(YEARS))

		YEARS_USE = intersect(YEARS, YEARS_pred)
		for (yr in 1:length(YEARS_USE))
		{
		  yr_pred = which(YEARS_pred == YEARS_USE[yr])
		  yr_mod = which(YEARS == YEARS_USE[yr])

		  Data.projection <- matrix(NA, nrow=length(unique(Predict_data$LONG))*length(unique(Predict_data$LAT)), ncol=(ncol(DESIGN.pred[[yr_pred]])))
		  Data.projection <- cbind(Data.projection, expand.grid(LONG=seq(range(Predict_data$LONG)[1],range(Predict_data$LONG)[2],length.out=length(unique(Predict_data$LONG))), LAT=seq(range(Predict_data$LAT)[1],range(Predict_data$LAT)[2],length.out=length(unique(Predict_data$LAT)))))
		  Data.projection[, 1:(ncol(Data.projection)-2)] <- DESIGN.pred[[yr_pred]]				# replace the above values in the Projection matrix data

		  Factors <- which(rownames(result$summary.fix)%in%colnames(DESIGN.pred[[yr_pred]]))

		  Data.projection <- Data.projection[,1:(ncol(Data.projection)-2)]

		  link_y_hat_bin[[yr]] = inla.mesh.project(mesh_proj,LATENT_FIELD[,yr_mod]) + matrix(as.matrix(Data.projection)%*%as.numeric(result$summary.fix[Factors,'mean']),nrow=length(unique(Predict_data$LONG)), ncol=length(unique(Predict_data$LAT)))
		}

		### Do mcmc sampling
		if (Mcmc_prediction ==TRUE)
		{
		  set.seed(Mcmc_seed)
		  s1000 <- inla.posterior.sample(Mcmc_sample, result)

		  mesh_proj = inla.mesh.projector(mesh, xlim=range(Predict_data$LONG), ylim=range(Predict_data$LAT), dims=c(length(unique(Predict_data$LONG)),length(unique(Predict_data$LAT))))

		  predictions <- function(x)
		  {
		    which_field_mesh <- grep("field", names(x$latent[,1]))
		    field_mesh = matrix(as.numeric(x$latent[,1])[which_field_mesh], mesh$n, length(YEARS))

		    link_y_hat <- matrix(NA, nrow=length(unique(Predict_data$LONG))*length(unique(Predict_data$LAT)), ncol=length(YEARS_USE))
		    for (yr in 1:length(YEARS_USE))
		    {
		      yr_pred = which(YEARS_pred == YEARS_USE[yr])
		      yr_mod = which(YEARS == YEARS_USE[yr])
		      Data.projection <- matrix(NA, nrow=length(unique(Predict_data$LONG))*length(unique(Predict_data$LAT)), ncol=(ncol(DESIGN.pred[[yr_pred]])))
		      Data.projection <- cbind(Data.projection, expand.grid(LONG=seq(range(Predict_data$LONG)[1],range(Predict_data$LONG)[2],length.out=length(unique(Predict_data$LONG))), LAT=seq(range(Predict_data$LAT)[1],range(Predict_data$LAT)[2],length.out=length(unique(Predict_data$LAT)))))
		      Data.projection[, 1:(ncol(Data.projection)-2)] <- DESIGN.pred[[yr_pred]]# match the (X.Y) data in Projection matrix data to the one present in pred.dat

		      Factors <- which(rownames(result$summary.fix)%in%colnames(DESIGN.pred[[yr_pred]]))

		      Data.projection <- Data.projection[,1:(ncol(Data.projection)-2)]

		      # extract all coefficients
		      all_coef <- x$latent[,1][(1+which_field_mesh[length(which_field_mesh)]):length(x$latent[,1])]
		      names_coef <- names(all_coef)
		      # determine which ones are present in the current model
		      which_coeff <- sapply(1:length(names(DESIGN.pred[[yr_pred]])), function(y) which(names(DESIGN.pred[[yr_pred]])[y] == names_coef))

		      link_y_hat[,yr] = c(inla.mesh.project(mesh_proj, field_mesh[,yr_mod]) + matrix(as.matrix(Data.projection)%*%as.numeric(all_coef[which_coeff]),nrow=length(unique(Predict_data$LONG)), ncol=length(unique(Predict_data$LAT))))
		    }

		    return(link_y_hat)
		  }

		  Posteriors <- lapply(s1000, predictions)

		  save(Posteriors, file=paste0(getwd(), "/", NAME,"/Posterior.sample.dmp"))
		  rm(s1000, Posteriors)

		}




		######################################################
		############### Run the positive model ###############
		######################################################

		NAME_pos <- paste("Pos",paste(Which_region, collapse="_"),paste(Which_species, collapse="_"),sep="_")

		# the data for the positive model
		Pos_dat <- subset(Data_work, subset=c(Pres>0))

		if (all(c("AI") %in% Which_region)) Pos_dat$LONGITUDE[which(Pos_dat$LONGITUDE>0)] <- Pos_dat$LONGITUDE[which(Pos_dat$LONGITUDE>0)]-360
		Coord <- cbind(Pos_dat$LONGITUDE, Pos_dat$LATITUDE)					# the coordinates
		bnd = inla.nonconvex.hull(Coord, convex=200000)		# I need to get the boundary limit of the region on interest

		if (all(sort(c("EBS", "SLP", "AI", "GOA"))== sort(Which_region))) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(sort(c("AI", "GOA"))== sort(Which_region))) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(sort(c("EBS", "GOA"))== sort(Which_region))) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(sort(c("EBS", "SLP"))== sort(Which_region))) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(c("EBS") == Which_region)) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(c("SLP") == Which_region)) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(1000000,2000000), cutoff = 75000, offset=c(50000,100000))
		if (all(c("GOA") == Which_region)) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(3000000,5000000), cutoff = 50000, offset=c(100000,100000))
		if (all(c("AI") == Which_region)) mesh_pos = inla.mesh.2d(loc=Coord, boundary = bnd, max.edge=c(1000000,3000000), cutoff = 30000, offset=c(50000,100000))

		spde=inla.spde2.matern(mesh_pos, alpha=3/2)	# exponential decaying spatial correlation

		### The design matrix

		N.bin 	<- 	dim(Pos_dat)[1]
		## for the positive data only
		dat.pos <- Pos_dat
		X.1.pos <- X.1.all[which(Data_work$Pres>0),]
		Covar.names <- colnames(X.1.pos)										# the list of Env_covariate for the estimation model
		XX.list <- as.list(X.1.pos)

		### Setting the spatio-temporal model

		###### Specifies how to allocate the ID for each mesh_pos
		if (Spacetime=="NO") field.indices = inla.spde.make.index("field", n.spde=mesh_pos$n)
		if (Spacetime=="AR") field.indices = inla.spde.make.index("field", n.spde=mesh_pos$n, n.group=length(YEARS))
		if (Spacetime=="IID") field.indices = inla.spde.make.index("field", n.spde=mesh_pos$n, n.repl=length(YEARS))

		field.group <- rep(1:length(YEARS), as.numeric(table(factor(Pos_dat$YEAR, levels=YEARS))))
		field.repl <- rep(1:length(YEARS), as.numeric(table(factor(Pos_dat$YEAR, levels=YEARS))))
		if (Spacetime=="AR") A <- inla.spde.make.A(mesh_pos, loc=Coord, group=field.group)
		if (Spacetime=="IID") A <- inla.spde.make.A(mesh_pos, loc=Coord, repl=field.repl)
		if (Spacetime=="NO") A <- inla.spde.make.A(mesh_pos, loc=Coord)
		if(Spacetime != "NADA")
		{
		  ## creating the list of object to be used to define covariates that I want to keep for the estimation model
		  effect.list <- list()
		  effect.list[[1]] <- c(field.indices, list(Intercept=1))
		  for (i in 1:ncol(X.1)) effect.list[[i+1]] <- XX.list[[i]]
		  names(effect.list) <- c("1", Covar.names)

		  ## creating the list of predictors (for the estimation model) that match the one specified above
		  A.list <- list()
		  A.list[[1]] <- A
		  for (i in 1:ncol(X.1)) A.list[[i+1]] <- 1

		  ## model frame specification
		  stack <- inla.stack(data=list(resp=(Pos_dat$CPUE)), A=A.list, effects=effect.list, tag='est')
		}

		## model formula specification --> this is the main part to test different models. For example, assuming different spatial structure every year, or serial autocorrelation
		formula0 = as.formula(paste0("log(resp) ~ ", paste(Covar.names, collapse="+")))
		formula1 = as.formula(paste0("log(resp) ~ -1 + Intercept +",  paste(Covar.names, collapse="+"), "+ f(field, model=spde)"))
		formula2 = as.formula(paste0("log(resp) ~ -1 + Intercept +",  paste(Covar.names, collapse="+"), "+ f(field, model=spde, group=field.group, control.group=list(model='ar1'))"))
		formula3 = as.formula(paste0("log(resp) ~ -1 + Intercept +",  paste(Covar.names, collapse="+"), "+ f(field, model=spde, replicate=field.repl, control.group=list(model='exchangeable'))"))

		# Run model
		if (Overwrite == TRUE | file.exists(NAME_pos)==FALSE)
		{
		  unlink(x=paste0(getwd(), "/", NAME_pos), recursive = TRUE, force = TRUE)
		  if (Spacetime=="NADA") result_pos = inla(formula0, family="gaussian", data=cbind(resp=dat$y, X.1), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), keep=FALSE, working.directory= NAME_pos)
		  if (Spacetime=="NO") result_pos = inla(formula1, family="gaussian", data=inla.stack.data(stack, spde=spde), verbose=TRUE, control.predictor=list(A=inla.stack.A(stack), compute=TRUE), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE), keep=FALSE, working.directory= NAME_pos)
		  if (Spacetime=="AR") result_pos = inla(formula2, family="gaussian", data=inla.stack.data(stack, spde=spde), verbose=TRUE, control.predictor=list(A=inla.stack.A(stack), compute=TRUE), control.inla= list(int.strategy = "eb"), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE, openmp.strategy="large"), keep=FALSE, num.threads=2, working.directory= NAME_pos)
		  if (Spacetime=="IID") result_pos = inla(formula3, family="gaussian", data=inla.stack.data(stack, spde=spde), verbose=TRUE, control.predictor=list(A=inla.stack.A(stack), compute=TRUE), control.compute = list(dic=TRUE, cpo=TRUE, config=TRUE, openmp.strategy="large"), keep=FALSE, num.threads=2, working.directory= NAME_pos)
		}

		# if you have already run the model
		if (Overwrite == FALSE)
		{
		  result_pos <- inla.collect.results(paste(NAME_pos, "/results.files", sep=""))
		}

		summary(result_pos)
		spde.result_pos = inla.spde.result(result_pos, "field", spde, do.transf=TRUE)

		# Calculate the estimated expected surface (for the expected value of parameters, replace by NA all values that are outside the range of the Predict_data)
		mesh_proj = inla.mesh.projector(mesh_pos, xlim=range(Predict_data$LONG), ylim=range(Predict_data$LAT), dims=c(length(unique(Predict_data$LONG)),length(unique(Predict_data$LAT))))

		link_y_hat_pos <- list()
		LATENT_FIELD = matrix(result_pos$summary.ran$field$mean, mesh_pos$n, length(YEARS))

		YEARS_USE = intersect(YEARS, YEARS_pred)
		for (yr in 1:length(YEARS_USE))
		{
		  yr_pred = which(YEARS_pred == YEARS_USE[yr])
		  yr_mod = which(YEARS == YEARS_USE[yr])
		  Data.projection <- matrix(NA, nrow=length(unique(Predict_data$LONG))*length(unique(Predict_data$LAT)), ncol=(ncol(DESIGN.pred[[yr_pred]])))
		  Data.projection <- cbind(Data.projection, expand.grid(LONG=seq(range(Predict_data$LONG)[1],range(Predict_data$LONG)[2],length.out=length(unique(Predict_data$LONG))), LAT=seq(range(Predict_data$LAT)[1],range(Predict_data$LAT)[2],length.out=length(unique(Predict_data$LAT)))))
		  Data.projection[, 1:(ncol(Data.projection)-2)] <- DESIGN.pred[[yr_pred]]				# replace the above values in the Projection matrix data

		  Factors <- which(rownames(result_pos$summary.fix)%in%colnames(DESIGN.pred[[yr_pred]]))

		  Data.projection <- Data.projection[,1:(ncol(Data.projection)-2)]

		  link_y_hat_pos[[yr]] = inla.mesh.project(mesh_proj,LATENT_FIELD[,yr_mod]) +  matrix(as.matrix(Data.projection)%*%as.numeric(result_pos$summary.fix[Factors,'mean']),nrow=length(unique(Predict_data$LONG)), ncol=length(unique(Predict_data$LAT)))
		}


		### Do mcmc sampling
		if (Mcmc_prediction ==TRUE)
		{
		  set.seed(Mcmc_seed)
		  s1000 <- inla.posterior.sample(Mcmc_sample, result_pos)

		  mesh_proj = inla.mesh.projector(mesh_pos, xlim=range(Predict_data$LONG), ylim=range(Predict_data$LAT), dims=c(length(unique(Predict_data$LONG)),length(unique(Predict_data$LAT))))

		  predictions <- function(x)
		  {
		    which_field_mesh <- grep("field", names(x$latent[,1]))
		    field_mesh = matrix(as.numeric(x$latent[,1])[which_field_mesh], mesh_pos$n, length(YEARS))

		    link_y_hat <- matrix(NA, nrow=length(unique(Predict_data$LONG))*length(unique(Predict_data$LAT)), ncol=length(YEARS_USE))
		    for (yr in 1:length(YEARS_USE))
		    {
		      yr_pred = which(YEARS_pred == YEARS_USE[yr])
		      yr_mod = which(YEARS == YEARS_USE[yr])
		      Data.projection <- matrix(NA, nrow=length(unique(Predict_data$LONG))*length(unique(Predict_data$LAT)), ncol=(ncol(DESIGN.pred[[yr_pred]])))
		      Data.projection <- cbind(Data.projection, expand.grid(LONG=seq(range(Predict_data$LONG)[1],range(Predict_data$LONG)[2],length.out=length(unique(Predict_data$LONG))), LAT=seq(range(Predict_data$LAT)[1],range(Predict_data$LAT)[2],length.out=length(unique(Predict_data$LAT)))))
		      Data.projection[, 1:(ncol(Data.projection)-2)] <- DESIGN.pred[[yr_pred]]# match the (X.Y) data in Projection matrix data to the one present in pred.dat

		      Factors <- which(rownames(result_pos$summary.fix)%in%colnames(DESIGN.pred[[yr_pred]]))

		      Data.projection <- Data.projection[,1:(ncol(Data.projection)-2)]

		      # extract all coefficients
		      all_coef <- x$latent[,1][(1+which_field_mesh[length(which_field_mesh)]):length(x$latent[,1])]
		      names_coef <- names(all_coef)
		      # determine which ones are present in the current model
		      which_coeff <- sapply(1:length(names(DESIGN.pred[[yr_pred]])), function(y) which(names(DESIGN.pred[[yr_pred]])[y] == names_coef))

		      link_y_hat[,yr] = c(inla.mesh.project(mesh_proj, field_mesh[,yr_mod]) + matrix(as.matrix(Data.projection)%*%as.numeric(all_coef[which_coeff]),nrow=length(unique(Predict_data$LONG)), ncol=length(unique(Predict_data$LAT))))
		    }

		    return(link_y_hat)
		  }

		  Posteriors <- lapply(s1000, predictions)
		  SD_obs_mcmc <- sapply(1:Mcmc_sample, function(x) as.numeric(sqrt(1/exp(s1000[[x]]$hyperpar[1]))))		# to get posterior sample of observation error SD. INLA hyperprior is on X=log(precision). therefore, to get the SD = sqrt(1/exp(X))

		  save(Posteriors, file=paste0(getwd(), "/", NAME_pos,"/Posterior.sample.dmp"))
		  save(SD_obs_mcmc, file=paste0(getwd(), "/", NAME_pos,"/Posterior.sample.obserror.dmp"))
		  rm(s1000, Posteriors)
		}


		###############################################################
		############### Combine the models for prediction #############
		###############################################################

		if (Mcmc_prediction==TRUE)
		{
		  load(paste0(getwd(), "/", NAME,"/Posterior.sample.dmp"))
		  link_y_hat_bin_mcmc <- Posteriors
		  load(paste0(getwd(), "/", NAME_pos,"/Posterior.sample.dmp"))
		  link_y_hat_pos_mcmc <- Posteriors
		  load(paste0(getwd(), "/", NAME_pos,"/Posterior.sample.obserror.dmp"))
		  SD_obs_mcmc <- SD_obs_mcmc

		  SDM_result <- list()
		  for (samp in 1:Mcmc_sample)
		  {
		    y_hat_pred_mcmc <- list()
		    X <-  rep(unique(Predict_data$LONG), length(unique(Predict_data$LAT)))
		    Y <- rep(unique(Predict_data$LAT), each=length(unique(Predict_data$LONG)))
		    for (yr in 1:length(YEARS_USE))
		    {
		      mcmc_predict <- data.frame(X,Y)
		      ## this is to calculate the predictive distribution for each mcmc
		      if (Abundance_or_catch == "Abundance") mcmc_predict$predicted <- inv.logit(link_y_hat_bin_mcmc[[samp]][,yr])*exp(link_y_hat_pos_mcmc[[samp]][,yr])
		      ## this is to calculate the sampled catch from each mcmc
		      if (Abundance_or_catch == "Catch") mcmc_predict$predicted <- rbinom(n=length(link_y_hat_bin_mcmc[[samp]][,yr]),size=1,inv.logit(link_y_hat_bin_mcmc[[samp]][,yr]))*exp(rnorm(length(link_y_hat_pos_mcmc[[samp]][,yr]),link_y_hat_pos_mcmc[[samp]][,yr],SD_obs_mcmc[samp]))
		      y_hat_pred_mcmc[[yr]] <- mcmc_predict
		    }

		    # Keeping only the prediction falling within a certain region
		    which_columns <- match(Which_region,names(Predict_data))
		    ifelse(length(which_columns)>1, To_keep <- which(apply(Predict_data[,which_columns], 1, sum)>0), To_keep <- which(Predict_data[,which_columns]>0))
		    ifelse(length(which_columns)>1, Outside_survey <- which(apply(Predict_data[,which_columns], 1, sum)==0), Outside_survey <- which(Predict_data[,which_columns]==0))

		    temp <- matrix(NA, nrow=nrow(y_hat_pred_mcmc[[1]]), ncol=length(YEARS_USE))
		    for (yr in 1:length(YEARS_USE))
		    {
		      new_dat <- y_hat_pred_mcmc[[yr]]
		      # Land <- which(is.na(new_dat[,3]) == TRUE)
		      new_dat[Outside_survey,3] <- NA
		      # new_dat[Land,3] <- NA
		      temp[,yr] <- new_dat[,3]
		    }

		    SDM_result[[samp]] <- temp
		  }
		}


		if (Mcmc_prediction==FALSE)
		{
		  y_hat_pred <- list()
		  for (yr in 1:length(YEARS_USE))
		  {
		    y_hat_pred[[yr]] <- inv.logit(link_y_hat_bin[[yr]])*exp(link_y_hat_pos[[yr]])
		  }

		  # Keeping only the prediction falling within a certain region
		  which_columns <- match(Which_region,names(Predict_data))
		  ifelse(length(which_columns)>1, To_keep <- which(apply(Predict_data[,which_columns], 1, sum)>0), To_keep <- which(Predict_data[,which_columns]>0))
		  ifelse(length(which_columns)>1, Outside_survey <- which(apply(Predict_data[,which_columns], 1, sum)==0), Outside_survey <- which(Predict_data[,which_columns]==0))

		  SDM_result <- list()
		  temp <- matrix(NA, nrow=length(y_hat_pred[[1]]), ncol=length(YEARS_USE))
		  for (yr in 1:length(YEARS_USE))
		  {
		    new_dat <- y_hat_pred[[yr]]
		    # Land <- which(is.na(new_dat) == TRUE)
		    new_dat[Outside_survey] <- NA
		    # new_dat[Land] <- NA
		    temp[,yr] <- new_dat
		  }

		  SDM_result[[1]] <- temp

		}

	### Object to return

		OUTPUT <- list()
		OUTPUT$YEARS_USE <- YEARS_USE
		OUTPUT$SDM_result <- SDM_result

		return(OUTPUT)
	}

