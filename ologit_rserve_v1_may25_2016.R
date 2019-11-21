library(VGAM)

Rplink <- function(PHENO, GENO, CLUSTER, COVAR)
{
    f1 <- function(s){    # s denotes the genotype vector at the tested SNP

    subid <- seq(1, length(PHENO))


       ##important: Need to remove subjects with missing data before analysis
       ##Some R libraries will crash if there are any missing data (NAs) in the dataset to be analyzed:
       ##olddat <- cbind(PHENO, s,  subid)
	   olddat <- cbind(PHENO,s,COVAR,subid)  
       ## Find rows containing missing data:
       badid <- unique(na.action(na.omit(olddat)))
	
        ## Remove rows containing missing data from dataset:
        if(length(badid) >= 1){
            PHENO  <- PHENO[-c(badid)]
            s      <- s[-c(badid)]
            COVAR <- COVAR[-c(badid),]
            subid  <- subid[-c(badid)]
        }
       COVAR_mat <- as.matrix(COVAR) # formatting change needed for some libraries

        ## Fit the regression model.

         ##Adam: I have fit this into a 'try' function so that if there is an error, it is
         ##handled gracefully (by default, if there is a single NA result SNP, all supplied
         ##SNPs (around 200 SNPs) are given NA values as well, regardless of if they should or shouldnt

	 ## The vglm function comes as part of the VAGM library and allows for proportional hazards regression

gee.fit <- try(
                       vglm(PHENO ~s + COVAR_mat,family=propodds ),
                                silent = TRUE)
	
        if(class(gee.fit) !=  "try-error")
        {
                #Get the Beta, SE, z stat
                #coefres <- summary(gee.fit)$coefficients[1,]
		 r <- coef(summary(gee.fit))["s",]
        } else {
                r <- c(-1,-1,-1,-1) #If the model didn't fit, just supply impossible values
        }
        ## Have the function return the results
        c(length(r), r)


    }  # f1()

    apply(GENO,2,f1)

} # Rplink()
		

