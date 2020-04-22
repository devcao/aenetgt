#'Generates individual testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se The testing sensitivity used for individual testing.
#' @param Sp The testing specificity used for individual testing.
#' @param cj This is an extraneous argument and is here only in order that all four assay.gen functions take the same arguments. Default is \code{cj=1} and an error will be returns if \code{cj!=1}.
#' @return a list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates individual level testing and stores the 
#' testing responses in accordance to the data structure 
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure 
#' see McMahan et al. (2017).
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 100
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- .95 # set assay sensitivity
#' Sp <- .97 # set assay specificity
#' # subject individuals to individual testing
#' assay.data <- individual.assay.gen(Y.true,Se,Sp,cj=1)
individual.assay.gen <-function(Y.true,Se,Sp,cj=1){
if(cj!=1 || max(length(Se),length(Sp))>1) 
{
	print("cj must be 1 for individual testing and Se and Sp must each be of length 1")
	
} else {
	
N <- length(Y.true)
	
Y<-matrix(-99,nrow=N, ncol=3) 		# initialize output matrix Y
Z<-matrix(-99,nrow=N,ncol=cj+3) 	# initialize output matrix Z
	
Z[,1] <- rbinom(N,1,(Se*Y.true+(1-Sp)*(1-Y.true)))
Z[,2] <- 1
Z[,3] <- 1
Z[,4] <- 1:N

Y[,2] <- 1
Y[,3] <- 1:N
return(list("Z"=Z, "Y"=Y))
}

}

#'Generates master pool testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se The master pool testing sensitivity.
#' @param Sp The master pool testing specificity.
#' @param cj The size of the master pools (Note: The number of individuals \code{length(Y.true)} should be 
#'      evenly divisible by \code{cj}, this is only for decoding purposes; i.e., 
#'      the regression methods do not require this condition).
#' @return a list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates Initial pool testing and stores the 
#' testing responses in accordance to the data structure 
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure 
#' see McMahan et al. (2017).
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 100
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- .95 # set master pool assay sensitivity
#' Sp <- .97 # set master poot assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to master pool testing
#' assay.data <- masterpool.assay.gen(Y.true,Se,Sp,cj)
masterpool.assay.gen <-function(Y.true,Se,Sp,cj){

if(length(Se)>1 | length(Sp)>1) 
{
	
	print("error: Se or Sp has length greater than 1")
	return(NA)
	
} else {
		
	N<-length(Y.true)
	Jmax<-N/cj 
	J<-1
	
	Y<-matrix(-99,nrow=N, ncol=3) 
	Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 
	
	
	for(j in 1:(N/cj)){
	prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
	Z[J,1]<-rbinom(n=1,size=1,prob=prob)
	Z[J,2]<-cj
	Z[J,3]<-1
	Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
	Y[((j-1)*cj+1):(j*cj),2]<-1
	Y[((j-1)*cj+1):(j*cj),3]<-J
	J<-J+1
	}
	
	J<-J-1
	Z<-Z[1:J,]
	
	return(list("Z"=Z, "Y"=Y))
	
	}
}

#'Generates Dorfman testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing sensitivity for the master pools and the second entry is the 
#'      test sensitivity for individual testing.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for the master pools and the second entry is the 
#'      test specificity for individual testing.
#' @param cj The size of the master pools (Note: The number of individuals \code{length(Y.true)} should be 
#'      evenly divisible by \code{cj}. This is only for decoding purposes; i.e., 
#'      the regression methods do not require this condition).
#' @return A list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates Dorfman decoding and stores the 
#' testing responses in accordance to the data structure 
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure 
#' see McMahan et al. (2017).
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 100
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95,.92) # set master pool and individual assay sensitivity
#' Sp <- c(.97,.98) # set master pool and individual assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to Dorfman testing
#' assay.data <- dorfman.assay.gen(Y.true,Se,Sp,cj)
dorfman.assay.gen <-function(Y.true,Se,Sp,cj){
N<-length(Y.true)
Jmax<-N+N/cj 
J<-1

Y<-matrix(-99,nrow=N, ncol=4) 
Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 
D<-numeric(N)


for(j in 1:(N/cj)){
prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se[1],1-Sp[1])
Z[J,1]<-rbinom(n=1,size=1,prob=prob)
Z[J,2]<-cj
Z[J,3]<-1
Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
Y[((j-1)*cj+1):(j*cj),2]<-1
Y[((j-1)*cj+1):(j*cj),3]<-J
J<-J+1
if(Z[J-1,1]==1){
for(k in ((j-1)*cj+1):(j*cj)){
prob<-ifelse(Y.true[k]>0,Se[2],1-Sp[2])
Z[J,1]<- rbinom(n=1,size=1,prob=prob)
Z[J,2]<-1
Z[J,3]<-2
Z[J,4]<-k
Y[k,2]<-2
Y[k,4]<-J
D[k] <- Z[J,1] # store individual assay results
J<-J+1
}
} else { D[((j-1)*cj+1):(j*cj)]<-0 } # store individual assay results
}

J<-J-1
Z<-Z[1:J,]

return(list("Z"=Z, "Y"=Y, "D" = D))
}

#'Generates array testing data.
#'
#' @param Y.true The true disease statuses of the individuals.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing sensitivity for the row/column pools and the second entry is the 
#'      test sensitivity for individual testing.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for the row/column pools and the second entry is the 
#'      test specificity for individual testing.
#' @param cj Row and column pool sizes to be used (Note: The number of individuals 
#'      should be evenly divisible by \code{cj*cj}. This is only for decoding 
#'      purposes; i.e., the regression methods do not require this condition)
#' @return A list containing objects \code{Z} and \code{Y}.
#'
#' This function simulates array decoding and stores the 
#' testing responses in accordance to the data structure 
#' required to fit the group testing regression model presented
#' in Gregory et al. (2018+). For the specifics of this structure 
#' see McMahan et al. (2017).
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 100
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95,.92) # set row/col and individual assay sensitivity
#' Sp <- c(.97,.98) # set row/col and individual assay specificity
#' cj <- 4 # set dimension of arrays 
#' # subject individuals to array testing
#' assay.data <- array.assay.gen(Y.true,Se,Sp,cj)
array.assay.gen <- function(Y.true, Se, Sp, cj){
N<-length(Y.true) 					# get number of individuals
Jmax<-2*N/cj + N					# specify maximum number of assays
J<-1							# initialize assay index
AT<-N/(cj^2)						# number of arrays (so N must be divisible by cj^2)

Y<-matrix(-99,nrow=N, ncol=5) 			# initialize output matrix Y
Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 			# initialize output matrix Z

Y.A<-array(-99,c(cj,cj,AT))				# initialize array object to contain true individual statuses in the arrays
ID.A<-array(-99,c(cj,cj,AT))				# initialize array object to contain individual indices in the array tests
ind<-1
for(i in 1:AT){
for(j in 1:cj){
for(m in 1:cj){
Y.A[m,j,i]<-Y.true[ind]				# populate Y.A with the true statuses for the individuals in each array
ID.A[m,j,i]<-ind					# populate ID.A with the IDs of the individuals in each array
ind<-ind+1
}}}

D <- numeric(N)					#*** initialize vector to be populated with the individual diagnoses

for(s in 1:AT){					# begin loop through arrays

array.yk<-Y.A[,,s]					# pull the true statuses from array s
array.id<-ID.A[,,s]					# pull the individual IDs from array s

a<-rep(0,nrow(array.yk))				# initialize vector in which to store the row assays on array s
b<-rep(0,ncol(array.yk))				# initialize vector in which to store the column assays on array s

for(i in 1:cj){					# carry out row assays on array s
   prob<- ifelse(sum(array.yk[i,])>0, Se[1], 1-Sp[1]) # compute probability of a positive assay on row i of array s
   g<- rbinom(n=1,size=1,prob=prob)			# generate random assay on row i of array s
   a[i]<-g						# store randomly generated assay on row i of array s in position i of the vector a
   Z[J,1]<-g 						# store this also in row J of Z, column 1
   Z[J,2]<-cj 						# store the pool size of this assay in row J of Z, column 2
   Z[J,3]<-1						# store a 1 in row J of Z, column 3, to be used later to reference Se[1], Sp[1]
   Z[J,4:(cj+3)]<-array.id[i,]			# store the indices of the individuals in row i of array s in row J of Z, in the remaining cj columns
   Y[array.id[i,],2]<-2				# store a 2 in the rows of Y corresponding to the individuals in row i of array s; they are tested twice (once in a row, once in a column)
   Y[array.id[i,],3]<-J				# store J, the index of the assay, in the rows of Y corresponding to the individuals in row i of array s
   J<-J+1						# increment assay index
}
for(j in 1:cj){					# carry out column assays on array s
   prob<- ifelse(sum(array.yk[,j])>0, Se[1], 1-Sp[1]) # compute probability of a positive assay on column j of array s
   g<- rbinom(n=1,size=1,prob=prob)			# generate random assay on column j of array s
   b[j]<-g						# store randomly generated assay on column j of array s in position j of the vector b
   Z[J,1]<-g 							 
   Z[J,2]<-cj 								
   Z[J,3]<-1							
   Z[J,4:(cj+3)]<-array.id[,j]			
   Y[array.id[,j],4]<-J					
   J<-J+1						
}

if(sum(a)>0 & sum(b)>0){				# if at least one row and at least one column of array s tested positive...
array.yk1<-as.matrix(array.yk[(a==1),(b==1)]) 	# pull true statuses of individuals at intersections of rows and columns with positive assays
array.id1<-as.matrix(array.id[(a==1),(b==1)]) 	# pull individual IDs at intersections of rows and columns with positive assays
for(i in 1:nrow(array.yk1)){				# carry out individual assays on individuals at intersections of rows and columns with positive assays.
for(j in 1:ncol(array.yk1)){			
   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])	# get probability of positive assay on individual i,j in intersection
   g<- rbinom(n=1,size=1,prob=prob)			# generate a random assay result for this individual using this probability
   Z[J,1]<-g 						# store this random assay result in row J of Z, column 1
   Z[J,2]<-1 						# store a 1 in row J of Z, column 1, to indicate that this assay was done on 1 individual
   Z[J,3]<-2						# store a 2 in row J of Z, column 2, to be used later to reference Se[2] and Sp[2]
   Z[J,4]<-array.id1[i,j]				# store in row J of Z, column 4, the ID of the individual to whom this assay corresponds
   Y[array.id1[i,j],2]<-3				# store a 3 in the row of Y corresponding to this individual to indicate that this is the third test in which the individual was assayed
   Y[array.id1[i,j],5]<-J				# store J, the index of the assay, in the row of Y corresponding to this individual
   D[array.id1[i,j]] <- g 				#*** store final diagnosis for this individual in the row of D corresponding to this individual
   J<-J+1						# increment assay index J
}}}

if(sum(a)>0 & sum(b)==0){				# if no columns but at least one row of array s tested positive...
array.yk1<-as.matrix(array.yk[(a==1),])		# pull true statuses of individuals in rows of array s which tested positive
array.id1<-as.matrix(array.id[(a==1),])		# pull IDs of individuals in rows of array s which tested positive
for(i in 1:nrow(array.yk1)){				# carry out individual assays on individuals in rows of array s which tested positive
for(j in 1:ncol(array.yk1)){			
   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])# get probability of a positive assay for this individual
   g<- rbinom(n=1,size=1,prob=prob)			# generate random assay result for this individual using this probability
   Z[J,1]<-g 						# store this random assay result in row J of Z, column 1
   Z[J,2]<-1 						# store a 1 in row J of Z, column 2, to indicate that this assay was carried out on 1 individual
   Z[J,3]<-2						# store a 2 in row J of Z, column 3, to be used later to reference Se[2] and Sp[2]
   Z[J,4]<-array.id1[i,j]				# store in row J of Z, column 4, the ID of the individual to whom this assay corresponds
   Y[array.id1[i,j],2]<-3				# store a 3 in the row of Y corresponding to this individual to indicate that this is the third test in which the individual was assayed 
   Y[array.id1[i,j],5]<-J				# store J, the index of the assay, in the row of Y corresponding to this individual
   D[array.id1[i,j]] <- g				#*** store final diagnosis for this individual in the row of D corresponding to this individual
   J<-J+1						# increment assay index J
}}}

if(sum(a)==0 & sum(b)>0){				# if no rows but at least one column of array s tested positive...
array.yk1<-as.matrix(array.yk[,(b==1)])
array.id1<-as.matrix(array.id[,(b==1)])
for(i in 1:nrow(array.yk1)){
for(j in 1:ncol(array.yk1)){
   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
   g<- rbinom(n=1,size=1,prob=prob)
   Z[J,1]<-g 
   Z[J,2]<-1 
   Z[J,3]<-2
   Z[J,4]<-array.id1[i,j]
   Y[array.id1[i,j],2]<-3
   Y[array.id1[i,j],5]<-J
   D[array.id1[i,j]] <- g				#*** store final diagnosis for this individual in the row of D corresponding to this individual
   J<-J+1
}}}

} 

J<-J-1
Z<-Z[1:J,]

return(list("Z"=Z, "Y"=Y, "D"=D))
}

#' Pulls individual diagnoses from group testing data if available
#'
#' @param \code{Z} output from one of the functions \code{individial.assay.gen()},\code{masterpool.assay.gen()},\code{dorfman.assay.gen()}, or \code{array.assay.gen()}.
#' @param \code{Y} output from one of the functions \code{individial.assay.gen()},\code{masterpool.assay.gen()},\code{dorfman.assay.gen()}, or \code{array.assay.gen()}.
#' @return a vector of 0s and 1s which are the individual diagnoses, NULL if \code{Z} and \code{Y} come from \code{masterpool.assay.gen()}.
#'
#' This function pulls the individual diagnoses from
#' matrices Z and Y when individual diagnoses are available.
#' So for masterpool testing, NULL will be returned.
pull.diagnoses <- function(Z,Y)
{
	if( (sum(Z[,2] > 1) == nrow(Z)) & sum(Z[,1]==1) > 0  ) # check to see if masterpool testing has been done.
	{	# all tests on pools of more than one individual? and all tests not negative?
		
		print("no individual tests available")
		return(NULL)
		
	} 
	
	# Only individuals who are tested individually can be diagnosed as positive.  
	# This is true in Dorfman and Array testing.
	
	Z.ind.assays <- Z[which(Z[,2]==1),]
	pos.ind <- Z.ind.assays[Z.ind.assays[,1]==1,4]
	D <- numeric(nrow(Y))
	D[pos.ind] <- 1
	
	return(D)
	
}

#' Computes conditional expectations of individual disease statuses.
#'
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param X Design matrix with first column a column of 1s.
#' @param b Parameter values at which to compute the conditional expectations.
#' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
#' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
#' @return The vector of conditional expectations.
#'
#' This function computes the conditional expectations of each individual
#' disease status, conditional on the observed assay data and the disease 
#' statuses of all other individuals.  This function is used in the EM algorithm
#' performed by the functions \code{mlegt}, \code{enetgt}, \code{enetgt.grid}, and 
#' \code{enetgt.grid.0}.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 100
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95,.92) # set master pool and individual assay sensitivity
#' Sp <- c(.97,.98) # set master pool and individual assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to Dorfman testing
#' assay.data <- dorfman.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' b <- data$b
#' EY <- EY.exact(Z,Y,X,b,Se,Sp)
EY.exact<-function(Z,Y,X,b,Se,Sp){

n<-dim(Y)[1]
p<-logit(b,X)
EY<-rep(-99,n)
id.go<-(1:n)[Y[,2]==1] # individuals tested in a group only
id.rt<-(1:n)[Y[,2]>1]  # individuals who are retested after being tested in a pool

##################################################################
# Handles MPT only or negative master pools under Dorfman testing
while(length(id.go)>0){
pid<-Y[id.go[1],3]  			# assay index for group in which individual was tested
Zj<-Z[pid,1]					# the group assay result
cj<-Z[pid,2]					# the size of the group
Se.p<-Se[Z[pid,3]]				# the group testing sensitivity
Sp.p<-Sp[Z[pid,3]]				# the group testing specificity
grp<-Z[pid,(4:(3+cj))]			# the individuals in the group
pg<-p[grp]						# p(x'b) for individuals in the group
PZt0<-prod(1-pg)				# probability of \tilde Z = 0
PZt1<-1-PZt0					# probability of \tilde Z = 1
PZ0<-(1-Se.p)*PZt1+Sp.p*PZt0	# probability of Z = 0
PZ1<-Se.p*PZt1+(1-Sp.p)*PZt0	# probability of Z = 1
A<-pg*Se.p
B<-pg*(1-Se.p)
# Completes the E-Step
EY[grp]<-Zj*A/PZ1 + (1-Zj)*B/PZ0
id.go<-id.go[(id.go %in% grp)==FALSE] 
}

###########################################################
# Handles positive master pools under Dorfman testing 

while(length(id.rt)>0){

tid<-Y[id.rt[1],c(3,4)]
pid<-tid[Z[tid,2]>1]
cj<-Z[pid,2]
Se.p<-Se[Z[pid,3]]
Sp.p<-Sp[Z[pid,3]]
grp<-Z[pid,(4:(3+cj))]
pg<-p[grp]

iid<-apply(Y[grp,c(3,4)],1,sum)-pid
Se.i<-Se[Z[iid,3]]
Sp.i<-Sp[Z[iid,3]]
Yij<-Z[iid,1]

TS<-expand.grid(rep(list(0:1),cj))
DS<-matrix(Yij,nrow=(2^cj),ncol=cj,byrow=TRUE)

PTS<-t(t(TS)*pg)+t(t(1-TS)*(1-pg))
PTO<- t(t(TS*DS)*Se.i) + t(t(TS*(1-DS))*(1-Se.i)) + t(t((1-TS)*DS)*(1-Sp.i)) + t(t((1-TS)*(1-DS))*Sp.i)
TSP<-apply(TS,1,max)
PMO<-Se.p*TSP+(1-Sp.p)*(1-TSP)

denom<-sum(apply(cbind(PMO,PTO,PTS),1,prod))
num<-rep(-99,cj)
for(i in 1:cj){
num[i]<-sum(apply(cbind(Se.p,PTO[TS[,i]==1,],PTS[TS[,i]==1,]),1,prod))
}
EY[grp]<-num/denom
id.rt<-id.rt[(id.rt %in% grp)==FALSE] 
}
return(EY)
}

#' Approximates the conditional expectations of individual disease statuses with Gibbs sampling.
#'
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param X Design matrix with first column a column of 1s.
#' @param b Parameter values at which to compute the conditional expectations.
#' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
#' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
#' @param GI The length of the Gibbs sampling Markov chain.
#' @return The vector of conditional expectations.
#'
#' This function uses a Gibbs sampler to appriximate the conditional expectation of 
#' each individual's disease status, conditional on the observed assay data and the disease 
#' statuses of all other individuals. This function is used in the EM algorithm
#' performed by the functions \code{mlegt}, \code{enetgt}, \code{enetgt.grid}, and 
#' \code{enetgt.grid.0} under array testing.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 100
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95,.92) # set master pool and individual assay sensitivity
#' Sp <- c(.97,.98) # set master pool and individual assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to array testing
#' assay.data <- array.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' b <- data$b
#' EY <- EY.approx(Z,Y,X,b,Se,Sp,GI=5000)
EY.approx<-function(Z,Y,X,b,Se,Sp,GI=5000){
n<-dim(Y)[1]
na<-length(Se)
# The E-Step
p<-logit(b,X)
Y[,1] <- rbinom(n,1,p) # this is the initial set of Y values for the Gibbs
W <- EYgibbs(n,p,Y,Z,Se,Sp,na,GI)
EY<-W/GI
return(EY)
}
#' Approximates the conditional covariance between all pairs of disease statuses with Gibbs sampling.
#'
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param X Design matrix with first column a column of 1s.
#' @param b Parameter values at which to compute the conditional expectations.
#' @param Se A vector of testing sensitivities of length \code{max(Z[,3])}.
#' @param Sp A vector of testing specificities of length \code{max(Z[,3])}.
#' @param EY The vector of conditional expectations of individual disease statuses.
#' @param GI The length of the Gibbs sampling Markov chain.
#' @return An approximation to the matrix of conditional covariances between all pairs of disease statuses.
#'
#' This function uses a Gibbs sampler to approximate the conditional covariance of 
#' all pairs of disease statuses, conditional on the observed assay data and the disease 
#' statuses of all other individuals.  This is used to implement Louis's method for computing 
#' the observed information matrix associated with the likelihood function, from which
#' Wald-type standard errors can be computed for the maximum likelihood estimators. 
#' For more details see Gregory et al. (2018+).
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model0(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95,.92) # set master pool and individual assay sensitivity
#' Sp <- c(.97,.98) # set master pool and individual assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to dorfman testing
#' assay.data <- dorfman.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' b.mle <- mlegt(X, Y, Z, Se, Sp, delta = .01)$b.mle # compute mle
#' EY <- EY.exact(Z,Y,X,b.mle,Se,Sp)
#' px <- as.numeric(logit(b.mle,X))
#' CovYiYj <- CovYiYj.approx(Z,Y,X,b.mle,Se,Sp,EY)
#' # use Louis' method to get the observed information matrix
#' Hess <- t(X) %*% ( - diag(px * (1 - px))  +  CovYiYj  ) %*% X	
#' b.cov.est <-  solve( - Hess ) 
#' b.mle.se <- sqrt(diag(b.cov.est)) # estimated standard errors of mles	
CovYiYj.approx <- function (Z, Y, X, b, Se, Sp, EY, GI = 50000)
{
	
	n <- nrow(Y)
	# W has in row i in the first column the number of the individuals involved in any group tests with individual i
	# including individual i itself.
	# The subsequent columns contain the IDs of those individuals
    W <- matrix(-99,n,n)
	n.cols.needed <- 2
	group.assays <- which(Z[,2]!=1)
	for(i in 1:n) # make it work with any mixture of testing regimes...
	{
		assays.i <- Y[i,3:(3+Y[i,2]-1)] # assays in which individual i participated
		group.assays.i <- intersect(group.assays, assays.i)  # group assays in which individual i participated
		if(length(group.assays.i)==0)
		{
		
			W[i,1] <- 1
			W[i,2] <- i
			
		} else {
			
			inds.with.i <- unique(as.numeric(Z[group.assays.i,-c(1:3)])) # individuals which participated with individual i in any group
            inds.with.i <- inds.with.i[inds.with.i>0] # remove any -99 values
			inds.with.i.but.leq.i <- inds.with.i[inds.with.i <= i]

			W[i,1:(length(inds.with.i.but.leq.i)+1)] <- c(length(inds.with.i.but.leq.i), inds.with.i.but.leq.i)	
			n.cols.needed <- max(n.cols.needed,length(inds.with.i.but.leq.i)+1)
			
		}
		
	}
	
	W <- W[,1:n.cols.needed,drop=FALSE]
	
    na <- length(Se)
    p <- logit(b, X)
    Y[,1] <- rbinom(n, 1, p)
    CovYiYj.lower.tri <- CovYiYjgibbs(n, p, Y, Z, W, Se, Sp, EY, na, GI = 50000)
    CovYiYj <- (CovYiYj.lower.tri + t(CovYiYj.lower.tri) - diag(diag(CovYiYj.lower.tri)))
    
    return(CovYiYj)  
      
}

#' Computes the mle on group testing data.
#'
#' @param X Design matrix with first column a column of 1s. 
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param binit Parameter value at which to initialize the EM-algorithm. The default is \code{b=1}, for which an initial value is chosen internally.
#' @param delta Convergence criterion.
#' @param E.approx Logical.  If \code{TRUE} then E-step done with \code{EY.approx()}. If \code{FALSE}, then E-step done with \code{EY.exact()}.
#' @param get.SEs Logical.  If \code{TRUE} then estimated standard errors for the maximum likelihood estimators are returned.
#' @return The maximum likelihood estimator.
#' 
#' This function implements an EM-algorithm to find the maximum likelihood estimator based on the observed data \code{X}, \code{Y}, \code{Z}, and the sensitivities and specificities in \code{Se}, \code{Sp}.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model0(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95,.92) # set master pool and individual assay sensitivity
#' Sp <- c(.97,.98) # set master pool and individual assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to dorfman testing
#' assay.data <- dorfman.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' mlegt.out <- mlegt(X, Y, Z, Se, Sp, delta = .01) # compute mle
mlegt <- function( X, Y, Z, Se, Sp, binit = 1, delta = 1e-3, E.approx = FALSE, get.SEs = FALSE)
{

	p <- ncol(X) - 1 # The first column of X is ones
	
	if(length(binit) == 1)
	{ 
		
		b1 <- c(-2,rep(0,p))

	} else {

		b1 <- binit
	}
	
	get.EY <- eval(parse(text = ifelse(E.approx, "EY.approx","EY.exact")))

	max.diff <- 1
	iter <- 1
	while(  max.diff > delta )
	{
		b0 <- b1
		EY <- get.EY(Z,Y,X,b0,Se,Sp)
		b1 <- logistic_enet( EY, X, 0, rep(1,ncol(X)), 0 , b0, delta )$b
		max.diff <- max(abs(b1 - b0))
		iter <- iter + 1
		if(is.na(max.diff)) break;
	}
	
	b.mle <- as.numeric(b1)
		
	if(get.SEs & is.na(b.mle[1])==FALSE)
	{
		# Get observed information matrix by Louis' Method:
		px <- as.numeric(logit(b.mle,X))
		CovYiYj <- CovYiYj.approx(Z,Y,X,b.mle,Se,Sp,EY)

		Hess <- t(X) %*% (  - diag(px * (1 - px))  +  CovYiYj  ) %*% X
	
		b.cov.est <-  solve( - Hess ) 
		b.mle.se <- sqrt(diag(b.cov.est))
		
	} else {
		
		b.mle.se <- rep( NA , p + 1 )
	}
	
	conv <- ifelse(is.na(max.diff),FALSE,TRUE)

	output <- list(	b.mle = b.mle,
					b.mle.se = b.mle.se,
					delta = delta,
					Se = Se,
					Sp = Sp,
					conv = conv,
					iter = iter
					)
					
	return(output)

}

#' Computes the elastic net estimator with weighted l1 norm on group testing data.
#'
#' @param X Design matrix with first column a column of 1s.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param lambda Tuning parameter.
#' @param theta Ridge versus lasso penalty mixer.
#' @param weights Vector of weights to be used in weighting the l1 penalty. Default is \code{weights=1}, which causes equal weights to be used for each coefficient.
#' @param binit Initial values for EM-algorithm.
#' @param delta Convergence criterion.
#' @param E.approx Logical.  If \code{TRUE} then E-step done with \code{EY.approx()}. If \code{FALSE}, then E-step done with \code{EY.exact()}.
#' @param get.SEs Logical.  If \code{TRUE} then estimated standard errors for the estimates are returned.
#' @return The elastic net estimator with weighted l1 norm under the choices of \code{lambda}, \code{theta}, and \code{weights}.
#' 
#' This function implements a penalized EM-algorithm to find the elastic net estimator with weighted l1 norm based on the observed data \code{X}, \code{Y}, \code{Z}, and the sensitivities and specificities in \code{Se}, \code{Sp}.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95) # set individual assay sensitivity
#' Sp <- c(.97) # set individual assay specificity
#' cj <- 1 # set size of master pools
#' # subject individuals to individual testing
#' assay.data <- individual.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' # compute the elastic net estimator with weights = 1
#' enetgt.out <- enetgt(X, Y, Z, Se, Sp, lambda=10, theta=.5, weights = 1) 
#' # compute adaptive elastic net with weights from the elastic net estimator
#' a.enetgt.out <- enetgt(X, Y, Z, Se, Sp, lambda=.5, theta=.5, 
#'                        weights = 1/abs(enetgt.out$b.enet[-1])) 
enetgt <- function( X, Y, Z, Se, Sp, lambda, theta, weights = 1, binit = 1, delta = 1e-3, E.approx = FALSE, get.SEs = FALSE)
{
	
	p <- ncol(X) - 1 # The first column of X is ones
	
	if(length(binit) == 1)
	{ 
		
		b1 <- c(-2,rep(0,p))

	} else {

		b1 <- binit
	}
	
	if(length(weights)==1)
	{
			weights <- rep(1,p)	
	} 
	
	get.EY <- eval(parse(text = ifelse(E.approx, "EY.approx","EY.exact")))
	
	max.diff <- 1
	iter <- 1
		
	while(  max.diff > delta )
	{
		b0 <- b1
		EY <- get.EY(Z,Y,X,b0,Se,Sp)
		b1 <- logistic_enet( EY, X, lambda, weights, theta , b0, delta )$b
		max.diff <- max(abs(b1 - b0))
		iter <- iter + 1
		if(is.na(max.diff)) break;
	}
		
	b.enet <- b1
	
	if(get.SEs)
	{
		# Get observed information matrix by Louis' Method:
		px <- as.numeric(logit(b.enet,X))
		CovYiYj <- CovYiYj.approx(Z,Y,X,b.enet,Se,Sp,EY)

		Hess <- t(X) %*% (  - diag(px * (1 - px))  +  CovYiYj  ) %*% X
	
		Sigma <- lambda * diag(c(0,  theta * weights / abs(b.enet[-1]) + (1-theta) ))
		nz <- which(diag(Sigma)!=Inf) # get nonzero elements of b.enet
		b.nz.enet.cov.est <- solve( - Hess[nz,nz] + Sigma[nz,nz] ) %*% ( - Hess[nz,nz] ) %*% solve( - Hess[nz,nz] + Sigma[nz,nz] )
		b.nz.enet.se <- sqrt(abs(diag(b.nz.enet.cov.est)))*sign(diag(b.nz.enet.cov.est))		
		
		b.enet.se <- numeric(p+1)
		b.enet.se[nz] <- b.nz.enet.se
	
	} else {
		
		b.enet.se <- rep( NA , p + 1 )
	}
		
	conv <- ifelse(is.na(max.diff),FALSE,TRUE)
		
	output <- list(	b.enet = b.enet,
					b.enet.se = b.enet.se,
					lambda = lambda,
					theta = theta,
					weights = weights,
					delta = delta,
					Se = Se,
					Sp = Sp,
					conv = conv,
					iter = iter # iterations to the enet estimator, not to aenet estimator
					)
					
	return(output)

}

#' Computes the elastic net estimators with weighted l1 norm on group testing data over a grid of lambda and theta values.
#'
#' @param X Design matrix with first column a column of 1s.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param n.lambda Number of lambda values for which to compute the estimator.
#' @param n.theta Number of theta values for which to compute the estimator.
#' @param weights Vector of weights to be used in weighting the l1 penalty. Default is \code{weights=1}, which causes equal weights to be used for each coefficient.
#' @param delta Convergence criterion.
#' @param E.approx Logical.  If \code{TRUE} then E-step done with \code{EY.approx()}. If \code{FALSE}, then E-step done with \code{EY.exact()}.
#' @param verbose Logical. If \code{TRUE} then progress is reported after computation of the estimator at each tuning parameter pair.
#' @param get.SEs Logical.  If \code{TRUE} then estimated standard errors for the estimates are returned.
#' @param ridge.include Logical.  If \code{TRUE} then the sequence of theta values begins with 0 to include the ridge estimator.
#' @return A list of output which includes an array \code{B.ENET} of solutions over the grid of lambda and theta values.
#'
#' This function implements a penalized EM-algorithm to find the elastic net estimator with weighted l1 norm based on the observed data \code{X},
#' \code{Y}, \code{Z}, and the sensitivities and specificities in \code{Se}, \code{Sp} over a grid of lambda and theta values.  Only the size of 
#' the grid must be given, and the function chooses its own set of lambda and theta values.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95) # set individual assay sensitivity
#' Sp <- c(.97) # set individual assay specificity
#' cj <- 1 # set size of master pools
#' # subject individuals to individual testing
#' assay.data <- individual.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' n.lambda <- 10
#' n.theta <- 2
#' # compute the elastic net estimator with weights = 1 over grid of lambda and theta values
#' enetgt.grid.out <- enetgt.grid(X, Y, Z, Se, Sp, n.lambda, n.theta, weights = 1) 
enetgt.grid <-function(X, Y, Z, Se, Sp, n.lambda = 5, n.theta = 3, weights = 1, delta = 0.001, E.approx = FALSE, verbose = FALSE, get.SEs = FALSE, ridge.include = FALSE) 
{
		
	if (length(weights) == 1) {
		
       weights <- rep(1, ncol(X) - 1)
    } 
    
    # keep covariates for which the weight is finite
    
    keep <- which(weights<Inf)
    weights.keep <- weights[keep]
    X.keep <- as.matrix(X[,c(1,1+keep)])
                      	
    n <- nrow(X.keep)
    p.keep <- ncol(X.keep) - 1
   	get.EY <- eval(parse(text = ifelse(E.approx, "EY.approx", "EY.exact")))

    b.keep <- c(-2, rep(0.01, p.keep))
    max.diff <- 1
    
    if(ridge.include)
    {
    	theta.seq <- c(0, 1/2^((n.theta - 2):0))
   	    theta.temp <- theta.seq[2]
   	    
    } else {
    	
    	theta.seq <- c(1/2^((n.theta - 1):0))
   	    theta.temp <- theta.seq[1]
    }

    
    D <- pull.diagnoses(Z, Y)
    
    if (length(D) == 0) 
    {
    	iter <- 1
        while ( (sum(b.keep != 0) > 1) & (iter < 25)) 
        {
    
            EY <- get.EY(Z, Y, X.keep, b.keep, Se, Sp)
            lambda <- max(abs(t(EY - mean(EY) * (1 - mean(EY))) %*% X.keep[, -1])/(theta.temp * weights.keep))
            b.keep <- enetgt(X.keep, Y, Z, Se, Sp, lambda, theta.temp, weights.keep, b.keep, delta = 0.01, E.approx)$b.enet
            iter <- iter + 1
        
        }
        
    }	else {
    	
        lambda <- max(abs(t(D - mean(D) * (1 - mean(D))) %*% X.keep[, -1])/theta.temp * weights.keep)
        b.keep <- enetgt(X.keep, Y, Z, Se, Sp, lambda, theta.temp, weights.keep, b.keep, delta = 0.01, E.approx)$b.enet
        
    }
    
    b.enet.keep <- b.keep
    
    lambda.max <- lambda
    lambda.min <- ifelse(n >= p.keep, 0.01, 0.001 * lambda)
    lambda.seq <- exp(log(lambda.min) + (n.lambda:1)/n.lambda * ((log(lambda.max) - log(lambda.min))))
    
    B.ENET.keep <- B.ENET.SE.keep <- array(NA, dim = c(n.lambda, p.keep + 1, n.theta))
    ITER <- matrix(NA, n.lambda, n.theta)
    
    for (t in 1:n.lambda) 
    	for (s in 1:n.theta) 
    	{
    		
	       	b.keep <- b.enet.keep
	        enetgt.out.keep <- enetgt(X.keep, Y, Z, Se, Sp, lambda.seq[t], theta.seq[s], weights.keep, b.keep, delta, E.approx, get.SEs=get.SEs)
	        
	        b.enet.keep <- enetgt.out.keep$b.enet
	        b.enet.se.keep <- enetgt.out.keep$b.enet.se
	         
	        B.ENET.keep[t, , s] <- b.enet.keep
	        B.ENET.SE.keep[t, , s] <- b.enet.se.keep	        
	        
	        ITER[t, s] <- enetgt.out.keep$iter
	        
	        if (verbose == TRUE) {print(paste("solution found for lambda, theta pair (", t, ",", s, ")"))}
	        
    	}
    
    p <- ncol(X) - 1
    B.ENET <- B.ENET.SE <- array(0, dim = c(n.lambda, p + 1, n.theta))
	B.ENET[,c(1,1+keep),] <- B.ENET.keep
	B.ENET.SE[,c(1,1+keep),] <- B.ENET.SE.keep	
	
    output <- list(	B.ENET = B.ENET, 
    				B.ENET.SE = B.ENET.SE,
    				lambda.seq = lambda.seq, 
        			theta.seq = theta.seq, 
        			delta = delta, 
        			E.approx = E.approx, 
        			Se = Se, 
        			Sp = Sp, 
        			n = n, 
        			p = p, 
        			ITER = ITER)
    
    return(output)
    
}

# #' Computes the elastic net estimators with weighted l1 norm on group testing data over a user-provided grid of lambda and theta values.
# #'
# #' @param X Design matrix with first column a column of 1s.
# #' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
# #' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
# #' @param Se A vector of testing sensitivities, where the first element is the
# #'      testing specificity for pools and the second entry is the 
# #'      test specificity for individual testing, if applicable.
# #' @param Sp A vector of testing specificities, where the first element is the
# #'      testing specificity for pools and the second entry is the 
# #'      test specificity for individual testing, if applicable.
# #' @param lambda.seq Sequence of lambda values at which the solution path will be computed.
# #' @param theta.seq Sequence of theta values at which the solution path will be computed.
# #' @param weights Vector of weights to be used in weighting the l1 penalty. Default is \code{weights=1}, which causes equal weights to be used for each coefficient.
# #' @param B.INIT Initial values for EM-algorithm.
# #' @param delta Convergence criterion.
# #' @param E.approx Logical.  If \code{TRUE} then E-step done with \code{EY.approx()}. If \code{FALSE}, then E-step done with \code{EY.exact()}.
# #' @param verbose Logical. If \code{TRUE} then progress is reported after computation of the estimator at each tuning parameter pair.
# #' @param get.SEs Logical.  If \code{TRUE} then estimated standard errors for the estimates are returned.
# #' @return A list of output which includes an array \code{B.ENET} of solutions over the grid of lambda and theta values.
# #'
# #' This function implements a penalized EM-algorithm to find the elastic net estimator with weighted l1 norm based on the observed data \code{X},
# #' \code{Y}, \code{Z}, and the sensitivities and specificities in \code{Se}, \code{Sp} over the grid of lambda and theta values specified in 
# #' \code{lambda.seq} and \code{theta.seq}.
enetgt.grid.0 <- function( X, Y, Z, Se, Sp, lambda.seq, theta.seq, weights = 1, B.INIT, delta = 1e-3, E.approx = FALSE, verbose = FALSE, get.SEs=FALSE)
{
	
	if (length(weights) == 1) {
		
       weights <- rep(1, ncol(X) - 1)
    } 
    
    # keep covariates for which the weight is finite
    
    keep <- which(weights<Inf)
    weights.keep <- weights[keep]
    X.keep <- X[,c(1,1+keep)]
                  
    if(ncol(X.keep)==1)####
    {
    	
    	print("no covariates in model: all weights equal to infinity")
      
	    output <- NA
      
	    return(output)
    	    	
    } else {

		n <- nrow(X.keep)
		p.keep <- ncol(X.keep) - 1
		
		get.EY <- eval(parse(text = ifelse(E.approx, "EY.approx","EY.exact")))
		
		n.lambda <- length(lambda.seq)
		n.theta <- length(theta.seq)
		
		B.ENET.keep <- array(NA,dim=c(n.lambda,p.keep+1,n.theta))
		B.INIT.keep <- B.INIT[,c(1,1+keep),]
		
		ITER <- matrix(NA,n.lambda,n.theta)
		
		for( t in 1:n.lambda)
			for(s in 1:n.theta)
			{
			
				b.keep <- B.INIT.keep[t,,s]
				
				enetgt.out.keep <- enetgt( X.keep, Y, Z, Se, Sp, lambda.seq[t], theta.seq[s], weights.keep, b.keep, delta, E.approx, get.SEs=get.SEs)
									
				b.enet.keep <- enetgt.out.keep$b.enet
										
				B.ENET.keep[t,,s] <- b.enet.keep
	
				ITER[t,s] <- enetgt.out.keep$iter
				
				if(verbose == TRUE)
				{
					
					print(paste("solution found for lambda, theta pair (", t,",",s,")"))
					
				}
			
			}
			
			
	    p <- ncol(X) - 1
	    B.ENET <-  array(0, dim = c(n.lambda, p + 1, n.theta))
		B.ENET[,c(1,1+keep),] <- B.ENET.keep
			
		output <- list(	B.ENET = B.ENET,
						lambda.seq = lambda.seq,
						theta.seq = theta.seq,
						delta = delta,
						E.approx = E.approx,
						Se = Se,
						Sp = Sp,
						n = n,
						p = p,
						ITER = ITER
						)
		
	  }
					
	return(output)

}


####################################################################################################################################
####################################################################################################################################
####### CROSSVALIDATION FUNCTIONS
####################################################################################################################################
####################################################################################################################################

#' Splits individual testing data into crossvalidation data sets.
#'
#' @param X Design matrix with first column a column of 1s.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param K The number of crossvalidation folds.
#' @return List of length \code{K} of the crossvalidation training and testing data sets.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95) # set individual assay sensitivity
#' Sp <- c(.97) # set individual assay specificity
#' cj <- 1 # set size of master pools
#' # subject individuals to individual testing
#' assay.data <- individual.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' individual.cv.fold.data <- get.individual.cv.fold.data(X,Y,Z,K=5)
get.individual.cv.fold.data <- function(X,Y,Z,K)
{ # this is like masterpool testing with a pool size of 1
	
	Z.pools <- Z # make Z.pools the same as Z, individuals are tested in a pool size of 1
	cj.vec <- Z.pools[,2]
	num.pools <- nrow(Z.pools) # 

	num.pools.in.folds <- c( rep(floor(num.pools/K),K-1), num.pools - floor(num.pools/K)*(K-1) ) # put extra pools in the last fold
		
	individual.cv.fold.data <- list()
	
	Z.rows.before.fold <- numeric(0)
	ind.before.fold <- numeric(0)
		
	for(k in 1:K)
	{
		
		
		pools.fold <- (c(0,cumsum(num.pools.in.folds))[k]+1):cumsum(num.pools.in.folds)[k]
		ind.fold <-  unique(sort(as.vector(Z.pools[pools.fold,-c(1:3)])))
		ind.fold <- ind.fold[ind.fold!=-99]
							
		Z.rows.fold <- unique(sort(as.vector(Y[ind.fold,-c(1:2)])))
		Z.rows.fold <- Z.rows.fold[Z.rows.fold!=-99]
							
		Z.rows.after.fold <- (1:nrow(Z))[-c(Z.rows.before.fold,Z.rows.fold)]
		ind.after.fold	<- (1:nrow(Y))[-c(ind.before.fold,ind.fold)]	
		
		# build training data
				
		X.train.fold <- X[-ind.fold,]
		
		temp <- Y
		temp[ind.after.fold,-c(1,2)] <- Y[ind.after.fold,-c(1,2)] - length(Z.rows.fold)	# adjust pool IDs for fold
		Y.train.fold <- temp[-ind.fold,]
		
		temp <- Z
		temp[Z.rows.after.fold,-c(1:3)] <- Z[Z.rows.after.fold,-c(1:3)] - length(ind.fold) # adjust pool IDs for fold
		Z.train.fold <- temp[-Z.rows.fold,]
				
		# build testing data
				
		X.test.fold <- X[ind.fold,]
		
		temp <- Y
		temp[ind.fold,-c(1,2)] <- Y[ind.fold,-c(1,2)] - length(Z.rows.before.fold)
		Y.test.fold <- temp[ind.fold,]
		
		temp <- Z
		temp[Z.rows.fold,-c(1:3)] <- Z[Z.rows.fold,-c(1:3)] - length(ind.before.fold)
		Z.test.fold <- temp[Z.rows.fold,]
		
		Z.rows.before.fold <- c(Z.rows.before.fold, Z.rows.fold)
		ind.before.fold <- c(ind.before.fold,ind.fold)
	
		individual.cv.fold.data[[length(individual.cv.fold.data) + 1]] <- list( X.train = X.train.fold,
																				Y.train = Y.train.fold,
																				Z.train = Z.train.fold,
																				X.test = X.test.fold,
																				Y.test = Y.test.fold,
																				Z.test = Z.test.fold,																													num.pools.removed = length(pools.fold),
																				num.pools = num.pools
																				)
	
	}
	
	return(individual.cv.fold.data)
	
}

#' Splits masterpool testing data into crossvalidation data sets
#'
#' @param X Design matrix with first column a column of 1s.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param K The number of crossvalidation folds.
#' @return List of length \code{K} of the crossvalidation training and testing data sets.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.95) # set master pool assay sensitivity
#' Sp <- c(.97) # set master pool assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to master pool testing
#' assay.data <- masterpool.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' masterpool.cv.fold.data <- get.masterpool.cv.fold.data(X,Y,Z,K=5)
get.masterpool.cv.fold.data <- function(X,Y,Z,K)
{
	
	Z.pools <- Z[Z[,2]>1,] # Z.pools will be the same as Z for master pool testing
	cj.vec <- Z.pools[,2]
	num.pools <- nrow(Z.pools) # same as nrow(Z) for master pool testing

	num.pools.in.folds <- c( rep(floor(num.pools/K),K-1), num.pools - floor(num.pools/K)*(K-1) ) # put extra pools in the last fold
		
	masterpool.cv.fold.data <- list()
	
	Z.rows.before.fold <- numeric(0)
	ind.before.fold <- numeric(0)
	
	for(k in 1:K)
	{
		
		pools.fold <- (c(0,cumsum(num.pools.in.folds))[k]+1):cumsum(num.pools.in.folds)[k]
		ind.fold <-  unique(sort(as.vector(Z.pools[pools.fold,-c(1:3)])))
		ind.fold <- ind.fold[ind.fold!=-99]
							
		Z.rows.fold <- unique(sort(as.vector(Y[ind.fold,-c(1:2)])))
		Z.rows.fold <- Z.rows.fold[Z.rows.fold!=-99]
							
		Z.rows.after.fold <- (1:nrow(Z))[-c(Z.rows.before.fold,Z.rows.fold)]
		ind.after.fold	<- (1:nrow(Y))[-c(ind.before.fold,ind.fold)]	
		
		# build training data
				
		X.train.fold <- X[-ind.fold,]
		
		temp <- Y
		temp[ind.after.fold,-c(1,2)] <- Y[ind.after.fold,-c(1,2)] - length(Z.rows.fold)	# adjust pool IDs for fold
		Y.train.fold <- temp[-ind.fold,]
		
		temp <- Z
		temp[Z.rows.after.fold,-c(1:3)] <- Z[Z.rows.after.fold,-c(1:3)] - length(ind.fold) # adjust pool IDs for fold
		Z.train.fold <- temp[-Z.rows.fold,]
				
		# build testing data
				
		X.test.fold <- X[ind.fold,]
		
		temp <- Y
		temp[ind.fold,-c(1,2)] <- Y[ind.fold,-c(1,2)] - length(Z.rows.before.fold)
		Y.test.fold <- temp[ind.fold,]
		
		temp <- Z
		temp[Z.rows.fold,-c(1:3)] <- Z[Z.rows.fold,-c(1:3)] - length(ind.before.fold)
		Z.test.fold <- temp[Z.rows.fold,]
		
		Z.rows.before.fold <- c(Z.rows.before.fold, Z.rows.fold)
		ind.before.fold <- c(ind.before.fold,ind.fold)
	
		masterpool.cv.fold.data[[length(masterpool.cv.fold.data) + 1]] <- list( X.train = X.train.fold,
																				Y.train = Y.train.fold,
																				Z.train = Z.train.fold,
																				X.test = X.test.fold,
																				Y.test = Y.test.fold,
																				Z.test = Z.test.fold,																		
																				num.pools.removed = length(pools.fold),
																				num.pools = num.pools
																				)
	
	}
	
	return(masterpool.cv.fold.data)
	
}

#' Splits Dorfman testing data into crossvalidation data sets.
#'
#' @param X Design matrix with first column a column of 1s.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param K The number of crossvalidation folds.
#' @return List of length \code{K} of the crossvalidation training and testing data sets.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.92,.95) # set master pool and individual assay sensitivity
#' Sp <- c(.97,.98) # set master pool and individual assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to Dorfman testing
#' assay.data <- dorfman.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' dorfman.cv.fold.data <- get.dorfman.cv.fold.data(X,Y,Z,K=5)
get.dorfman.cv.fold.data <- function(X,Y,Z,K)
{
			
	Z.pools <- Z[Z[,2]>1,]
	cj.vec <- Z.pools[,2]
	num.pools <- nrow(Z.pools)

	num.pools.in.folds <- c( rep(floor(num.pools/K),K-1), num.pools - floor(num.pools/K)*(K-1) ) # put extra pools in the last fold
		
	dorfman.cv.fold.data <- list()
	
	Z.rows.before.fold <- numeric(0)
	ind.before.fold <- numeric(0)

	for(k in 1:K)
	{
		
		pools.fold <- (c(0,cumsum(num.pools.in.folds))[k]+1):cumsum(num.pools.in.folds)[k]
		ind.fold <-  unique(sort(as.vector(Z.pools[pools.fold,-c(1:3)])))
		ind.fold <- ind.fold[ind.fold!=-99]
							
		Z.rows.fold <- unique(sort(as.vector(Y[ind.fold,-c(1:2)])))
		Z.rows.fold <- Z.rows.fold[Z.rows.fold!=-99]
							
		Z.rows.after.fold <- (1:nrow(Z))[-c(Z.rows.before.fold,Z.rows.fold)]
		ind.after.fold	<- (1:nrow(Y))[-c(ind.before.fold,ind.fold)]	
		
		# build training data
				
		X.train.fold <- X[-ind.fold,]
		
		temp <- Y
		temp[ind.after.fold,-c(1,2)] <- Y[ind.after.fold,-c(1,2)] - length(Z.rows.fold)	# adjust pool IDs for fold
		Y.train.fold <- temp[-ind.fold,]
		
		temp <- Z
		temp[Z.rows.after.fold,-c(1:3)] <- Z[Z.rows.after.fold,-c(1:3)] - length(ind.fold) # adjust pool IDs for fold
		Z.train.fold <- temp[-Z.rows.fold,]
				
		# build testing data
				
		X.test.fold <- X[ind.fold,]
		
		temp <- Y
		temp[ind.fold,-c(1,2)] <- Y[ind.fold,-c(1,2)] - length(Z.rows.before.fold)
		Y.test.fold <- temp[ind.fold,]
		
		temp <- Z
		temp[Z.rows.fold,-c(1:3)] <- Z[Z.rows.fold,-c(1:3)] - length(ind.before.fold)
		Z.test.fold <- temp[Z.rows.fold,]
		
		Z.rows.before.fold <- c(Z.rows.before.fold, Z.rows.fold)
		ind.before.fold <- c(ind.before.fold,ind.fold)
			
		dorfman.cv.fold.data[[length(dorfman.cv.fold.data) + 1]] <- list( 	X.train = X.train.fold,
																			Y.train = Y.train.fold,
																			Z.train = Z.train.fold,
																			X.test = X.test.fold,
																			Y.test = Y.test.fold,
																			Z.test = Z.test.fold,
																			num.pools.removed = length(pools.fold),
																			num.pools = num.pools
																			)
	}
	
	return(dorfman.cv.fold.data)
	
}


#' Splits array testing data into crossvalidation data sets.
#'
#' @param X Design matrix with first column a column of 1s.
#' @param Y Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param K The number of crossvalidation folds; \code{K} may be changed internally if the number of arrays is not divisible by \code{K}.
#' @return List of length \code{K} of the crossvalidation training and testing data sets.
#'
#' @examples
#' # generate individual covariate values and disease statuses
#' N <- 160
#' data <- model1(N)
#' X <- data$X
#' Y.true <- data$Yi
#' Se <- c(.92,.95) # set row/col and individual assay sensitivity
#' Sp <- c(.97,.98) # set row/col and individual assay specificity
#' cj <- 4 # set size of master pools
#' # subject individuals to array testing
#' assay.data <- array.assay.gen(Y.true,Se,Sp,cj)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' array.cv.fold.data <- get.array.cv.fold.data(X,Y,Z,K=5)
get.array.cv.fold.data <- function(X,Y,Z,K)
{
	
	pools.before.fold <- numeric(0)
	individuals.before.fold <- numeric(0)
	
	n <- nrow(X)
	cj <- dim(Z)[2]-3
	n.arrays <- floor(n / cj^2)
	
	
	if(n.arrays %% K != 0) # if number of arrays not divisible by K take arrays away
	{
	
		K.new <- ifelse(sum(n.arrays %% K:7 == 0)==0,Inf , min((K:7)[(n.arrays %% K:7 == 0)]))
				
		if(K.new == Inf){
			
			n.keep <- floor(n / (K * cj^2)) * K*cj^2
			
			X.keep <- X[1:n.keep,]
			Z.keep <- Z[  Z[,4]<=n.keep , ]
			Y.keep <- Y[1:n.keep,]			
			
			X <- X.keep
			Z <- Z.keep
			Y <- Y.keep
			
			n.arrays <- n.keep / cj^2
			
		} else {
			
			K <- K.new
			
			print(paste("K changed to ",K," so that the number of arrays (",n.arrays,") is divisible by K",sep=""))
		
		}
	}
			
	array.cv.fold.data <- list()
	
	for(k in 1:K)
	{
		
		arrays.fold.k <- ((n.arrays/K)*(k-1) + 1):((n.arrays/K)*k)
		fold.k <-  ((min(arrays.fold.k) - 1 )* cj^2 + 1):(max(arrays.fold.k) * cj^2 )
	
		individuals.after.fold <- (1:nrow(Y))[-c(individuals.before.fold,fold.k)]
	
		pools.fold.k <- (1:nrow(Z))[Z[,4] %in% fold.k]
		pools.after.fold <- (1:nrow(Z))[-c(pools.before.fold,pools.fold.k)]
		
		# build training data
				
		X.train.fold <- X[-fold.k,]
		
		Y.temp <- Y
		Y.temp[individuals.after.fold,3:5] <- Y[individuals.after.fold,3:5] - length(pools.fold.k)	# adjust pool IDs for fold
		Y.train.fold <- Y.temp[-fold.k,]
		
		Z.temp <- Z
		Z.temp[pools.after.fold,4:(3 + cj)] <- Z[pools.after.fold,4:(3 + cj)] - (n.arrays/K)*cj^2 # adjust pool IDs for fold
		Z.train.fold <- Z.temp[-pools.fold.k,]
		
		# build testing data
		
		X.test.fold <- X[fold.k,]
							
		Y.temp <- Y
		Y.temp[fold.k,3:ncol(Y)] <- Y[fold.k,3:ncol(Y)] - length(pools.before.fold)
		Y.test.fold <- Y.temp[fold.k,]
		
		Z.temp <- Z
		Z.temp[pools.fold.k,4:(3 + cj)] <- Z[pools.fold.k,4:(3 + cj)] - (k-1)*(n.arrays/K)*cj^2
		Z.test.fold <- Z.temp[pools.fold.k,]

		# augment some vectors of indices
		
		pools.before.fold <- c(pools.before.fold,pools.fold.k)
		individuals.before.fold <- c(individuals.before.fold,fold.k)
	
	
		array.cv.fold.data[[length(array.cv.fold.data) + 1]] <- list( X.train = X.train.fold,
																			Y.train = Y.train.fold,
																			Z.train = Z.train.fold,
																			X.test = X.test.fold,
																			Y.test = Y.test.fold,
																			Z.test = Z.test.fold,
																			num.pools.removed = length(arrays.fold.k),
																			num.pools = n.arrays
																			)
	
	}
	
	return(array.cv.fold.data)
	
}

# #' Splits data which are a mix of individual and Dorfman group testing data into crossvalidation data sets
# #'
# #' @param X.dorf Design matrix fir individuals tested under Dorfman regime, with first column a column of 1s.
# #' @param Y.dorf Group testing output from a Dorfman testing regime.
# #' @param Z.dorf Group testing output from an Dorfman testing regime.
# #' @param X.ind Design matrix fir individuals tested under individual regime, with first column a column of 1s.
# #' @param Y.ind Group testing output from a individual testing regime.
# #' @param Z.ind Group testing output from an individual testing regime.
# #' @param K =  the number of crossvalidation folds.
# #' @return List of length \code{K} of the crossvalidation training and testing data sets.
# #'
# #' Splits data into crossvalidation data sets where \code{X.dorf}, \code{Y.dorf}, and \code{Z.dorf} come from a Dorfman testing regime and 
# #' \code{X.ind}, \code{Y.ind} and \code{Z.ind} come from an individual testing regime.
get.dorfman.individual.cv.fold.data <- function(X.dorf,Y.dorf,Z.dorf,X.ind,Y.ind,Z.ind,K)
{
		
	dorfman.cv.fold.data <- get.dorfman.cv.fold.data(X.dorf,Y.dorf,Z.dorf,K)
	individual.cv.fold.data <- get.individual.cv.fold.data(X.ind,Y.ind,Z.ind,K)	
	
	dorfman.individual.cv.fold.data <- list()
	
	for(k in 1:K)
	{
		
		dorf.fold <- dorfman.cv.fold.data[[k]]
		ind.fold <- individual.cv.fold.data[[k]]

		num.pools.removed <- dorf.fold$num.pools.removed + ind.fold$num.pools.removed
		num.pools <- dorf.fold$num.pools + ind.fold$num.pools

		all.fold.train <- dorfman.individual.combine(	dorf.fold$X.train,dorf.fold$Y.train,dorf.fold$Z.train,
														ind.fold$X.train,ind.fold$Y.train,ind.fold$Z.train)
														
		all.fold.test <- dorfman.individual.combine(	dorf.fold$X.test,dorf.fold$Y.test,dorf.fold$Z.test,
														ind.fold$X.test,ind.fold$Y.test,ind.fold$Z.test)												
														
		dorfman.individual.cv.fold.data[[length(dorfman.individual.cv.fold.data) + 1]] <- list( X.train = all.fold.train$X,
																								Y.train = all.fold.train$Y,
																								Z.train = all.fold.train$Z,
																								X.test = all.fold.test$X,
																								Y.test = all.fold.test$Y,
																								Z.test = all.fold.test$Z,
																								num.pools.removed = num.pools.removed,
																								num.pools = num.pools
																								)
		
	}	
	
	return(dorfman.individual.cv.fold.data)
	
}
#' Peforms crossvalidation to select a combination of tuning parameters.
#'
#' @param cv.fold.data Output from any of the functions \code{get.xx.cv.fold.data}, where \code{xx} is \code{individual}, \code{masterpool}, \code{dorfman}, or \code{array}.
#' @param regime testing regime which must be equal to one of \code{individual}, \code{masterpool}, \code{dorfman}, and \code{array}.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param lambda.seq Sequence of lambda values at which the weighted elastic net estimator will be computed.
#' @param theta.seq Sequence of theta values at which the weighted elastic net estimator will be computed.
#' @param weights Vector of weights to be used in weighting the l1 penalty. Default is \code{weights=1}, which causes equal weights to be used for each coefficient.
#' @param B.INIT Array of values to be used to initialize the EM-algorithm at each of the tuning paramter combinations specified in \code{lambda.seq} and \code{theta.seq}.
#' @param delta Convergence criterion.
#' @param E.approx Logical.  If \code{TRUE} then E-step done with \code{EY.approx()}. If \code{FALSE}, then E-step done with \code{EY.exact()}.
#' @param verbose Logical. If \code{TRUE} then progress is reported after computation of the estimator at each tuning parameter pair.
#' @param plot Logical. If \code{TRUE} then a plot is produced showing the mean of the crossvalidation estimates of the deviance over the null deviance.
#' @return A list which includes the chosen values of the tuning parameters according to crossvalidation.
#'
#' @examples
# generate covariate values and disease statuses for 200 individuals from model0:
#' data <- model0(200)
#' X <- data$X
#' Y.true <- data$Y
#' # subject individuals to individual testing
#' Se <- c(.94) # individual testing sensitivities
#' Sp <- c(.95) # individual testing specificities
#' assay.data <- individual.assay.gen(Y.true,Se,Sp,cj=1)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' # compute the mle on the individual testing data:
#' mlegt.out <- mlegt(X,Y,Z,Se,Sp,delta=.01)
#' b.mle <- mlegt.out$b.mle
#' # compute adaptive elastic net estimator over a grid of tuning parameter values
#' n.lambda <- 8
#' n.theta <- 2
#' enetgt.grid.out <- enetgt.grid(X,Y,Z,Se,Sp,n.lambda,n.theta,weights = 1/abs(b.mle[-1]),delta=.01)
#' # make a choice of the tuning parameter using 3-fold crossvalidation:
#' cv.fold.data <- get.individual.cv.fold.data(X,Y,Z,K=3)
#' cv.enetgt.grid.out <- cv.enetgt.grid(cv.fold.data,"individual",Se,Sp,enetgt.grid.out$lambda.seq,
#'						enetgt.grid.out$theta.seq,weights=1/abs(b.mle[-1]),
#'						B.INIT=enetgt.grid.out$B.ENET,delta=.01)
#' b.aenet.cv <- enetgt.grid.out$B.ENET[cv.enetgt.grid.out$cv.ind[1],,cv.enetgt.grid.out$cv.ind[2]]
#' b.aenet.cv
cv.enetgt.grid <- function( cv.fold.data,regime, Se, Sp, lambda.seq, theta.seq, weights, B.INIT, delta = 1e-2, E.approx = 
FALSE,verbose=FALSE,plot=FALSE)
{
		
	neg.ll <- eval(parse(text = paste("neg.ll",regime,sep=".")))
	
	n.lambda <- length(lambda.seq)
	n.theta <- length(theta.seq)
	K <- length(cv.fold.data)
	
	null.ll.CV.folds.enet <- matrix(NA,K,n.theta)
	ll.CV.folds.enet <- dev.ratio.CV.folds.enet <- array(NA, dim = c(n.lambda,n.theta,K))
	CV.FOLDS.B.ENET <- array(NA, dim = c(n.lambda,ncol(cv.fold.data[[1]]$X.train),n.theta,K))
	
	for(k in 1:K)
	{
		
		fold.k <- cv.fold.data[[k]]
		
		cv.lambda.scale <- (fold.k$num.pools - fold.k$num.pools.removed)/fold.k$num.pools
		
		enetgt.grid.0.out.fold.k <- enetgt.grid.0( 	fold.k$X.train, 
													fold.k$Y.train, 
													fold.k$Z.train,  
													Se, 
													Sp, 
													lambda.seq = cv.lambda.scale*lambda.seq, 
													theta.seq = theta.seq,
													weights = weights,
													B.INIT = B.INIT,
												 	delta = delta,
												 	E.approx = E.approx,
												 	verbose = verbose)
												 	
		# browser()
		
		for(j in 1:n.theta)
		{
			null.ll.CV.folds.enet[k,j] <- neg.ll(enetgt.grid.0.out.fold.k$B.ENET[1,,j],fold.k$Z.test,fold.k$X.test,Se,Sp)
		}				
		
		for(i in 1:n.lambda)
			for(j in 1:n.theta)
			{
				
				b.enet <- enetgt.grid.0.out.fold.k$B.ENET[i,,j]

				CV.FOLDS.B.ENET[i,,j,k] <- b.enet
				ll.CV.folds.enet[i,j,k] <- neg.ll(b.enet,fold.k$Z.test,fold.k$X.test,Se,Sp)
				dev.ratio.CV.folds.enet[i,j,k] <- ll.CV.folds.enet[i,j,k] / null.ll.CV.folds.enet[k,j]
	
			}
			
		print(paste("fold ",k," complete",sep=""))
	
	}
	
	# CV dev ratio analysis:
	
	dev.ratio.CV.enet <- apply(dev.ratio.CV.folds.enet,c(1,2),mean)
	dev.ratio.CV.enet.se <- apply(dev.ratio.CV.folds.enet,c(1,2),sd)/sqrt(K)
	se.at.min <- dev.ratio.CV.enet.se[which.min(dev.ratio.CV.enet[,2]),2]
	
	cv.ind <- which(dev.ratio.CV.enet == min(dev.ratio.CV.enet , na.rm=TRUE), arr.ind=TRUE)[1,]
	cv.ind.se <- dev.ratio.CV.enet.se[cv.ind[1],cv.ind[2]]
	cv.ind.lambda1se <- c(min(which(dev.ratio.CV.enet[,cv.ind[2]] <= (dev.ratio.CV.enet[cv.ind[1],cv.ind[2]] + cv.ind.se))),cv.ind[2])
	
	cv.lambda <- lambda.seq[cv.ind[1]]
	cv.theta <- theta.seq[cv.ind[2]]
	cv.lambda1se <- lambda.seq[cv.ind.lambda1se[1]]
	
	if(plot == TRUE)
	{	
		
		plot(NA,xlim=range(log(lambda.seq)),ylim=range(dev.ratio.CV.enet[,cv.ind[2]]+se.at.min, dev.ratio.CV.enet[,cv.ind[2]] - se.at.min),xlab="log(Lambda)",ylab="dev / nulldev")
		
		points(dev.ratio.CV.enet[,cv.ind[2]]~log(lambda.seq),col="red",pch=19)
		lines(c(dev.ratio.CV.enet[,cv.ind[2]]+se.at.min)~log(lambda.seq),col="red",lty=3)	
		lines(c(dev.ratio.CV.enet[,cv.ind[2]]-se.at.min)~log(lambda.seq),col="red",lty=3)	
	
		abline(v=log(cv.lambda))
		abline(v=log(cv.lambda1se))
		
	}	
		
	output <- list( CV.FOLDS.B.ENET = CV.FOLDS.B.ENET,
					cv.ind = cv.ind,
					cv.lambda = cv.lambda,
					cv.theta = cv.theta,
					cv.ind.lambda1se = cv.ind.lambda1se,
					cv.lambda1se = cv.lambda1se)

	return(output)
	
}

#' Chooses tuning parameters via AIC, BIC, and ERIC.
#'
#' @param enetgt.grid.out List object returned by the function \code{enetgt.grid}.
#' @param Z Group testing output from one of the functions \code{individual.assay.gen}, \code{masterpool.assay.gen}, \code{dorfman.assay.gen}, or \code{array.assay.gen}.
#' @param X Design matrix with first column a column of 1s.
#' @param Se A vector of testing sensitivities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param Sp A vector of testing specificities, where the first element is the
#'      testing specificity for pools and the second entry is the 
#'      test specificity for individual testing, if applicable.
#' @param regime Testing regime which must be equal to one of \code{individual}, \code{masterpool}, \code{dorfman}, \code{array}, and \code{dorfman.individual}.
#' @return A list object including the tuning parameter combinations selected by AIC, BIC, and ERIC.
#'
#' Computes the AIC, BIC and ERIC criteria over a grid of lambda, theta pairs coming from \code{enetgt.grid.out}, 
#' which is an object returned by \code{enetgt.grid()}.
#'
#' @examples
# generate covariate values and disease statuses for 200 individuals from model0:
#' data <- model0(200)
#' X <- data$X
#' Y.true <- data$Y
#' # subject individuals to individual testing
#' Se <- c(.94) # individual testing sensitivities
#' Sp <- c(.95) # individual testing specificities
#' assay.data <- individual.assay.gen(Y.true,Se,Sp,cj=1)
#' Z <- assay.data$Z
#' Y <- assay.data$Y
#' # compute the mle on the individual testing data:
#' mlegt.out <- mlegt(X,Y,Z,Se,Sp,delta=.01)
#' b.mle <- mlegt.out$b.mle
#' # compute adaptive elastic net estimator over a grid of tuning parameter values
#' n.lambda <- 8
#' n.theta <- 2
#' enetgt.grid.out <- enetgt.grid(X,Y,Z,Se,Sp,n.lambda,n.theta,weights = 1/abs(b.mle[-1]),delta=.01)
#' # make choices of the tuning parameters according to the aic, bic, and eric criteria
#' aic.bic.eric.enetgt.grid.out <- aic.bic.eric.enetgt.grid(enetgt.grid.out,Z,X,Se,Sp,"individual")
#' b.aenet.aic <- aic.bic.eric.enetgt.grid.out$b.enet.aic
#' b.aenet.bic <- aic.bic.eric.enetgt.grid.out$b.enet.bic
#' b.aenet.eric <- aic.bic.eric.enetgt.grid.out$b.enet.eric
#' b.aenet.aic 
#' b.aenet.bic
#' b.aenet.eric
aic.bic.eric.enetgt.grid <- function(enetgt.grid.out,Z,X,Se,Sp,regime)
{
		
	neg.ll <- eval(parse(text = paste("neg.ll",regime,sep=".")))
		
	lambda.seq <- enetgt.grid.out$lambda.seq
	theta.seq <-  enetgt.grid.out$theta.seq

	n.lambda <- length(lambda.seq)
	n.theta <- length(theta.seq)
	n <- enetgt.grid.out$n
	
	aic <- bic <- eric <- DF.HAT <- matrix(NA,n.lambda,n.theta)
		
	for(i in 1:n.lambda)
		for(j in 1:n.theta)
		{
			
			b <- enetgt.grid.out$B.ENET[i,,j]
			px <- logit(b,X)			
			W <- diag(px * (1 - px))			
			A <- which(b!=0)

			if(sum(A) > 1) 
			{
				
				XtWX <- t(X[,A[-1]]) %*% W %*% X[,A[-1]] # remove intercept column from X		
				df.hat <- 1 + sum(diag( solve(XtWX + lambda.seq[i]*(1-theta.seq[j])) %*%  XtWX ))
				
			}	else { # if there is only the intercept
				
				df.hat <- 1 
				
			}	
			
			minus2.ll <- 2 * neg.ll(b,Z,X,Se,Sp)	
			aic[i,j] <- minus2.ll + 2 * df.hat
			bic[i,j] <- minus2.ll + df.hat * log(n)
			eric[i,j] <- minus2.ll + 1 * df.hat * log(n/lambda.seq[i]) 	
			DF.HAT[i,j] <- df.hat 			
			
		}
				
	if( sum(is.na(aic))==prod(dim(aic)))
	{
		
		aic.ind <- NA
		aic.lambda <- NA
		aic.theta <- NA
		b.enet.aic <- rep(NA,ncol(X))
		
	} else {
		
		aic.ind <- which(aic == min(aic, na.rm=TRUE), arr.ind=TRUE)[1,]
		aic.lambda <- lambda.seq[aic.ind[1]]
		aic.theta <- theta.seq[aic.ind[2]]	
		b.enet.aic <- enetgt.grid.out$B.ENET[aic.ind[1],,aic.ind[2]]
	
	}
			
	if( sum(is.na(bic))==prod(dim(bic)))
	{
		
		bic.ind <- NA
		bic.lambda <- NA
		bic.theta <- NA
		b.enet.bic <- rep(NA,ncol(X))
		
	} else {
		
		
		bic.ind <- which(bic == min(bic, na.rm=TRUE), arr.ind=TRUE)[1,]
		bic.lambda <- lambda.seq[bic.ind[1]]
		bic.theta <- theta.seq[bic.ind[2]]	
		b.enet.bic <- enetgt.grid.out$B.ENET[bic.ind[1],,bic.ind[2]]
		
	}
		
	if( sum(is.na(eric))==prod(dim(eric)))
	{
		
		eric.ind <- NA
		eric.lambda <- NA
		eric.theta <- NA
		b.enet.eric <- rep(NA,ncol(X))
		
	} else {
			
		eric.ind <- which(eric == min(eric, na.rm=TRUE), arr.ind=TRUE)[1,]
		eric.lambda <- lambda.seq[eric.ind[1]]
		eric.theta <- theta.seq[eric.ind[2]]	
		b.enet.eric <- enetgt.grid.out$B.ENET[eric.ind[1],,eric.ind[2]]
	
	}
	
	output <- list( aic.ind = aic.ind,
					aic.lambda = aic.lambda,
					aic.theta = aic.theta,
					b.enet.aic = b.enet.aic,
					bic.ind = bic.ind,
					bic.lambda = bic.lambda,
					bic.theta = bic.theta,
					b.enet.bic = b.enet.bic,
					eric.ind = eric.ind,
					eric.lambda = eric.lambda,
					eric.theta = eric.theta,
					b.enet.eric = b.enet.eric,
					DF.HAT = DF.HAT)
					
	return(output)
	
}
# #' Combines Dorfman and individual testing data.
# #'
# #' @param X.dorf Design matrix fir individuals tested under Dorfman regime, with first column a column of 1s 
# #' @param Y.dorf Group testing output from a Dorfman testing regime
# #' @param Z.dorf Group testing output from an Dorfman testing regime
# #' @param X.ind Design matrix fir individuals tested under individual regime, with first column a column of 1s 
# #' @param Y.ind Group testing output from a individual testing regime
# #' @param Z.ind Group testing output from an individual testing regime
# #' @return List with items \code{X}, \code{Y}, \code{Z}
# #'
# #' Combines data from individuals tested according to a dorfman
# #' assay scheme and individuals who are individually tested.
dorfman.individual.combine <- function(X.dorf,Y.dorf,Z.dorf,X.ind,Y.ind,Z.ind)
{
	
		
	X.dorf.ind <- rbind(X.dorf,X.ind)
	
	temp <- cbind(Y.ind,rep(-99,nrow(X.ind)))
	temp[,3] <- Y.ind[,3] + max(Y.dorf[,-c(1,2)])
	
	Y.dorf.ind <- rbind(Y.dorf,temp)
	
	temp <- cbind(Z.ind,matrix(-99,nrow(Z.ind),ncol(Z.dorf)-ncol(Z.ind)))
	temp[,4] <- Z.ind[,4] + max(Z.dorf[,-c(1:3)])
	
	Z.dorf.ind <- rbind(Z.dorf,temp)
	
	
	output <- list( X = X.dorf.ind,
					Y = Y.dorf.ind,
					Z = Z.dorf.ind )
	
	return(output)

}
# #' Evaluates the negative log-likelihood for individual testing data.
# #'
# #' @param b Parameter vector.
# #' @param X Design matrix with first column a column of 1s.
# #' @param Z Group testing output from \code{individual.assay.gen}.
# #' @param Se Vector of testing sensitivities of length \code{max(Z[,3])}.
# #' @param Sp Vector of testing specificities of length \code{max(Z[,3])}.
# #' @return The value of the negative log-likelihood for individual testing data.
neg.ll.individual <- function(b,Z,X,Se,Sp)  # this is the exact log-likelihood for individual testing
{

	L <- 0
		
	for(i in 1:nrow(X))
	{
		
		Zi <- Z[i,1]
		test <- Z[i,3]
		ind.i <- (1 - Sp[test])^Zi * Sp[test]^(1 - Zi) * ( 1 - logit(b,X[i,])) + Se[test]^Zi * (1-Se[test])^(1 - Zi) *  logit(b,X[i,]) 
		
		L <- L - log(ind.i)
		
	}
		
	return(L)
	
}	
# #' Evaluates the negative log-likelihood for master pool testing data.
# #'
# #' @param b Parameter vector.
# #' @param X Design matrix with first column a column of 1s.
# #' @param Z Group testing output from \code{masterpool.assay.gen}.
# #' @param Se Vector of testing sensitivities of length \code{max(Z[,3])}.
# #' @param Sp Vector of testing specificities of length \code{max(Z[,3])}.
# #' @return The value of the negative log-likelihood for master pool testing data.
neg.ll.masterpool <- function(b,Z,X,Se,Sp)  # this is the exact log-likelihood for master pool testing.
{
	
	L <- 0
	
	for(j in 1:nrow(Z))
	{
		
		c.j <- Z[j,2]
		ind.j <- Z[j,4:(c.j+3)]
		
		Zj <- Z[j,1]
		
		test <- Z[j,3]
		
		grp.j <- (1 - Sp[test])^Zj * Sp[test]^(1 - Zj) * prod( 1 - logit(b,X[ind.j,]) ) + Se[test]^Zj * (1-Se[test])^(1 - Zj) * ( 1 - prod( 1 - logit(b,X[ind.j,]))) 
	
		L <- L - log(grp.j)
	
	}
	
	return(L)
	
}	

#####################################################################
# R function: A helper function for neg.ll.dorfman().  Computes the 
#             joint probability 
#
# P( Z_j , \{Y_i = y_i, i \in \mathcal{P}_j \}, \{ \tilde Y_i = \tilde y_i, i \in \mathcal{P}_j \} )
#
#             for the purpose of summing over all 
#
#       \{ \tilde Y_i = \tilde y_i, i \in \mathcal{P}_j \}
#
#            This function gets applied to the rows of a matrix in neg.ll.dorfman().
# #' A helper function for \code{neg.ll.dorfman()}
# #'
# #' @param tYji Vector of individual disease statuses
# #' @param Yji Vector of individual assays
# #' @param Se.ind Sensitivy of individual tests
# #' @param Sp.ind Specificity of individual tests
# #' @param Se.pool Sensitivy of pool tests
# #' @param Sp.pool Specificity of pool tests
# #' @param Xmat Design matrix with first column a column of 1s 
# #' @param b parameter vector
# #' @return A quantity needed for the function \code{neg.ll.Dorfman()}
pzj1YjitYji <- function(tYji,Yji,Se.ind,Sp.ind,Se.pool,Sp.pool,Xmat,b)
{
	
	prod1 <- prod(  (Sp.ind^(1-Yji) * (1 - Sp.ind)^ Yji )^(1 - tYji) * ( (Se.ind^ Yji * (1 - Se.ind)^(1-Yji)) )^ tYji )
	
	prod2 <- prod( logit(b,Xmat)^ tYji * ( 1 - logit(b,Xmat))^(1-tYji))
	
	res <- Se.pool ^ max(tYji) * (1 - Sp.pool) ^ (1-max(tYji)) * prod1 * prod2
	
	return(res)
	
}
# #' Evaluates the negative log-likelihood for Dorfman testing data.
# #'
# #' @param b Parameter vector.
# #' @param X Design matrix with first column a column of 1s.
# #' @param Z Group testing output from \code{dorfman.assay.gen}.
# #' @param Se Vector of testing sensitivities of length \code{max(Z[,3])}.
# #' @param Sp Vector of testing specificities of length \code{max(Z[,3])}.
# #' @return The value of the negative log-likelihood for Dorfman testing data.
neg.ll.dorfman <- function(b, Z, X, Se, Sp)
{

	Pneg <- which(Z[,1]==0 & Z[,2] > 1 )
	Ppos <- which(Z[,1]==1 & Z[,2] > 1 )
		
	L <- 0
		
	for(j in Pneg)	
	{
		
		c.j <- Z[j,2]
		
		ind.j <- Z[j,4:(c.j+3)]
		
		pz0 <- Sp[Z[j,3]] * prod( 1 - logit(b,X[ind.j,]) ) + ( 1 - Se[Z[j,3]]) * ( 1 - prod(  1 - logit(b,X[ind.j,])  ))
		
		L <- L - log(pz0)	
		
	}		
	
	for(j in Ppos)
	{
		
		c.j <- Z[j,2]
		
		ind.j <- Z[j,4:(c.j+3)]
		
		Yji <- Z[(Z[,2]==1 & Z[,4] %in% ind.j),1]
		
		Se.ind <- Se[Z[(Z[,2]==1 & Z[,4] %in% ind.j),3]]
		Sp.ind <- Sp[Z[(Z[,2]==1 & Z[,4] %in% ind.j),3]]
		
		Se.pool <- Se[Z[j,3]]
		Sp.pool <- Sp[Z[j,3]]
		
		Yj.grid <- expand.grid(rep(list(0:1),c.j))
		
		pz1Yji <- sum(apply(X = Yj.grid, MARGIN=1, FUN = pzj1YjitYji, Yji = Yji, Se.ind = Se.ind,
							Sp.ind = Sp.ind ,Se.pool = Se.pool,Sp.pool=Sp.pool, Xmat = X[ind.j,], b = b ))
		
		L <- L - log(pz1Yji)
				
	}

	return(L)
	
}
# #' Gets a Monte Carlo approximation to the negative log-likelihood for array testing data.
# #'
# #' @param b Parameter vector.
# #' @param X Design matrix with first column a column of 1s.
# #' @param Z Group testing output from \code{array.assay.gen}.
# #' @param Se Vector of testing sensitivities of length \code{max(Z[,3])}.
# #' @param Sp Vector of testing specificities of length \code{max(Z[,3])}.
# #' @param B Number of Monte Carlo draws.
# #' @return Approximate value of the negative log-likelihood for array testing data.
neg.ll.array <- function(b, Z, X, Se, Sp, B = 10000)
{

	n <- nrow(X)
	
	cj <- ncol(Z[,-c(1:3)]) # extract the array dimension, which we assume to be the same for all arrays
	J <- floor(n / cj^2)

	cjsq <- cj^2

	ll <- numeric()
	
	for(j in 1:J)
	{
		
		ind.j <- ((j-1)*cjsq+1):(j*cjsq)
		which.rows.of.Z.j <- (1:nrow(Z))[Z[,4] %in% ind.j]	
		Zj <-Z[which.rows.of.Z.j,]
	
		Zjr <- Zj[1:cj,1]
		Zjc <- Zj[(cj+1):(2*cj),1]
		Yji <- Zj[-c(1:(2*cj)),1]
		whichjretest <- Zj[-c(1:(2*cj)),4] - (j-1)*cj^2
	
		pxji <- logit(b,X[ind.j,])	
		
		ll[j] <- llj_array(Zjr, Zjc, Yji, whichjretest, pxji,  Se , Sp ,B = B)$llj
		
	 	# print(j)
	 	
	}
		
	neg.ll.approx <- - sum(ll)

	# output <- list( neg.ll.approx = neg.ll.approx,
					# ll = ll)	

	return(neg.ll.approx)
	
	
}
# #' Computes the negative log-likelihood when data from \code{dorfman.individual.combine()}
# #' 
# #' @param b parameter vector
# #' @param X Design matrix coming from \code{dorfman.individual.combine()}
# #' @param Z Group testing output coming from \code{dorfman.individual.combine()}
# #' @param Se Vector of testing sensitivities of length 2
# #' @param Sp Vector of testing specificities of length 2
# #' @return The value of the negative log-likelihood
# #'
# #' This function computes the negative log-likelihood for a set of Dorfman data and a set of individual-testing data
# #' together.  Needed for Iowa data. Note: it is assumed that the individual testing sensitivity and specificity 
# #' are the second elements of Se and Sp, respectively.
neg.ll.dorfman.individual  <- function(b,Z,X,Se,Sp)
{
	# need dorfman assays to come in the first rows and individual assay data afterwards
	max.id.dorf <- max(Z[Z[,2]==4,4:ncol(Z)])
	Z.dorf <- Z[Z[,4]<= max.id.dorf,]
	X.dorf <- X[1:max.id.dorf,]
	Z.ind <- Z[Z[,4]> max.id.dorf,]
	X.ind <- X[-c(1:max.id.dorf),]
	
	neg.ll <- neg.ll.dorfman(b,Z.dorf,X.dorf,Se,Sp) + neg.ll.individual(b,Z.ind,X.ind,Se,Sp)
	
	return(neg.ll)
	
}
#' Compute probabilities based on the logit link.
#'
#' @param b Parameter vector.
#' @param X Design matrix with first column a column of 1s.
#' @return The fitted probability according to the logit link function.
logit<-function(b,X){
p<-as.numeric(exp(X%*%b)/(1+exp(X%*%b)))
return(p)
}

# #' Computes probabilities based on the logit link
# #'
# #' @param b
# #' @param X
# #' @returns The fitted probability according to the logit link function
# Lam.ii<-function(bk,X){
# p<-logit(bk,X)
# res<-p*(1-p)
# return(as.vector(res))
# }
# Z.i<-function(y,X,bk){
# Xb<- X %*% bk
# p<-logit(bk,X)
# Z<-Xb+(y-p)/(p*(1-p))
# return(as.vector(Z))
# }

#' Generates data from model0.
#'
#' @param n Number of individuals for whom to generate data.
#' @return List with design matrix \code{X}, disease statuses \code{Yi}, model identifier \code{model}, and parameter vector \code{b}.
#'
#' This is a small model which is not really interesting but on which everything runs fast, because of its smallness.
#'
#' @examples
#' N <- 100
#' data <- model0(N)
model0 <- function(n)
{

	p <- 3
	b <- c(-1,3,1.5,0) # first entry is the intercept
	
	C <- matrix(.75,p,p) + diag(rep(.25,p))
	
	X <- cbind(1,scale(matrix(rnorm(p*n),nrow=n)) %*% chol(C))
	px <- logit(b,X)
	Yi <- rbinom(n,1,px)

	output <- list(	X = X,
					Yi = Yi,
					model = "model0",
					b = b)
					
	return(output)
	
}

#' Generates data from model1.
#'
#' @param n Number of individuals for whom to generate data.
#' @return List with design matrix \code{X}, disease statuses \code{Yi}, model identifier \code{model}, and parameter vector \code{b}.
#'
#' @examples
#' N <- 100
#' data <- model1(N)
model1 <- function(n)
{
    p <- 10
    b <- c(-4, 2/3, 1/3, 1, 0, 0, 0, 0, 0, 0, 0)
    C <- matrix(0.5, p, p) + diag(rep(0.5, p))
    X <- cbind(1, scale(matrix(rnorm(p * n), nrow = n)) %*% chol(C))
    px <- logit(b, X)
    Yi <- rbinom(n, 1, px)
    output <- list(X = X, Yi = Yi, model = "model1", b = b)
    return(output)
}

#' Generates data from model2.
#'
#' @param n Number of individuals for whom to generate data.
#' @return List with design matrix \code{X}, disease statuses \code{Yi}, model identifier \code{model}, and parameter vector \code{b}.
#'
#' @examples
#' N <- 100
#' data <- model2(N)
model2 <- function(n) 
{
    p <- 18
    b <- c(-4, -1, -0.75,-0.5, 0.5,0.75,1,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0)
    C <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
    X <- cbind(1, scale(matrix(rnorm(p * n), nrow = n)) %*% chol(C))
    px <- logit(b, X)
    Yi <- rbinom(n, 1, px)
    output <- list(X = X, Yi = Yi, model = "model2", b = b)
    return(output)
}
#' Generates data from model3.
#'
#' @param n Number of individuals for whom to generate data.
#' @return List with design matrix \code{X}, disease statuses \code{Yi}, model identifier \code{model}, and parameter vector \code{b}.
#'
#' @examples
#' N <- 100
#' data <- model3(N)
model3 <- function(n)
{

	p <- 24
    b <- c(-3, rep(c(0.5, 0, 0), 8))
    C <- diag(rep(1, 8)) %x% 0.5^abs(outer(1:3, 1:3, FUN = "-"))
    X <- cbind(1, scale(matrix(rnorm(p * n), nrow = n)) %*% chol(C))
    px <- logit(b, X)
    Yi <- rbinom(n, 1, px)
    output <- list(X = X, Yi = Yi, model = "model3", b = b)
    return(output)
	
}



logistic_enet = function(EY, X, lambda, weights, theta, b0, delta){
  
  return( logistic_enet(EY, X, lambda, weights, theta, b0, 
                        delta) )
  
}
