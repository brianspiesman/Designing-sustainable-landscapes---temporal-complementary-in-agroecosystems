rep = 1
T = 3000							#number of iterations (generations)
RC = 16							#Number rows and columns in square matrix

disp = .2                                  	#proportion of population dispersing from each cell in each timestep
dispPgood = 0.1						#pred dispersal rate out of "good" habitat
dispPbad = seq(0,1,by=.1)				#pred dispersal rate out of "bad" habitat

Kgood = 100  #runif(1,1.8,2)				#prey K in good habitat
Kbad = 10  #runif(1,.7,.9)				#prey K in bad habitat
Mgood = 0.1							#mortality rate in good habitat
Mbad = 0.1							#mortality rate in bad habitat

g = 1.5                                   	#prey intrinsic growth rate
b = 0.7							#predator conversion efficiency
a = 0.6							#predator attack rate
h = 0.1							#predator handling time
mP = 0.1							#mortality

disperse_i <-function(move_i, i){
	move <- i + move_i
		if(move > RC) {move = move - RC}
		if(move < 1)  {move = move + RC}
	move
}

disperse_j <-function(move_j, j){
	move <- j + move_j
		if(move > RC) {move = move - RC}
		if(move < 1)  {move = move + RC}
	move
}

Nabund = matrix(0,rep,length(dispPbad))
Pabund = matrix(0,rep,length(dispPbad))
prH = matrix(0,rep,length(dispPbad))
HabA = matrix(0,rep,length(dispPbad))
HabB = matrix(0,rep,length(dispPbad))
VarHab = matrix(0,rep,length(dispPbad))
TotHab = matrix(0,rep,length(dispPbad))
VarLocHab = matrix(0,rep,length(dispPbad))
StatTotalHab = matrix(0,rep,length(dispPbad))

	pb = winProgressBar(title = "progress bar", min = 0, max = rep, width = 300)				#Setup progress bar
for (r in 1:rep){
	Sys.sleep(0.0001);setWinProgressBar(pb, r, title=paste(round(r/rep*100, 1),"% complete"))		#Start progress bar

##Generate landscape
	Landscape = matrix(1,RC,RC)				#Landscape of all 1's
	H = runif(1,0,1)						#Proportion of landscape that can be hospitable
	Habitat = (RC*RC)*H					#Number of matrix entries that can be hospitable
	NoHabitat = matrix(0,RC,RC)				#Matrix of inhospitable (0) and hospitable (1) locations
	while(sum(NoHabitat) < Habitat){
		x = round(runif(1,1,(RC*RC)))
		NoHabitat[x] = 1
	}

	Lands1 = Landscape*NoHabitat
	H = runif(1,0,1)						#Proportion of landscape that can be hospitable
	Habitat = (RC*RC)*H					#Number of matrix entries that can be hospitable
	NoHabitat = matrix(0,RC,RC)				#Matrix of inhospitable (0) and hospitable (1) locations
	while(sum(NoHabitat) < Habitat){
		x = round(runif(1,1,(RC*RC)))
		NoHabitat[x] = 1
	}
	Lands2 = Landscape*NoHabitat
	StaticMap = Lands1+Lands2
	StaticMap = StaticMap/StaticMap; StaticMap[is.nan(StaticMap)] = 0
	LandChange = 1						#Landscape changes if = 1

	for(d in 1:length(dispPbad)){
	N = array(1,dim=c(RC,RC,T))				#Array of prey population sizes in each landscape cell in each time step
	P = array(1,dim=c(RC,RC,T))				#Array of predator population sizes in each landscape cell in each time step
		for (t in 2:T){
			if(LandChange == 1){
				if(t %% 2 == 0) {Lands = Lands2}
				else{Lands = Lands1}
			}
			else{Lands = Lands1}

################Population dynamics and dispersal
##Prey dispersal
			temp_disp = matrix(0,RC,RC) 
      		for (i in 1:RC){
            		for (j in 1:RC){
					move_i = rpois(1, lambda = 3)+1                 #dispersal distance sampled from poisson distribution
               			move_j = rpois(1, lambda = 3)+1               	#dispersal distance sampled from poisson distribution
						if(move_i > 15){move_i=15}
						if(move_j > 15){move_j=15}
               			target_i = disperse_i(move_i,i)
                			target_j = disperse_j(move_j,j)
               			temp_disp[target_i, target_j] = ((disp*N[i,j,(t-1)]) + temp_disp[target_i,target_j])        #add number dispersing to proper cell of temp_disp matrix
               			temp_disp[i,j] = (temp_disp[i,j] - (disp*N[i,j,(t-1)]))                                     #subtract individuals from the cell in which they just moved
           			}
			}
			N[,,t] = temp_disp[,] + N[,,(t-1)]

##Prey pop dynamics
			for (i in 1:RC){
          			for (j in 1:RC){
                 			if (Lands[i,j] == 1) {K = Kgood; m = Mgood}					#K (carrying capacity) when in landscape habitat type 1
                  		else{K = Kbad; m = Mbad}								#K (carrying capacity) when in landscape habitat type 0      		
               			N[i,j,t] = (N[i,j,t]*exp(g*(1-(N[i,j,t]/K)))) - (m*N[i,j,t]) 		#ricker model
						if (N[i,j,t]<0.0000001){N[i,j,t]=0}						#trap - makes local pop size zero if below 0.0000001
				}
			}

##Predator dispersal
			temp_disp = matrix(0,RC,RC) 
     			for (i in 1:RC){
           			for (j in 1:RC){
					if((b*a*N[i,j,t])/(1 + (a*h*N[i,j,t]*P[i,j,t])) > 1){
						dispersal=dispPgood
					}
					else{
						dispersal = dispPbad[d]
					}
					move_i = rpois(1, lambda = 3)+1                 #dispersal distance sampled from poisson distribution
               			move_j = rpois(1, lambda = 3)+1               	#dispersal distance sampled from poisson distribution
						if(move_i > 15){move_i=15}
						if(move_j > 15){move_j=15}
                			target_i = disperse_i(move_i,i)
                			target_j = disperse_j(move_j,j)
					
               			temp_disp[target_i, target_j] = ((dispersal*P[i,j,(t-1)]) + temp_disp[target_i,target_j])        #add number dispersing to proper cell of temp_disp matrix
               			temp_disp[i,j] = (temp_disp[i,j] - (dispersal*P[i,j,(t-1)]))                                     #subtract individuals from the cell in which they just moved
           			}
			}
			P[,,t] = temp_disp[,] + P[,,(t-1)]

##Pred pop dynamics
			for (i in 1:RC){
          			for (j in 1:RC){
               			N[i,j,t] = (N[i,j,t]) - (a*N[i,j,t]*P[i,j,t])/(1+(a*h*N[i,j,t]*P[i,j,t]))
						if (N[i,j,t]<0.0000001){N[i,j,t]=0}						#trap - makes local pop size zero if below 0.0000001
					P[i,j,t] = (b*a*N[i,j,t]*P[i,j,t])/(1 + (a*h*N[i,j,t]*P[i,j,t])) - (mP*P[i,j,t])
						if (P[i,j,t]<0.0000001){P[i,j,t]=0}						#trap - makes local pop size zero if below 0.0000001

				}
			}
		}
		Nequalibabund = array(0,dim=c(RC,RC))
		for (i in 1:RC){
    			for (j in 1:RC){
            		Nequalibabund[i,j] = mean(N[i,j,(T-200):T])		#mean abundance in each cell over the final 4 timesteps
    			}
		}
		Nabundmap = Nequalibabund

		Nlandabund = seq(1:T)					#total abundance in the landscape at each timestep
		for (t in 1:T){
			Nlandabund[t] = sum(N[,,t])
		}

		Pequalibabund = array(0,dim=c(RC,RC))
		for (i in 1:RC){
    			for (j in 1:RC){
            		Pequalibabund[i,j] = mean(P[i,j,(T-200):T])		#mean abundance in each cell over the final 4 timesteps
    			}
		}
		Pabundmap = Pequalibabund

		Plandabund = seq(1:T)					#total abundance in the landscape at each timestep
		for (t in 1:T){
			Plandabund[t] = sum(P[,,t])
		}

		Nabund[r,d] = sum(Nequalibabund)
		Pabund[r,d] = sum(Pequalibabund)
		prH[r,d] = H
		HabA[r,d] = sum(Lands1)
		HabB[r,d] = sum(Lands2)
		TotHab[r,d] = sum(HabA[r],HabB[r])
		VarHab[r,d] = var(c(HabA[r],HabB[r]))
		StatTotalHab[r,d] = sum(StaticMap)
		VarLocMat = matrix(0,RC,RC)
		for (i in 1:RC){
			for (j in 1:RC){
				VarLocMat[i,j] = var(c(Lands1[i,j],Lands2[i,j]))
			}
		}
		VarLocHab[r,d] = mean(VarLocMat)
	}
}
close(pb)		#Turn off progress bar

difPabund = matrix(0,rep,10)
for (i in 1:rep){
	for (j in 1:10){
		difPabund[i,j] = Pabund[i,j+1] - Pabund[i,1]
	}
}
difNabund = matrix(0,rep,10)
for (i in 1:rep){
	for (j in 1:10){
		difNabund[i,j] = Nabund[i,j+1] - Nabund[i,1]
	}
}

gen = seq(1:T)
mx = max(c(Nlandabund,Plandabund))
mn = min(c(Nlandabund,Plandabund))
plot(Nlandabund  ~ gen, type="o", xlab = "Generation", ylab = "Total abundance in landscape", col="green")
points(Plandabund ~ gen, type="o",col="blue")
