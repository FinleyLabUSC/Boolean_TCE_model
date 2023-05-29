print(getwd())
setwd("Projects/TCE"); 
print(getwd())
# Set this before run:
initState = "acute"	# TCF1hi Tex naive

TCR = NFkB = NFATC2 = AP1 = IRF4 = IFNg = AKT = mTOR = glycolysis = NFATC1.med = NFATC1.lo = NFATC1 = NR4A1 = 
  FOXO1 = AP1.DNA = BLIMP1 = BCL6 = TCF1 = BATF = IL2 = IL12 = IL21 = PD1 = BATF.IRF4 = STAT3 = 0.0
startState <- sort(c(TCR,NFkB,NFATC2,AP1,IRF4,IFNg,AKT,mTOR,glycolysis,NFATC1.med,NFATC1.lo,NFATC1,NR4A1,FOXO1,AP1.DNA,BLIMP1,
                     BCL6,TCF1,BATF,IL2,IL12,IL21,PD1,BATF.IRF4,STAT3))
names(startState) <- sort(c("TCR","NFkB","NFATC2","AP1","IRF4","IFNg","AKT","mTOR","glycolysis","NFATC1.med","NFATC1.lo",
                            "NFATC1","NR4A1","FOXO1","AP1.DNA","BLIMP1","BCL6","TCF1","BATF","IL2","IL12","IL21","PD1","BATF.IRF4","STAT3"))


prevS2 = prevS1 = NULL
prevS2 = prevS1 = S = startState

if (initState == "acute") {
	startState[c("TCR", "IL2", "IL12", "IL21")] <- 1
}

texS <- rep(0, length(startState))
names(texS) <- names(startState)
texS[c("PD1","NFATC1","NR4A1","BLIMP1")] <- 1

print(getwd())
# allSS = allS = goodIndx = NULL

allStates = S
# TIME = 0

for (N in 1:10000) {
	TIME = 0 
	S <- startState
	allStates <- startState
	
	while(!identical(prevS2, S) & TIME < 100) {
	  
	  allStates <- rbind(allStates, S) #collecting all states of genes at each TIME step per simulation 
	  
	  write.csv(allStates,sprintf("50sims_testing/allStates10k_sim%d.csv", N),sep="\t", col.names=TRUE, row.names=TRUE, append = FALSE) #saving all iterations per simulation to sep files

		TIME <- TIME + 1
		
		if (TIME == 2) prevS1 = S else if (TIME > 2) {
			prevS2 <- prevS1
			prevS1 <- S
		}
		
		indx <- sample(1:22, 22) #index for 22 boolean rules
		
		for (i in indx) {
			nextCMD <- rev(readLines("TexSimCmds.R", n=indx[i]))[1] #reading in rules by line 
			eval(parse(text=nextCMD)) 
		}
	}
	
	print(paste0(N, "   ", TIME, "   ", Sys.time()))	#printing end time per simulation 
	

	allStates <- rbind(allStates, S)
	write.csv(allStates,sprintf("50sims_testing/allStates10k_sim%d.csv", N),sep="\t", col.names=TRUE, row.names=TRUE, append = FALSE)


}


#OLD CODE 
# write.csv(allStates,sprintf("50sims_testing/allStates10k_sim%d.csv", N),sep="\t", col.names=TRUE, row.names=TRUE, append = FALSE)

# allStates <- rbind(allStates, S) 
# write.csv(allStates,sprintf("allStates10k_sim%d.csv", N),sep="\t", col.names=TRUE, row.names=TRUE, append = FALSE)
# goodSS <- 0
# if (identical(texS, S)) {
# 	goodSS <- 1
# 	goodIndx <- rbind(goodIndx, indx)
# }

# allSS <- c(allSS, goodSS)
# print('printing TIME')
# print(paste0(N, "   ", TIME, "   ", Sys.time()))	# (paste(S, collapse=" : "))


# initS <- rbind(initS, S)
# filename <- sprintf("10sims_testing/allStates10k_sim%d.csv", N)
# write.table(initS, file = filename,  sep="\t", append = FALSE)
# write.csv(initS,sprintf("10sims_testing/csvallStates10k_sim%d.csv", N),sep="\t", col.names=TRUE, row.names=TRUE, append = FALSE)
