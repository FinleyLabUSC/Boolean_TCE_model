print(getwd())
setwd("Projects/TCE"); #settign the working directory where this code will be run 
print(getwd()) #check to make sure you are in the right directory 
# Set this before run:
initState = "acute"	# TCF1hi Tex naive

TCR = NFkB = NFATC2 = AP1 = IRF4 = IFNg = AKT = mTOR = glycolysis = NFATC1.med = NFATC1.lo = NFATC1 = NR4A1 = 
  FOXO1 = AP1.DNA = BLIMP1 = BCL6 = TCF1 = BATF = IL2 = IL12 = IL21 = PD1 = BATF.IRF4 = STAT3 = 0.0  #starting all genes OFF (0)
startState <- sort(c(TCR,NFkB,NFATC2,AP1,IRF4,IFNg,AKT,mTOR,glycolysis,NFATC1.med,NFATC1.lo,NFATC1,NR4A1,FOXO1,AP1.DNA,BLIMP1,
                     BCL6,TCF1,BATF,IL2,IL12,IL21,PD1,BATF.IRF4,STAT3)) #creates a start state variable for each gene 
names(startState) <- sort(c("TCR","NFkB","NFATC2","AP1","IRF4","IFNg","AKT","mTOR","glycolysis","NFATC1.med","NFATC1.lo",
                            "NFATC1","NR4A1","FOXO1","AP1.DNA","BLIMP1","BCL6","TCF1","BATF","IL2","IL12","IL21","PD1","BATF.IRF4","STAT3")) #labeling the starting state genes 


prevS2 = prevS1 = NULL
prevS2 = prevS1 = S = startState

if (initState == "acute") {
  # startState[c("TCR", "BCL6", "TCF1")] <- 1
	startState[c("TCR", "IL2", "IL12", "IL21")] <- 1
	# startState[c("TCR", "BCL6", "BLIMP1", "TCF1")] <- 1
}

# texS <- rep(0, length(startState))
# names(texS) <- names(startState)
# texS[c("PD1","NFATC1","NR4A1","BLIMP1", "TCR")] <- 1

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
	  write.csv(allStates,sprintf("NFATC1_mod/allStates10k_PD1block_sim%d.csv", N),sep="\t", col.names=TRUE, row.names=TRUE, append = FALSE) #saving all iterations per simulation to sep files
	  
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
	write.csv(allStates,sprintf("NFATC1_mod/allStates10k_PD1block_sim%d.csv", N),sep="\t", col.names=TRUE, row.names=TRUE, append = FALSE) # saving all iterations per simulation to sep files

}
