Sunwghts = function(Ei, Ri, Li, Ui, Bi1, Bi2) {
	stnd = 1 - Ei 
	Einew = Ei + stnd; Rinew = Ri + stnd
	Linew = Li + stnd; Uinew = Ui + stnd
	Bi1new = Bi1 + stnd; Bi2new = Bi2 + stnd
	diff = Linew - Rinew
	vals = c(0, diff, Uinew) 
	vals = sort(unique(vals))
	if (vals[length(vals)] == Inf) {
		uj = vals[1:(length(vals)-2)] 
	} else {
		uj = vals[1:(length(vals)-1)]}
	n = length(Ei); m = length(uj)
	alphaij = matrix(0, nrow=n, ncol=m); alphaijstar = matrix(1, nrow=n, ncol=m)
	betaij = matrix(0, nrow=n, ncol=m); betaijstar = matrix(1, nrow=n, ncol=m)

	for (i in 1:n) {
		for (j in 1:m) {
			xij = seq(1, Rinew[i], by=(Rinew[i]-1)/99)
			if (Linew[i] - Rinew[i] <= uj[j] & uj[j] <= Uinew[i]) { 
				alphaij[i,j] = 1
				if (Linew[i] == Uinew[i]) {
					alphaijstar[i,j] = 1/100 }
				else {
					alphaijstar[i, j] = sum( (xij + uj[j] >= Linew[i]-Rinew[i]) & (xij + uj[j] <= Uinew[i] - 1) )/100 }
			}
			if (Bi1new[i] - Rinew[i] <= uj[j] & uj[j] <= Bi2new[i]) { 
				betaij[i, j] = 1
				betaijstar[i, j] = sum(Bi1new[i] - uj[j] <= xij & xij <= Bi2new[i] - uj[j])/100
				}
		}
	}
	fnew = rep(1/length(uj), length(uj))
	tol = 10^(-5)
	counter = 0
	epsilon = 1
	muij = matrix(0, nrow=n, ncol=m); nuij = matrix(0, nrow=n, ncol=m)
	pij = rep(1, length(uj))
	while (epsilon > tol) {
		fold = fnew
		for (i in 1:n) {
			for (j in 1:m) {
				muij[i,j] = (alphaij[i,j]*alphaijstar[i,j]*fold[j])/sum(alphaij[i,]*alphaijstar[i,]*fold)
				nuij[i,j] = ((1-betaij[i,j]*betaijstar[i,j])*fold[j])/sum(betaij[i,]*betaijstar[i,]*fold)
			}
		}
		M = sum(muij+nuij)
		for (j in 1:m) {
			fnew[j] = (1/M)*(sum(muij[,j]+nuij[,j]))
		}
		epsilon = sum(abs(fnew - fold))
		counter = counter + 1
	}
	return(list(uj, fnew, counter))			
}
