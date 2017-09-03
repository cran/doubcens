DGLwghts = function(X_L, X_R, Z_L, Z_R) {
	N = length(X_L) 
	x_val = c()
	t_val = c()
	for (i in 1:N) {
		x_val = c(x_val, X_L[i], X_R[i])
		t_val = c(t_val, Z_L[i] - X_L[i], Z_L[i] - X_R[i], 
					Z_R[i] - X_L[i], Z_R[i] - X_R[i])
	}
	x_val = sort(unique(x_val[x_val > 0]))
	t_val = sort(unique(t_val[t_val > -Inf]))
	alpha = matrix(0, nrow=length(x_val)*length(t_val), ncol=N)
	combos = expand.grid(1:length(x_val), 1:length(t_val))
	jth_index = combos[,1]; kth_index = combos[,2]
	for (i in 1:N) {
		for (j in 1:length(combos$Var1)) {
			if (X_L[i] <= x_val[combos[j,1]] & x_val[combos[j,1]] <= X_R[i] &
			Z_L[i] <= x_val[combos[j,1]] + t_val[combos[j,2]] & 
			x_val[combos[j,1]] + t_val[combos[j,2]] <= Z_R[i]) 
			{
				alpha[j,i]=1 
			}
		}
	}
	w_new = rep(1/length(x_val), length(x_val))
	f_new = rep(1/length(t_val), length(t_val))
	mu = matrix(0, nrow=length(alpha[,1]), ncol=length(alpha[1,]))
	for (i in 1:N) {
		w_group = w_new[jth_index]
		f_group = f_new[kth_index]
		for (j in 1:length(combos$Var1)) {
			mu[j,i] = (alpha[j,i]*w_new[combos[j,1]]*f_new[combos[j,2]]) / 
			(sum(alpha[,i]*w_group*f_group))
		}
	}
	eps = 10^{-5}
	diff_onset = w_new
	diff_failure = f_new
	tol = sum(abs(diff_onset)) + sum(abs(diff_failure))
	counter = 0 
	while (tol > eps) {
		w_old = w_new
		f_old = f_new
		for (j in 1:length(x_val)) {
			w_new[j] = sum(mu[combos$Var1==j,]) / N
		}
		for (k in 1:length(t_val)) {
			f_new[k] = sum(mu[combos$Var2==k,]) / N
		}
		diff_onset = w_old - w_new
		diff_failure = f_old - f_new
		tol = sum(abs(diff_onset)) + sum(abs(diff_failure))
		for (i in 1:N) {
			w_group = w_new[jth_index]
			f_group = f_new[kth_index]
			for (j in 1:length(combos$Var1)) {
				mu[j,i] = (alpha[j,i]*w_new[combos[j,1]]*f_new[combos[j,2]]) / 
				(sum(alpha[,i]*w_group*f_group))
			}
		}
		counter = counter + 1
	}
	return(list(x_val, w_new, t_val, f_new, counter))
}
