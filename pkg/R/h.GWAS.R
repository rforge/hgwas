h.GWAS <-
function(y, X, Z, family = 'normal', link = 'logit', X.disp = NULL, 
	ac.mat = NULL, rho = 0, phi.start = 1, lambda.start = .1, alpha.start = 1,
	conv.crit = 1e-5, max.iter = 200, plotting = TRUE) {

## rho = .911 is the default value for smoothed DHGLM
## rho cannot be 1, otherwise chol() would not work on ac.mat
	
N <- nrow(X)
p <- ncol(X)
k <- ncol(Z)
if (is.null(X.disp)) {
	X.disp <- rep(1, k)
	pd <- 1
} else {
	pd <- ncol(X.disp)
}

# ----- autocorrelation matrix of the variance of the random effects	
if (is.null(ac.mat)) ac.mat <- rho**(toeplitz(1:k) - 1)
L <- t(chol(ac.mat))
cat('chol(ac.mat) ready.\n')

# ----- interlinked GLMs

if (family == 'normal') {

# ----- create augmented data frame 
Y <- c(y, rep(0, k))
zero <- matrix(0, k, p)
Aug.Data <- data.frame(rbind(cbind(X, Z), cbind(zero, diag(k))))
# ----- initialize the parameter values
alpha0 <- alpha.start
phi0 <- phi.start
lambda <- rep(lambda.start, k)
CONV <- FALSE
niter <- 1

if (plotting) pdf('IWLS_NNN.pdf', width = 21, height = 14)

# ----- start IWLS
while (!CONV & niter < max.iter) {
	# ----- estimate the mean: fixed and random effects
	lm1 <- lm(Y ~ . - 1, data = Aug.Data, weights = c(rep(1/phi0, N), 1/lambda))
	q <- as.numeric(hatvalues(lm1))
	dev <- (as.numeric(residuals(lm1)))^2
	dev_q <- dev/(1 - q)
	# ----- estimate dispersion parameter of y
	gl1 <- glm(as.numeric(dev_q[1:N]) ~ 1, family = Gamma(link = log), weights = as.numeric((1 - q[1:N])/2))
	phi <- as.numeric(exp(coef(gl1)))
	# ----- estimate the marker dispersion: fixed and random effects
	z <- c(as.numeric(dev_q[(N + 1):(N + k)]), rep(0, k))
	z[1:k] <- log(lambda) + (z[1:k] - lambda)/lambda ## log link is applied!
	zero2 <- matrix(0, k, pd)
	Aug.Data2 <- data.frame(rbind(cbind(X.disp, L), cbind(zero2, diag(k))))
	lm2 <- lm(z ~ . - 1, data = Aug.Data2, weights = c((1 - q[(N + 1):(N + k)])/2, 1/rep(alpha0, k))) ## mind the weights!
	q2 <- as.numeric(hatvalues(lm2))
	dev2 <- (as.numeric(residuals(lm2)))^2
	dev2_q <- dev2/(1 - q2)
	lambda <- exp(as.numeric(lm2$fitted.values)[1:k])
	# ----- estimate the dispersion of the marker dispersion random effects
	gl2 <- glm(as.numeric(dev2_q[(k + 1):(k + k)]) ~ 1, family = Gamma(link = log), weights = as.numeric((1 - q2[(k + 1):(k + k)])/2))
	alpha <- as.numeric(exp(coef(gl2)))
	# ----- criteria & output
	CONV <- max(abs(c(alpha, phi) - c(alpha0, phi0))) < conv.crit
	cat('-------- iteration', niter, '--------\n')
	cat('max(lambda) =', max(lambda), '\n')
	cat('min(lambda) =', min(lambda), '\n')
	cat('alpha =', alpha, '\n')
	cat('|alpha - alpha0| =', abs(alpha - alpha0), '\n')
	cat('|phi - phi0| =', abs(phi - phi0), '\n')
	cat('-----------------------------\n')
	cat('\n')
	# ----- update dispersion parameters
	alpha0 <- alpha
	phi0 <- phi
	niter <- niter + 1 
	# ----- plots for current iteration
	if (plotting) {
		par(mfrow = c(2, 1))
		plot(1:k, coef(lm1)[(p + 1):(p + k)], xlab = 'Marker index', ylab = 'Marker effects', type = 'h', lwd = .6, col = 4, bty = 'n')
		plot(1:k, lambda, xlab = 'Marker index', ylab = 'Marker dispersion', type = 'h', col = 2, bty = 'n')
	}
}

if (plotting) dev.off()

}

if (family == 'binary') {

# ----- create augmented data frame 
Y <- c(y, rep(0, k))
mu <- rep(mean(y), N)
zero <- matrix(0, k, p)
Aug.Data <- data.frame(rbind(cbind(X, Z), cbind(zero, diag(k))))
# ----- initialize the parameter values
alpha0 <- alpha.start
phi0 <- 1 ## fix the binomial dispersion parameter to be 1
lamda <- rep(lambda.start, k)
CONV <- FALSE
niter <- 1

if (plotting) pdf('IWLS_BNN.pdf', width = 21, height = 14)

# ----- start IWLS
while (!CONV & niter < max.iter) {
	V.mu <- mu*(1-mu)
	# ----- link functions, linear predictors and derivatives
	if(link=='logit') {eta <- log(mu/(1-mu));deta.dmu <- 1/V.mu}
	if(link=='probit') {eta <- log(mu/(1-mu));deta.dmu <- 1/dnorm(qnorm(mu))}
	if(link=='cloglog') {eta <- cloglog(mu);deta.dmu <- -1/(1-mu)/log(1-mu)}
	# ----- weight matrix 
	Gamma_M.inv <- c(rep(1/phi0, N), 1/lamda)
	W_M0 <- 1/deta.dmu^2/V.mu
	W_M1 <- rep(1, k)
	W_Ma <- c(W_M0,W_M1)
	Sigma_M.inv <- W_Ma*Gamma_M.inv
	# ----- working response
	Y[1:N] <- eta + (y - mu)*deta.dmu
	# ----- estimate the mean: fixed and random effects 
	lm1 <- lm(Y ~ . - 1, data = Aug.Data, weights = Sigma_M.inv)
	q <- hatvalues(lm1)
	dev <- (residuals(lm1))^2
	eta_hat <- predict(lm1)[1:N]
	mu_hat <- exp(eta_hat)/(1 + exp(eta_hat))
	dev_q <- dev/(1 - q)
	phi <- 1 ## fix the binomial dispersion parameter to be 1
	# ----- estimate the marker dispersion: fixed and random effects
	z <- c(as.numeric(dev_q[(N + 1):(N + k)]), rep(0, k))
	z[1:k] <- log(lamda) + (z[1:k] - lamda)/lamda
	zero2 <- matrix(0, k, pd)
	Aug.Data2 <- data.frame(rbind(cbind(X.disp, L), cbind(zero2, diag(k))))
	lm2 <- lm(z ~ . - 1, data = Aug.Data2, weights = c(as.numeric(1 - q[(N + 1):(N + k)])/2, 1/rep(alpha0, k)))
	q2 <- hatvalues(lm2)
	dev2 <- (residuals(lm2))^2
	dev2_q <- dev2/(1 - q2)
	lambda <- exp(as.numeric(lm2$fitted.values)[1:k])
	# ----- estimate the dispersion of the marker dispersion random effects
	gl2 <- glm(as.numeric(dev2_q[(k + 1):(k + k)]) ~ 1, family = Gamma(link = log), weights = as.numeric((1 - q2[(k + 1):(k + k)])/2))
	alpha <- as.numeric(exp(coef(gl2)))
	# ----- criteria & output
	CONV <- max(abs(c(alpha, phi) - c(alpha0, phi0))) < conv.crit
	cat('-------- iteration', niter, '--------\n')
	cat('max(lambda) =', max(lambda), '\n')
	cat('min(lambda) =', min(lambda), '\n')
	cat('alpha =', alpha, '\n')
	cat('|alpha - alpha0| =', abs(alpha - alpha0), '\n')
	cat('|phi - phi0| =', abs(phi - phi0), '\n')
	cat('-----------------------------\n')
	cat('\n')
	# ----- update dispersion parameters
	if (niter%%5 == 0) alpha0 <- alpha # ----- update alpha ONLY EVERY 5 ITERATIONS
	phi0 <- phi
	niter <- niter + 1 
	# ----- plots for current iteration
	if (plotting) {
		par(mfrow = c(2, 1))
		plot(1:k, coef(lm1)[(p + 1):(p + k)], xlab = 'Marker index', ylab = 'Marker effects', type = 'h', lwd = .6, col = 4, bty = 'n')
		plot(1:k, lambda, xlab = 'Marker index', ylab = 'Marker dispersion', type = 'h', col = 2, bty = 'n')
	}
}

if (plotting) dev.off()

}

result <- list(phi = phi, alpha = alpha, lambda = lambda, 
	beta = coef(lm1)[1:p], gamma = as.numeric(lm2$coefficients[1:pd]),
	v = coef(lm1)[(p + 1):(p + k)], b = as.numeric(lm2$coefficients[(pd + 1):(pd + k)]), niter = niter)
return(result)

# ----- end of function dhglm()	
}

