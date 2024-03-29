require(dirmult)


Loglik <- function(Y, X, b, model) {
	# Compute the log likelihood, constant part discarded
	# The likelihood is scaled. Be careful when computing AIC BIC
	# Thu Apr 14 12:04:40 2022 ------------------------------
	## add the coerce into matrix
  X <- as.matrix(X)
	g <- exp(X %*% t(b))	# n * q
	gs <- rowSums(g)
	ys <- rowSums(Y)
	if (model == "dirmult") {
		res <- 	sum(lgamma(gs) - lgamma(ys + gs) +
						rowSums(lgamma(Y + g) - lgamma(g)))
	}
  if (model == "mult") {
		res <- 	sum(rowSums(Y * log(g)) - ys * log(gs))
	}
  if (model == "dir") {
		res <- 	sum(lgamma(gs) + rowSums((g-1) * log(Y) - lgamma(g)))
	}

  res / nrow(X)
}

Score <- function(Y, X, b, model) {
	# Compute the Score function at b
  X <- as.matrix(X)
	S <- 0
	g <- exp(X %*% t(b))	# n * q
	gs <- rowSums(g)
	ys <- rowSums(Y)
	
	# Thu Apr 14 09:54:24 2022 ------------------------------
	# so far just trust the equation
	# ?digamma the first derivative of a Gamma function
	## digamma(x) = d/dx{ ln \Gamma(x)} = \Gamma'(x) / \Gamma(x)
	
	if (model == "dirmult") {
		S <-  t((digamma(gs) - digamma(ys + gs) +
					digamma(Y + g) - digamma(g)) * g) %*% X
	}

	if (model == "mult") {
		S <- t((Y / g - ys / gs) * g) %*% X
	}

	if (model == "dir") {
		S <-  t((digamma(gs) - digamma(g) + log(Y)) * g) %*% X
	}
	
	S / nrow(X)
}

Hessian <- function(Y, X, b, model){
	#	Compute the diagonal of the hessian matrix at b

  X <- as.matrix(X)
	H <- 0
	g <- exp(X %*% t(b))	# n * q
	gs <- rowSums(g)
	ys <- rowSums(Y)
	if (model == "dirmult") {
			H <- t((trigamma(gs) - trigamma(ys + gs) +
						trigamma(Y + g) - trigamma(g)) * g^2 +
					(digamma(gs) - digamma(ys+gs) +
						digamma(Y + g) - digamma(g)) * g) %*% X^2
	}
  if (model == "mult") {
			H <- t((-Y / g^2 + ys / gs^2) * g^2  +
					(Y / g - ys / gs) * g) %*% X^2
	}
	if (model == "dir") {
			H <- t((trigamma(gs) - trigamma(g)) * g^2  +
					(digamma(gs)  - digamma(g) + log(Y)) * g) %*% X^2
	}
	#	Divided by the sample size
	H / nrow(X)
}


