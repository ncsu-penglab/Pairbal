###################################################################################
#	@author Hayden Brochu
#
#	This contains functions needed to identify microbial reference frames and to
# perform feature reduction via 'Pairbal' analysis.
###################################################################################

# Load packages
.cran_packages <- c("parallel", "pbmcapply", "tidyverse")
sapply(.cran_packages, require, character.only = T)

###################################################################################

# Accessory function to check threads.
#	nthreads: number of cores to use.
#	quietly: Whether or not status messages should be reported to user. Default: FALSE
checkThreads <- function(nthreads, quietly = F) {
	if (is.na(nthreads)) {
		nthreads <- min(4, detectCores())
		if (!quietly) message("Defaulting to ", nthreads, " threads.")
	} else if (nthreads > detectCores()) {
		nthreads <- min(4, detectCores())
		if (!quietly) message("Requested too many cores. Defaulting to ", nthreads, ".")
	} else if (nthreads < 1 | nthreads %% 1 != 0) {
		nthreads <- min(4, detectCores())
		if (!quietly) message("Requested cores is not a positive integer. Defaulting to ", nthreads, ".")
	}
	nthreads # Return the number of threads to be used.
}

# Accessory function to check counts.
#	cnts: count matrix with rows = microbial features, cols = samples. Data.frame. Matrix will be coerced.
#	quietly: Whether or not status messages should be reported to user. Default: FALSE
checkCounts <- function(cnts, nthreads, p, quietly = F) {
	if (class(cnts) != "data.frame") {
		if (class(cnts) == "matrix") {
			if (!quietly) message("Counts are a matrix. Coercing to data.frame.")
			cnts <- as.data.frame(cnts)
		} else {
			stop("Counts not a data.frame or matrix. Cannot be processed.")
		}
	}
	cnts # Return count matrix.
}

# Accessory function to check pseudocount.
#	cnts: count matrix with rows = microbial features, cols = samples. Data.frame. Matrix will be coerced.
#	p: pseudocount. Only added if zeroes detected in the count matrix
#	quietly: Whether or not status messages should be reported to user. Default: FALSE
checkPseudocount <- function(cnts, p, quietly = F) {
	# Check pseudocount provided
	if (all(cnts > 0)) {
		if (!quietly) message("All counts are non-zero. Proceeding with input count matrix as is.")
	} else {
		if (is.na(p)) {
			if (!quietly) message("Zero values detected and pseudocount not provided. Defaulting to pseudocount of 1.")
			cnts <- cnts + 1
		} else {
			if (!quietly) message("Zero values detected. Adding user-defined pseudocount of ", p)
			cnts <- cnts + p
		}
	}
	cnts # Return count matrix.
}

# Function to compute isometric log-ratio of x, y. These may be vectors (indicating they only contain 1 ASV) or matrices (contain >1 ASV).
#	x: Denominator of ILR
#	y: Numerator of ILR
computeILR <- function(x, y) {
	if(is.vector(x)) x <- t(as.matrix(x))
	n.x <- nrow(x) # Number of ASVs in denominator of ILR
	x.logmean <- colSums(log(x)) / n.x
	if(is.vector(y)) y <- t(as.matrix(y))
	if (!all(colnames(x) == colnames(y))) stop("Sample names not consistent between numerator and denominator features.")
	n.y <- nrow(y) # Number of ASVs in numerator of ILR
	y.logmean <- colSums(log(y)) / n.y
	y.logmean - x.logmean # return ILR
}

# Recursive function used for steps in MC simulation
# lrvars: log-ratio covariance matrix
# cands: candidate ASVs to add to reference frame
# All other variables associated with current reference frame to which we are proposing a change
stepMC <- function(lrvars, cands, it, V, s, x, y, A, B, C, Dp, Ep, Dn, En, s.init, n) {
	# Compute positive transitions
	if (length(Dp) > 0) { # i.e. there are candidates we can add
		ap <- 1 / ((s + 1) * (n - s - 1)) - 1 / (s * (n - s))
		bp <- 1 / (2 * s^2) - 1 / (2 * (s + 1)^2)
		cp <- 1 / (2 * (n - s)^2) - 1 / (2 * (n - s - 1)^2)
		dp <- -n / ((s + 1)^2 * (n - s - 1))
		ep <- n / ((s + 1) * (n - s - 1)^2)
		tp <- ap * A + bp * B + cp * C + dp * Dp + ep * Ep
	} else { # edge case. Invariant set contains all candidate ASVs, so we cannot add to it.
		tp <- NULL
	}

	# Compute negative transitions
	if (length(Dn) > 1) { # i.e. there are candidates we can remove
		an <- 1 / ((s - 1) * (n - s + 1)) - 1 / (s * (n - s))
		bn <- 1 / (2 * s^2) - 1 / (2 * (s - 1)^2)
		cn <- 1 / (2 * (n - s)^2) - 1 / (2 * (n - s + 1)^2)
		dn <- n / ((s - 1)^2 * (n - s + 1))
		en <- -n / ((s - 1) * (n - s + 1)^2)
		tn <- an * A + bn * B + cn * C + dn * Dn + en * En
	} else { # edge case. Invariant set contains only one candidate ASV, so we cannot remove it.
		tn <- NULL
	}

	## !! Note that either tn or tp will be non-NULL, since edge cases cannot occur at the same time

	# Adjust transition vectors by setting positive terms to zero, then generate step probability distribution
	# If all transitions are positive, then we've reached a local minimum and will stop.
	t <- c(tp, tn)
	t <- t[t < 0]
	if (length(t) > 0) {
		pstep <- t / sum(t)
		pstep <- cumsum(pstep[pstep > 0])
	} else {
		ret <- list(ASVs = y, var = as.numeric(V), transitions = it, s.init = s.init)
		return(ret)
	}

	# Randomly choose our step and retrieve the ASV
	z <- pstep[pstep > runif(1)] %>% names %>% .[1]

	# Compute updated variance and variance terms.
	if (z %in% y) { # removal of ASV from invariant set
		x <- c(x, z)
		y <- y[y != z]
		s <- length(y)
		A <- A + Dn[z] - En[z]
		B <- B - 2 * Dn[z]
		C <- C + 2 * En[z]
		cands.num <- c(names(Dp), z)
		Dp <- c(Dp, sum(lrvars[z, y])) - unlist(lrvars[z, cands.num]); names(Dp)[length(Dp)] <- z
		Ep <- c(Ep, sum(lrvars[z, x])) + unlist(lrvars[z, cands.num]); names(Ep)[length(Ep)] <- z
		Dn <- Dn[y] - unlist(lrvars[z, y])
		En <- En[y] + unlist(lrvars[z, y])
	} else { # addition of ASV to invariant set
		x <- x[x != z]
		y <- c(y, z)
		s <- length(y)
		A <- A - Dp[z] + Ep[z]
		B <- B + 2 * Dp[z]
		C <- C - 2 * Ep[z]
		cands.num <- names(Dp)[names(Dp) != z]
		Dp <- Dp[cands.num] + unlist(lrvars[z, cands.num])
		Ep <- Ep[cands.num] - unlist(lrvars[z, cands.num])
		Dn <- c(Dn, sum(lrvars[z, y])) + unlist(lrvars[z, y]); names(Dn)[length(Dn)] <- z
		En <- c(En, sum(lrvars[z, x])) - unlist(lrvars[z, y]); names(En)[length(En)] <- z
	}
	V <- A / (s * (n - s)) - B / (2 * s^2) - C / (2 * (n - s)^2)

	# Recursively try to step to new variance state
	it <- it + 1
	stepMC(lrvars, cands, it, V, s, x, y, A, B, C, Dp, Ep, Dn, En, s.init, n)
}

# Function to run the MC simulation
runMC <- function(lrvars, cands, smin.it, smax.it, seed.it) {
	# Initial state. We will keep track of how many transitions we make
	it <- 0

	# Constants needed for computations
	asvs.all <- rownames(lrvars)
	n <- length(asvs.all)

	# Get random set of ASVs
	set.seed(seed.it)
	s <- sample(smin.it:smax.it, size = 1)
	s.init <- s
	y <- sample(cands, size = s)

	# Get numerator ASVs
	x <- asvs.all[!asvs.all %in% y]

	# Compute initial variance
	A <- sum(lrvars[x, y])
	B <- sum(lrvars[y, y])
	C <- sum(lrvars[x, x])
	V <- A / (s * (n - s)) - B / (2 * s^2) - C / (2 * (n - s)^2)

	# Compute variance terms needed for variance transitions
	cands.num <- cands[cands %in% x]
	Dp <- rowSums(lrvars[cands.num, y])
	Ep <- rowSums(lrvars[cands.num, x])
	Dn <- rowSums(lrvars[y, y])
	En <- rowSums(lrvars[y, x])

	# Recursively try to step to new variance state
	stepMC(lrvars, cands, it, V, s, x, y, A, B, C, Dp, Ep, Dn, En, s.init, n)
}

# Function to compute the log-ratio covariance matrix
computeLRcovar <- function(cnts, p = NA, nthreads = NA, quietly = F) {
	# Check number of cores / threads requested
	nthreads <- checkThreads(nthreads, quietly = quietly)

	# Check count matrix
	cnts <- checkCounts(cnts, quietly = quietly)

	# Check pseudocount if needed
	if (min(cnts) == 0) {
		cnts <- checkPseudocount(cnts, p, quietly = quietly)
	}

	# Compute log-ratio covariance matrix
	lcnts <- cnts %>% log
	lcnts.lst <- split(t(lcnts), rep(1:nrow(lcnts), each = ncol(lcnts)))
	lrvars <- pbmclapply(lcnts.lst, function(x) apply(t(x - t(lcnts)), 1, function(y) var(y)), mc.cores = nthreads)
	lrvars <- do.call(cbind, lrvars)
	lrvars %>% as.data.frame %>% setNames(rownames(.))
}

# Wrapper function for MCMC procedure
# Required input: count matrix
InvariantSetMC <- function(cnts, min.pct = 0.9, min.cnt = 5, smin = 0.1, smax = 0.9, mc.it = NA, nthreads = NA, p = NA, lrvars = NA, seed = 100, quietly = F) {
	#################
	# 0) Check input

	# Check number of cores / threads requested
	nthreads <- checkThreads(nthreads, quietly = quietly)

	# Check count matrix
	cnts <- checkCounts(cnts, quietly = quietly)

	# Check pseudocount provided
	cnts <- checkPseudocount(cnts, p, quietly = quietly)

	# Check MC iterations
	if (is.na(mc.it)) {
		mc.it <- 250
	} else if (mc.it < 1 | mc.it %% 1 != 0) {
		stop("Number of MC iterations must be a positive integer. User-supplied value: ", mc.it)
	}

	# Check filtering cutoffs
	if (min.cnt < 0) stop(min.cnt, " is not a valid count cutoff. Must be >= 0")
	if (min.pct < 0 | min.pct > 1) stop(min.pct, " is not a valid percentage cutoff. Must be [0, 1]")

	# Get candidate ASVs based on filtering thresholds
	cands <- rownames(cnts[rowSums(cnts > min.cnt) >= min.pct * ncol(cnts),])

	# Check that size range is compatible with number of candidates
	k <- length(cands)
	if (!quietly) message("Identified ", k, " candidate ASVs for the invariant set")
	if (k < 2) stop("Insufficient ASVs after filtering. Consider relaxing min.cnt and/or min.pct.")
	if (smin > 1 | smin < 0) {
		stop("Minimum size of ", smin, " is out of range. Needs to be [0, 1]")
	} else if (smax > 1 | smax < 0) {
		stop("Maximum size of ", smax, " is out of range. Needs to be [0, 1]")
	} else if (smax < smin) {
		stop("Maximum size of ", smax, " is less than minimum size of ", smin, ". smax >= smax")
	}

	# Get minimum and maximum initial invariant set sizes using smin/smax quantile thresholds
	smin.it <- quantile(1:k, smin) %>% round %>% as.numeric
	smax.it <- quantile(1:k, smax) %>% round %>% as.numeric

	#################
	# 1) Compute variances for all pairwise balances

	if (!all(is.na(lrvars))) {
		if (!quietly) message("Log-ratio covariance matrix supplied. Checking that it is compatible with count matrix..")
		if (nrow(lrvars) == nrow(cnts) & ncol(lrvars) == nrow(cnts) & all(rownames(lrvars) == rownames(cnts)) & all(colnames(lrvars) == rownames(cnts))) {
			if (!quietly) message("Covariance matrix is compatible")
			if (class(lrvars) == "matrix") {
				if (!quietly) message("Covariance matrix is a matrix. Coercing to data.frame.")
				lrvars <- lrvars %>% as.data.frame
			}
		} else {
			if (!quietly) message("Covariance matrix is not compatible. Will re-generate now on the fly.")
			lrvars <- NA
		}
	}
	if (all(is.na(lrvars))) {
		if (!quietly) message("Generating log-ratio covariance matrix..")
		lrvars <- computeLRcovar(cnts, nthreads = nthreads, quietly = T)
	}

	#################
	# 2) Run MCMC procedure

	if (!quietly) message("Running MCMC procedure..")
	minima <- seed:(mc.it + seed - 1) %>% as.list %>% pbmclapply(function(seed.it) runMC(lrvars, cands, smin.it, smax.it, seed.it), mc.cores = nthreads)
	minima
}

# Compute invariant normalized balances
#	cnts: count matrix with rows = microbial features, cols = samples. Data.frame. Matrix will be coerced.
#	invariantASVs: invariant ASVs identified by computeInvariantSet function
#	topASVs: Numerator ASVs (cannot be in the invariant set). These can be specified (i.e. results from Pairbal) else all non-invariant set ASVs are assumed
#	p: pseudocount. Only added if zeroes detected in the count matrix. Default: 1
#	nthreads: number of processes to use
#	clr: whether the centered log ratio should be computed. This will override invariantASVs and topASVs parameters. Default: FALSE
#	quietly: Whether or not status messages should be reported to user. Default: FALSE
computeInvariantBal <- function(cnts, invariantASVs = NA, topASVs = NA, p = 1, nthreads = NA, clr = F, quietly = F) {
	### 0) Check input

	# Check number of cores / threads requested
	nthreads <- checkThreads(nthreads, quietly = quietly)

	# Check count matrix
	cnts <- checkCounts(cnts, quietly = quietly)

	# Check pseudocount provided
	cnts <- checkPseudocount(cnts, p, quietly = quietly)

	if (clr) { # CLR set to true. Compute centered log ratio.
		invariantASVs <- rownames(cnts)
		topASVs <- rownames(cnts)
	} else {
		# Check top and invariant ASVs
		if (any(is.na(invariantASVs))) stop("Invariant ASVs not specified or at least one is a missing value. Set clr = T if you intend to compute centered log ratios.")
		if (!all(invariantASVs %in% rownames(cnts))) stop("Some invariant ASVs not in count matrix rownames.")
		if (any(is.na(topASVs))) {
			message("Numerator ASVs not specified. Computing balance values for all non-invariant ASVs")
			topASVs <- rownames(cnts)[!rownames(cnts) %in% invariantASVs]
		} else { # numerator ASVs specified
			if (!all(topASVs %in% rownames(cnts))) stop("Some top ASVs not in count matrix rownames.")
		}
		if (length(topASVs) == 0) stop("Number of numerator ASVs is zero. Check invariant set?")
	}

	### 1) Compute balances

	ilrs <- topASVs %>% as.list %>% pbmclapply(function(x) computeILR(cnts[invariantASVs,], cnts[x,]), mc.cores = nthreads)
	do.call(cbind, ilrs) %>% as.data.frame %>% setNames(topASVs)
}

# Function to compute log-ratio for all pairwise combinations of microbial features (ASVs or OTUs).
#	cnts: count matrix with rows = microbial features, cols = samples. Data.frame. Matrix will be coerced.
#	nthreads: number of cores to use
#	p: pseudocount. Only added if zeroes detected in the count matrix
computePWB <- function(cnts, nthreads = NA, p = NA) {
	### 0) Check input
	# Check number of cores / threads requested
	nthreads <- checkThreads(nthreads)

	# Check count matrix
	cnts <- checkCounts(cnts)

	# Check pseudocount provided
	cnts <- checkPseudocount(cnts, p)

	### 1) Compute all pairwise balances
	lcnts <- cnts %>% log
	lcnts.lst <- split(t(lcnts), rep(1:nrow(lcnts), each = ncol(lcnts)))
	names(lcnts.lst) <- rownames(lcnts)
	lcnts.lst <- lcnts.lst[1:(length(lcnts.lst) - 1)] # Drop last entry since it will only produce duplicate pairwise balances
	lrs <- pbmclapply(as.list(names(lcnts.lst)), function(x) {
		xx <- lcnts.lst[[x]] # numerator values
		y <- lcnts[as.logical(cumsum(rownames(lcnts) == x)),] %>% .[-1,] # denominator values, subset keeping only ASVs that appear after ASV in our numerator
		xx <- t(xx - t(y)) %>% as.data.frame # compute log-ratios
		xx$numerator <- x
		xx$denominator <- rownames(xx)
		rownames(xx) <- NULL
		xx
	}, mc.cores = nthreads)
	bind_rows(lrs) # Return as a single data.frame
}

# Accessory function to test the contrast on a vector of values
#	x: balance values
#	contrast.vec: logical vector of same length as x indicating the bivariate groups (TRUE or FALSE)
#	test.type: type of statistical test to use. Defaults to two-sided wilcoxon rank sum test ("wilcox"), can also choose two-sided t-test ("t")
checkContrast <- function(x, contrast.vec, test.type = "wilcox") {
	if (test.type == "wilcox") {
		test.out <- wilcox.test(x[contrast.vec], x[!contrast.vec], "two.sided", exact = F)
	} else if (test.type == "t"){ # wilcox.test
		test.out <- t.test(x[contrast.vec], x[!contrast.vec], "two.sided")
	} else {
		message("Invalid test specified. Specify either \"wilcox\" or \"t\"")
		return(NA)
	}
	test.out$p.value
}

# Filter pairwise balances using a log fold change cutoff and p-value threshold. Three steps:
#	1) Compute contrasts for each balance.
#	2) Filter balances using log fold change cutoff
#	3) Filter balances using p-value cutoff and sort by p-value.
#
#	balances: Data.frame with all log-ratios. This is the output from computePWB.
#	meta: metadata with samples as rownames. Data.frame, if matrix it will be coerced to data.frame.
#	bivar: bivariate condition for the contrast. Must be a column name in the metadata provided.
#	group.oi: Group of interest that will have log fold change computed relative to other group. Must be in bivar. Defaults to randomly chosen group.
#	min.abs.lfc: minimum absolute log fold change, defaults to 1. Must be >= 0
#	p.cutoff: p-value cutoff, defaults to 0.05. Must be (0,1]
#	test.type: type of statistical test to use. Defaults to two-sided wilcoxon rank sum test ("wilcox"), can also choose two-sided t-test ("t")
filterPWB <- function(balances, meta, bivar, group.oi = NA, min.abs.lfc = NA, p.cutoff = NA, test.type = NA) {
	### 0) Check input.

	# Check test type
	if (is.na(test.type)) {
		test.type <- "wilcox"
		message("Test type not provided. Using default: ", test.type)
	} else if (!test.type %in% c("wilcox", "t")) {
		stop("Invalid test specified. Specify either \"wilcox\" or \"t\"")
	}

	# Check balances
	if (class(balances) != "data.frame") {
		if (class(balances) == "matrix") {
			message("Balances are supplied as a matrix. Coercing to data.frame.")
			balances <- data.frame(balances, check.names = F, stringsAsFactors = F)
		} else {
			stop("Balances not a data.frame or matrix. Make sure this is the output from computePWB.")
		}
	}
	if (!all(colnames(balances)[(ncol(balances) - 1):ncol(balances)] == c("numerator", "denominator"))) {
		stop("Balance input does not match expected format from computePWB output.\n",
			"\tColumn names should be sample names followed by numerator and denominator.")
	}
	samples <- colnames(balances)[1:(ncol(balances) - 2)]

	# Check metadata and bivariate condition provided
	if (class(meta) != "data.frame") {
		if (class(meta) == "matrix") {
			message("Input metadata is a matrix. Coercing to data.frame.")
			meta <- data.frame(meta, check.names = F, stringsAsFactors = F)
		} else {
			stop("Metadata is not a data.frame or matrix. Cannot be processed.")
		}
	}
	if (length(samples) != nrow(meta)) stop("Sample names in metadata (rownames) are not consistent with the those in the balance list.")
	if (!all(sort(samples) == sort(rownames(meta)))) stop("Sample names in metadata (rownames) are not consistent with the those in the balance list.")
	if (!bivar %in% colnames(meta)) stop("Bivariate condition not found in metadata.")
	factors <- unique(meta[,bivar])
	num.uniq <- length(factors)
	if (num.uniq != 2) stop("Bivariate condition does not have two factors. Found ", num.uniq, " factors: ", paste0(factors, collapse = ", "))

	# Check p-value cutoff provided
	if (is.na(p.cutoff)) {
		p.cutoff <- 0.05
		message("p-value cutoff not provided. Using cutoff of ", p.cutoff)
	}
	if (p.cutoff <= 0 | p.cutoff > 1) stop("p-value cutoff (p.cutoff) must be (0,1].")

	# Check absolute log fold change provided
	if (is.na(min.abs.lfc)) {
		min.abs.lfc <- 1
		message("Absolute log fold change cutoff not provided. Using cutoff of ", min.abs.lfc)
	} else if (min.abs.lfc < 0) {
		stop("Absolute log fold change cutoff (min.abs.lfc) must be greater than zero.")
	}

	### 1) Compute contrasts for each balance.

	# Make contrast.
	contrast.vec <- meta[samples, bivar]
	if (!is.factor(contrast.vec)) contrast.vec <- factor(contrast.vec)
	grps <- levels(contrast.vec)
	if (is.na(group.oi)) {
		group.oi <- grps[1]
	} else {
		if (!group.oi %in% grps) stop("Group of interest (", group.oi, ") not found in bivariate groups: ", paste0(grps, collapse = ", "))
	}
	contrast.vec <- contrast.vec == group.oi
	message("Contrast successfully parsed: ", sum(contrast.vec), " ", group.oi, " vs. ", sum(!contrast.vec), " ", grps[grps != group.oi])

	### 2) Filter balances using log fold change cutoff

	# Compute fold-change for each balance, then filter using lfc cutoff
	balances$lfc <- rowMeans(balances[, samples[contrast.vec]]) - rowMeans(balances[, samples[!contrast.vec]])
	rowsToKeep <- abs(balances$lfc) >= min.abs.lfc
	message(sum(rowsToKeep), " of ", nrow(balances), " balances remain after applying minimum fold-change cutoff of ", min.abs.lfc)
	balances <- balances[rowsToKeep,]

	### 3) Filter balances using p-value cutoff and sort by p-value.

	# Compute p-values using user-defined test (t-test or wilcoxon rank-sum)
	balances$p <- apply(balances[,samples], 1, checkContrast, contrast.vec, test.type)
	rowsToKeep <- balances$p < p.cutoff
	message(sum(rowsToKeep), " of ", nrow(balances), " balances remain after applying p-value cutoff of ", p.cutoff)

	# Sort balances by p-value.
	balances[order(balances$p),]
}

# Function to retrieve dimensionally reduced set of ASVs.
#	balances: List of filtered pairwise balances. Expected format is output from filterPWB function.
#	ret.all: logical. If true, will return all ASVs with both initial and greedy adjusted balance frequencies. Default: FALSE
getTopFeatures <- function(balances, ret.all = F) {
	# Check balances
	if (class(balances) != "data.frame") {
		if (class(balances) == "matrix") {
			message("Balances are supplied as a matrix. Coercing to data.frame.")
			balances <- data.frame(balances, check.names = F, stringsAsFactors = F)
		} else {
			stop("Balances not a data.frame or matrix. Make sure this is the output from computePWB.")
		}
	}
	if (!all(colnames(balances)[(ncol(balances) - 3):ncol(balances)] == c("numerator", "denominator", "lfc", "p"))) {
		stop("Balance input does not match expected format from filterPWB output.\n",
			"\tColumn names should be sample names followed by numerator, denominator, lfc, then p.")
	}

	# Check ret.all parameter
	if (!is.logical(ret.all)) stop(ret.all, " is not a logical value. Must be set to TRUE or FALSE. Default: FALSE")

	# Get ASV frequencies (coming from both numerator and denominator)
	asv.dat <- data.frame(table(c(balances$numerator, balances$denominator)), stringsAsFactors = F)
	colnames(asv.dat) <- c("ASV", "freq")
	asv.dat <- asv.dat[order(asv.dat$freq, decreasing = T),]

	###########
	# 2) Re-assign balances based on ASV frequencies (greedy approach, which removes ASVs that only appear in balances
	#	with more frequent ASVs).

	# Using greedy approach, go through the balances and "dock" balance counts from the less frequent ASVs.
	tmp <- balances
	asv.dat.g <- asv.dat
	asv.dat.g$freq.g <- asv.dat.g$freq
	i <- 1
	while(nrow(tmp) > 0) {
		asv.dat.g <- asv.dat.g[order(asv.dat.g$freq.g, decreasing = T),]
		asv.u <- asv.dat.g[i, "ASV"]
		rowsToRemove <- tmp$numerator == asv.u | tmp$denominator == asv.u
		asvs.d <- unlist(tmp[rowsToRemove, c("numerator", "denominator")])
		asvs.d <- asvs.d[asvs.d != asv.u]
		tmp <- tmp[!rowsToRemove,]
		asvsToAdjust <- asv.dat.g$ASV %in% asvs.d
		asv.dat.g[asvsToAdjust, "freq.g"] <- asv.dat.g[asvsToAdjust, "freq.g"] - 1
		i <- i + 1
	}

	# Remove ASVs with greedy adjusted frequencies of zero. If ret.all = TRUE, skip this.
	if (!ret.all) asv.dat.g <- asv.dat.g[asv.dat.g$freq.g > 0,]

	# Order ASVs by their final frequencies. Order may have changed slightly.
	asv.dat.g <- asv.dat.g[order(asv.dat.g$freq.g, decreasing = T),]

	# Infer direction of ASVs based on their consensus effect on treatment group (use meandiff)
	asv.dat.g$direction <- apply(asv.dat.g, 1, function(x) {
		up <- sum((balances$numerator == x["ASV"] & balances$lfc > 0) | (balances$denominator == x["ASV"] & balances$lfc < 0))
		ifelse(up > as.numeric(x["freq"]) / 2, "positive", "negative")
	})

	# Convert ASV column to character
	asv.dat.g$ASV <- asv.dat.g$ASV %>% as.character

	# Return final feature matrix.
	asv.dat.g
}
