#setwd('P:/paper/CoxExact/simu3/test');
#source('ExactCox.txt');

ExactCox = function(time, status, trt, hr = 1, alternative = 'two.sided', conf.int = FALSE, conf.level = 0.95) {
	levels = unique(trt);
	stopifnot(length(levels) == 2);
	trt = (trt == levels[2]);
	levels2 = unique(status);
	stopifnot(length(levels2) == 2 & min(levels2) >= 0 & max(levels2) == 1);

AlphaInverse = function(alpha0, x, y, m, n, d, margin) {
	prop1 = m * margin / (m * margin + n);
	q0 = qbinom(alpha0, d-(x+y), prop1) + x;
	if (alpha0 <= 0) return(0)
	if (alpha0 >= 1) {
		result = pbinom(q0, d, prop1);
		return(result);
	}
	tp11 = pbinom(q0-x-1, d-(x+y), prop1);
	tp12 = pbinom(q0-x, d-(x+y), prop1);
	tp21 = pbinom(q0-1, d, prop1);
	tp22 = pbinom(q0, d, prop1);
	if (tp12 > 0) {
		slope = (tp22 - tp21) / (tp12 - tp11);
		result = slope * (alpha0 - tp11) + tp21;
		return(result);
	}
	else return(0)
}

AlphaInverse1 = function(alpha0, x, y, m, n, d, margin) {
	prop1 = m * margin / (m * margin + n);
	q0 = qbinom(alpha0, d-(x+y), prop1) + x;
	if (alpha0 <= 0) return(0)
	if (alpha0 >= 1) {
		result = pbinom(q0, d, prop1);
		return(result);
	}
	tp11 = pbinom(q0-x-1, d-(x+y), prop1);
	tp12 = pbinom(q0-x, d-(x+y), prop1);
	di = x + y;
	Di1 = d - (x + y);
	tp21 = sum(BiasedUrn::dFNCHypergeo(0:di, m, n, di, margin) * pbinom(q0-(0:di)-1, Di1, prop1));
	tp22 = sum(BiasedUrn::dFNCHypergeo(0:di, m, n, di, margin) * pbinom(q0-(0:di), Di1, prop1));
	if (tp12 > 0) {
		slope = (tp22 - tp21) / (tp12 - tp11);
		result = slope * (alpha0 - tp11) + tp21;
		return(result);
	}
	else return(0)
}

CoxExactTest = function(X, margin = 1, left = TRUE) {
	if (NROW(X) < 1) {
		pval = list();
		pval$left = 1;
		return(pval);
	}
	if (!left) {
		names(X) = c('atRisk2', 'atRisk1', 'numEvent2', 'numEvent1');
		margin = 1 / margin;
	}
	numEvent = X$numEvent1 + X$numEvent2;
	nTime = dim(X)[1];
	nTotalEvent = sum(numEvent);
	nCurrentEvent = 0;
	result = list();
	pval = 1;
	atRisk1 = X$atRisk1;
	atRisk2 = X$atRisk2;
	e2 = X$numEvent2;
	e1 = X$numEvent1;
	for (i in nTime:1) {
		nCurrentEvent = nCurrentEvent + numEvent[i];
		if (numEvent[i] <= 1) pval = AlphaInverse(pval, e2[i], e1[i], atRisk2[i], atRisk1[i], nCurrentEvent, margin)
		else pval = AlphaInverse1(pval, e2[i], e1[i], atRisk2[i], atRisk1[i], nCurrentEvent, margin)
		if (pval > 1) pval = 1
	}
	result$left = pval;
	return(result);
}


CoxExactCi = function(X, alpha = 0.025) {
	pvalFun = function(hazRatio) {
		if (alpha > 0.5) {
			pval = CoxExactTest(X, hazRatio, left = TRUE);
			return(1 - alpha - pval$left)
		}
		else {
			pval = CoxExactTest(X, hazRatio, left = FALSE);
			return(pval$left - alpha)
		}
	}

	if (sum(X$numEvent2) < 0.5 & alpha <= 0.5) return(0)
	if (sum(X$numEvent1) < 0.5 & alpha >= 0.5) return(Inf)
	fit = uniroot(pvalFun, interval = c(1e-8, 1e+8), tol = 1e-8, extendInt = 'upX', trace=1);
	return(fit$root);
}

# Calculate the at risk numbers (m_i, n_i)
ConvertSurvData = function(survData) {
	result = survData[order(survData$time), ];
	N = dim(result)[1];

	atRisk1 = rep(NA, N);
	atRisk2 = atRisk1;
	numEvent1 = result$status & (result$trt == 0);
	numEvent2 = result$status & (result$trt == 1);
	atRisk1[1] = sum(result$trt == 0);
	atRisk2[1] = sum(result$trt == 1);
	trt = result$trt;
	time = result$time;
	if (N < 2) return(result)
	for (i in 2:N) {
		tmp1 = c(0, 1);
		atRisk1[i] = atRisk1[i - 1] - 1 + trt[i - 1];
		atRisk2[i] = atRisk2[i - 1] - trt[i - 1];
		if (time[i] == time[i - 1]) {
			numEvent1[i] = numEvent1[i - 1] + numEvent1[i];
			numEvent2[i] = numEvent2[i - 1] + numEvent2[i];
		}
	}
	tmp2 = rep(TRUE, N);
	tmp2[1:(N-1)] = (result$time[1:(N-1)] != result$time[2:N]);
	tmp3 = rep(TRUE, N);
	tmp3[2:N] = (result$time[1:(N-1)] != result$time[2:N]);
	result4 = data.frame(cbind(atRisk1[tmp3], atRisk2[tmp3], numEvent1[tmp2], numEvent2[tmp2]));
	names(result4) = c('atRisk1', 'atRisk2', 'numEvent1', 'numEvent2');
	result5 = result4[result4$numEvent1 + result4$numEvent2 > 1e-8 &
		result4$atRisk1 > 1e-8 & result4$atRisk2 > 1e-8, ];
	return(result5);
}

	X1 = data.frame(cbind(time, status, trt));
	LifeTable = ConvertSurvData(X1);
	rval = list();

	if (alternative == 'greater') {
		cet = CoxExactTest(LifeTable, margin = hr, left = FALSE);
		rval$p.value = cet$left;
	}
	else if (alternative == 'less') {
		cet = CoxExactTest(LifeTable, margin = hr, left = TRUE);
		rval$p.value = cet$left;
	}
	else if (alternative == 'two.sided') {
		cet1 = CoxExactTest(LifeTable, margin = hr, left = TRUE);
		cet2 = CoxExactTest(LifeTable, margin = hr, left = FALSE);
		cet = min(cet1$left, cet2$left) * 2;
		rval$p.value = min(cet, 1);
	}
	else stop('Arguments for alternative can only be: greater, less, two.sided.')

	if (conf.int) {
		lower_ci = CoxExactCi(LifeTable, (1 - conf.level) / 2);
		upper_ci = CoxExactCi(LifeTable, (1 + conf.level) / 2);
		rval$conf.int = c(lower_ci, upper_ci);
	}
	rval$alternative = alternative;
	call = match.call();
	class(rval) <- "ExactCox"
	return(rval)
}
