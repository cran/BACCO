"betahat.fun" <-
function (xold, Ainv, d, give.variance = FALSE, func = regressor.basis) 
{
    H <- regressor.multi(xold, func = func)
    H <- as.matrix(H)
    colnames(H)[1] <- "constant"
    out <- solve(quad.form(Ainv, H), crossprod(crossprod(Ainv, 
        H), d))
    out <- as.vector(out)
    names(out) <- colnames(H)
    if (give.variance) {
        jj.sigma <- drop(sigmahatsquared(H, Ainv, d))
        return(list(betahat = out, sigmahatsquared = jj.sigma, 
            variance = jj.sigma * solve(quad.form(Ainv, H))))
    }
    else {
        return(out)
    }
}
"betahat.fun.A" <-
function (xold, A, d, give.variance = FALSE, func = regressor.basis) 
{
    H <- regressor.multi(xold, func = func)
    H <- as.matrix(H)
    colnames(H)[1] <- "constant"
    out <- solve(quad.form.inv(A, H), crossprod(H, solve(A, d)))
    out <- as.vector(out)
    names(out) <- rownames(H)
    if (give.variance) {
        jj.sigma <- drop(sigmahatsquared.A(H, A, d))
        return(list(betahat = out, sigmahatsquared = jj.sigma, 
            variance = jj.sigma * solve(quad.form.inv(A, H))))
    }
    else {
        return(out)
    }
}
"corr" <-
function (x1, x2, scales = NULL, pos.def.matrix = NULL, power = 2, 
    coords = "cartesian", spherical.distance.function = NULL) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales, nrow = length(scales))
    }
    x1 <- as.vector(unlist(x1))
    x2 <- as.vector(unlist(x2))
    if (power == 2) {
        corrscale <- function(x) {
            exp(-abs(x))
        }
    }
    else {
        corrscale <- function(x) {
            exp(-(abs(x))^(power/2))
        }
    }
    m <- switch(coords, cartesian = x1 - x2, spherical = spherical.distance.function(x1, 
        x2), stop("coords must be either Cartesian or Spherical"))
    return(corrscale(quad.form(pos.def.matrix, m)))
}
"corr.matrix" <-
function (xold, yold = NULL, use.neatversion = TRUE, distance.function = corr, 
    ...) 
{
    if (is.null(yold)) {
        yold <- xold
    }
    if (use.neatversion) {
        out <- apply(xold, 1, function(y) {
            apply(yold, 1, function(x) {
                distance.function(x, y, ...)
            })
        })
        return(out)
    }
    else {
        n <- nrow(xold)
        m <- nrow(yold)
        A <- matrix(NA, m, n)
        for (i in 1:n) {
            for (j in 1:m) {
                A[j, i] <- distance.function(xold[i, ], yold[j, 
                  ], ...)
            }
        }
        colnames(A) <- rownames(xold)
        rownames(A) <- rownames(yold)
        return(A)
    }
}
"estimator" <-
function (val, A, d, scales = NULL, pos.def.matrix = NULL, power = 2) 
{
    d.missing.estimated <- d + NA
    for (i in 1:nrow(val)) {
        val.oneshort <- val[-i, ]
        val.missing <- val[i, ]
        d.oneshort <- d[-i]
        d.missing <- d[i]
        A.oneshort <- A[-i, -i]
        Ainv.oneshort <- solve(A.oneshort)
        d.missing.estimated[i] <- interpolant(val.missing, d.oneshort, 
            val.oneshort, Ainv = Ainv.oneshort, scales = scales, 
            pos.def.matrix = pos.def.matrix, power = power, g = TRUE)$mstar.star
    }
    return(d.missing.estimated)
}
"interpolant" <-
function (x, d, xold, Ainv = NULL, A = NULL, use.Ainv = TRUE, 
    scales = NULL, pos.def.matrix = NULL, func = regressor.basis, 
    power = 2, give.full.list = FALSE) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix (used to calculate tx)")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    if (is.null(A)) {
        A <- corr.matrix(xold, pos.def.matrix = pos.def.matrix, 
            power = power)
    }
    if (is.null(Ainv) & use.Ainv) {
        Ainv <- solve(A)
    }
    tx <- apply(xold, 1, corr, x2 = x, pos.def.matrix = pos.def.matrix, 
        power = power)
    hx <- unlist(func(x))
    H <- regressor.multi(xold, func = func)
    if (use.Ainv) {
        n <- nrow(Ainv)
        jj <- betahat.fun(xold, Ainv, d, give.variance = TRUE, 
            func = func)
        betahat <- jj$betahat
        beta.var <- jj$variance
        sigmahat.square <- jj$sigmahatsquared
        beta.marginal.sd <- sqrt(diag(beta.var))
        prior <- crossprod(hx,betahat)
        mstar.star <- prior + crossprod(crossprod(Ainv, 
            tx), d - H %*% betahat)
        cstar.x.x <- corr(x, x, pos.def.matrix = pos.def.matrix, 
            power = power) - quad.form(Ainv, tx)
        cstar.star <- cstar.x.x + quad.form.inv(quad.form(Ainv, 
            H), hx - crossprod(H, crossprod(Ainv, tx)))
    }
    else {
        n <- nrow(A)
        jj <- betahat.fun.A(xold, A, d, give.variance = TRUE, 
            func = func)
        betahat <- jj$betahat
        sigmahat.square <- jj$sigmahatsquared
        beta.var <- jj$variance
        beta.marginal.sd <- sqrt(diag(beta.var))
        prior <- drop(crossprod(hx,betahat))
        mstar.star <- drop(prior + crossprod(solve(A, 
            tx), d - H %*% betahat))
        cstar.x.x <- corr(x, x, pos.def.matrix = pos.def.matrix, 
            power = power) - quad.form.inv(A, tx)
        cstar.star <- cstar.x.x + quad.form.inv(quad.form.inv(A, 
            H), hx - crossprod(H, solve(A, tx)))
    }
    if (give.full.list) {
        return(list(betahat = betahat, prior = prior,
                    beta.var = beta.var, beta.marginal.sd = beta.marginal.sd, 
            sigmahat.square = sigmahat.square, mstar.star = mstar.star, 
            cstar = cstar.x.x, cstar.star = cstar.star, Z = sqrt(abs(sigmahat.square * 
                cstar.star))))
    }
    else {
        return(as.vector(mstar.star))
    }
}
"interpolant.quick" <-
function (x, d, xold, Ainv, scales = NULL, pos.def.matrix = NULL, 
    func = regressor.basis, give.Z = FALSE, power = 2) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    betahat <- betahat.fun(xold, Ainv, d, func = func)
    H <- regressor.multi(xold, func = func)
    mstar.star <- rep(NA, nrow(x))
    prior <- rep(NA,nrow(x))

    if (give.Z) {
        Z <- rep(NA, nrow(x))
        sigmahat.square <- sigmahatsquared(H, Ainv, d)
    }
    for (i in 1:nrow(x)) {
        hx <- func(x[i, ])
        tx <- apply(xold, 1, corr, x2 = x[i, ], pos.def.matrix = pos.def.matrix, 
            power = power)
        prior[i] <- crossprod(hx,betahat)
        mstar.star[i] <- prior[i] + crossprod(crossprod(Ainv, 
            tx), (d - H %*% betahat))
        if (give.Z) {
            cstar.x.x <- 1 - quad.form(Ainv, tx)
            cstar.star <- cstar.x.x + quad.form.inv(quad.form(Ainv, 
                H), hx - crossprod(H, crossprod(Ainv, tx)))
            Z[i] <- sqrt(abs(sigmahat.square * cstar.star))
        }
    }
    if (give.Z) {
        return(list(mstar.star = mstar.star, Z = Z, prior = prior))
    }
    else {
        return(mstar.star)
    }
}

"var.conditional" <-
function (x, xold, d, A, Ainv, scales = NULL, pos.def.matrix = NULL, 
    func = regressor.basis, power = 2) {
    if (is.null(scales) & is.null(pos.def.matrix)) {
      stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
      stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
      pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    
    H <- regressor.multi(xold, func = func)
    hx <- regressor.multi(x,func=func)
    txvec <- 
      apply(x,1,function(y){apply(xold,1,function(x){corr(x,y,pos.def.matrix=pos.def.matrix, power=power)})})
    bit1 <- quad.form(Ainv,txvec)
    jj <- hx - crossprod(txvec,Ainv) %*% H
    bit2 <- quad.form.inv(quad.form(Ainv,H),t(jj))
    
    cstar <- corr.matrix(xold=x,pos.def.matrix=pos.def.matrix) - bit1 + bit2
    cstar <- as.matrix(cstar) 
    
    rownames(cstar) <- rownames(x)
    colnames(cstar) <- rownames(x)
    
    return(
           cstar*
           sigmahatsquared(H=regressor.multi(xold, func = func),
                           Ainv=Ainv, d=d)
           )
}

"cond.sample" <- function (n=1, x, xold, d, A, Ainv, scales = NULL, pos.def.matrix = NULL, 
    func = regressor.basis, power = 2) {
  mstar <- 
    interpolant.quick(x=x, d=d, xold=xold, Ainv=Ainv, scales=scales,
                      pos.def.matrix=pos.def.matrix,func=func,
                      give.Z=FALSE)

jj.sigma <- var.conditional(x, xold, d, A, Ainv, scales = scales, pos.def.matrix = pos.def.matrix, 
               func = func, power = power)

  random.bit <- 
    rmvt(n=n,sigma=jj.sigma, df=length(d))

  out <- sweep(random.bit,2,mstar,"+")
  colnames(out) <- rownames(x)
  return(out)
}

"latin.hypercube" <-
function (n, d, normalize = FALSE) 
{
    if (normalize) {
        f <- function(...) {
            sample(0:(n - 1))/(n - 1)
        }
    }
    else {
        f <- function(...) {
            (sample(1:n) - 0.5)/n
        }
    }
    out <- sapply(1:d, f)
    colnames(out) <- letters[1:d]
    return(out)
}
"makeinputfiles" <-
function (number.of.runs = 100, gaussian = TRUE, directoryname = "~/goldstein/genie-cgoldstein/", 
    filename = "QWERTYgoin", expert.estimates, area.outside = 0.05) 
{
    pad <- function(y, len, padchar = "0") {
        paste(paste(rep(padchar, len - nchar(y)), collapse = ""), 
            y, collapse = "", sep = "")
    }
    a <- expert.estimates
    minimum.values <- a$low
    maximum.values <- a$high
    lh.normalized.uniformdist <- latin.hypercube(number.of.runs, 
        nrow(a))
    lh.real <- lh.normalized.uniformdist
    lh.real[] <- NA
    if (gaussian) {
        for (i in 1:ncol(lh.real)) {
            a <- minimum.values[i]
            b <- maximum.values[i]
            if (a > 0) {
                lh.real[, i] <- exp(qnorm(p = lh.normalized.uniformdist[, 
                  i], mean = (log(a) + log(b))/2, sd = (log(b) - 
                  log(a))/(2 * qnorm(1 - area.outside/2))))
            }
            else {
                lh.real[, i] <- qunif(p = lh.normalized.uniformdist[, 
                  i], min = 0, max = b)
            }
        }
        lh.real[is.nan(lh.real)] <- 0
    }
    else {
        lh.normalized <- lh.normalized.uniformdist
        unnormalize <- function(normalized.vector) {
            minimum.values + normalized.vector * (maximum.values - 
                minimum.values)
        }
        lh.real <- t(apply(lh.normalized, 1, unnormalize))
    }
    write.table(x = lh.real, file = paste(directoryname, "list_of_inputfiles.txt", 
        sep = ""), quote = FALSE, row.names = TRUE, col.names = FALSE)
    for (i in 1:number.of.runs) {
        fullfilename <- paste(directoryname, filename, ".", i - 
            1, sep = "")
        f <- function(string, append = TRUE) {
            write(string, file = fullfilename, append = append)
        }
        f("400001 400000 400000  4000  400", append = FALSE)
        f("n")
        f("100   5")
        f("20.0")
        f("20.0")
        f("0.90")
        for (j in 1:8) {
            f(lh.real[i, j])
        }
        f("0.0")
        f("0.0")
        for (j in 9:11) {
            f(lh.real[i, j])
        }
        f("0.0")
        for (j in 12:12) {
            f(lh.real[i, j])
        }
        f("1")
        for (j in 13:14) {
            f(lh.real[i, j])
        }
        f("0.0")
        f("0.0")
        f("0.0")
        for (j in 16:17) {
            f(lh.real[i, j])
        }
        f(paste("tmp/", pad(i - 1, 3), sep = ""))
        f("tmp/tmp.avg")
    }
    return(0)
}
"model" <-
function (x) 
{
    if (1 == 1) {
        centrepoint <- x
        x[] <- 4
        return(sum((x - centrepoint)^2))
    }
    else {
        return(as.real(x[1] < x[2]))
    }
}
"optimal.scales" <-
function (val, scales.start, d, use.like = TRUE, give.answers = FALSE, 
    ...) 
{
    if (use.like) {
        objective.fun <- function(scales, val, d) {
            -scales.likelihood(scales = exp(scales), xold = val, 
                d = d)
        }
    }
    else {
        objective.fun <- function(scales, val, d) {
            A <- corr.matrix(val, scales = exp(scales))
            error <- abs(d - estimator(val, A, d, scales = exp(scales)))
            return(sum(error^2))
        }
    }
    jj <- optim(par = log(scales.start), objective.fun, val = val, 
        d = d, ...)
    if (give.answers) {
        return(jj)
    }
    else {
        return(exp(jj$par))
    }
}
"pad" <-
function (x, len, padchar = "0", strict = TRUE) 
{
    n <- nchar(x)
    if (nchar(padchar) > 1) {
        stop("padchar must be a single character")
    }
    if (n > len) {
        if (strict) {
            stop("input arg too long")
        }
        else {
            return(substr(as.character(x), n - len + 1, n))
        }
    }
    return(paste(paste(rep(padchar, len - n), collapse = ""), 
        x, sep = "", collapse = ""))
}
"prior.B" <-
function (H, Ainv, B0 = NULL) 
{
    if (is.null(B0)) {
        return(quad.form(Ainv, H))
    }
    else {
        return(B0 + quad.form(Ainv, H))
    }
}
"prior.b" <-
function (H, Ainv, d, b0 = NULL, B0 = NULL) 
{
    B <- prior.B(H, Ainv, B0)
    if (is.null(b0)) {
        b <- crossprod(solve(B), crossprod(H, Ainv)) %*% d
    }
    else {
        b <- solve(B) %*% (B0 %*% b0 + crossprod(H, Ainv) %*% 
            d)
    }
    return(b)
}
"quad.form" <-
function (M, x, chol = FALSE) 
{
    if (chol == FALSE) {
        return(drop(crossprod(crossprod(M, x), x)))
    }
    else {
        jj <- crossprod(M, x)
        return(drop(crossprod(jj, jj)))
    }
}
"quad.form.inv" <-
function (M, x) 
{
    drop(crossprod(x, solve(M, x)))
}
"regressor.basis" <-
function (x) 
{
    x <- c(1, x)
    names(x)[1] <- "const"
    return(x)
}
"regressor.multi" <-
function (x.df, func = regressor.basis) 
{
  out <- t(apply(x.df, 1, func))
  if(nrow(out) == 1){
    return(t(out))
  } else {
    return(out)
  }
}
"s.chi" <-
function (H, Ainv, d, s0 = 0, fast.but.opaque = TRUE) 
{
    if (fast.but.opaque) {
        out <- s0 + quad.form(Ainv - quad.form(quad.form(solve(quad.form(Ainv, 
            H)), t(H)), Ainv), d)
    }
    else {
        out <- s0 + t(d) %*% (Ainv - Ainv %*% H %*% solve(t(H) %*% 
            Ainv %*% H) %*% t(H) %*% Ainv) %*% d
    }
    return(out)
}
"sample.from.exp.est" <-
function (number.of.runs, expert.estimates, gaussian = TRUE, 
    area.outside = 0.05) 
{
    a <- expert.estimates
    minimum.values <- a$low
    maximum.values <- a$high
    lh.normalized.uniformdist <- latin.hypercube(number.of.runs, 
        16)
    lh.real <- lh.normalized.uniformdist
    lh.real[] <- NA
    if (gaussian) {
        for (i in 1:ncol(lh.real)) {
            a <- minimum.values[i]
            b <- maximum.values[i]
            lh.real[, i] <- exp(qnorm(p = lh.normalized.uniformdist[, 
                i], mean = (log(a) + log(b))/2, sd = (log(b) - 
                log(a))/(2 * qnorm(1 - area.outside/2))))
        }
        lh.real[is.nan(lh.real)] <- 0
    }
    else {
        lh.normalized <- lh.normalized.uniformdist
        unnormalize <- function(normalized.vector) {
            minimum.values + normalized.vector * (maximum.values - 
                minimum.values)
        }
        lh.real <- t(apply(lh.normalized, 1, unnormalize))
    }
    return(lh.real)
}
"scales.likelihood" <-
function (pos.def.matrix = NULL, scales = NULL, power = 2, xold,
          use.Ainv = TRUE, d, func = regressor.basis) 
{
    if (is.null(scales) & is.null(pos.def.matrix)) {
        stop("need either scales or a pos.definite.matrix")
    }
    if (!is.null(scales) & !is.null(pos.def.matrix)) {
        stop("scales *and* pos.def.matrix supplied.  corr() needs one only.")
    }
    if (is.null(pos.def.matrix)) {
        pos.def.matrix <- diag(scales,nrow=length(scales))
    }
    H <- regressor.multi(xold, func = func)
    q <- ncol(H) - 1
    n <- nrow(H)
    A <- corr.matrix(xold, pos.def.matrix = pos.def.matrix, power =
                     power)
    bit2 <- 1/sqrt(det(A))

    if(use.Ainv){
      Ainv <- solve(A)
      bit1 <- sigmahatsquared(H, Ainv, d)^(-(n - q)/2)
      bit3 <- 1/sqrt(det(quad.form(Ainv, H)))
    } else {
      bit1 <- sigmahatsquared.A(H, A, d)^(-(n - q)/2)
      bit3 <- 1/sqrt(det(quad.form.inv(A, H)))
    }
    return(drop(bit1 * bit2 * bit3))
}
"sigmahatsquared" <-
function (H, Ainv, d) 
{
    n <- nrow(Ainv)
    q <- ncol(H) - 1
    H <- as.matrix(H)
    out <- quad.form(Ainv - quad.form.inv(quad.form(Ainv, H), 
        crossprod(H, Ainv)), d)/(n - q - 2)
    return(out)
}
"sigmahatsquared.A" <-
function (H, A, d) 
{
    n <- nrow(A)
    q <- ncol(H) - 1
    (quad.form.inv(A, d) - quad.form(quad.form.inv(quad.form.inv(A, 
        H), t(solve(A, H))), d))/(n - q - 2)
}
