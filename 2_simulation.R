mvqtl_simulation <- function(qtl,
                             var_covar,
                             seed) {


  # set up the environment
  {
    # misc functions
    {
      cat <- function() {}

      rint <- function(v, c = 3/8, ...) {
        qnorm(p = (rank(v) - c)/(length(v) - 2*c + 1), ...)
      }

      tryNA <- function(expr) {
        suppressMessages(
          suppressWarnings(
            tryCatch(expr = expr,
                     error = function(e) NA,
                     finally = NA)))
      }

      tryNULL <- function(expr) {
        suppressMessages(
          suppressWarnings(
            tryCatch(expr = expr,
                     error = function(e) NULL,
                     finally = NULL)))
      }

      try0 <- function(expr) {
        suppressMessages(
          suppressWarnings(
            tryCatch(expr = expr,
                     error = function(e) 0,
                     finally = 0)))
      }
    }

    # LM testing function
    {
      lm_test <- function(resp, locus) {
        LL1 <- logLik(lm(formula = resp ~ locus))
        LL0 <- logLik(lm(formula = resp ~ 1))
        return(2*(LL1 - LL0))
      }
    }

    # Cao testing functions
    {
      addiVar <- function(gt, depvar, covariates) {
        temp <- cbind(depvar, covariates)
        p <- ncol(temp)
        lm1 <- lm(depvar~., data=as.data.frame(temp))
        betainit <- as.numeric(coef(lm1))
        sigmainit <- summary(lm1)$sigma
        designmat <- cbind(rep(1, length(gt)),covariates)
        loglik <- function(param) {
          sigmasq <- param[p+1]+gt*param[p+2]
          res <- depvar-as.matrix(designmat)%*%matrix(param[1:p], ncol=1)
          sum(log(sigmasq))+t(res)%*%diag(1/sigmasq)%*%res
        }
        ui <- rbind(c(rep(0, p), 1, 0),c(rep(0, p), 1, 2))
        ci <- c(0, 0)
        mle <- constrOptim(theta = c(betainit, sigmainit^2, sigmainit^2*0.1),
                           f = loglik, ui=ui, ci=ci, method="Nelder-Mead")
        -0.5*(length(gt)*log(2*pi)+loglik(mle$par))
      }

      LRT_mean_var <- function(gt,
                               depvar,
                               covariates = NULL,
                               meanModel = "genotypic",
                               varModel = "genotypic",
                               test = c("mean","var","both"),
                               ...)
      {
        ## check if there is missing data
        dataAll <- cbind(gt, depvar, covariates)
        if(sum(is.na(dataAll)) > 0)  stop("missing value found")

        ## check the number of different genotypes
        total.rHomo <- sum(gt==2)
        total.rare <- sum(gt==2)+sum(gt==1)
        if (total.rare==0) stop("only one genetic group")
        if (total.rHomo==0) {
          meanModel <- "dominant"
          # cat("only two genotypes, using dominant model \n")
        }


        data1 <- as.data.frame(cbind(depvar, covariates))

        if(meanModel=="genotypic") {
          hetero <- rHomo <- rep(0, length(gt))
          hetero[gt==1] <- 1
          rHomo[gt==2] <- 1
          data2 <- as.data.frame(cbind(depvar, hetero, rHomo, covariates))
        }

        if(meanModel=="additive")  data2 <- as.data.frame(cbind(depvar, gt, covariates))

        if(meanModel=="dominant") {
          gt[gt==2] <- 1
          data2 <- as.data.frame(cbind(depvar, gt, covariates))
        }

        gtf <<- as.factor(gt)

        if (varModel=="additive") {
          if (meanModel!="additive") stop("variance model can only be additive when mean model is additive")
          if (test=="mean") {
            m.null <- addiVar(gt=gt, depvar=depvar, covariates=covariates)
            m.alter <- addiVar(gt=gt, depvar=depvar, covariates=cbind(gt, covariates))
            dfreedom <- 1
          }

          if (test=="var") {
            m.null <- as.numeric(logLik(lm(depvar~., data=data2)))
            m.alter <- addiVar(gt=gt, depvar=depvar, covariates=cbind(gt, covariates))
            dfreedom <- 1
          }

          if (test=="both") {
            m.null <- as.numeric(logLik(lm(depvar~., data=data1)))
            m.alter <- addiVar(gt=gt, depvar=depvar, covariates=cbind(gt, covariates))
            dfreedom <- 2
          }
        } else {
          if (test=="mean") {
            gls.null <- gls(depvar~., data=data1, weights=varIdent(form = ~ 1|gtf), method="ML", control=glsControl(opt="optim"), ...)
            gls.alter <- gls(depvar~., data=data2, weights=varIdent(form = ~ 1|gtf), method="ML", control=glsControl(opt="optim"), ...)
            m.null <- as.numeric(logLik(gls.null))
            m.alter <- as.numeric(logLik(gls.alter))
            if (meanModel=="genotypic") dfreedom <- 2  else  dfreedom <- 1
          }

          if (test=="var") {
            gls.null <- lm(depvar~., data=data2)
            gls.alter <- gls(depvar~., data=data2, weights=varIdent(form = ~ 1|gtf), method="ML", control=glsControl(opt="optim"), ...)
            m.null <- as.numeric(logLik(gls.null))
            m.alter <- as.numeric(logLik(gls.alter))
            if (meanModel=="dominant") dfreedom <- 1  else  dfreedom <- 2
          }

          if (test=="both") {
            gls.null <- lm(depvar~., data=data1)
            gls.alter <- gls(depvar~., data=data2, weights=varIdent(form = ~ 1|gtf), method="ML", control=glsControl(opt="optim"), ...)
            m.null <- as.numeric(logLik(gls.null))
            m.alter <- as.numeric(logLik(gls.alter))
            if (meanModel=="genotypic") dfreedom <- 4
            if (meanModel=="additive") dfreedom <- 3
            if (meanModel=="dominant") dfreedom <- 2
          }
        }

        LRT_Stat <- 2*(m.alter-m.null)

        return(LRT_Stat)
      }
    }

    # DGLM testing functions
    {
      dglm_ml <- function(...) { -0.5*(dglm::dglm(..., method = 'ml')$m2loglik) }

      mean_test <- function(resp, mean_locus, var_locus, covar) {
        if (is.null(covar)) {
          LLj <- dglm_ml(formula = resp ~ mean_locus, dformula = ~ var_locus)
          LLv <- dglm_ml(formula = resp ~ 1,          dformula = ~ var_locus)
        } else {
          LLj <- dglm_ml(formula = resp ~ mean_locus, dformula = ~ var_locus + covar)
          LLv <- dglm_ml(formula = resp ~ 1,          dformula = ~ var_locus + covar)
        }
        return(2*(LLj - LLv))
      }

      var_test <- function(resp, mean_locus, var_locus, covar) {
        if (is.null(covar)) {
          LLj <- dglm_ml(formula = resp ~ mean_locus, dformula = ~ var_locus)
          LLm <- dglm_ml(formula = resp ~ mean_locus, dformula = ~ 1)
        } else {
          LLj <- dglm_ml(formula = resp ~ mean_locus, dformula = ~ var_locus + covar)
          LLm <- dglm_ml(formula = resp ~ mean_locus, dformula = ~ covar)
        }
        return(2*(LLj - LLm))
      }

      joint_test <- function(resp, locus, covar) {
        if (is.null(covar)) {
          LLj <- dglm_ml(formula = resp ~ locus, dformula = ~ locus)
          LL0 <- dglm_ml(formula = resp ~ 1,     dformula = ~ 1)
        } else {
          LLj <- dglm_ml(formula = resp ~ locus, dformula = ~ locus + covar)
          LL0 <- dglm_ml(formula = resp ~ 1,     dformula = ~ covar)
        }
        return(2*(LLj - LL0))
      }
    }
  }

  set.seed(seed)

  LM__asymp_ps    <- LM__rint_ps    <- LM__residperm_ps    <- LM__locusperm_ps    <- rep(NA, num_sims_per_split)
  Caom__asymp_ps  <- Caom__rint_ps  <- Caom__residperm_ps  <- Caom__locusperm_ps <-  rep(NA, num_sims_per_split)
  DGLMm_no_covar__asymp_ps <- DGLMm_no_covar__rint_ps <- DGLMm_no_covar__residperm_ps <- DGLMm_no_covar__locusperm_ps <- rep(NA, num_sims_per_split)
  DGLMm_ye_covar__asymp_ps <- DGLMm_ye_covar__rint_ps <- DGLMm_ye_covar__residperm_ps <- DGLMm_ye_covar__locusperm_ps <- rep(NA, num_sims_per_split)

  Lev__asymp_ps   <- Lev__rint_ps   <- Lev__residperm_ps   <- Lev__locusperm_ps   <- rep(NA, num_sims_per_split)
  Caov__asymp_ps  <- Caov__rint_ps  <- Caov__residperm_ps  <- Caov__locusperm_ps <-  rep(NA, num_sims_per_split)
  DGLMv_no_covar__asymp_ps <- DGLMv_no_covar__rint_ps <- DGLMv_no_covar__residperm_ps <- DGLMv_no_covar__locusperm_ps <- rep(NA, num_sims_per_split)
  DGLMv_ye_covar__asymp_ps <- DGLMv_ye_covar__rint_ps <- DGLMv_ye_covar__residperm_ps <- DGLMv_ye_covar__locusperm_ps <- rep(NA, num_sims_per_split)

  Caoj__asymp_ps  <- Caoj__rint_ps  <- Caoj__residperm_ps  <- Caoj__locusperm_ps  <- rep(NA, num_sims_per_split)
  DGLMj_no_covar__asymp_ps <- DGLMj_no_covar__rint_ps <- DGLMj_no_covar__residperm_ps <- DGLMj_no_covar__locusperm_ps <- rep(NA, num_sims_per_split)
  DGLMj_ye_covar__asymp_ps <- DGLMj_ye_covar__rint_ps <- DGLMj_ye_covar__residperm_ps <- DGLMj_ye_covar__locusperm_ps <- rep(NA, num_sims_per_split)

  LM_residperm_num_perms    <- LM_locusperm_num_perms    <- rep(NA, num_sims_per_split)
  Caom_residperm_num_perms  <- Caom_locusperm_num_perms <-  rep(NA, num_sims_per_split)
  DGLMm_no_covar_residperm_num_perms <- DGLMm_no_covar_locusperm_num_perms <- rep(NA, num_sims_per_split)
  DGLMm_ye_covar_residperm_num_perms <- DGLMm_ye_covar_locusperm_num_perms <- rep(NA, num_sims_per_split)

  Lev_residperm_num_perms   <- Lev_locusperm_num_perms   <- rep(NA, num_sims_per_split)
  Caov_residperm_num_perms  <- Caov_locusperm_num_perms <-  rep(NA, num_sims_per_split)
  DGLMv_no_covar_residperm_num_perms <- DGLMv_no_covar_locusperm_num_perms <- rep(NA, num_sims_per_split)
  DGLMv_ye_covar_residperm_num_perms <- DGLMv_ye_covar_locusperm_num_perms <- rep(NA, num_sims_per_split)

  Caoj_residperm_num_perms  <- Caoj_locusperm_num_perms  <- rep(NA, num_sims_per_split)
  DGLMj_no_covar_residperm_num_perms <- DGLMj_no_covar_locusperm_num_perms <- rep(NA, num_sims_per_split)
  DGLMj_ye_covar_residperm_num_perms <- DGLMj_ye_covar_locusperm_num_perms <- rep(NA, num_sims_per_split)

  for (sim_idx in 1:num_sims_per_split) {

    message('Starting sim ', sim_idx, ' of ', num_sims_per_split)
    locus <- sample(x = c(-1, 0, 1), size = num_obs, replace = TRUE, prob = c(0.25, 0.5, 0.25))
    covar <- rep(x = 1:5, times = num_obs/5)

    # SIMULATION

    # set up residual variance of phenotype
    {
      if (qtl %in% c('none', 'mqtl')) {
        a <- 0
        d <- 0
      }

      if (qtl == 'vqtl') {
        a <- vqtl_effect
        d <- 0
      }

      if (qtl == 'mvqtl') {
        a <- mvqtl_var_effect
        d <- 0
      }

      locus_linpred <- locus*a

      if (var_covar) {
        covar_effect <- seq(from = -0.4, to = 0.4, length.out = 5)[covar]
      } else {
        covar_effect <- 0
      }

      vars <- exp(2*(locus_linpred + covar_effect))
      message('vars:', head(vars))
    }

    # set up expected value of phenotype
    {
      mean_resid_var <- mean(vars)

      if (qtl %in% c('none', 'vqtl')) {
        a <- 0
        d <- 0
      }

      if (qtl == 'mqtl') {
        a <- sqrt(mqtl_pve*mean_resid_var/(1-mqtl_pve))
        d <- 0
      }

      if (qtl == 'mvqtl') {
        a <- sqrt(mvqtl_mean_pve*mean_resid_var/(1-mvqtl_mean_pve))
        d <- 0
      }

      e_y <- (locus==(-1))*-a+(locus==0)*d + (locus==1)*a
      message('Es:', head(e_y))
    }

    # simulate phenotype
    y <- rnorm(n = num_obs,
               mean = e_y,
               sd = sqrt(vars))

    message('y:', head(y))

    locus <- factor(locus)
    locus_df <- length(unique(locus)) - 1

    covar <- factor(covar)


    # TESTING

    # mQTL tests
    if (TRUE) {

      # compute test stats
      LM_LRS <- tryNA(lm_test(resp = y, locus = locus))
      LM_rint_LRS <- tryNA(lm_test(resp = rint(y), locus = locus))

      Caom_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = y, test = "mean"))
      Caom_rint_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = rint(y), test = "mean"))

      DGLMm_no_covar_LRS <- tryNA(mean_test(resp = y, mean_locus = locus, var_locus = locus, covar = NULL))
      DGLMm_rint_no_covar_LRS <- tryNA(mean_test(resp = rint(y), mean_locus = locus, var_locus = locus, covar = NULL))

      DGLMm_ye_covar_LRS <- tryNA(mean_test(resp = y, mean_locus = locus, var_locus = locus, covar = covar))
      DGLMm_rint_ye_covar_LRS <- tryNA(mean_test(resp = rint(y), mean_locus = locus, var_locus = locus, covar = covar))

      # perms for empirical evaluation
      if (num_perms > 0) {

        message('starting mean perms')

        lm_resids <- residuals(lm(y ~ 1))
        lm_hats <- y - lm_resids

        cao_resids <- residuals(dglm(formula = y ~ 1, dformula = ~ locus))
        cao_hats <- y - cao_resids

        dglm_no_covar_resids <- cao_resids # no difference, why calculate?
        dglm_no_covar_hats <- y - dglm_no_covar_resids

        dglm_ye_covar_resids <- residuals(dglm(formula = y ~ 1, dformula = ~ locus + covar))
        dglm_ye_covar_hats <- y - dglm_ye_covar_resids

        LM_LRSs_residperm <- Caom_LRSs_residperm <- DGLMm_no_covar_LRSs_residperm <- DGLMm_ye_covar_LRSs_residperm <- rep(NA, num_perms)
        LM_LRSs_locusperm <- Caom_LRSs_locusperm <- DGLMm_no_covar_LRSs_locusperm <- DGLMm_ye_covar_LRSs_locusperm <- rep(NA, num_perms)

        for (perm_idx in 1:num_perms) {

          lm_PB <- lm_hats + sample(lm_resids)
          cao_PB <- cao_hats + sample(cao_resids)
          dglm_no_covar_PB <- dglm_no_covar_hats + sample(dglm_no_covar_resids)
          dglm_ye_covar_PB <- dglm_ye_covar_hats + sample(dglm_ye_covar_resids)

          LM_LRSs_residperm[perm_idx] <- tryNA(lm_test(resp = lm_PB, locus = locus))
          Caom_LRSs_residperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = cao_PB, test = "mean"))
          DGLMm_no_covar_LRSs_residperm[perm_idx] <- tryNA(mean_test(resp = dglm_no_covar_PB, mean_locus = locus, var_locus = locus, covar = NULL))
          DGLMm_ye_covar_LRSs_residperm[perm_idx] <- tryNA(mean_test(resp = dglm_ye_covar_PB, mean_locus = locus, var_locus = locus, covar = covar))

          LM_LRSs_locusperm[perm_idx] <- tryNA(lm_test(resp = y, locus = sample(locus)))
          Caom_LRSs_locusperm[perm_idx] <- tryNA(mean_test(resp = y, mean_locus = sample(locus), var_locus = locus, covar = NULL))
          DGLMm_no_covar_LRSs_locusperm[perm_idx] <- tryNA(mean_test(resp = y, mean_locus = sample(locus), var_locus = locus, covar = NULL))
          DGLMm_ye_covar_LRSs_locusperm[perm_idx] <- tryNA(mean_test(resp = y, mean_locus = sample(locus), var_locus = locus, covar = covar))

        }
      }

      # evaluate test stats
      LM__asymp_ps[sim_idx] <- pchisq(q = LM_LRS, df = 2, lower.tail = FALSE)
      LM__rint_ps[sim_idx] <- pchisq(q = LM_rint_LRS, df = 2, lower.tail = FALSE)
      if (num_perms > 0) {
        LM__residperm_ps[sim_idx] <- mean(x = LM_LRSs_residperm >= LM_LRS, na.rm = TRUE)
        LM__locusperm_ps[sim_idx] <- mean(x = LM_LRSs_locusperm >= LM_LRS, na.rm = TRUE)
        LM_residperm_num_perms[sim_idx] <- sum(!is.na(LM_LRSs_residperm))
        LM_locusperm_num_perms[sim_idx] <- sum(!is.na(LM_LRSs_locusperm))
      }

      Caom__asymp_ps[sim_idx] <- pchisq(q = Caom_LRS, df = locus_df, lower.tail = FALSE)
      Caom__rint_ps[sim_idx] <- pchisq(q = Caom_rint_LRS, df = locus_df, lower.tail = FALSE)
      if (num_perms > 0) {
        Caom__residperm_ps[sim_idx] <- mean(x = Caom_LRSs_residperm >= Caom_LRS, na.rm = TRUE)
        Caom_residperm_num_perms[sim_idx] <- sum(!is.na(Caom_LRSs_residperm))

        Caom__locusperm_ps[sim_idx] <- mean(x = Caom_LRSs_locusperm >= Caom_LRS, na.rm = TRUE)
        Caom_locusperm_num_perms[sim_idx] <- sum(!is.na(Caom_LRSs_locusperm))
      }

      DGLMm_no_covar__asymp_ps[sim_idx] <- pchisq(q = DGLMm_no_covar_LRS, df = locus_df, lower.tail = FALSE)
      DGLMm_no_covar__rint_ps[sim_idx] <- pchisq(q = DGLMm_rint_no_covar_LRS, df = locus_df, lower.tail = FALSE)
      DGLMm_ye_covar__asymp_ps[sim_idx] <- pchisq(q = DGLMm_ye_covar_LRS, df = locus_df, lower.tail = FALSE)
      DGLMm_ye_covar__rint_ps[sim_idx] <- pchisq(q = DGLMm_rint_ye_covar_LRS, df = locus_df, lower.tail = FALSE)

      if (num_perms > 0) {
        DGLMm_no_covar__residperm_ps[sim_idx] <- mean(x = DGLMm_no_covar_LRSs_residperm >= DGLMm_no_covar_LRS, na.rm = TRUE)
        DGLMm_no_covar_residperm_num_perms[sim_idx] <- sum(!is.na(DGLMm_no_covar_LRSs_residperm))

        DGLMm_no_covar__locusperm_ps[sim_idx] <- mean(x = DGLMm_no_covar_LRSs_locusperm >= DGLMm_no_covar_LRS, na.rm = TRUE)
        DGLMm_no_covar_locusperm_num_perms[sim_idx] <- sum(!is.na(DGLMm_no_covar_LRSs_locusperm))

        DGLMm_ye_covar__residperm_ps[sim_idx] <- mean(x = DGLMm_ye_covar_LRSs_residperm >= DGLMm_ye_covar_LRS, na.rm = TRUE)
        DGLMm_ye_covar_residperm_num_perms[sim_idx] <- sum(!is.na(DGLMm_ye_covar_LRSs_residperm))

        DGLMm_ye_covar__locusperm_ps[sim_idx] <- mean(x = DGLMm_ye_covar_LRSs_locusperm >= DGLMm_ye_covar_LRS, na.rm = TRUE)
        DGLMm_ye_covar_locusperm_num_perms[sim_idx] <- sum(!is.na(DGLMm_ye_covar_LRSs_locusperm))

      }
    }

    # vQTL tests
    if (TRUE) { #any(locus_var_effect_add, !locus_mean_effect_add)) {

      # compute test stats
      Lev_stat <-  tryNA(leveneTest(y = y, group = locus)$`F value`[1])
      Lev_rint_stat <- tryNA(leveneTest(y = rint(y), group = locus)$`F value`[1])

      Caov_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = y, test = "var"))
      Caov_rint_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = rint(y), test = "var"))

      DGLMv_no_covar_LRS <- tryNA(var_test(resp = y, mean_locus = locus, var_locus = locus, covar = NULL))
      DGLMv_no_covar_rint_LRS <- tryNA(var_test(resp = rint(y), mean_locus = locus, var_locus = locus, covar = NULL))

      DGLMv_ye_covar_LRS <- tryNA(var_test(resp = y, mean_locus = locus, var_locus = locus, covar = covar))
      DGLMv_ye_covar_rint_LRS <- tryNA(var_test(resp = rint(y), mean_locus = locus, var_locus = locus, covar = covar))

      # perms for empirical evaluation
      if (num_perms > 0) {

        message('starting var perms')

        lm_resids <- residuals(lm(y ~ locus))
        lm_hats <- y - lm_resids

        cao_resids <- lm_resids #(no difference, why calculate?)
        cao_hats <- y - cao_resids

        dglm_no_covar_resids <- lm_resids #(no difference, why calculate?)
        dglm_no_covar_hats <- y - dglm_no_covar_resids

        dglm_ye_covar_resids <- residuals(dglm(formula = y ~ locus, dformula = ~ covar))
        dglm_ye_covar_hats <- y - dglm_ye_covar_resids

        Lev_stats_residperm <- Caov_LRSs_residperm <- DGLMv_no_covar_LRSs_residperm <- DGLMv_ye_covar_LRSs_residperm <- rep(NA, num_perms)
        Lev_stats_locusperm <- Caov_LRSs_locusperm <- DGLMv_no_covar_LRSs_locusperm <- DGLMv_ye_covar_LRSs_locusperm <- rep(NA, num_perms)
        for (perm_idx in 1:num_perms) {

          lm_PB <- lm_hats + sample(lm_resids)
          cao_PB <- cao_hats + sample(cao_resids)
          dglm_no_covar_PB <- dglm_no_covar_hats + sample(dglm_no_covar_resids)
          dglm_ye_covar_PB <- dglm_ye_covar_hats + sample(dglm_ye_covar_resids)

          Lev_stats_residperm[perm_idx] <- tryNA(leveneTest(y = lm_PB, group = locus)$`F value`[1])
          Caov_LRSs_residperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = cao_PB, test = "var"))
          DGLMv_no_covar_LRSs_residperm[perm_idx] <- tryNA(var_test(resp = dglm_no_covar_PB, mean_locus = locus, var_locus = locus, covar = NULL))
          DGLMv_ye_covar_LRSs_residperm[perm_idx] <- tryNA(var_test(resp = dglm_ye_covar_PB, mean_locus = locus, var_locus = locus, covar = covar))

          Lev_stats_locusperm[perm_idx] <- tryNA(leveneTest(y = y, group = sample(locus))$`F value`[1])
          Caov_LRSs_locusperm[perm_idx] <- tryNA(var_test(resp = y, mean_locus = locus, var_locus = sample(locus), covar = NULL))
          DGLMv_no_covar_LRSs_locusperm[perm_idx] <- tryNA(var_test(resp = y, mean_locus = locus, var_locus = sample(locus), covar = NULL))
          DGLMv_ye_covar_LRSs_locusperm[perm_idx] <- tryNA(var_test(resp = y, mean_locus = locus, var_locus = sample(locus), covar = covar))

        }
      }

      # evaluate test stats
      Lev__asymp_ps[sim_idx] <- pf(q = Lev_stat, df1 = locus_df, df2 = num_obs - locus_df, lower.tail = FALSE)
      Lev__rint_ps[sim_idx] <- pf(q = Lev_rint_stat, df1 = locus_df, df2 = num_obs - locus_df, lower.tail = FALSE)
      if (num_perms > 0) {
        Lev__residperm_ps[sim_idx] <- mean(x = Lev_stats_residperm >= Lev_stat, na.rm = TRUE)
        Lev__locusperm_ps[sim_idx] <- mean(x = Lev_stats_locusperm >= Lev_stat, na.rm = TRUE)
        Lev_residperm_num_perms[sim_idx] <- sum(!is.na(Lev_stats_residperm))
        Lev_locusperm_num_perms[sim_idx] <- sum(!is.na(Lev_stats_locusperm))
      }

      Caov__asymp_ps[sim_idx] <- pchisq(q = Caov_LRS, df = locus_df, lower.tail = FALSE)
      Caov__rint_ps[sim_idx] <- pchisq(q = Caov_rint_LRS, df = locus_df, lower.tail = FALSE)
      if (num_perms > 0) {
        Caov__residperm_ps[sim_idx] <- mean(x = Caov_LRSs_residperm >= Caov_LRS, na.rm = TRUE)
        Caov__locusperm_ps[sim_idx] <- mean(x = Caov_LRSs_locusperm >= Caov_LRS, na.rm = TRUE)
        Caov_residperm_num_perms[sim_idx] <- sum(!is.na(Caov_LRSs_residperm))
        Caov_locusperm_num_perms[sim_idx] <- sum(!is.na(Caov_LRSs_locusperm))
      }

      DGLMv_no_covar__asymp_ps[sim_idx] <- pchisq(q = DGLMv_no_covar_LRS, df = locus_df, lower.tail = FALSE)
      DGLMv_no_covar__rint_ps[sim_idx] <- pchisq(q = DGLMv_no_covar_rint_LRS, df = locus_df, lower.tail = FALSE)

      DGLMv_ye_covar__asymp_ps[sim_idx] <- pchisq(q = DGLMv_ye_covar_LRS, df = locus_df, lower.tail = FALSE)
      DGLMv_ye_covar__rint_ps[sim_idx] <- pchisq(q = DGLMv_ye_covar_rint_LRS, df = locus_df, lower.tail = FALSE)

      if (num_perms > 0) {

        DGLMv_no_covar__residperm_ps[sim_idx] <- mean(x = DGLMv_no_covar_LRSs_residperm >= DGLMv_no_covar_LRS, na.rm = TRUE)
        DGLMv_no_covar_residperm_num_perms[sim_idx] <- sum(!is.na(DGLMv_no_covar_LRSs_residperm))

        DGLMv_no_covar__locusperm_ps[sim_idx] <- mean(x = DGLMv_no_covar_LRSs_locusperm >= DGLMv_no_covar_LRS, na.rm = TRUE)
        DGLMv_no_covar_locusperm_num_perms[sim_idx] <- sum(!is.na(DGLMv_no_covar_LRSs_locusperm))

        DGLMv_ye_covar__residperm_ps[sim_idx] <- mean(x = DGLMv_ye_covar_LRSs_residperm >= DGLMv_ye_covar_LRS, na.rm = TRUE)
        DGLMv_ye_covar_residperm_num_perms[sim_idx] <- sum(!is.na(DGLMv_ye_covar_LRSs_residperm))

        DGLMv_ye_covar__locusperm_ps[sim_idx] <- mean(x = DGLMv_ye_covar_LRSs_locusperm >= DGLMv_ye_covar_LRS, na.rm = TRUE)
        DGLMv_ye_covar_locusperm_num_perms[sim_idx] <- sum(!is.na(DGLMv_ye_covar_LRSs_locusperm))
      }

    }

    # joint tests
    if (TRUE) {
      # compute test stats
      Caoj_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = y, test = "both"))
      Caoj_rint_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = rint(y), test = "both"))

      DGLMj_no_covar_LRS <- tryNA(joint_test(resp = y, locus = locus, covar = NULL))
      DGLMj_no_covar_rint_LRS <- tryNA(joint_test(resp = rint(y), locus = locus, covar = NULL))

      DGLMj_ye_covar_LRS <- tryNA(joint_test(resp = y, locus = locus, covar = covar))
      DGLMj_ye_covar_rint_LRS <- tryNA(joint_test(resp = rint(y), locus = locus, covar = covar))

      # perms for empirical evaluation
      if (num_perms > 0) {

        message('starting joint perms')

        lm_resids <- residuals(lm(y ~ 1))
        lm_hats <- y - lm_resids

        cao_resids <- lm_resids #(no difference, why calculate?)
        cao_hats <- y - cao_resids

        dglm_no_covar_resids <- lm_resids # no difference, why calculate?
        dglm_no_covar_hats <- y - lm_resids

        dglm_ye_covar_resids <- residuals(dglm(formula = y ~ 1, dformula = ~ covar))
        dglm_ye_covar_hats <- y - dglm_ye_covar_resids


        Caoj_LRSs_residperm <- DGLMj_no_covar_LRSs_residperm <- DGLMj_ye_covar_LRSs_residperm <- rep(NA, num_perms)
        Caoj_LRSs_locusperm <- DGLMj_no_covar_LRSs_locusperm <- DGLMj_ye_covar_LRSs_locusperm <- rep(NA, num_perms)
        for (perm_idx in 1:num_perms) {

          cao_PB <- cao_hats + sample(cao_resids)
          dglm_no_covar_PB <- dglm_no_covar_hats + sample(dglm_no_covar_resids)
          dglm_ye_covar_PB <- dglm_ye_covar_hats + sample(dglm_ye_covar_resids)

          Caoj_LRSs_residperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = cao_PB, test = "both"))
          DGLMj_no_covar_LRSs_residperm[perm_idx] <- tryNA(joint_test(resp = dglm_no_covar_PB, locus = locus, covar = NULL))
          DGLMj_ye_covar_LRSs_residperm[perm_idx] <- tryNA(joint_test(resp = dglm_ye_covar_PB, locus = locus, covar = covar))

          Caoj_LRSs_locusperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(sample(locus)) - 1, depvar = y, test = "both"))
          DGLMj_no_covar_LRSs_locusperm[perm_idx] <- tryNA(joint_test(resp = y, locus = sample(locus), covar = NULL))
          DGLMj_ye_covar_LRSs_locusperm[perm_idx] <- tryNA(joint_test(resp = y, locus = sample(locus), covar = covar))
        }
      }

      # evaluate test stats
      Caoj__asymp_ps[sim_idx] <- pchisq(q = Caoj_LRS, df = 2*locus_df, lower.tail = FALSE)
      Caoj__rint_ps[sim_idx] <- pchisq(q = Caoj_rint_LRS, df = 2*locus_df, lower.tail = FALSE)
      if (num_perms > 0) {
        Caoj__residperm_ps[sim_idx] <- mean(x = Caoj_LRSs_residperm >= Caoj_LRS, na.rm = TRUE)
        Caoj__locusperm_ps[sim_idx] <- mean(x = Caoj_LRSs_locusperm >= Caoj_LRS, na.rm = TRUE)
        Caoj_residperm_num_perms[sim_idx] <- sum(!is.na(Caoj_LRSs_residperm))
        Caoj_locusperm_num_perms[sim_idx] <- sum(!is.na(Caoj_LRSs_locusperm))
      }

      DGLMj_no_covar__asymp_ps[sim_idx] <- pchisq(q = DGLMj_no_covar_LRS, df = 2*locus_df, lower.tail = FALSE)
      DGLMj_no_covar__rint_ps[sim_idx] <- pchisq(q = DGLMj_no_covar_rint_LRS, df = 2*locus_df, lower.tail = FALSE)

      DGLMj_ye_covar__asymp_ps[sim_idx] <- pchisq(q = DGLMj_ye_covar_LRS, df = 2*locus_df, lower.tail = FALSE)
      DGLMj_ye_covar__rint_ps[sim_idx] <- pchisq(q = DGLMj_ye_covar_rint_LRS, df = 2*locus_df, lower.tail = FALSE)

      if (num_perms > 0) {

        DGLMj_no_covar__residperm_ps[sim_idx] <- mean(x = DGLMj_no_covar_LRSs_residperm >= DGLMj_no_covar_LRS, na.rm = TRUE)
        DGLMj_no_covar_residperm_num_perms[sim_idx] <- sum(!is.na(DGLMj_no_covar_LRSs_residperm))

        DGLMj_no_covar__locusperm_ps[sim_idx] <- mean(x = DGLMj_no_covar_LRSs_locusperm >= DGLMj_no_covar_LRS, na.rm = TRUE)
        DGLMj_no_covar_locusperm_num_perms[sim_idx] <- sum(!is.na(DGLMj_no_covar_LRSs_locusperm))

        DGLMj_ye_covar__residperm_ps[sim_idx] <- mean(x = DGLMj_ye_covar_LRSs_residperm >= DGLMj_ye_covar_LRS, na.rm = TRUE)
        DGLMj_ye_covar_residperm_num_perms[sim_idx] <- sum(!is.na(DGLMj_ye_covar_LRSs_residperm))

        DGLMj_ye_covar__locusperm_ps[sim_idx] <- mean(x = DGLMj_ye_covar_LRSs_locusperm >= DGLMj_ye_covar_LRS, na.rm = TRUE)
        DGLMj_ye_covar_locusperm_num_perms[sim_idx] <- sum(!is.na(DGLMj_ye_covar_LRSs_locusperm))

      }
    }

  }

  data_frame(seed,

             var_covar,
             qtl,

             LM__asymp_ps,
             LM__rint_ps,
             LM__residperm_ps,
             LM__locusperm_ps,

             Caom__asymp_ps,
             Caom__rint_ps,
             Caom__residperm_ps,
             Caom__locusperm_ps,

             DGLMm_no_covar__asymp_ps,
             DGLMm_no_covar__rint_ps,
             DGLMm_no_covar__residperm_ps,
             DGLMm_no_covar__locusperm_ps,

             DGLMm_ye_covar__asymp_ps,
             DGLMm_ye_covar__rint_ps,
             DGLMm_ye_covar__residperm_ps,
             DGLMm_ye_covar__locusperm_ps,

             Lev__asymp_ps,
             Lev__rint_ps,
             Lev__residperm_ps,
             Lev__locusperm_ps,

             Caov__asymp_ps,
             Caov__rint_ps,
             Caov__residperm_ps,
             Caov__locusperm_ps,

             DGLMv_no_covar__asymp_ps,
             DGLMv_no_covar__rint_ps,
             DGLMv_no_covar__residperm_ps,
             DGLMv_no_covar__locusperm_ps,

             DGLMv_ye_covar__asymp_ps,
             DGLMv_ye_covar__rint_ps,
             DGLMv_ye_covar__residperm_ps,
             DGLMv_ye_covar__locusperm_ps,

             Caoj__asymp_ps,
             Caoj__rint_ps,
             Caoj__residperm_ps,
             Caoj__locusperm_ps,

             DGLMj_no_covar__asymp_ps,
             DGLMj_no_covar__rint_ps,
             DGLMj_no_covar__residperm_ps,
             DGLMj_no_covar__locusperm_ps,

             DGLMj_ye_covar__asymp_ps,
             DGLMj_ye_covar__rint_ps,
             DGLMj_ye_covar__residperm_ps,
             DGLMj_ye_covar__locusperm_ps,

             LM_residperm_num_perms,
             LM_locusperm_num_perms,

             Caom_residperm_num_perms,
             Caom_locusperm_num_perms,

             DGLMm_no_covar_residperm_num_perms,
             DGLMm_no_covar_locusperm_num_perms,

             DGLMm_ye_covar_residperm_num_perms,
             DGLMm_ye_covar_locusperm_num_perms,

             Lev_residperm_num_perms,
             Lev_locusperm_num_perms,

             Caov_residperm_num_perms,
             Caov_locusperm_num_perms,

             DGLMv_no_covar_residperm_num_perms,
             DGLMv_no_covar_locusperm_num_perms,

             DGLMv_ye_covar_residperm_num_perms,
             DGLMv_ye_covar_locusperm_num_perms,

             Caoj_residperm_num_perms,
             Caoj_locusperm_num_perms,

             DGLMj_no_covar_residperm_num_perms,
             DGLMj_no_covar_locusperm_num_perms,

             DGLMj_ye_covar_residperm_num_perms,
             DGLMj_ye_covar_locusperm_num_perms)


}
