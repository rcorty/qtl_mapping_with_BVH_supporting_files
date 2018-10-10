mvqtl_simulation <- function(job_idx) {

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

    # read simulated data and subset
    {

      message('Working directory is...')
      print(getwd())

      first_sim_idx <- (job_idx - 1)*num_sims_per_job + 1
      last_sim_idx <- job_idx*num_sims_per_job
      sim_data <- read_rds(path = paste0(ifelse(test = Sys.info()[1] == 'Darwin',
                                                yes = '',
                                                no = '../'),
                                         'simulation_studies/',
                                         '2_data/',
                                         'simulated_data_nsims=', num_sims, '.gz')) %>%
        filter(between(x = sim_idx, left = first_sim_idx, right = last_sim_idx)) %>%
        mutate(locus = factor(locus),
               covar = factor(covar))

      locus_df <- 2
      num_obs <- 300
    }

    # define containers for results
    {
      result_skeleton <- sim_data %>%
        select(sim_idx, bvh, qtl) %>%
        distinct()

      LM_asymp_ps    <- LM_rint_ps    <- LM_residperm_ps    <- LM_locusperm_ps    <- rep(NA, 8*num_sims_per_job)
      Caom_asymp_ps  <- Caom_rint_ps  <- Caom_residperm_ps  <- Caom_locusperm_ps <-  rep(NA, 8*num_sims_per_job)
      DGLMm_asymp_ps <- DGLMm_rint_ps <- DGLMm_residperm_ps <- DGLMm_locusperm_ps <- rep(NA, 8*num_sims_per_job)

      LM_residperm_num_perms    <- LM_locusperm_num_perms    <- rep(NA, 8*num_sims_per_job)
      Caom_residperm_num_perms  <- Caom_locusperm_num_perms <-  rep(NA, 8*num_sims_per_job)
      DGLMm_residperm_num_perms <- DGLMm_locusperm_num_perms <- rep(NA, 8*num_sims_per_job)

      Lev_asymp_ps   <- Lev_rint_ps   <- Lev_residperm_ps   <- Lev_locusperm_ps   <- rep(NA, 8*num_sims_per_job)
      Caov_asymp_ps  <- Caov_rint_ps  <- Caov_residperm_ps  <- Caov_locusperm_ps <-  rep(NA, 8*num_sims_per_job)
      DGLMv_asymp_ps <- DGLMv_rint_ps <- DGLMv_residperm_ps <- DGLMv_locusperm_ps <- rep(NA, 8*num_sims_per_job)

      Lev_residperm_num_perms   <- Lev_locusperm_num_perms   <- rep(NA, 8*num_sims_per_job)
      Caov_residperm_num_perms  <- Caov_locusperm_num_perms <-  rep(NA, 8*num_sims_per_job)
      DGLMv_residperm_num_perms <- DGLMv_locusperm_num_perms <- rep(NA, 8*num_sims_per_job)

      Caoj_asymp_ps  <- Caoj_rint_ps  <- Caoj_residperm_ps  <- Caoj_locusperm_ps  <- rep(NA, 8*num_sims_per_job)
      DGLMj_asymp_ps <- DGLMj_rint_ps <- DGLMj_residperm_ps <- DGLMj_locusperm_ps <- rep(NA, 8*num_sims_per_job)

      Caoj_residperm_num_perms  <- Caoj_locusperm_num_perms  <- rep(NA, 8*num_sims_per_job)
      DGLMj_residperm_num_perms <- DGLMj_locusperm_num_perms <- rep(NA, 8*num_sims_per_job)
    }
  }

  set.seed(seed = job_idx)

  for (sim_num in 1:num_sims_per_job) {

    message('---- ---- ---- Starting sim ', sim_num, ' of ', num_sims_per_job)
    sim_num_offset <- 8*(sim_num - 1)

    for (this_bvh in c(FALSE, TRUE)) {

      message('---- ---- Starting BVH=', this_bvh)
      bvh_offset <- 4*this_bvh

      for (this_qtl in c('none', 'mqtl', 'vqtl', 'mvqtl')) {

        # browser()

        message('---- Starting QTL=', this_qtl)
        qtl_offset <- which(x = this_qtl == c('none', 'mqtl', 'vqtl', 'mvqtl'))

        result_idx <- sim_num_offset + bvh_offset + qtl_offset

        message('result index:', result_idx)

        this_sim_idx <- sim_num + first_sim_idx - 1
        this_data <- sim_data %>%
          filter(sim_idx == this_sim_idx,
                 bvh == this_bvh,
                 qtl == this_qtl)

        # browser()

        attach(this_data)

        # mQTL tests
        {
          # compute test stats
          LM_LRS <- tryNA(lm_test(resp = phen, locus = locus))
          LM_rint_LRS <- tryNA(lm_test(resp = rint(phen), locus = locus))

          Caom_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = phen, test = "mean"))
          Caom_rint_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = rint(phen), test = "mean"))

          DGLMm_LRS <- tryNA(mean_test(resp = phen, mean_locus = locus, var_locus = locus, covar = covar))
          DGLMm_rint_LRS <- tryNA(mean_test(resp = rint(phen), mean_locus = locus, var_locus = locus, covar = covar))


          # perms for empirical evaluation
          if (num_perms > 0) {

            message('starting mean perms')

            lm_resids <- residuals(lm(phen ~ 1))
            lm_hats <- phen - lm_resids

            cao_resids <- residuals(dglm(formula = phen ~ 1, dformula = ~ locus))
            cao_hats <- phen - cao_resids

            dglm_resids <- residuals(dglm(formula = phen ~ 1, dformula = ~ locus + covar))
            dglm_hats <- phen - dglm_resids

            LM_LRSs_residperm <- Caom_LRSs_residperm <- DGLMm_LRSs_residperm <- rep(NA, num_perms)
            LM_LRSs_locusperm <- Caom_LRSs_locusperm <- DGLMm_LRSs_locusperm <- rep(NA, num_perms)

            for (perm_idx in 1:num_perms) {

              lm_PB <- lm_hats + sample(lm_resids)
              cao_PB <- cao_hats + sample(cao_resids)
              dglm_PB <- dglm_hats + sample(dglm_resids)

              LM_LRSs_residperm[perm_idx] <- tryNA(lm_test(resp = lm_PB, locus = locus))
              Caom_LRSs_residperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = cao_PB, test = "mean"))
              DGLMm_LRSs_residperm[perm_idx] <- tryNA(mean_test(resp = dglm_PB, mean_locus = locus, var_locus = locus, covar = covar))

              LM_LRSs_locusperm[perm_idx] <- tryNA(lm_test(resp = phen, locus = sample(locus)))
              Caom_LRSs_locusperm[perm_idx] <- tryNA(mean_test(resp = phen, mean_locus = sample(locus), var_locus = locus, covar = NULL))
              DGLMm_LRSs_locusperm[perm_idx] <- tryNA(mean_test(resp = phen, mean_locus = sample(locus), var_locus = locus, covar = covar))
            }
          }

          # evaluate test stats
          LM_asymp_ps[result_idx] <- pchisq(q = LM_LRS, df = locus_df, lower.tail = FALSE)
          LM_rint_ps[result_idx] <- pchisq(q = LM_rint_LRS, df = locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            LM_residperm_ps[result_idx] <- mean(x = LM_LRSs_residperm >= LM_LRS, na.rm = TRUE)
            LM_locusperm_ps[result_idx] <- mean(x = LM_LRSs_locusperm >= LM_LRS, na.rm = TRUE)
            LM_residperm_num_perms[result_idx] <- sum(!is.na(LM_LRSs_residperm))
            LM_locusperm_num_perms[result_idx] <- sum(!is.na(LM_LRSs_locusperm))
          }

          Caom_asymp_ps[result_idx] <- pchisq(q = Caom_LRS, df = locus_df, lower.tail = FALSE)
          Caom_rint_ps[result_idx] <- pchisq(q = Caom_rint_LRS, df = locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            Caom_residperm_ps[result_idx] <- mean(x = Caom_LRSs_residperm >= Caom_LRS, na.rm = TRUE)
            Caom_locusperm_ps[result_idx] <- mean(x = Caom_LRSs_locusperm >= Caom_LRS, na.rm = TRUE)
            Caom_residperm_num_perms[result_idx] <- sum(!is.na(Caom_LRSs_residperm))
            Caom_locusperm_num_perms[result_idx] <- sum(!is.na(Caom_LRSs_locusperm))
          }

          DGLMm_asymp_ps[result_idx] <- pchisq(q = DGLMm_LRS, df = locus_df, lower.tail = FALSE)
          DGLMm_rint_ps[result_idx] <- pchisq(q = DGLMm_rint_LRS, df = locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            DGLMm_residperm_ps[result_idx] <- mean(x = DGLMm_LRSs_residperm >= DGLMm_LRS, na.rm = TRUE)
            DGLMm_locusperm_ps[result_idx] <- mean(x = DGLMm_LRSs_locusperm >= DGLMm_LRS, na.rm = TRUE)
            DGLMm_residperm_num_perms[result_idx] <- sum(!is.na(DGLMm_LRSs_residperm))
            DGLMm_locusperm_num_perms[result_idx] <- sum(!is.na(DGLMm_LRSs_locusperm))
          }
          }

        # vQTL tests
        {

          # compute test stats
          Lev_stat <-  tryNA(leveneTest(y = phen, group = locus)$`F value`[1])
          Lev_rint_stat <- tryNA(leveneTest(y = rint(phen), group = locus)$`F value`[1])

          Caov_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = phen, test = "var"))
          Caov_rint_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = rint(phen), test = "var"))

          DGLMv_LRS <- tryNA(var_test(resp = phen, mean_locus = locus, var_locus = locus, covar = covar))
          DGLMv_rint_LRS <- tryNA(var_test(resp = rint(phen), mean_locus = locus, var_locus = locus, covar = covar))

          # perms for empirical evaluation
          if (num_perms > 0) {

            message('starting var perms')

            cao_resids <- residuals(lm(phen ~ locus))
            cao_hats <- phen - cao_resids

            dglm_resids <- residuals(dglm(formula = phen ~ locus, dformula = ~ covar))
            dglm_hats <- phen - dglm_resids

            Lev_stats_residperm <- Caov_LRSs_residperm <- DGLMv_LRSs_residperm <- rep(NA, num_perms)
            Lev_stats_locusperm <- Caov_LRSs_locusperm <- DGLMv_LRSs_locusperm <- rep(NA, num_perms)
            for (perm_idx in 1:num_perms) {

              lm_PB <- lm_hats + sample(lm_resids)
              cao_PB <- cao_hats + sample(cao_resids)
              dglm_PB <- dglm_hats + sample(dglm_resids)

              Lev_stats_residperm[perm_idx] <- tryNA(leveneTest(y = lm_PB, group = locus)$`F value`[1])
              Caov_LRSs_residperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = cao_PB, test = "var"))
              DGLMv_LRSs_residperm[perm_idx] <- tryNA(var_test(resp = dglm_PB, mean_locus = locus, var_locus = locus, covar = covar))

              Lev_stats_locusperm[perm_idx] <- tryNA(leveneTest(y = phen, group = sample(locus))$`F value`[1])
              Caov_LRSs_locusperm[perm_idx] <- tryNA(var_test(resp = phen, mean_locus = locus, var_locus = sample(locus), covar = NULL))
              DGLMv_LRSs_locusperm[perm_idx] <- tryNA(var_test(resp = phen, mean_locus = locus, var_locus = sample(locus), covar = covar))

            }
          }

          # evaluate test stats
          Lev_asymp_ps[result_idx] <- pf(q = Lev_stat, df1 = locus_df, df2 = num_obs - locus_df, lower.tail = FALSE)
          Lev_rint_ps[result_idx] <- pf(q = Lev_rint_stat, df1 = locus_df, df2 = num_obs - locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            Lev_residperm_ps[result_idx] <- mean(x = Lev_stats_residperm >= Lev_stat, na.rm = TRUE)
            Lev_locusperm_ps[result_idx] <- mean(x = Lev_stats_locusperm >= Lev_stat, na.rm = TRUE)
            Lev_residperm_num_perms[result_idx] <- sum(!is.na(Lev_stats_residperm))
            Lev_locusperm_num_perms[result_idx] <- sum(!is.na(Lev_stats_locusperm))
          }

          Caov_asymp_ps[result_idx] <- pchisq(q = Caov_LRS, df = locus_df, lower.tail = FALSE)
          Caov_rint_ps[result_idx] <- pchisq(q = Caov_rint_LRS, df = locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            Caov_residperm_ps[result_idx] <- mean(x = Caov_LRSs_residperm >= Caov_LRS, na.rm = TRUE)
            Caov_locusperm_ps[result_idx] <- mean(x = Caov_LRSs_locusperm >= Caov_LRS, na.rm = TRUE)
            Caov_residperm_num_perms[result_idx] <- sum(!is.na(Caov_LRSs_residperm))
            Caov_locusperm_num_perms[result_idx] <- sum(!is.na(Caov_LRSs_locusperm))
          }

          DGLMv_asymp_ps[result_idx] <- pchisq(q = DGLMv_LRS, df = locus_df, lower.tail = FALSE)
          DGLMv_rint_ps[result_idx] <- pchisq(q = DGLMv_rint_LRS, df = locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            DGLMv_residperm_ps[result_idx] <- mean(x = DGLMv_LRSs_residperm >= DGLMv_LRS, na.rm = TRUE)
            DGLMv_locusperm_ps[result_idx] <- mean(x = DGLMv_LRSs_locusperm >= DGLMv_LRS, na.rm = TRUE)
            DGLMv_residperm_num_perms[result_idx] <- sum(!is.na(DGLMv_LRSs_residperm))
            DGLMv_locusperm_num_perms[result_idx] <- sum(!is.na(DGLMv_LRSs_locusperm))
          }

        }

        # joint tests
        {
          # compute test stats
          Caoj_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = phen, test = "both"))
          Caoj_rint_LRS <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = rint(phen), test = "both"))

          DGLMj_LRS <- tryNA(joint_test(resp = phen, locus = locus, covar = covar))
          DGLMj_rint_LRS <- tryNA(joint_test(resp = rint(phen), locus = locus, covar = covar))

          # perms for empirical evaluation
          if (num_perms > 0) {

            message('starting joint perms')

            cao_resids <- residuals(lm(phen ~ 1))
            cao_hats <- phen - cao_resids

            dglm_resids <- residuals(dglm(formula = phen ~ 1, dformula = ~ covar))
            dglm_hats <- phen - dglm_resids

            Caoj_LRSs_residperm <- DGLMj_LRSs_residperm <- rep(NA, num_perms)
            Caoj_LRSs_locusperm <- DGLMj_LRSs_locusperm <- rep(NA, num_perms)
            for (perm_idx in 1:num_perms) {

              cao_PB <- cao_hats + sample(cao_resids)
              dglm_PB <- dglm_hats + sample(dglm_resids)

              Caoj_LRSs_residperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(locus) - 1, depvar = cao_PB, test = "both"))
              DGLMj_LRSs_residperm[perm_idx] <- tryNA(joint_test(resp = dglm_PB, locus = locus, covar = covar))

              Caoj_LRSs_locusperm[perm_idx] <- tryNA(LRT_mean_var(gt = as.numeric(sample(locus)) - 1, depvar = phen, test = "both"))
              DGLMj_LRSs_locusperm[perm_idx] <- tryNA(joint_test(resp = phen, locus = sample(locus), covar = covar))
            }
          }

          # evaluate test stats
          Caoj_asymp_ps[result_idx] <- pchisq(q = Caoj_LRS, df = 2*locus_df, lower.tail = FALSE)
          Caoj_rint_ps[result_idx] <- pchisq(q = Caoj_rint_LRS, df = 2*locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            Caoj_residperm_ps[result_idx] <- mean(x = Caoj_LRSs_residperm >= Caoj_LRS, na.rm = TRUE)
            Caoj_locusperm_ps[result_idx] <- mean(x = Caoj_LRSs_locusperm >= Caoj_LRS, na.rm = TRUE)
            Caoj_residperm_num_perms[result_idx] <- sum(!is.na(Caoj_LRSs_residperm))
            Caoj_locusperm_num_perms[result_idx] <- sum(!is.na(Caoj_LRSs_locusperm))
          }

          DGLMj_asymp_ps[result_idx] <- pchisq(q = DGLMj_LRS, df = 2*locus_df, lower.tail = FALSE)
          DGLMj_rint_ps[result_idx] <- pchisq(q = DGLMj_rint_LRS, df = 2*locus_df, lower.tail = FALSE)
          if (num_perms > 0) {
            DGLMj_residperm_ps[result_idx] <- mean(x = DGLMj_LRSs_residperm >= DGLMj_LRS, na.rm = TRUE)
            DGLMj_locusperm_ps[result_idx] <- mean(x = DGLMj_LRSs_locusperm >= DGLMj_LRS, na.rm = TRUE)
            DGLMj_residperm_num_perms[result_idx] <- sum(!is.na(DGLMj_LRSs_residperm))
            DGLMj_locusperm_num_perms[result_idx] <- sum(!is.na(DGLMj_LRSs_locusperm))
          }
        }

        detach(name = 'this_data')
      }
    }
  }

  message('Done w all sims.')

  result_skeleton %>%
    mutate(job_idx = job_idx,

           LM_asymp_ps = LM_asymp_ps,
           LM_rint_ps  = LM_rint_ps,
           LM_residperm_ps = LM_residperm_ps,
           LM_locusperm_ps = LM_locusperm_ps,

           Caom_asymp_ps = Caom_asymp_ps,
           Caom_rint_ps = Caom_rint_ps,
           Caom_residperm_ps = Caom_residperm_ps,
           Caom_locusperm_ps = Caom_locusperm_ps,

           DGLMm_asymp_ps = DGLMm_asymp_ps,
           DGLMm_rint_ps = DGLMm_rint_ps,
           DGLMm_residperm_ps = DGLMm_residperm_ps,
           DGLMm_locusperm_ps = DGLMm_locusperm_ps,

           Lev_asymp_ps = Lev_asymp_ps,
           Lev_rint_ps = Lev_rint_ps,
           Lev_residperm_ps = Lev_residperm_ps,
           Lev_locusperm_ps = Lev_locusperm_ps,

           Caov_asymp_ps = Caov_asymp_ps,
           Caov_rint_ps = Caov_rint_ps,
           Caov_residperm_ps = Caov_residperm_ps,
           Caov_locusperm_ps = Caov_locusperm_ps,

           DGLMv_asymp_ps = DGLMv_asymp_ps,
           DGLMv_rint_ps = DGLMv_rint_ps,
           DGLMv_residperm_ps = DGLMv_residperm_ps,
           DGLMv_locusperm_ps = DGLMv_locusperm_ps,

           Caoj_asymp_ps = Caoj_asymp_ps,
           Caoj_rint_ps = Caoj_rint_ps,
           Caoj_residperm_ps = Caoj_residperm_ps,
           Caoj_locusperm_ps = Caoj_locusperm_ps,

           DGLMj_asymp_ps = DGLMj_asymp_ps,
           DGLMj_rint_ps = DGLMj_rint_ps,
           DGLMj_residperm_ps = DGLMj_residperm_ps,
           DGLMj_locusperm_ps = DGLMj_locusperm_ps,

           LM_residperm_num_perms = LM_residperm_num_perms,
           LM_locusperm_num_perms = LM_locusperm_num_perms,

           Caom_residperm_num_perms = Caom_residperm_num_perms,
           Caom_locusperm_num_perms = Caom_locusperm_num_perms,

           DGLMm_residperm_num_perms = DGLMm_residperm_num_perms,
           DGLMm_locusperm_num_perms = DGLMm_locusperm_num_perms,

           Lev_residperm_num_perms = Lev_residperm_num_perms,
           Lev_locusperm_num_perms = Lev_locusperm_num_perms,

           Caov_residperm_num_perms = Caov_residperm_num_perms,
           Caov_locusperm_num_perms = Caov_locusperm_num_perms,

           DGLMv_residperm_num_perms = DGLMv_residperm_num_perms,
           DGLMv_locusperm_num_perms = DGLMv_locusperm_num_perms,

           Caoj_residperm_num_perms = Caoj_residperm_num_perms,
           Caoj_locusperm_num_perms = Caoj_locusperm_num_perms,

           DGLMj_residperm_num_perms = DGLMj_residperm_num_perms,
           DGLMj_locusperm_num_perms = DGLMj_locusperm_num_perms) %>%
    select(job_idx, sim_idx, bvh, qtl, LM_asymp_ps:DGLMj_locusperm_num_perms)

}

