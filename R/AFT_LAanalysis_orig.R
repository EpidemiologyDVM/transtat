library(transtat)

# Data frame for LA household data
LAdata <- read.csv('LAdata_2019-12.csv', header = TRUE)

# Data summary
dim(LAdata)
sum(LAdata$index)
by(LAdata, LAdata$hhid, function(h) sum(h$index))
sum(LAdata$sickyn)


# households missing data
sum(is.na(LAdata$male))
sum(is.na(LAdata$age))
sum(is.na(LAdata$proph))
sum(with(LAdata, is.na(male) | is.na(age) | is.na(proph)))
hh_NA <- with(LAdata, unique(hhid[is.na(male) | is.na(age) | is.na(proph)]))
with(subset(LAdata, hhid %in% hh_NA), sum(index))
with(subset(LAdata, hhid %in% hh_NA), sum(sickyn))

# Generate pairwise data 
pairdata <- function(incpd=2, latpd=0, infpd=6) {
  # Sum of incubation and infectious periods must be less than 10 days.

  # infection time
  LAdata$Eday <- with(LAdata, ifelse(!is.na(onset_day), onset_day - incpd, NA))

  # onset of infectiousness
  LAdata$Iday <- with(LAdata, ifelse(!is.na(onset_day), Eday + latpd, NA))

  # end of infectiousness
  LAdata$Rday <- with(LAdata, ifelse(!is.na(onset_day), Iday + infpd, NA))

  # infection time of index case
  LAdata$CTday <- sapply(
    LAdata$hhid, 
    function(h) min(LAdata$onset_day[LAdata$hhid == h], na.rm = TRUE) - incpd
  )

  # age category
  LAdata$adult <- ifelse(LAdata$age < 18, 0, 1)

  # individuals with no prophylaxis
  # CTday + 14 is the last time at which an individual could be infected
  # and have symptom onset within 14 days of the index symptom onset.
  LAext_proph0 <- with(
    subset(LAdata, index == 0 & !is.na(proph) & proph == 0),
    data.frame(
      id_sus = id,
      hhid = hhid,
      start = CTday,
      stop = ifelse(sickyn, onset_day - incpd, CTday + 14),
      infset = sickyn,
      external = 1,
      male_sus = male,
      adult_sus = adult,
      proph_sus = 0
    )
  )
  
  # person-time prior to (possible) prophylaxis
  LAext_proph1a <- with(
    subset(LAdata, index == 0 & (is.na(proph) | proph == 1)),
    data.frame(
      id_sus = id,
      hhid = hhid,
      start = CTday,
      stop = CTday + incpd,
      infset = 0,
      external = 1,
      male_sus = male,
      adult_sus = adult,
      proph_sus = 0
    )
  )

  # person-time after (possible) prophylaxis
  LAext_proph1b <- with(
    subset(LAdata, index == 0 & (is.na(proph) | proph == 1)),
    data.frame(
      id_sus = id,
      hhid = hhid,
      start = CTday + incpd,
      stop = ifelse(sickyn, onset_day - incpd, CTday + 14),
      infset = 0,
      external = 1,
      male_sus = male,
      adult_sus = adult,
      proph_sus = proph
    )
  )

  # combine external row data and add infectiousness covariates (all zero)
  LAext <- rbind(LAext_proph0, LAext_proph1a, LAext_proph1b)
  LAext <- subset(LAext, start < stop)
  LAext$id_inf <- 0
  LAext$male_inf <- 0
  LAext$adult_inf <- 0
  LAext$proph_inf <- 0

  # Generate pairwise data within households
  hhlist <- split(LAdata, LAdata$hhid)
  hhpairs <- function(hh) {
    infecteds <- with(hh, id[sickyn == 1])
    ipairs <- function(i) {
      # get covariates for i
      hhid <- hh$hhid[hh$id == i]
      male_inf <- hh$male[hh$id == i]
      adult_inf <- hh$adult[hh$id == i]
      proph_inf <- hh$proph[hh$id == i]

      Iday_i <- hh$Iday[hh$id == i]
      riskset <- subset(hh, is.na(Eday) | Eday > Iday_i)
      if (nrow(riskset) > 0) {
        riskset$hhid <- hhid

        # infectiousness covariates
        riskset$id_inf <- i
        riskset$male_inf <- male_inf
        riskset$adult_inf <- adult_inf
        riskset$proph_inf <- proph_inf

        # susceptibility covariates
        riskset$id_sus <- riskset$id
        riskset$male_sus <- riskset$male
        riskset$adult_sus <- riskset$adult
        riskset$proph_sus <- riskset$proph

        # start, stop, and infset indicator
        riskset$start <- 0
        riskset$stop <- with(
          riskset, ifelse(is.na(Eday), infpd, pmin(infpd, Eday - Iday_i))
        )
        riskset$infset <- with(
          riskset, ifelse((is.na(Eday) | Eday - Iday_i > infpd), 0, 1)
        )
        riskset$external <- 0
        return(riskset[, names(LAext)])
      }
    }
    hhpairlist <- lapply(infecteds, ipairs)
    hhpairlist[["make.row.names"]] <- FALSE
    return(do.call(rbind, hhpairlist))
  }

  pairs <- lapply(hhlist, hhpairs)
  pairs[["make.row.names"]]  = FALSE
  LAint <- do.call(rbind, pairs)
  return(rbind(do.call(rbind, pairs), LAext))
}

# Fit full model with complete data and calculate AIC
pairdat <- pairdata()
#hh_NAproph <- with(pairdat, unique(hhid[is.na(proph_inf) | is.na(proph_sus)]))
pairdat_complete <- subset(pairdat, !(hhid %in% hh_NA))

# Fit full model for all nine combinations of dist and xdist
treg_exp <- transreg(
  Surv(start, stop, infset) ~
    male_inf + male_sus + adult_inf + adult_sus + proph_inf + proph_sus
    + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "exponential",
  optim_method = "BFGS"
)
for (dist in c("exponential", "weibull", "loglogistic")) {
  for (xdist in c("exponential", "weibull", "loglogistic")) {
    print(paste("dist =", dist))
    print(paste("xdist =", xdist))
    try({
      treg <- transreg(
        Surv(start, stop, infset) ~
          male_inf + male_sus + adult_inf + adult_sus + proph_inf + proph_sus
          + ext(external),
        sus = "id_sus",
        data = pairdat_complete,
        dist = dist,
        xdist = xdist,
        optim_method = "BFGS",
        init = coef(treg_exp)
      )
    print(AIC(treg))
    })
  }
}

# Best combination is exponential internal with log-logistic external
# Build model using backwards elimination to find best AIC
# full model; AIC 207.66
treg1 <- transreg(
  Surv(start, stop, infset) ~
    male_inf + male_sus + adult_inf + adult_sus + proph_inf + proph_sus
    + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg_exp)
)
AIC(treg1)
summary(treg1)

# remove proph_inf (p = 0.904); get AIC 206.13
treg2 <- transreg(
  Surv(start, stop, infset) ~
    male_inf + male_sus + adult_inf + adult_sus + proph_sus
    + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg1)[-match("proph_inf", names(coef(treg1)))]
)
AIC(treg2)
summary(treg2)

# remove male_inf (p = 0.544); get AIC 204.75
treg3 <- transreg(
  Surv(start, stop, infset) ~
    male_sus + adult_inf + adult_sus + proph_sus
    + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg2)[-match("male_inf", names(coef(treg2)))]
)
AIC(treg3)
summary(treg3)

# remove male_sus (p = 0.381); get AIC 203.59
treg4 <- transreg(
  Surv(start, stop, infset) ~
    adult_inf + adult_sus + proph_sus
    + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg3)[-match("male_sus", names(coef(treg3)))]
)
AIC(treg4)
summary(treg4)

# remove adult_inf (p = 0.244); get AIC 203.87
treg5 <- transreg(
  Surv(start, stop, infset) ~
    adult_sus + proph_sus
    + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)[-match("adult_inf", names(coef(treg4)))]
)
AIC(treg5)
summary(treg5)

# remove adult_sus (p = 0.153); get AIC 204.28
treg6 <- transreg(
  Surv(start, stop, infset) ~
    proph_sus
    + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg5)[-match("adult_sus", names(coef(treg5)))]
)
AIC(treg6)
summary(treg6)

# null model; AIC 206.98
treg7 <- transreg(
  Surv(start, stop, infset) ~
    ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg6)[-match("proph_sus", names(coef(treg6)))]
)
AIC(treg7)
summary(treg7)

# Final model is treg4 (adult_inf, adult_sus, proph_sus)
# Check interaction for susceptibility covariates
# adult_sus (p = 0.894 and AIC = 205.57)
treg4_intx_adult <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus * ext(external) + proph_sus,
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)
)
AIC(treg4_intx_adult)
summary(treg4_intx_adult)

# proph_sus (p = .868 and AIC = 202.371)
treg4_intx_proph <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + proph_sus * ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)
)
AIC(treg4_intx_proph)
summary(treg4_intx_proph)

# Get p-value for joint effect of proph_sus and proph_sus:external
treg4_null_proph <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + ext(external),
  sus = "id_sus",
  data = pairdat_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)[-match(c("proph_sus"), names(coef(treg4)))]
)
1 - pchisq(2 * (logLik(treg4) - logLik(treg4_null_proph)), 2)

# Summary of final model (Table 4)
summary(treg4, conf.type = "lr")
exp(coef(treg4))
exp(confint(treg4, type = "lr"))

# Household SAR with and without prophylaxis
hSAR <- function(model, covars, infpd=6) {
  lnrate <- as.numeric(covars %*% coef(model))
  sd_lnrate <- sqrt(as.numeric(covars %*% vcov(model) %*% covars))

  sar <- 1 - exp(-exp(lnrate) * infpd)
  sar_ci <- 1 - exp(-exp(lnrate + c(-1, 1) * qnorm(.975) * sd_lnrate) * infpd)

  return(list(sar = sar, sar_ci = sar_ci))
}

sarCC0 <- hSAR(treg4, c(1, 0, 0, 0, 0, 0))
sarCC1 <- hSAR(treg4, c(1, 0, 0, 1, 0, 0))

sarCA0 <- hSAR(treg4, c(1, 0, 1, 0, 0, 0))
sarCA1 <- hSAR(treg4, c(1, 0, 1, 1, 0, 0))

sarAC0 <- hSAR(treg4, c(1, 1, 0, 0, 0, 0))
sarAC1 <- hSAR(treg4, c(1, 1, 0, 1, 0, 0))

sarAA0 <- hSAR(treg4, c(1, 1, 1, 0, 0, 0))
sarAA1 <- hSAR(treg4, c(1, 1, 1, 1, 0, 0))

sar0_exp <- hSAR(treg6_exp, c(1, 0, 0))
sar1_exp <- hSAR(treg6_exp, c(1, 1, 0))

# Fit model without external data
treg4_int <- transreg(
  Surv(start, stop, infset) ~ adult_inf + adult_sus + proph_sus,
  sus = "id_sus",
  data = subset(pairdat_complete, external == 0),
  dist = "exponential",
  optim_method = "BFGS",
  init = coef(treg4)[-c(5, 6)]
)
summary(treg4_int, conf.type = "lr")

treg4_int_weib <- transreg(
  Surv(start, stop, infset) ~ adult_inf + adult_sus + proph_sus,
  sus = "id_sus",
  data = subset(pairdat_complete, external == 0),
  dist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)[-c(5, 6)]
)
summary(treg4_int_weib, conf.type = "lr")

treg4_int_llog <- transreg(
  Surv(start, stop, infset) ~ adult_inf + adult_sus + proph_sus,
  sus = "id_sus",
  data = subset(pairdat_complete, external == 0),
  dist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)[-c(5, 6)]
)
summary(treg4_int_llog, conf.type = "lr")

# Sensitivity analysis
# Latent period 1 day
lp1 <- pairdata(latpd = 1)
lp1_complete <- subset(lp1, !(hhid %in% hh_NA))

treg4_lp1 <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + proph_sus + ext(external),
  sus = "id_sus",
  data = lp1_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)
)
summary(treg4_lp1, conf.type = "lr")
exp(coef(treg4_lp1))

treg4_lp1_int <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + proph_sus,
  sus = "id_sus",
  data = subset(lp1_complete, external == 0),
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)[1:4]
)
summary(treg4_lp1_int, conf.type = "lr")


# Infectious period 5 days
ip5 <- pairdata(infpd = 5)
ip5_complete <- subset(ip5, !(hhid %in% hh_NA))

treg4_ip5 <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + proph_sus + ext(external),
  sus = "id_sus",
  data = ip5_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)
)
summary(treg4_ip5, conf.type = "lr")
exp(coef(treg4_ip5))

treg4_ip5_int <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + proph_sus,
  sus = "id_sus",
  data = subset(ip5_complete, external == 0),
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)[1:4]
)
summary(treg4_ip5_int, conf.type = "lr")

# Infectious period 7 days
ip7 <- pairdata(infpd = 7)
ip7_complete <- subset(ip7, !(hhid %in% hh_NA))

treg4_ip7 <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + proph_sus + ext(external),
  sus = "id_sus",
  data = ip7_complete,
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)
)
summary(treg4_ip7, conf.type = "lr")
exp(coef(treg4_ip7))

treg4_ip7_int <- transreg(
  Surv(start, stop, infset) ~ 
    adult_inf + adult_sus + proph_sus,
  sus = "id_sus",
  data = subset(ip7_complete, external == 0),
  dist = "exponential",
  xdist = "loglogistic",
  optim_method = "BFGS",
  init = coef(treg4)[1:4]
)
summary(treg4_ip7_int, conf.type = "lr")



