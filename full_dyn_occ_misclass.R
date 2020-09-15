# Here we create a Dynamic occupancy model with misclassification using
# disaggregated individual-level data.This model estimates misclassification
# probabilities with NEON Niwot Ridge carabid data 2015-2018

library(dclone) #jags.parfit()
library(MCMCvis)
library(dplyr)
library(patchwork)
library(ggmcmc)
library(readr)
library(ggplot2)
options(digits = 3)

# Load jags input
source("source/jags_input.R")
jagsinput <- return_jags_input("full")
list2env(jagsinput, .GlobalEnv)
rm(jagsinput)

# Plot discrepancies for samples identified to the species level
para_new %>%
  filter(!is.na(exp_sciname)) %>% #taxonRank == "species",
  count(para_morph_combo, exp_sciname) %>%
  arrange(n) %>%
  mutate(discrepancy = para_morph_combo != exp_sciname) %>%
  ggplot(aes(x = para_morph_combo,
             y = exp_sciname,
             color = discrepancy)) +
  geom_point(aes(size = n)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab("Species ID from parataxonomist") +
  ylab("Species ID from expert taxonomist") +
  scale_color_manual(values = c("black", "red"))

# JAGS model --------------------------------------------------------------
# Run model in JAGS. 
jags_d <- list(nsite = dim(L)[1],
               nsurv = dim(L)[2], 
               nyear = dim(L)[3], 
               K_exp = dim(alpha)[1], 
               K_para = dim(alpha)[2],
               L = L, 
               alpha = alpha,
               Ltot = sum(L), 
               site = y_df$plotID_idx, #needs to be numeric
               year = as.numeric(as.factor(y_df$col_year)), #needs to be numeric starting from 1
               # if the individual was labeled by the expert, true ID is known
               k = y_df$expertID_idx,
               # for all individuals, we get paratxonomist IDs
               y = y_df$parataxID_idx,
               z = z.dat,
               R = diag(rep(1, 4)))
JAGSinits <- function(){ list(z = z.init) }
nc <- 4
ni <- 20000
cl <- makeCluster(nc)
jm <- jags.parfit(cl = cl,
                  data = jags_d,
                  params = c("logit_psi1", "logit_phi", "logit_gamma",
                             "psi",
                             "lambda", "Theta",
                             "eps_site", "eps_spec", "Tau_spec", "Tau_site",
                             "log_growth", "turnover"),
                  model = "full_dyn_occ_misclass_JAGS.txt",
                  n.chains = nc,
                  n.adapt = 2000,
                  n.update = 2000,
                  thin = ni/1000,
                  n.iter = ni) 

dir.create("output", showWarnings = FALSE)
jm_summ <- MCMCsummary(jm)
saveRDS(jm, "output/full_jm.rds")
saveRDS(jm_summ, "output/full_jmsumm.rds")

# View JAGS jmput --------------------------------------------------------

jm <- readRDS("output/full_jm.rds")
jm_summ <- readRDS("output/full_jmsumm.rds")

# Did model converge?
hist(jm_summ$Rhat, breaks=40)
range(jm_summ$Rhat, na.rm=TRUE)

# Look at high Rhat values
jm_summ <- rownames_to_column(jm_summ)
jm_summ %>% filter(Rhat > 1.1)

MCMCsummary(jm, params = paste0('Theta\\[48,66\\]'), ISB=F, Rhat=T) 
MCMCsummary(jm, params = paste0('Theta\\[48,48\\]'), ISB=F, Rhat=T) 

MCMCsummary(jm, params = paste0('Theta\\[46,67\\]'), ISB=F, Rhat=T)
MCMCsummary(jm, params = paste0('Theta\\[46,46\\]'), ISB=F, Rhat=T)

# Look at raw numbers
# theta
theta_summ <- MCMCsummary(jm, params = 'Theta')
saveRDS(theta_summ, "output/full_theta_summ.rds")
theta_summ <- readRDS("output/full_theta_summ.rds")
hist(theta_summ$Rhat)
hist(theta_summ$mean)
range(theta_summ$Rhat, na.rm=TRUE)
theta_summ <- rownames_to_column(theta_summ)

# # Look at some theta posteriors with high Rhat values
MCMCtrace(jm, params = paste0('Theta\\[1,1\\]'), type = 'both', ind = F, pdf=F, ISB=F, Rhat=T) 
theta_prior <- MCMCpack::rdirichlet(1000, c(100,rep(1.5,dim(alpha)[2]-1) ))
plot(density(theta_prior[, 1]))


theta_df <- data.frame(expert_index = rep(1:dim(alpha)[1], dim(alpha)[2]) ,
                       paramorph_index = rep(1:dim(alpha)[2], each = dim(alpha)[1]),
                       expert_sciname = rep(rownames, dim(alpha)[2]),
                       para_morph = rep(colnames, each = dim(alpha)[1]))  %>% 
    mutate(theta_median = theta_summ$"50%")
#make para_morph a factor to force plotting in order
theta_df$para_morph = factor(theta_df$para_morph, levels=colnames) 

# Plot heatmap of theta values
png("figures/thetaconfusion.png", width=850, height=600)
print(ggplot(theta_df, aes(x=para_morph, y=expert_sciname, fill= theta_median)) + 
  geom_tile() +
  scale_fill_viridis_c("Posterior\nmedian") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Parataxonomist ID") + ylab("Expert Taxonomist ID") +
  scale_y_discrete(limits = rev(levels(as.factor(theta_df$expert_sciname)))))
dev.off()

# psi - occupancy prob.
(psi_summary <- MCMCsummary(jm, params = 'psi', round=2))
MCMCtrace(jm, params = 'psi', type = 'trace', ind = F, pdf=F)

# Plot all species' occupancy through seasons
png("figures/occthroughtime.png", width=850, height=600)
print(rownames_to_column(psi_summary) %>%
  as_tibble %>%
  separate("rowname", into = c("site", "species", "year"), sep = ",") %>%
  mutate_at(c('site', 'species', 'year'), readr::parse_number) %>%
  mutate(site = dimnames(z.init)[[1]][site], 
         species = dimnames(z.init)[[2]][species], 
         year = dimnames(z.init)[[3]][year]) %>%
  ggplot(aes(year, mean, group = site)) + 
  geom_line() + 
  facet_wrap(~species) + 
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), color = NA, alpha = .05) +
  xlab("Year") + 
  ylab("Occupancy probability") + 
  theme_minimal())
dev.off()


# Visualize relationships among site and species level parameters ---------

get_rho_df <- function(level) {
  stopifnot(level %in% c("site", "spec"))
  jm %>%
    lapply(function(x) {
      cols <- paste0("Tau_", level) %>%
        grepl(colnames(x))
      as.data.frame(x[, cols]) %>%
        mutate(iter = 1:n()) %>%
        tidyr::pivot_longer(cols = -iter) %>%
        separate(name, into = c("row", "col"), sep = ",") %>%
        mutate_at(c("row", "col"), readr::parse_number)
    }) %>%
    bind_rows(.id="chain") %>%
    unite("iter", c("chain", "iter"), sep = "_") %>%
    split(.$iter) %>%
    lapply(FUN = function(d) {
      Tau <- reshape2::acast(data = select(d, -iter), row~col, value.var="value")
      Sigma <- solve(Tau)
      Rho <- cov2cor(Sigma)
      reshape2::melt(Rho, varnames = c("row", "col"), value.name = "rho") %>%
        as_tibble %>%
        filter(row > col)
    }) %>%
    bind_rows(.id = "iter")
}

rho_site <- get_rho_df(level = "site")
rho_spec <- get_rho_df(level = "spec")

rho_site %>%
  ggplot(aes(rho)) + 
  geom_density() + 
  facet_grid(row~col)

rho_spec %>%
  ggplot(aes(rho)) + 
  geom_density() + 
  facet_grid(row~col)

eps_d <- ggmcmc::ggs(jm, family = "eps_") %>%
  mutate(level = ifelse(grepl("site", Parameter), "site", "spec")) %>%
  separate(Parameter, into = c("row", "col"), sep = ",") %>%
  mutate_at(c('row', 'col'), parse_number) %>%
  unite("iter", c("Iteration" , "Chain")) %>%
  group_by(row, col, level) %>%
  summarize(med = median(value)) %>%
  ungroup %>%
  mutate(param = c("e_psi1", "e_phi", "e_gamma", "e_lambda")[col])

diag_hists <- function(lev, column, xmin=-5, xmax=5, binwidth=.25, 
                       xlabel) {
  stopifnot(lev %in% eps_d$level)
  df <- eps_d %>%
    filter(level == lev, 
           col == column)
  df %>%
    ggplot(aes(x=med)) + 
      geom_histogram(binwidth = binwidth, fill = "#1F968BFF") + 
      xlab(xlabel) + 
      ylab("") + 
      xlim(xmin, xmax) + 
      theme(axis.text.y = element_blank())
}

diag_hists(lev = "site", column = 1, 
           xlabel = as.expression(bquote(epsilon^(psi[1]))))
diag_hists(lev = "spec", column = 1, 
           xlabel = as.expression(bquote(epsilon^(psi[1]))))

lower_diag <- function(d1, d2, lev, xmin=-5, xmax=5) {
  df <- eps_d %>%
    select(-col) %>%
    filter(level == lev) %>%
    pivot_wider(names_from = param, values_from = med)
  df$x <- df[[d1]]
  df$y <- df[[d2]]
  
  df %>%
    ggplot(aes(x, y)) + 
    geom_point(color = "#482677FF") + 
    xlab("") + 
    ylab("") + 
    xlim(xmin, xmax) + 
    ylim(xmin, xmax)
}

lower_diag("e_psi1", "e_lambda", lev = "site")
lower_diag("e_lambda", "e_psi1", lev = "site")


upper_diag <- function(d, r, c, xlabel) {
  stopifnot(r > c)
  d %>%
    filter(row == r, col == c) %>%
    ggplot(aes(rho)) + 
      geom_density(fill = "#482677FF", alpha = .7) + 
      xlab(xlabel) + 
      ylab("") + 
      xlim(-1, 1) + 
      theme(axis.text.y = element_blank())
}

upper_diag(rho_spec, r = 4, c = 3, 
           xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi)))))
upper_diag(rho_site, r = 4, c = 3, 
           xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi)))))

# Cobble together all of the panels (site level)
site_xmin <- -1.7
site_xmax <- 1.2

wrap_plots(
  list(
    # Column 1
    diag_hists(lev = "site", column = 1, 
               xlabel = as.expression(bquote(epsilon^(psi[1]))), 
               xmin = site_xmin, xmax = site_xmax), 
    lower_diag("e_psi1", "e_phi", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_psi1", "e_gamma", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_psi1", "e_lambda", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    # Column 2
    upper_diag(rho_site, r = 2, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(phi), epsilon^(psi[1]))))),
    diag_hists(lev = "site", column = 2, 
               xlabel = as.expression(bquote(epsilon^(phi))), 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_phi", "e_gamma", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_phi", "e_lambda", lev = "site", 
               xmin = site_xmin, xmax = site_xmax),
    # Column 3
    upper_diag(rho_site, r = 3, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(psi[1]))))),
    upper_diag(rho_site, r = 3, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(phi))))),
    diag_hists(lev = "site", column = 3, 
               xlabel = as.expression(bquote(epsilon^(gamma))), 
               xmin = site_xmin, xmax = site_xmax),
    lower_diag("e_gamma", "e_lambda", lev = "site", 
               xmin = site_xmin, xmax = site_xmax), 
    # Column 4
    upper_diag(rho_site, r = 4, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(psi[1]))))),
    upper_diag(rho_site, r = 4, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi))))),
    upper_diag(rho_site, r = 4, c = 3, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(gamma))))),
    diag_hists(lev = "site", column = 4, 
               xlabel = as.expression(bquote(epsilon^(lambda))), 
               xmin = site_xmin, xmax = site_xmax)
    ), 
  nrow = 4, byrow = FALSE)
#dir.create("fig", showWarnings = FALSE)
#ggsave("fig/site-level-ranefs.pdf", width = 7, height = 5)
ggsave("figures/site-level-ranefs.png", width = 7, height = 5)

# Cobble together all of the panels (species level)
spec_xmin <- -4.7
spec_xmax <- 1.8

wrap_plots(
  list(
    # Column 1
    diag_hists(lev = "spec", column = 1, 
               xlabel = as.expression(bquote(epsilon^(psi[1]))), 
               xmin = spec_xmin, xmax = spec_xmax), 
    lower_diag("e_psi1", "e_phi", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_psi1", "e_gamma", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_psi1", "e_lambda", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    # Column 2
    upper_diag(rho_spec, r = 2, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(phi), epsilon^(psi[1]))))),
    diag_hists(lev = "spec", column = 2, 
               xlabel = as.expression(bquote(epsilon^(phi))), 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_phi", "e_gamma", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_phi", "e_lambda", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax),
    # Column 3
    upper_diag(rho_spec, r = 3, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(psi[1]))))),
    upper_diag(rho_spec, r = 3, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(gamma), epsilon^(phi))))),
    diag_hists(lev = "spec", column = 3, 
               xlabel = as.expression(bquote(epsilon^(gamma))), 
               xmin = spec_xmin, xmax = spec_xmax),
    lower_diag("e_gamma", "e_lambda", lev = "spec", 
               xmin = spec_xmin, xmax = spec_xmax), 
    # Column 4
    upper_diag(rho_spec, r = 4, c = 1, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(psi[1]))))),
    upper_diag(rho_spec, r = 4, c = 2, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(phi))))),
    upper_diag(rho_spec, r = 4, c = 3, 
               xlabel = as.expression(bquote(R(epsilon^(lambda), epsilon^(gamma))))),
    diag_hists(lev = "spec", column = 4, 
               xlabel = as.expression(bquote(epsilon^(lambda))), 
               xmin = spec_xmin, xmax = spec_xmax)
  ), 
  nrow = 4, byrow = FALSE)
#ggsave("fig/spec-level-ranefs.pdf", width = 7, height = 5)
ggsave("figures/spec-level-ranefs.png", width = 7, height = 5)


# phi - survival probability
logit_phi_summary <- MCMCsummary(jm, params = 'logit_phi', round=2)
MCMCtrace(jm, params = 'logit_phi', type = 'trace', ind = F, pdf=F)

# gamma - colonization probability
logit_gamma_summary <- MCMCsummary(jm, params = 'logit_gamma', round=2)
MCMCtrace(jm, params = 'logit_gamma', type = "trace", ind = F, pdf=F)

# lambda - expected abundance, given occupancy
lambda_summary <- MCMCsummary(jm, params = 'lambda', round=2)
MCMCsummary(jm, params = 'lambda', round=2) 

# Growth
log_growth_summary <- MCMCsummary(jm, params = 'log_growth', round=2)
MCMCsummary(jm, params = 'log_growth', round=2) #takes 5ish min

# Turnover
turnover_summary <- MCMCsummary(jm, params = 'turnover', round=2)
MCMCsummary(jm, params = 'turnover', round=2) #takes 5ish min

# Z
MCMCsummary(jm, params = 'Z', round=2)
MCMCtrace(jm, params = 'Z', type = 'density', ind = F, pdf=F)


# Visualize persistence and colonization
rownames_to_column(logit_phi_summary) %>%
      as_tibble %>%
      separate("rowname", into = c("site", "species"), sep = ",") %>%
      mutate_at(c('site', 'species'), readr::parse_number) %>%
      mutate(site = dimnames(z.init)[[1]][site], 
             species = dimnames(z.init)[[2]][species],
             mean = inv.logit(mean),
             param="Persistence") %>%
      rbind( rownames_to_column(logit_gamma_summary) %>%
                   as_tibble %>%
                   separate("rowname", into = c("site", "species"), sep = ",") %>%
                   mutate_at(c('site', 'species'), readr::parse_number) %>%
                   mutate(site = dimnames(z.init)[[1]][site], 
                          species = dimnames(z.init)[[2]][species],
                          mean = inv.logit(mean),
                          param="Colonization") ) %>%
      ggplot(aes(site, mean, col=param)) + 
      geom_pointrange(aes(ymin = inv.logit(`2.5%`), ymax = inv.logit(`97.5%`)), alpha=0.5) +
      facet_wrap(~species) + 
      xlab("Site") + 
      ylab("Rate") + 
      theme_minimal()

# Visualize growth and turnover
rownames_to_column(log_growth_summary) %>%
  as_tibble %>%
  separate("rowname", into = c("site", "species", "year"), sep = ",") %>%
  mutate_at(c('site', 'species', 'year'), readr::parse_number) %>%
  rename(perc2.5 = "2.5%",
         perc97.5 = '97.5%') %>%
  mutate(site = dimnames(z.init)[[1]][site], 
         species = dimnames(z.init)[[2]][species],
         year = dimnames(z.init)[[3]][year],
         mean = exp(mean),
         perc2.5 = exp(perc2.5),
         perc97.5 = exp(perc97.5),
         param="Growth") %>%
  # rbind( rownames_to_column(turnover_summary) %>%
  #          as_tibble %>%
  #          separate("rowname", into = c("site", "species", "year"), sep = ",") %>%
  #          mutate_at(c('site', 'species', 'year'), readr::parse_number) %>%
  #          rename(perc2.5 = "2.5%",
  #                 perc97.5 = '97.5%') %>%
  #          mutate(site = dimnames(z.init)[[1]][site], 
  #                 species = dimnames(z.init)[[2]][species], 
  #                 year = dimnames(z.init)[[3]][year],
  #                 param="Turnover" ) ) %>%
  ggplot(aes(year, mean, group = site, col=param)) + 
  geom_line() + 
  facet_wrap(~species) + 
  geom_ribbon(aes(ymin = perc2.5, ymax = perc97.5), color = NA, alpha = .05) +
  geom_hline(yintercept = 1) +
  xlab("Year") + 
  ylab("Rate") + 
  theme_minimal()


rownames_to_column(turnover_summary) %>%
           as_tibble %>%
           separate("rowname", into = c("site", "species", "year"), sep = ",") %>%
           mutate_at(c('site', 'species', 'year'), readr::parse_number) %>%
           rename(perc2.5 = "2.5%",
                  perc97.5 = '97.5%') %>%
           mutate(site = dimnames(z.init)[[1]][site],
                  species = dimnames(z.init)[[2]][species],
                  year = dimnames(z.init)[[3]][year],
                  param="Turnover" )  %>%
  ggplot(aes(year, mean, group = site, col=param)) + 
  geom_line() + 
  facet_wrap(~species) + 
  geom_ribbon(aes(ymin = perc2.5, ymax = perc97.5), color = NA, alpha = .05) +
  geom_hline(yintercept = 1) +
  xlab("Year") + 
  ylab("Rate") + 
  theme_minimal()

# Visualize encounter rate, lambda
rownames_to_column(lambda_summary) %>%
  as_tibble %>%
  separate("rowname", into = c("site", "species"), sep = ",") %>%
  mutate_at(c('site', 'species'), readr::parse_number) %>%
  mutate(site = dimnames(z.init)[[1]][site], 
         species = dimnames(z.init)[[2]][species]) %>%
  ggplot(aes(site, mean)) + 
  geom_pointrange(aes(ymin = (`2.5%`), ymax = (`97.5%`)), alpha=0.5) +
  facet_wrap(~species) + 
  xlab("Site") + 
  ylab("Encounter Rate") + 
  theme_minimal()
  