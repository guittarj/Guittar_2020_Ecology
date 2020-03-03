#analysis script for Guittar et al. 2020 in Ecology
#Manuscript title: Quantifying the roles of seed dispersal, filtering, and climate on regional patterns of grassland biodiversity
#Authors: John guittar, Deborah Goldberg, Kari Klanderud, Astrid Berge, Marta Boixaderes,  Eric Meineri, Joachim Topper, Vigdis Vandvik.
#script by John Guittar. For questions email guittarj@gmail.com

#### Setup ####

# set working direction
wd <- '[insert working directory]'
setwd(wd)

# load packages
source("custom_functions.r")
loadpax(pkg = c('grid','knitr','vegan','lme4','tidyverse','gridExtra','kableExtra','ggpubr','lubridate','cowplot','VGAM'))

# load abundance data
x1 <- read.csv("data\\transitions_PersistentSppOccurOnce.csv", stringsAsFactors = FALSE)
x2 <- read.csv("data\\transitions_PersistentSppOccurTwice.csv", stringsAsFactors = FALSE)
x3 <- read.csv("data\\transitions_PersistentSppOccurThrice.csv", stringsAsFactors = FALSE)
x4 <- read.csv("data\\transitions_PersistentSppOccurAlways.csv", stringsAsFactors = FALSE)

#load veg abundance data for all years
cover <- read.csv('data\\cover.csv', header = TRUE, row.names = 1)

#load adult survey metadata, site metadata
cover_meta <- read.csv('data\\cover_meta.csv', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
site_meta <- read.csv("data\\site_meta.csv", stringsAsFactors = FALSE)

# Load species 'dictionary'.
dict <- read.csv("data\\SpeciesCodeDictionary.csv", stringsAsFactors = FALSE)

#Load trait data
traits <- read.csv(file = 'data\\traits.csv', stringsAsFactors = FALSE)

# seed bank data was sampled at a greater area than seedling and seed rain data... so...
# here I randomly subsample the seed bank data (without replacement) down to the fraction (based on area)
# equal to what it should be
frac <- (4 * 0.25^2) / (0.64^2)
x1 <- subbank(x1, frac)
x2 <- subbank(x2, frac)
x3 <- subbank(x3, frac)
x4 <- subbank(x4, frac)

# Create version of data where seed rain and seed bank are merged into one 'seed' category
xS_sp1 <- mergeseeds(x1)
xS_sp2 <- mergeseeds(x2)
xS_sp3 <- mergeseeds(x3)
xS_sp4 <- mergeseeds(x4)

# create versions with and without unidentified 'sp' individuals
# i do this now because it is easy to forget later...
# and in some cases I want to know the total density of seedlings irrespective of species identity
xS1 <- filter(xS_sp1, sp != 'sp')
xS2 <- filter(xS_sp2, sp != 'sp')
xS3 <- filter(xS_sp3, sp != 'sp')
xS4 <- filter(xS_sp4, sp != 'sp')
all_cutoffs <- bind_rows(
  mutate(xS1, group = '1/4 local occurrences'),
  mutate(xS2, group = '2/4 local occurrences'),
  mutate(xS3, group = '3/4 local occurrences'),
  mutate(xS4, group = '4/4 local occurrences'))
  
#for most analyses we will use x3 (persistent = 3 or more occurences)
x_sp <- x3
x <- filter(x3, sp != 'sp')
xS_sp <- xS_sp3
xS <- xS3

# calculate putative climate origins for both versions of data (seed rain and bank separated, and combined)
x_origins <- add_origins(x, 'categorical')
xS_origins <- add_origins(xS, 'categorical')

#trait label names
trait_labs <- c(
  leaf.area   = expression('Log Leaf Area ('*mm^2*')'), 
  max.height  = expression('Log Max'~'Height (m)'), 
  seed.mass   = expression('Log'~'Seed'~'Mass'~'(mg)'), 
  sla         = expression('Log SLA ('*m^2*kg^-1*')'),
  buds        = expression('Bud'~'Number'),
  lat         = expression('Lateral'~'Spread'~'(%)'), 
  offs        = expression('Offspring'~'(%)'), 
  conper      = expression('Conn. Persistence'~'(%)')
)

# stage labels and function
stage_lvs <- c(
  'rain' = 'Seed rain',
  'bank' = 'Seed bank', 
  'seed' = 'All seeds',
  'germ' = 'Emerged',
  'est' = 'Established',
  'mature' = 'Mature vegetation')

#### GLMS ####
#Note I originally tried to do a zero-inflated negative binomial but it has ~3x df and has a much higher AIC. Thus zero inflated models aren't an inferior choice here.

#A good stack overflow post on glm QC: https://stats.stackexchange.com/questions/70558/diagnostic-plots-for-count-regression
#the strategy here is to do GLMs at three tiers of resolution: (1) species grouped into transient/persistent; (2) species grouped into persistent, same temp/precip, different temp/precip; (3) species grouped into persistent, same temp/precip, warmer/drier, cooler/wetter, unknown temp/precip.

#organize data for models
j <- xS_origins %>%
  left_join(site_meta[, c('site','temp','precip')], by = 'site') %>%
  mutate(
    temp = as.numeric(scale(temp, scale = FALSE)), 
    precip = as.numeric(scale(precip, scale = FALSE)),
    idT = factor(idT, levels = c('Local','Same temperature','Cooler','Warmer','Unknown')),
    idP = factor(idP, levels = c('Local','Same precipitation','Drier','Wetter','Unknown'))) %>%
  spread(stage, abun, fill = 0)

# remove double zeroes so GLMs work
tmp <- filter(j, !(seed == 0 & germ == 0))

#calculate how many data points occur for each category, useful for later stats
germ_n <- c(table(tmp$idT), table(tmp$idP))[c(1,2,7,3,4,8,9,5)]
paste('N =', sum(germ_n))

#transform seed abun data with yeo-johnson transformation from VGAM
#first find lambda value which maximizes normality
#pval is shapiro wilks test for non-normality, so bigger pvals are better
y <- tmp$seed
nn <- length(y)
ltry <- seq(-0.15, -0.25, len = 1000)  # Try these values of lambda
lltry <- length(ltry)
psi <- matrix(as.numeric(NA), nn, lltry)
for (ii in 1:lltry) psi[, ii] <- yeo.johnson(y, lambda = ltry[ii])
colnames(psi) <- c(paste0("lambda_", ltry))
seed_lambda <- reshape2::melt(lapply(apply(psi, 2, shapiro.test), function(x) x$p.value)) %>%
  transmute(lambda = as.numeric(gsub('lambda_', "", L1)), pval = value) %>%
  filter(pval == max(pval)) %>%
  pull(lambda)

paste("best Yeo-Johnson lambda for seedling emergence GLM =", round(seed_lambda,3))

#transform seed abundance data using best lambda
tmp$seed_trans <- yeo.johnson(tmp$seed, lambda = seed_lambda)

#run all the GLMs with different parameters
germNull <-   MASS::glm.nb(germ ~ seed_trans, data = tmp)
germClim <-   MASS::glm.nb(germ ~ seed_trans + temp + precip, data = tmp)
germS <-      MASS::glm.nb(germ ~ seed_trans + temp * id + precip * id, data = tmp)
germT <-      MASS::glm.nb(germ ~ seed_trans + temp + precip + idT, data = tmp)
germP <-      MASS::glm.nb(germ ~ seed_trans + temp + precip + idP, data = tmp)

### GLMs to model individuals established
tmp <- filter(j, germ > 0)

#calculate how many data points occur for each category, useful for later stats
est_n <- c(table(tmp$idT), table(tmp$idP))[c(1,2,7,3,4,8,9,5)]
paste('N =', sum(est_n))

#tranform emergence data to normalize it
#first calculate best alpha
y <- tmp$germ
nn <- length(y)
ltry <- seq(-50, 50, len = 1000)  # Try these values of lambda
lltry <- length(ltry)
psi <- matrix(as.numeric(NA), nn, lltry)
for (ii in 1:lltry) psi[, ii] <- yeo.johnson(y, lambda = ltry[ii])
colnames(psi) <- c(paste0("lambda_", ltry))
germ_lambda <- reshape2::melt(lapply(apply(psi, 2, shapiro.test), function(x) x$p.value)) %>%
  transmute(lambda = as.numeric(gsub('lambda_', "", L1)), pval = value) %>%
  filter(pval == max(pval)) %>%
  pull(lambda)

paste("best Yeo-Johnson lambda for emerging seedlings =", round(germ_lambda,3))
tmp$germ_trans <- yeo.johnson(tmp$germ, lambda = germ_lambda)

estNull <-    MASS::glm.nb(est ~ germ_trans, data = tmp)
estClim <-    MASS::glm.nb(est ~ germ_trans + temp + precip, data = tmp)
estS <-      MASS::glm.nb(est ~ germ_trans + temp * id + precip * id, data = tmp)
estT <-      MASS::glm.nb(est ~ germ_trans + temp + precip  + idT, data = tmp)
estP <-      MASS::glm.nb(est ~ germ_trans + temp + precip  + idP, data = tmp)

# store coefficients
coefs <- list()
coefs$germNull <- as.data.frame(summary(germNull)$coefficients)
coefs$germClim <- as.data.frame(summary(germClim)$coefficients)
coefs$germS <- as.data.frame(summary(germS)$coefficients)
coefs$germT <- as.data.frame(summary(germT)$coefficients)
coefs$germP <- as.data.frame(summary(germP)$coefficients)
coefs$estNull <- as.data.frame(summary(estNull)$coefficients)
coefs$estClim <- as.data.frame(summary(estClim)$coefficients)
coefs$estS <- as.data.frame(summary(estS)$coefficients)
coefs$estT <- as.data.frame(summary(estT)$coefficients)
coefs$estP <- as.data.frame(summary(estP)$coefficients)

#variable labels for table
myvars <- c(
  'seed_trans' = 'Seed no. (transformed)',
  'germ_trans' = 'Seedling no. (transformed)',
  'temp' = 'Local temp.', 
  'precip' = 'Local precip.',
  'idTransient' = 'Transient',
  'temp:idTransient' = 'Transient * Local temp.',
  'idTransient:precip' = 'Transient * Local precip.',
  'idTSame temperature' = 'Transients from similar temp.',
  'idTCooler' = 'Transients from cooler into warmer',
  'idTWarmer' = 'Transients from warmer into cooler',
  'idTUnknown' = 'Transients from unknown climates',
  'idPSame precipitation' = 'Transients from similar precip.',
  'idPDrier' = 'Transients from drier into wetter',
  'idPWetter' = 'Transients from wetter into drier',
  'idPUnknown' = 'Transients from unknown climates'
)

#process/clean coefficients for table 1 etc
mod_coefs <- lapply(coefs, function(x) data.frame(variable = row.names(x), x))
mod_coefs <- do.call(rbind.data.frame, mod_coefs)
mod_coefs$Model <- unlist(lapply(strsplit(row.names(mod_coefs), '\\.'), function(x) x[1]))
mod_coefs <- mod_coefs %>%
  filter(!variable %in% c('(Intercept)')) %>%
  transmute(
    Model, 
    Variable = myvars[match(variable, names(myvars))], 
    Estimate = z.value, 
    p.value = Pr...z..) %>%
  arrange(Model, Variable)

#### Main text figures ####

###
fig2_seedsByTemp <- function(){}

#calculate log10 seed density and seed species richness per m2
j <- xS %>%
  filter(stage == 'seed') %>%
  mutate(den = abun / (4 * 0.25^2)) %>%
  group_by(site, id) %>%
  summarise(Abundance = sum(den), Richness = length(unique(sp))) %>%
  gather(variable, val, -site, -id) %>%
  left_join(site_meta[, c('site','temp','precip')], by = 'site') %>%
  mutate(id = ifelse(id == 'Persistent', 'Locally-persistent', 'Locally-transient'))

#calculate regression statistics - temperature
#log10 transform abundances for stats because that is how I will show data
stats_temp <- j %>%
  mutate(val = ifelse(variable == 'Abundance', log10(val), val)) %>%
  group_by(variable, id) %>%
  do(model = summary(lm(val ~ temp, data = .))) %>%
  mutate(pval = model$coef[[2,4]],
         r2 = model$r.squared)

stats_temp

#calculate regression statistics - precipitation
stats_precip <- j %>%
  mutate(val = ifelse(variable == 'Abundance', log10(val), val)) %>%
  group_by(variable, id) %>%
  do(model = summary(lm(val ~ precip, data = .))) %>%
  mutate(pval = model$coef[[2,4]],
         r2 = model$r.squared)

stats_precip

#add significance to data
j <- j %>%
  left_join(stats_temp[, c('id', 'variable', 'pval')], by = c('id','variable')) %>%
  mutate(sig = pval < 0.05)

#plot
p <- ggplot(j, aes(x = temp, y = val, color = id)) +
  geom_point(aes(shape = id)) +
  stat_smooth(aes(lty = sig), method = 'lm', show.legend = FALSE, se = FALSE) +
  scale_color_manual(values = c('black','red'), name = '') +
  scale_shape_manual(values = c(16,1), name = '') +
  scale_x_continuous(breaks = c(6:10)) +
  theme_classic() +
  theme(
    legend.position = 'bottom',
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  labs(x = expression('Mean Summer Temperature ('*degree*'C)'))

myleg <- get_legend(p)

p <- p + theme(legend.position = 'none')

p1 <- p %+% filter(j, variable == 'Abundance') +
  scale_linetype_manual(values = c(2,1)) +
  scale_y_log10() +
  labs(y = expression('Seeds/'*m^2), tag = 'a')

p2 <- p %+% filter(j, variable == 'Richness') +
  labs(y = 'Species', tag = 'b')

ps <- plot_grid(p1, p2, align = 'hv')

fig2_seedsByTemp <- plot_grid(ps, myleg, ncol = 1, rel_heights = c(1, .12))

plot_it(fig2_seedsByTemp, name = 'fig2_seedsByTemp', dir = "images\\", width = 6, height = 3.25)

fig2_seedsByTemp

##
fig3_transitionProbs <- function(){}

#spread and rename species status (with an extra whitespace to help readability)
j <- xS %>%
  spread(stage, abun, fill = 0) %>%
  mutate(id = ifelse(id == 'Persistent', ' Locally-persistent ', ' Locally-transient '))

# calculate maximums to use for setting data windows
maxtran <- max(filter(j, id == ' Locally-transient ')$seed)
maxgerm <- max(filter(j, id == ' Locally-transient ')$germ)

p1 <- j %>%
  filter(seed <= maxtran) %>%
  ggplot(aes(x = seed + 1, y = germ + 1, color = id, fill = id)) +
    geom_abline(slope = 1, lty = 3) +
    geom_point(aes(x = germ + 1, y = seed + 1), alpha = 0) +
    geom_point(alpha = 0.1) +
    stat_smooth(aes(x = jitter(seed + 1)), alpha = 0.5, method = 'loess') +
    scale_color_manual(values = c('black','red')) +
    scale_fill_manual(values = c('black','red')) +
    scale_x_log10(breaks = c(1,10,100,1000)) +
    scale_y_log10(breaks = c(1,10,100,1000)) +
    theme_classic() +
    theme(legend.position = 'bottom', legend.title = element_blank()) +
    labs(x = 'Seeds + 1', y = 'Emerged Seedlings + 1', tag = 'a')

p2 <- j %>%
  filter(germ > 0) %>%
  filter(germ <= maxgerm) %>%
  ggplot(aes(x = germ + 1, y = est + 1, color = id, fill = id)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, lty = 3) +
    stat_smooth(method = 'loess') +
    scale_color_manual(values = c('black','red')) +
    scale_fill_manual(values = c('black','red')) +
    scale_x_log10(breaks = c(2,10)) +
    scale_y_log10() +
    theme_classic() +
    theme(legend.position = 'bottom', legend.title = element_blank()) +
    labs(x = 'Emerged seedlings + 1', y = 'Established seedlings + 1', tag = 'b')

fig3_transitionProbs <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = 'bottom')

plot_it(fig3_transitionProbs, name = 'fig3_transitionProbs', dir = "images\\", width = 6.5, height = 3.5)

fig3_transitionProbs

###
fig4_traitsBiplots <- function(){}

#calculate community means, using combined seed rain seed bank
j <- xS %>%
  left_join(traits, by = 'sp') %>%
  gather(trait, val, seed.mass, sla, leaf.area, max.height, conper, offs, lat,  buds) %>%
  mutate(trait = factor(trait_labs[match(trait, names(trait_labs))], levels = trait_labs)) %>%
  filter(stage == 'seed') %>%
  ungroup() %>%
  group_by(site, stage, trait, id) %>%
  summarise(cm = mean(val, na.rm = T)) %>%
  group_by(site, stage, trait) %>%
  spread(id, cm) %>%
  filter(!is.na(Transient + Persistent)) %>%
  mutate(diff = Persistent - Transient) %>%
  left_join(site_meta, by = c("site"))

#test for trend in difference with temp + precip, to see if trait-based filtering varies with climate
mods.lm <- j %>%
  group_by(trait, stage) %>%
  do(mod = summary(lm(diff ~ temp + precip, data = .))) %>%
  mutate(p.temp = mod$coefficients[[12]],
         t.temp = mod$coefficients[[11]],
         mod = NULL)

#now do t-tests to see if transient and persistent spp groups differ
mods <- j %>%
  group_by(trait, stage) %>%
  do(mod = t.test(.$Transient, .$Persistent, paired = TRUE)) %>%
  mutate(pval = mod$p.value) %>%
  left_join(mods.lm, by = c("trait", "stage"))

#add modeling statistics to dataframe j.
#factorize climate to match figure 1 plotting scheme
j <- j %>%
  left_join(mods[, c('trait','stage','pval')], by = c("trait", "stage")) %>%
  mutate(sig = ifelse(pval < 0.05, '*', NA)) %>% 
  ungroup() %>%
  mutate(sig = ifelse(is.na(sig), 'p > 0.05', 'p < 0.05')) %>%
  mutate(MST = factor(c(6,9,10.5)[match(temp.level, c(1:3))], levels = c(6,9,10.5)),
         MAP = factor(c(650,1300,2000,2900)[match(precip.level, c(1:4))], levels = c(650,1300,2000,2900)))

fig4_traitsBiplots <- ggplot(j, aes(x = Persistent, y = Transient)) +
  geom_abline(slope = 1, lty = 3) +
  geom_point(aes(color = MAP, fill = interaction(sig, MAP), shape = MST), size = 3, alpha = 0.8) +
  geom_point(aes(y = Persistent, x = Transient), alpha = 0) +
  scale_color_manual(values = c('#87CEEB', '#749ED5', '#5057B2', '#140D8F')) +
  scale_fill_manual(values = c('#87CEEB', 'transparent', '#749ED5', 'transparent', '#5057B2', 'transparent', '#140D8F', 'transparent')) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN) +
  scale_shape_manual(values = c(24,21,25)) +
  facet_wrap(~trait, scales = 'free', ncol = 4, labeller = label_parsed) + 
  theme_bw() +
  theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      legend.position = 'none') +
  labs(x = 'Trait mean of locally-persistent species', y = 'Trait mean of locally-transient species')

plot_it(fig4_traitsBiplots, name = 'fig4_traitsBiplots', dir = "images\\", width = 7, height = 3.75)

fig4_traitsBiplots

##
table1_emergenceMods <- function(){}

j <- mod_coefs %>%
  filter(Model %in% c('germNull','germClim','germS','germP','germT')) %>%
  mutate(Variable = as.character(Variable)) %>%
  complete(Model, Variable) %>%
  mutate(Variable = factor(Variable, levels = myvars[-11])) %>%
  arrange(Model, Variable) %>%
  mutate(
    Estimate = sprintf("%.2f", Estimate),
    Estimate = ifelse(p.value < 0.001, paste0("**", Estimate),
                      ifelse(p.value < 0.05, paste0("*", Estimate), Estimate)),
    p.value = NULL) %>%
  spread(Model, Estimate, fill = NA) %>%
  mutate(Variable = as.character(Variable)) %>%
  select(Variable, germNull, germClim, germS, germT, germP)

AICs <- c('$\\Delta$AIC', sprintf("%.2f", c(c(AIC(germNull), AIC(germClim), AIC(germS), AIC(germT), AIC(germP)) - AIC(germNull))))

j <- rbind(AICs, j)

options(knitr.kable.NA = '-')

table1_emergenceMods <- kable(j, format = 'latex', escape = F, booktabs = T, linesep = "", align = c('l','r','r','r','r','r','r'),
                           col.names = linebreak(c("", "Null model", "Site climate", "Site climate +\nSp. status", "Site climate +\nSp. pref. temp.", "Site climate +\nSp. pref. precip"), align = 'r')) %>%
  kable_styling() %>%
  group_rows("General predictors", 2, 4) %>%
  group_rows("Transient/Persistent predictors", 5, 7) %>%
  group_rows("Origin-based predictors", 8, 14) %>%
  footnote(threeparttable = TRUE, general = 'Standardized coefficients (z-scores) from different GLM models (columns) predicting numbers of emerged seedlings by species and site. In column headers, "Sp. status" refers to whether the species is locally-transient or locally-persistent, and "Sp. pref. temp./precip." refers to the nearest temperatures/precipitations at which we found the species to have a persistent adult population, which we used to infer the climate from which they likely dispersed. Data comprise all recorded seeds and seedlings that could be identified to species. N is equal to 692, the number of unique emerged seedling species-by-site combinations. Asterisks denote significance (*: p < 0.05, **: p < 0.001). Dashes denote predictors that were not included in a given model.', general_title = "Table 1.", title_format = "bold") %>%
  landscape()

kable_as_image(table1_emergenceMods, filename = "images\\table1_emergenceMods", file_format = 'png', density = 500)
kable_as_image(table1_emergenceMods, filename = "images\\table1_emergenceMods", file_format = "pdf")

#### Supplementary figures ####
###
figS1_vegNMDS <- function(){}

#filter so only control turfs are present
control_turfs <- cover %>%
  mutate(id = row.names(cover)) %>%
  gather(sp, abun, -id) %>%
  left_join(cover_meta[, c('id','siteID','TTtreat','Year')], by = 'id') %>%
  filter(!is.na(id) & TTtreat %in% c('TTC','TT1'))

control_turfs <- control_turfs %>%
  group_by(siteID, Year, sp) %>%
  summarise(abun = sum(abun)) %>%
  spread(sp, abun, fill = 0)

ord <- metaMDS(control_turfs[, c(3:ncol(control_turfs))])

figS1_vegNMDS <- data.frame(control_turfs[, c('siteID', 'Year')], scores(ord)) %>%
  left_join(cover_meta, by = c("siteID", "Year")) %>%
  mutate(
    MST = factor(c(6,9,10.5)[match(Temperature_level, c(1:3))], levels = c(6,9,10.5)),
    MAP = factor(c(650,1300,2000,2900)[match(Precipitation_level, c(1:4))], levels = c(650,1300,2000,2900))) %>%
  ggplot(aes(x = NMDS1, y = NMDS2, color = MAP, fill = MAP)) +
  geom_point(aes(shape = MST)) +
  stat_ellipse(aes(group = interaction(MAP, MST))) +
  scale_color_manual(values = c('#87CEEB', '#749ED5', '#5057B2', '#140D8F')) +
  scale_fill_manual(values = c('#87CEEB', '#749ED5', '#5057B2', '#140D8F')) +
  scale_shape_manual(values = c(24,21,25)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none')

plot_it(figS1_vegNMDS, name = 'figS1_vegNMDS', dir = "images\\", width = 3.5, height = 3)

figS1_vegNMDS

###
figS2_rarefactions <- function(){}

#filter so only control turfs are present
j <- cover %>%
  mutate(id = row.names(cover)) %>%
  gather(sp, abun, -id) %>%
  left_join(cover_meta[, c('id','siteID','TTtreat','Year')], by = 'id') %>%
  filter(!is.na(id) & TTtreat %in% c('TTC','TT1')) %>% 
  mutate(site = substring(siteID, 1, 3)) %>%
  filter(abun > 0) %>%
  mutate(turf = gsub('_.*', '', id), id = NULL) %>%
  left_join(x[, c('site','sp', 'id')], by = c("sp", "site")) %>%
  filter(!is.na(id)) %>%
  mutate(id = ifelse(id == 'Persistent', 'Locally-persistent', 'Locally-transient'))

#create empty list of 100 elements and perform rarefaction of data for each
df_all <- list()
for (i in 1:100) {
  df_all[[i]] <- rarefun(j)
}

#summarise rarefactions
#remove Ram turf 10 because it doesn't exist...
tmp <- bind_rows(df_all) %>%
  group_by(id, site, turf) %>%
  summarise(sd = sd(spp), spp = mean(spp)) %>%
  left_join(site_meta, by = 'site') %>%
  mutate(MST = factor(c(6,9,10.5)[match(temp.level, c(1:3))], levels = c(6,9,10.5)),
         MAP = factor(c(650,1300,2000,2900)[match(precip.level, c(1:4))], levels = c(650,1300,2000,2900))) %>%
  filter(!(site == 'Ram' & turf == 10))

#plot
figS2_rarefactions <- ggplot(tmp, aes(x = turf, y = spp, color = MAP)) + 
  geom_point(aes(shape = MST, fill = MAP)) +
  geom_line(aes(group = site)) +
  scale_x_continuous(breaks = seq_along(1:10)) +
  scale_color_manual(values = c('#87CEEB', '#749ED5', '#5057B2', '#140D8F')) +
  scale_fill_manual(values = c('#87CEEB', '#749ED5', '#5057B2', '#140D8F')) +
  scale_shape_manual(values = c(24,21,25)) +
  expand_limits(y = 0) +
  facet_wrap(~id) + 
  labs(y = 'Number of species observed', x = 'Number of plots surveyed') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none')

plot_it(figS2_rarefactions, name = 'figS2_rarefactions', dir = 'images\\', width = 6, height = 3.5)

figS2_rarefactions

###
figS3_statusHistogram <- function(){}

#filter so only control turfs are present
j <- cover %>%
  mutate(id = row.names(cover)) %>%
  gather(sp, abun, -id) %>%
  left_join(cover_meta[, c('id','siteID','TTtreat','Year')], by = 'id') %>%
  filter(!is.na(id) & TTtreat %in% c('TTC','TT1')) %>% 
  mutate(site = substring(siteID, 1, 3)) %>%
  filter(abun > 0) %>%
  mutate(turf = gsub('_.*', '', id), id = NULL) %>%
  left_join(x[, c('site','sp', 'id')], by = c("sp", "site")) %>%
  filter(!is.na(id)) %>%
  mutate(id = ifelse(id == 'Persistent', 'Locally-persistent', 'Locally-transient'))

figS3_statusHistogram <- j %>% 
  group_by(site, sp, id) %>% 
  summarise(abun = sum(abun)) %>% 
  group_by(site) %>% 
  mutate(abun = abun / sum(abun)) %>% 
  ggplot(aes(x = abun, fill = id)) + 
  geom_histogram(position = 'dodge', bins = 50, alpha = 0.7) + 
  scale_y_sqrt(breaks = c(10,50,100)) + 
  scale_x_sqrt(breaks = c(0.01,0.1,0.2,0.3,0.4)) + 
  labs(x = 'Relative abundance at site', y = 'Number of site-level populations') + 
  scale_fill_manual(values = c('black','red'), name = '') +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'bottom')

plot_it(figS3_statusHistogram, name = 'figS3_statusHistogram', dir = "images\\", width = 6.5, height = 4)

figS3_statusHistogram

###
figS4_seedsByTempByStage <- function(){}

#calculate log10 seed density and seed species richness per m2
j <- x %>%
  filter(stage %in% c('rain','bank')) %>%
  mutate(den = abun / (4 * 0.25^2)) %>%
  group_by(stage, site, id) %>%
  summarise(Abundance = sum(den), Richness = length(unique(sp))) %>%
  gather(variable, val, -stage, -site, -id) %>%
  left_join(site_meta[, c('site','temp','precip')], by = 'site') %>%
  mutate(id = ifelse(id == 'Persistent', 'Locally-persistent', 'Locally-transient'))

#calculate regression statistics - temperature
#log10 transform abundances for stats because that is how I will show data
stats_temp <- j %>%
  mutate(val = ifelse(variable == 'Abundance', log10(val), val)) %>%
  group_by(stage, variable, id) %>%
  do(model = summary(lm(val ~ temp, data = .))) %>%
  mutate(pval = model$coef[[2,4]],
         r2 = model$r.squared)

stats_temp

#calculate regression statistics - precipitation
stats_precip <- j %>%
  mutate(val = ifelse(variable == 'Abundance', log10(val), val)) %>%
  group_by(variable, id) %>%
  do(model = summary(lm(val ~ precip, data = .))) %>%
  mutate(pval = model$coef[[2,4]],
         r2 = model$r.squared)

stats_precip

#add significance to data
j <- j %>%
  left_join(stats_temp[, c('stage','id', 'variable', 'pval')], 
            by = c('stage','id','variable')) %>%
  mutate(sig = pval < 0.05)

#plot
p <- ggplot(j, aes(x = temp, y = val, color = id)) +
  geom_point(aes(shape = id)) +
  stat_smooth(aes(lty = sig), method = 'lm', show.legend = FALSE, se = FALSE) +
  scale_color_manual(values = c('black','red'), name = '') +
  scale_shape_manual(values = c(16,1), name = '') +
  scale_x_continuous(breaks = c(6:10)) +
  theme_classic() +
  theme(
    legend.position = 'bottom',
    strip.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "pt")) +
  labs(x = expression('Mean Summer Temperature ('*degree*'C)'))

myleg <- get_legend(p)

p <- p + theme(legend.position = 'none')

p1 <- p %+% filter(j, stage == 'rain' & variable == 'Abundance') +
  scale_linetype_manual(values = c(2,1)) +
  scale_y_log10() +
  labs(x = '', y = expression('Seeds/'*m^2), tag = 'a')

p2 <- p %+% filter(j, stage == 'rain' & variable == 'Richness') +
  labs(y = 'Species', x = '', tag = 'b') +
  scale_linetype_manual(values = c(2))

p3 <- p %+% filter(j, stage == 'bank' & variable == 'Abundance') +
  scale_linetype_manual(values = c(2,1)) +
  scale_y_log10() +
  labs(y = expression('Seeds/'*m^2), tag = 'c')

p4 <- p %+% filter(j, stage == 'bank' & variable == 'Richness') +
  labs(y = 'Species', tag = 'd') +
  scale_linetype_manual(values = c(1))

ps <- plot_grid(p1, p2, p3, p4, ncol = 2, align = 'hv')

figS4_seedsByTempByStage <- plot_grid(ps, myleg, ncol = 1, rel_heights = c(1, .12))

plot_it(figS4_seedsByTempByStage, name = 'figS4_seedsByTempByStage', dir = "images\\", width = 6, height = 6)

figS4_seedsByTempByStage


##
figS5_allCutoffs <- function(){}

#fix labels
j <- all_cutoffs %>%
  spread(stage, abun, fill = 0) %>%
  mutate(id = ifelse(id == 'Persistent', 'Locally-persistent', 'Locally-transient'))

"As the cutoff for locally-persistent became more stringent, and the cutoff for locally-transient (by definition) relaxed, the total number of locally-transient species-by-site combinations in the combined seed rain and seed bank rose from 119 (1989 seeds), to 149 (2549 seeds), to 167 (3665 seeds), to 205 (5007 seeds)"
j %>% group_by(group) %>% filter(id == 'Locally-transient') %>% summarise(rich = length(seed[seed > 0]), abun = sum(seed))

# calculate maximums to use for setting data windows
j1 <- j %>%
  filter(seed <= max(seed[id == 'Locally-transient'])) %>% 
  mutate(transition = 'Emergence')

j2 <- j %>%
  filter(germ > 0) %>%
  group_by(group) %>% 
  filter(germ <= max(germ[id == 'Locally-transient'])) %>% 
  mutate(transition = 'Establishment')

#emergence
p1 <- j1 %>%
  ggplot(aes(color = id, fill = id)) +
  geom_abline(slope = 1, lty = 3) +
  geom_point(aes(x = germ + 1, y = seed + 1), alpha = 0) +
  geom_point(aes(x = seed + 1, y = germ + 1), alpha = 0.1) +
  stat_smooth(aes(x = jitter(seed + 1), y = germ + 1), alpha = 0.5, method = 'loess') +
  scale_color_manual(values = c('black','red')) +
  scale_fill_manual(values = c('black','red')) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~group, ncol = 4, scales = 'free') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.position = 'bottom') +
  labs(x = 'Seeds + 1', y = 'Emerged Seedlings + 1')

#establishment
p2 <- j2 %>%
  ggplot(aes(x = germ + 1, y = est + 1, color = id, fill = id)) +
  geom_point(aes(x = est + 1, y = germ + 1), alpha = 0, data = filter(j2, est > 1)) +
  geom_point(alpha = 0.1) +
  geom_abline(slope = 1, lty = 3) +
  stat_smooth(method = 'loess') +
  scale_color_manual(values = c('black','red')) +
  scale_fill_manual(values = c('black','red')) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~group, ncol = 4, scales = 'free') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.position = 'none',
    legend.title = element_blank()) +
  labs(x = 'Emerged seedlings + 1', y = 'Established seedlings + 1')

#combine
figS5_allCutoffs <- ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = 'bottom')

#plot
plot_it(figS5_allCutoffs, name = 'figS5_allCutoffs', dir = "images\\", width = 6.5, height = 4.5)

figS5_allCutoffs


###
tableS1_siteStats <- function(){}

# first take average density across blocks
# then sum densities of seed rain and seed bank
# I want to keep sp's here for den measure
j <- rbind(
  x_sp,
  filter(xS_sp, stage == 'seed'),
  filter(x_sp, stage != 'mature') %>% mutate(stage = 'seeds.seedlings'),
  mutate(x_sp, stage = 'all'))

j_den <- j %>%
  filter(!stage %in% c('mature','all')) %>%
  group_by(stage, site) %>%
  mutate(den = abun / (4 * 0.25^2)) %>%
  summarise(den = sum(den)) %>%
  group_by(stage) %>%
  summarise(
    den_sd = sd(c(den, rep(0, 12 - length(den)))), 
    den = mean(c(den, rep(0, 12 - length(den)))))
  
j_den_trans <- j %>%
  filter(id == 'Transient') %>%
  filter(!stage %in% c('mature','all')) %>%
  group_by(stage, site) %>%
  mutate(den = abun / (4 * 0.25^2)) %>%
  summarise(den = sum(den)) %>%
  group_by(stage) %>%
  summarise(
    den_sd = sd(c(den, rep(0, 12 - length(den)))), 
    den = mean(c(den, rep(0, 12 - length(den)))))

j_site_richness <- j %>%
  filter(sp != 'sp') %>%
  group_by(site, stage) %>%
  summarise(rich = length(unique(sp))) %>%
  group_by(stage) %>%
  summarise(
    site_sd = sd(c(rich, rep(0, 12 - length(rich))), na.rm = T), 
    site_rich = mean(c(rich, rep(0, 12 - length(rich)))))

j_site_richness_trans <- j %>%
  filter(sp != 'sp') %>%
  filter(id == 'Transient') %>%
  group_by(site, stage) %>%
  summarise(rich = length(unique(sp))) %>%
  group_by(stage) %>%
  summarise(
    site_sd = sd(c(rich, rep(0, 12 - length(rich))), na.rm = T), 
    site_rich = mean(c(rich, rep(0, 12 - length(rich)))))
  
j_regional_richness <- j %>%
  filter(sp != 'sp') %>%
  group_by(stage) %>%
  summarise(rich = length(unique(sp)))

j_regional_richness_trans <- j %>%
  filter(sp != 'sp') %>%
  filter(id == 'Transient') %>%
  group_by(stage) %>%
  summarise(rich = length(unique(sp)))

j_all <- j_den %>%
  full_join(j_site_richness, by = 'stage') %>%
  full_join(j_regional_richness, by = 'stage')

j_trans <- j_den_trans %>%
  full_join(j_site_richness_trans, by = 'stage') %>%
  full_join(j_regional_richness_trans, by = 'stage')

#custom stage names
#stage name dictionary
stage_lvs2 <- c('Seed Rain' = 'rain', 'Seed Bank' = 'bank', 'All Seeds' = 'seed', 'Emerged Seedlings' = 'germ', 'Established Seedlings' = 'est', 'Seeds and Seedlings' = 'seeds.seedlings', 'Adults' = 'mature', 'All Stages' = 'all')

j <- rbind(mutate(j_all, group = 'all'), 
           mutate(j_trans, group = 'trans')) %>%
  filter(!(group == 'trans' & stage %in% c('mature','all'))) %>%
  transmute(
    group,
    Stage = factor(names(stage_lvs2)[match(stage, stage_lvs2)], names(stage_lvs2)),
    `Density (per sq. m.)` = ifelse(is.na(den), '', paste(round(den), '±', round(den_sd))),
    `Site richness` = paste(round(site_rich), '±', round(site_sd)),
    `Regional richness` = ifelse(group == 'all', rich, NA)) %>%
  arrange(group, Stage)
         
tmp <- j %>% select(-group)

tableS1_siteStats <- kable(tmp, format = 'latex', booktabs = T, escape = F, align = c('l', 'r','r','r')) %>%
    kable_styling() %>%
    row_spec(0, bold = T) %>%
    group_rows("All individuals", 1, 8) %>%
    group_rows("Individuals of locally-transient species only", 9, 14) %>%
    footnote(general = 'Densities and species richness values within and across sites for each life stage. The density of "all individuals" includes unidentified seeds and seedlings. Regional richness is not shown for locally-transient species because locally-transient/locally-persistent species status can vary by site, and thus cannot be summarized in this way across sites.', general_title = "Table S1.", threeparttable = T, title_format = "bold") %>%
  landscape()

kable_as_image(tableS1_siteStats, filename = "images\\tableS1_siteStats", file_format = 'png', density = 500)
kable_as_image(tableS1_siteStats, filename = "images\\tableS1_siteStats", file_format = "pdf")

###
tableS2_individuals <- function(){}

j <- rbind(x_origins, filter(xS_origins, stage == 'seed')) %>%
  mutate(id = ifelse(id == 'Persistent', 'Locally-persistent', 
                      ifelse(id == 'Transient', 'Locally-transient', id)))

j <- rbind(
  mutate(j, id = 'All individuals'),
  mutate(j),
  mutate(j, id = ifelse(idT == 'Unknown', 'Unknown temperature', idT)) %>% filter(id != 'Local'),
  mutate(j, id = ifelse(idP == 'Unknown', 'Unknown precipitation', idP)) %>% filter(id != 'Local')) %>%
  mutate(
    id = factor(id, levels = c('All individuals','Locally-persistent','Locally-transient','Same temperature','Cooler','Warmer','Unknown temperature','Same precipitation','Drier','Wetter','Unknown precipitation')),
    stage = stage_labels(stage)) %>%
  group_by(id, stage) %>%
  summarise(abun = round(sum(abun))) %>%
  group_by(stage) %>%
  mutate(percent = sprintf("%.1f", 100 * abun / abun[id == 'All individuals'], digits = 1)) %>%
  ungroup() %>%
  gather(var, val, abun, percent) %>%
  mutate(var = paste(stage, var, sep = "_"), stage = NULL) %>%
  spread(var, val, fill = 0)

j <- j[, c('id','Seed rain_abun','Seed rain_percent',
                'Seed bank_abun','Seed bank_percent',
                'All seeds_abun','All seeds_percent',
                'Emerged_abun','Emerged_percent',
                'Established_abun','Established_percent')]

colnames(j) <- c('Species status',rep(c('No.','%'), 5)) 
j$`Species status` <- as.character(j$`Species status`)
j$`Species status`[c(4,7,8,11)] <- c('Same', 'Unknown','Same','Unknown')

tableS2_individuals <- kable(j, format = 'latex', escape = T, booktabs = T, align = c('l','r','r','r','r','r','r','r','r','r','r'), linesep = "") %>%
  kable_styling() %>%
  add_indent(c(2:11)) %>%
  #add_indent(c(4:11)) %>%
  group_rows("Temperature of nearest persistent adult population", 4, 7) %>%
  group_rows("Precipitation of nearest persistent adult population", 8, 11) %>%
  row_spec(0, bold = T) %>%
  add_header_above(c(" " = 1, "Seed rain" = 2, "Seed bank" = 2, "All seeds" = 2, "Emerged" = 2, "Established" = 2), bold = T) %>%
  footnote(general = "Numbers of individuals recorded, grouped by locally-transient/locally-persistent species status (top three rows). Then, locally-transient individuals are further grouped by the putative climate preferences of their species, i.e., the temperatures/precipitations of the nearest sites at which persistent adult populations are known to occur, relative the climates of the local sites. Percentages are of all individuals (top row). Individuals that could not be identified to species were not included.", general_title = "Table S2.", threeparttable = T, title_format = "bold") %>%
  landscape()

kable_as_image(tableS2_individuals, filename = "images\\tableS2_individuals", density = 500)
kable_as_image(tableS2_individuals, filename = "images\\tableS2_individuals", file_format = "pdf")

###
tableS3_establishmentMods <- function(){}

j <- mod_coefs %>%
  filter(Model %in% c('estNull','estClim','estS','estP','estT')) %>%
  mutate(Variable = as.character(Variable)) %>%
  complete(Model, Variable) %>%
  mutate(Variable = factor(Variable, levels = myvars[-11])) %>%
  arrange(Model, Variable) %>%
  mutate(
    Estimate = sprintf("%.2f", Estimate),
    Estimate = ifelse(p.value < 0.001, paste0("**", Estimate),
                           ifelse(p.value < 0.05, paste0("*", Estimate), Estimate)),
    p.value = NULL) %>%
  spread(Model, Estimate, fill = NA) %>%
  mutate(Variable = as.character(Variable)) %>%
  select(Variable, estNull, estClim, estS, estT, estP)

AICs <- c('$\\Delta$AIC', sprintf("%.2f", c(c(AIC(estNull), AIC(estClim), AIC(estS), AIC(estT), AIC(estP)) - AIC(estNull))))

j <- rbind(AICs, j)

options(knitr.kable.NA = '-')

tableS3_establishmentMods <- kable(j, format = 'latex', escape = F, booktabs = T, linesep = "", align = c('l','r','r','r','r','r','r'),
                               col.names = linebreak(c("", "Null model", "Site climate", "Site climate +\nSp. status", "Site climate +\nSp. pref. temp.", "Site climate +\nSp. pref. precip"), align = 'r')) %>%
  kable_styling() %>%
  group_rows("General predictors", 2, 4) %>%
  group_rows("Transient/Persistent predictors", 5, 7) %>%
  group_rows("Origin-based predictors", 8, 14) %>%
  footnote(threeparttable = TRUE, general = 'Standardized coefficients (z-scores) from different GLM models (columns) predicting numbers of established seedlings by species and site. In column headers, "Sp. status" refers to whether the species is locally-transient or locally-persistent, and "Sp. pref. temp./precip." refers to the nearest temperatures/precipitations at which we found the species to have a persistent adult population, which we used to infer the climate from which they likely dispersed. The predictor "Seedling no. (transformed)" refers to the numbers of emerged seedlings of each species at each site, normalized with Yeo-Johnson transformations (refer to Methods). Data consisted of all recorded emerged/established seedlings that could be identified to species. N is equal to 692, the number of unique emerged seedling species-by-site combinations. Asterisks denote significance (*: p < 0.05, **: p < 0.001). Dashes denote predictors that were not included in a given model.', general_title = "Table S3.", title_format = "bold") %>%
  landscape()

kable_as_image(tableS3_establishmentMods, filename = 'images\\tableS3_establishmentMods', density = 500)
kable_as_image(tableS3_establishmentMods, filename = "images\\tableS3_establishmentMods", file_format = "pdf")

#### Misc. stats ####

"The study area comprises 12 semi-natural calcareous grassland sites in southern Norway that host at least 144 non-woody vascular plant species at the adult life stage..."
length(unique(x$sp[x$stage == 'mature']))

"...and at least 126 at the seed stage"
length(unique(xS$sp[xS$stage == 'seed']))

"Transient seeds occurred at all 12 grassland sites"
x %>% filter(id == 'Transient') %>% summarise(sites = length(unique(site)))

"How many seedlings did we record (including unidentified seedlings)? How many established?"
sum(x_sp$abun[x_sp$stage == 'germ'])
sum(x_sp$abun[x_sp$stage == 'est'])

"How many adult species were persistent? i.e., persistent adult richness?"
x %>% filter(stage == 'mature' & id == 'Persistent') %>% summarise(adult_rich = length(unique(sp)))

"How many adults were never seen as seeds or seedlings outside of thier local sites"
sum(!unique(x$sp[x$stage == 'mature' & x$id == 'Persistent']) %in% xS$sp[xS$stage %in% c('seed','seedling') & xS$id == 'Transient'])

"How many persistent adult species had no seed or seedlings observed anywhere"
sum(!unique(x$sp[x$stage == 'mature' & x$id == 'Persistent']) %in% xS$sp[xS$stage %in% c('seed','germ')])

"What were those species?"
unique(x$sp[x$stage == 'mature' & x$id == 'Persistent'])[!unique(x$sp[x$stage == 'mature' & x$id == 'Persistent']) %in% xS$sp[xS$stage %in% c('seed','germ')]]

"At how many sites did transient species emerge?" 
x %>% filter(id == 'Transient' & stage == 'germ') %>% group_by(site, stage) %>% summarise(abun = sum(abun), rich = length(unique(sp))) %>% ungroup() %>% mutate(all_trans = sum(abun), all_sp = sum(rich), all = sum(x$abun[x$stage == 'germ']), per = all_trans/all) %>% filter(abun > 0) %>% nrow()

"At how many sites did transient species establish?"
"How many transient species established at sites?"
x %>% filter(id == 'Transient' & stage == 'est') %>% group_by(site, stage) %>% summarise(abun = sum(abun), rich = length(unique(sp))) %>% ungroup() %>% mutate(all_trans = sum(abun), all_sp = sum(rich), all = sum(x$abun[x$stage == 'est']), per = all_trans/all) %>% filter(abun > 0) %>% mutate(total_sites = length(all))

"how many species dispersed into cooler/wetter sites"
xS_origins %>%
  filter(stage == 'seed') %>%
  group_by(idT) %>%
  summarise(abun = sum(abun), rich = length(unique(sp))) %>%
  ungroup() %>%
  mutate(percentage = abun / sum(abun))
xS_origins %>%
  filter(stage == 'seed') %>%
  group_by(idP) %>%
  summarise(abun = sum(abun), rich = length(unique(sp))) %>%
  ungroup() %>%
  mutate(percentage = abun / sum(abun))

#significant trend in species richness in adult vegetation
x %>%
  filter(stage == 'mature') %>%
  left_join(site_meta, by = 'site') %>%
  group_by(site, temp, precip) %>%
  summarise(rich = length(unique(sp[abun > 0]))) %>%
  ungroup() %>%
  do(
    modT = summary(lm(rich ~ temp, data = .)),
    modP = summary(lm(rich ~ precip, data = .))) %>%
  mutate(pvalT = modT[[1]]$coef[[2,4]], r2T = modT[[1]]$r.squared,
         pvalP = modP[[1]]$coef[[2,4]], r2P = modP[[1]]$r.squared)


# Transition probabilities, NOT GROUPED BY SPECIES.
xS_sp %>%
  filter(stage != 'mature') %>%
  group_by(site, stage) %>%
  mutate(den = abun / (4 * 0.25^2)) %>%
  summarise(den = sum(den)) %>%
  spread(stage, den, fill = 0) %>%
  ungroup() %>%
  summarise(
    Emergence = paste(sprintf("%.2f", mean(germ/seed)), "±", 
                      sprintf("%.2f", sd(germ/seed))),
    Establishment = paste(sprintf("%.2f", mean(est/germ)), "±", 
                          sprintf("%.2f", sd(est/germ))),
    `Seed-to-Establishment` = paste(sprintf("%.2f", mean(est/seed)), "±", 
                                    sprintf("%.2f", sd(est/seed)))) %>%
  gather(Transition, `Rate`, Emergence, Establishment, `Seed-to-Establishment`)
