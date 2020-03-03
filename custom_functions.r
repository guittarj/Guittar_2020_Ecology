# Custom functions for Guittar et al 2020 Ecology

add_origins <- function(x, method = c('numeric', 'categorical')) {
	
	if (method == 'numeric') {
    x1 <- merge(
		  add_origin_numeric(x, 'temp.level') %>% rename(idT = id),
		  add_origin_numeric(x, 'precip.level') %>% rename(idP = id)
	  )
  } else {
    x1 <- merge(
      add_origin_categorical(x, 'temp.level') %>% rename(idT = id),
      add_origin_categorical(x, 'precip.level') %>% rename(idP = id)
    )
  }

  x1 <- left_join(x, x1, by = c("site", "sp", "stage", "abun"))
  	
  return(x1)

}

add_origin_categorical <- function(x, y = c('temp.level','precip.level')) {
  # Add putative climate origins to list of site/species
   
  if(y == 'temp.level') clab <- c('Cooler','Warmer','Same temperature')
  if(y == 'precip.level') clab <- c('Drier','Wetter','Same precipitation')
  
  y <- c('site', y)
  
  tmp <- site_meta[, y]
  names(tmp)[2] <- 'lev'
  
  #add climate column
  x <- left_join(x, tmp, by = 'site')
  
  #create set of local communities
  #don't count species labeled as transients (id = 'Transient')
  loc <- filter(x, stage == 'mature' & id == 'Persistent' & abun > 0)
    
  # Here, if a species is present in both lower and higher climates, it is considered to be from the same climate. Kind of conservative but makes sense. 
  for (i in c(1:nrow(x))) {
    if (x$id[i] != 'Persistent') {
      if (x$sp[i] %in% filter(loc, lev < x$lev[i])$sp) x$id[i] <- clab[1]
      if (x$sp[i] %in% filter(loc, lev > x$lev[i])$sp) x$id[i] <- clab[2]
      if (x$sp[i] %in% filter(loc, lev == x$lev[i])$sp) x$id[i] <- clab[3]
      if (x$sp[i] %in% filter(loc, lev < x$lev[i])$sp &
          x$sp[i] %in% filter(loc, lev > x$lev[i])$sp) x$id[i] <- clab[3]
    }
  }
  
  # clean up remaining categories.
  x$lev <- NULL
  x$id[x$id == 'Persistent'] <- 'Local'
  x$id[x$id == 'Transient'] <- 'Unknown'
  
  return(x)
}

add_origin_numeric <- function(x, y = c('temp.level','precip.level')) {
  
  # Add putative climate origins to list of site/species   
  y <- c('site', y)
  
  tmp <- site_meta[, y]
  names(tmp)[2] <- 'lev'
  
  #add climate column
  x <- left_join(x, tmp, by = 'site')
  
  #create set of local communities
  #don't count species labeled as transients (id = 'Transient')
  loc <- filter(x, stage == 'mature' & id == 'Persistent' & abun > 0)
    
  # Here, if a species is present in both lower and higher climates, it is considered to be from the same climate. Kind of conservative but makes sense. 
  for (i in c(1:nrow(x))) {
    if (x$id[i] != 'Persistent' & x$sp[i] %in% loc$sp) {
      if (x$sp[i] %in% filter(loc, lev < x$lev[i])$sp) x$id[i] <- max(filter(loc, sp == x$sp[i])$lev - x$lev[i])
      if (x$sp[i] %in% filter(loc, lev > x$lev[i])$sp) x$id[i] <- min(filter(loc, sp == x$sp[i])$lev - x$lev[i])
      if (x$sp[i] %in% filter(loc, lev == x$lev[i])$sp) x$id[i] <- 0
      if (x$sp[i] %in% filter(loc, lev < x$lev[i])$sp &
          x$sp[i] %in% filter(loc, lev > x$lev[i])$sp) x$id[i] <- 0
    }
  }
  
  # clean up remaining categories.
  x$lev <- NULL
  x$id[x$id == 'Persistent'] <- 0
  x$id[x$id == 'Transient'] <- NA
  
  return(x)
}

loadpax <- function(pkg){
  # (1) checks package installation, (2) installs them if not, then (3) loads them
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

stat_sum_df <- function(fun, geom = 'crossbar', ...) {
  # For plotting
  stat_summary(fun.data = fun, geom = geom, width = 0.2, ...)
}


plot_it <- function(plot, name, dir = "", height, width, res = 300) {
  # custom function to plot a high quality image.
  #to alter quality, increase res
  
  png(filename = paste0(dir, name, '.png'), units = 'in',
      type = "cairo", width = width, height = height, res = res)
  
  #only works for ggplot. base plot() is dumb
  plot(plot)
  
  papersize <- ifelse(width < 7.5 & height < 10, 'a4', 
                ifelse(width < 10 & height < 7.5, 'a4r', 'special'))
  
  pdf(file = paste0(dir, name, '.pdf'), width = width, height = height, paper = papersize)
  
  #only works for ggplot. base plot() is dumb
  plot(plot)
  
  #close out any lingering devs...
  while(dev.cur() != 1) dev.off()
  
}

stage_labels <- function(s, lvs = stage_lvs) {
  #add labels for stages
  s <- as.vector(lvs[match(s, names(lvs))])
  s <- factor(s, levels = lvs[lvs %in% s])
  return(s)
}

subbank <- function(x, frac) {
  #subsample seed bank community down to a stated fraction
  
	set.seed(7)
	tmp <- x %>% filter(stage == 'bank')
	subsampled_bank <- tmp[0, ]

	for (i in unique(tmp$site)) {
	  tmp0 <- filter(tmp, site == i)
	  tmp0 <- tmp0[rep(c(1:nrow(tmp0)), times = tmp0$abun), ]
	  tmp0 <- tmp0[sample(c(1:nrow(tmp0)), size = round(nrow(tmp0) * frac)), ]
	  tmp0$abun <- 1
	  subsampled_bank <- rbind(subsampled_bank, tmp0)
	}
  
	subsampled_bank <- subsampled_bank %>%
	  group_by(site, stage, sp, id) %>%
	  summarise(abun = sum(abun))

	x <- rbind(as.data.frame(filter(x, stage != 'bank')), 
	           as.data.frame(subsampled_bank))

	return(x)

}

mergeseeds <- function(x) {
  #merge seed rain and seed bank
  
	x <- x %>%
	  mutate(stage = ifelse(stage %in% c('rain','bank'), 'seed', stage)) %>%
	  group_by(site, stage, sp, id) %>%
	  summarise(abun = sum(abun)) %>%
	  ungroup()

	return(x)

}

rarefun <- function(df) {
 #calculates cumulative numbers of species observed 
  
  df_tmp <- df %>% 
    group_by(site) %>%
    mutate(turf = as.numeric(factor(turf, levels = sample(unique(turf))))) %>%
    group_by(site, id) %>%
    arrange(id, turf, sp) %>% 
    mutate(spp = as.numeric(factor(sp, levels = unique(sp)))) %>%
    select(-sp, -TTtreat, -Year, -abun, -siteID) %>%
    ungroup() %>%
    complete(site, turf, id, fill = list(spp = 0)) %>%
    group_by(site, id) %>%
    arrange(id, turf, spp) %>% 
    mutate(spp = cummax(spp)) %>%
    group_by(id, site, turf) %>%
    summarise(spp = max(spp))
  
  return(df_tmp)
}

scaleFUN <- function(x) {
  #transformation function to customize decimal

    x <- sprintf("%.1f", x)
  
  return(x)

}
