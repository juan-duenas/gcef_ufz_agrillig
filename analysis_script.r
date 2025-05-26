# load packages
pkgs <- c("vegan", "MASS", "tidyverse", "phyloseq", "ggpubr", "mgcv", "ggordiplots")

vapply(pkgs, FUN = library, FUN.VALUE = logical(1L), logical.return = TRUE, character.only = TRUE)


# Abundance#####
#load qPCR table, filter and select desired columns

aa <- read_csv("bf1.csv")%>%
      filter(!is.na(treat) & !is.na(time))%>%
      select(c(1,3:8,16:17))%>%
      rename(id=ID, bac_mab=b_copies_per_g_soil, fun_mab=f_copies_per_g_soil,
             fun_mab.log=log_f_copies_per_g_soil, bac_mab.log=log_b_copies_per_g_soil)%>%
      mutate(field=case_when(field %in% "Ambient" ~ "Current",
                             field %in% "Future" ~ "Future"))
# There are duplicated records for time 0  (10 in total) run unique(aa$id)
#cast long and then wide for ease of plotting
aa.l <- aa %>% select(-B_F_ratio)%>%
               pivot_longer(cols = 5:8,
               names_to = c("taxa","index"),
               names_sep = "\\_",
               values_drop_na = TRUE)%>%
               pivot_wider(names_from = index, values_from = value, values_fill = 0)

# Field Effects -> cycle 0 
fe <- aa.l %>%
  filter(time == "0")

# plot
ggplot(fe, aes(x = field, y = mab, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 0",
    x = "Climate",
    y = "copies of DNA per mg soil",
  ) +
  theme_bw()

# Lab effects - after 4 cycles
l4 <- aa.l %>% filter(time == "4") 

#richness
ggplot(l4, aes(x = field, y = mab, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 4",
    x = "Climate",
    y = "Copies of DNA per gram of soil",
  ) +
  theme_bw()

#all cycles together
#richness
ggplot(aa.l, aes(x = time, y = mab, color = field)) +
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5) +  # Add data points as dots) +
  geom_smooth(method = "glm", formula = y~x, method.args = list(family = gaussian(link = "log")))+
  scale_color_manual(values = c("Current" = "#117733", "Future" = "#661100")) +
  facet_wrap(lab~taxa, scale="free")+
  labs(
    title ="All DR cylcles",
    x = "cycle",
    y = "Copies of DNA per mg of Soil",
  ) +
  theme_bw()

# Alpha diversity ####
# paths

#load images with phyloseq objects
load("DR_bac_raw_rar.RData")
load("DR_fun_raw_rar.RData")

#Estimate Diversity in various organisms
Rb <- estimate_richness(b.r, measures=c("Observed")) %>% mutate(id=sample_data(b.r)[,1]$id) # rarefied
Sb <- estimate_richness(b, measures=c("InvSimpson")) %>% mutate(id=sample_data(b)[,1]$id) #raw
Rf <- estimate_richness(f.r, measures=c("Observed")) %>% mutate(id=sample_data(f.r)[,1]$id) #rarefied
Sf <- estimate_richness(f, measures=c("InvSimpson")) %>% mutate(id=sample_data(f)[,1]$id) #raw - ignore warning -> it was like that from the start
tr <- sample_data(b)%>% data.frame(., row.names = NULL) 

#merge everything in a dataset
dr.d <- list(tr,Rb,Sb,Rf,Sf) %>% reduce(inner_join, by="id") %>% #join treatments with div estimate tables
         rename(r.bacteria="Observed.x", is.bacteria="InvSimpson.x", r.fungi="Observed.y", is.fungi="InvSimpson.y")%>%
        pivot_longer(cols = 5:8,
                     names_to = c("index","taxa"),
                     names_sep = "\\.",
                     values_drop_na = TRUE)%>%
        pivot_wider(names_from = index, values_from = value, values_fill = 0)

rm(Rb,Sb,Rf,Sf,tr)

# creating subsets according to hypothesis 

# Field Effects -> cycle 0
fe <- dr.d %>%
  filter(time == "0")

#richness
ggplot(fe, aes(x = field, y = r, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 0",
    x = "Climate",
    y = "Richness",
    ) +
  theme_bw()

#Inverse Simpson
ggplot(fe, aes(x = field, y = is, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 0",
    x = "Climate",
    y = "Inverse Simpson",
  ) +
  theme_bw()

# cycle 1 - lab effects - immediate
l1 <- dr.d %>% filter(time == "1") 
  
#richness
ggplot(l1, aes(x = field, y = r, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 1",
    x = "Climate",
    y = "Richness",
  ) +
  theme_bw()

#Inverse Simpson
ggplot(l1, aes(x = field, y = is, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 1",
    x = "Climate",
    y = "Inverse Simpson",
  ) +
  theme_bw()

# Lab effects - after 4 cycles
l4 <- dr.d %>% filter(time == "4") 

#richness
ggplot(l4, aes(x = field, y = r, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 4",
    x = "Climate",
    y = "Richness",
  ) +
  theme_bw()

#Inverse Simpson
ggplot(l4, aes(x = field, y = is, fill = field)) +
  geom_violin(alpha=0.3)+
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5, color = "black") +  # Add data points as dots) +
  scale_fill_manual(values = c("Current" = "#117733", "Future" = "#661100"), guide = "none" ) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               color = "#000000")+
  facet_wrap(~taxa, scale="free")+
  labs(
    title ="Cycle 4",
    x = "Climate",
    y = "Inverse Simpson",
  ) +
  theme_bw()


#all cycles together
#richness
ggplot(dr.d, aes(x = time, y = r, color = field)) +
      geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5) +  # Add data points as dots) +
      geom_smooth(method = "glm", formula = y~x, method.args = list(family = poisson(link = "log")))+
      scale_color_manual(values = c("Current" = "#117733", "Future" = "#661100")) +
      facet_wrap(lab~taxa, scale="free")+
      labs(
      title ="All DR cylcles",
       x = "cycle",
      y = "Richness",
      ) +
  theme_bw()

#Inverse Simpson
ggplot(dr.d, aes(x = time, y = is, color = field)) +
  geom_jitter(shape = 16, size = 2, width = 0.1, alpha = 0.5) +  # Add data points as dots) +
  geom_smooth(method = "glm", formula = y~x, method.args = list(family = gaussian(link = "identity")))+
  scale_color_manual(values = c("Current" = "#117733", "Future" = "#661100")) +
  facet_wrap(lab~taxa, scale="free")+
  labs(
    title ="All DR cylcles",
    x = "cycle",
    y = "Inverse Simpson",
  ) +
  theme_bw()

# Beta diversity ####
#load datasets
load("DR_bac_raw_rar.RData")
load("DR_fun_raw_rar.RData")

# remove rarefied datasets
rm(b.r, f.r)

# filter lab manip - include only samples that went through cycles
b.dr <- subset_samples(b, lab%in%c("cycle_drw"))
f.dr <- subset_samples(f, lab%in%c("cycle_drw"))

# load medatata and datasets for RDA
asvtb <- decostand(t(as(otu_table(b.dr), "matrix")),"pa")
trb <- as(sample_data(b.dr), "data.frame")%>% mutate(id2=rownames(.)) #%>% data.frame(row.names = NULL)
asvtf <- decostand(t(as(otu_table(f.dr), "matrix")), "pa")
trf <- as(sample_data(f.dr), "data.frame")%>% mutate(id2=rownames(.))# %>% data.frame(row.names = NULL)
#setting permutation schema
perm.t<-how(nperm=999)

# RDAs for testing siginificance of factors
rda.b <- rda(asvtb~trb$field * trb$time , trb)
round(RsquareAdj(rda.b)$r.squared,3)#total variance explained
#significance tests
(test <- anova.cca(rda.b, permutations = perm.t, by="term", model="reduced"))

rda.f <- rda(asvtf~trf$field * trf$time , trf)
round(RsquareAdj(rda.f)$r.squared,3)#total variance explained
#significance tests
(test <-rbind(test, anova.cca(rda.f, permutations = perm.t, by="term", model="reduced")))

##Plots
clrs=c("#117733", "#661100")
# bacteria
# #extract sample scores to plot with ggplot2
o.b <- gg_ordiplot (rda.b, groups = trb$field, hull = F, label = T,
                    spiders = F, ellipse = T, plot = F, choices = c(1, 2), scaling=1)

#extract points
points.b <-o.b$df_ord %>% mutate(id2=rownames(.)) %>% left_join(trb, by="id2") %>% 
          select(x,y,field,time) %>% rename(tr=field)%>% mutate(time=as.factor(time))        

# extract ellipses
elip.b <- o.b$df_ellipse
colnames(elip.b) <- c("tr", "x","y")

#plot with ggplots
rda.p.b <- ggplot()
(rda.p.b <- rda.p.b + 
    geom_vline(xintercept=0.0, color="Grey", linewidth=0.8, linetype=2)+
    geom_hline(yintercept=0.0, color="Grey", linewidth=0.8, linetype=2)+    
    geom_point(data=points.b, aes(x=x, y=y, shape = time, colour=tr), alpha=0.8,  size=2) +
    geom_path(data=elip.b, aes(x=x, y=y, colour=tr), size=1)+
    scale_colour_manual(values = clrs)+
    #scale_x_continuous(limits = c(min(points$x)-0.2,max(points$x)+0.2))+
    labs(title="Bacteria - RDA", x="Axis 1", y="Axis 2", shape="Cycle", colour="Climate")+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 12),
          legend.key = element_blank(),  #removes the box around each legend item
          legend.position = "bottom", #legend at the bottom
          legend.text = element_text(size=12),
          #legend.title = element_blank(),
          panel.border = element_rect(colour = "Black", fill = F),
          panel.grid = element_blank()))

# fungi
# #extract sample scores to plot with ggplot2
o.f <- gg_ordiplot (rda.f, groups = trf$field, hull = F, label = T,
                    spiders = F, ellipse = T, plot = F, choices = c(1, 2), scaling=1)

#extract points
points.f <-o.f$df_ord %>% mutate(id2=rownames(.)) %>% left_join(trf, by="id2") %>% 
          select(x,y,field,time) %>% rename(tr=field)%>% mutate(time=as.factor(time))  # extracts points

# extract ellipses
elip.f <- o.f$df_ellipse
colnames(elip.f) <- c("tr", "x","y")

#plot with ggplots
rda.p.f <- ggplot()
(rda.p.f <- rda.p.f + 
    geom_vline(xintercept=0.0, color="Grey", linewidth=0.8, linetype=2)+
    geom_hline(yintercept=0.0, color="Grey", linewidth=0.8, linetype=2)+    
    geom_point(data=points.f, aes(x=x, y=y, color=tr, shape=time), alpha=0.8,  size=2) +
    geom_path(data=elip.f, aes(x=x, y=y, colour=tr), size=1)+
    scale_colour_manual(values = clrs)+
    #scale_x_continuous(limits = c(min(points$x)-0.2,max(points$x)+0.2))+
    labs(title="Fungi - RDA", x="Axis 1", y="Axis 2", shape="Cycle", colour="Climate")+
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 12),
          legend.key = element_blank(),  #removes the box around each legend item
          legend.position = "bottom", #legend at the bottom
          legend.text = element_text(size=12),
          #legend.title = element_blank(),
          panel.border = element_rect(colour = "Black", fill = F),
          panel.grid = element_blank()))

ggpubr::ggarrange(rda.p.b, rda.p.f, common.legend = T)
