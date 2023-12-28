#start
#Figure 1
#Figure 1A

## take the proper subset of data
atrc <- fecundity %>% 
  filter(vial%in%c("35u","30p","4d","11c","28n","2b","71k","1a","75o","3c","70j","15h","16i","13g","Gno","Ax")) %>% 
  filter(!is.na(fecundity)) %>% 
  filter(fecundity!="NaN") %>% 
  droplevels()

## adjust levels to maintain plot order
atrc$vial <- factor(atrc$vial, levels = c("35u","30p","4d","11c","28n","1a","71k","2b","75o","3c","70j","15h","16i","13g","Gno","Ax"))

## build vectors for plot shading
plot_colors <- c("red","red","red","red","red","blue","blue","blue","blue","blue","blue","green","green","purple","magenta","black")
plot_lines <- c("solid","longdash","dotted","dotdash","dashed","solid","longdash","dotted","dotdash","dashed","twodash","solid","longdash","solid","solid","solid")

## make the plot
ggplot(atrc, aes(x = date3, y=fecundity, group = vial, col = vial, fill = vial, linetype = vial)) + 
  geom_smooth(aes(y = fecundity, linetype = vial), alpha = 0.05, size = .5) +
  scale_color_manual(values = plot_colors, name = NULL, labels = c(expression(italic(A.~aceti)),expression(italic(A.~pasteurianus)),expression(italic(A.~pomorum)~DmCS),expression(italic(A.~tropicalis)~DmCS),expression(italic(A.~tropicalis)~NBRC),expression(italic(L.~brevis)~DmCS),expression(italic(L.~brevis)~gravesensis),expression(italic(L.~plantarum)~DmCS),expression(italic(L.~plantarum)~WCFS1),expression(italic(L.~fructivorans)~DmCS),expression(italic(L.~fructivorans)~KCTC),expression(italic(E.~coli)),expression(italic(P.~putida)),expression(italic(B.~subtilis)),'5-species','Axenic')) + 
  scale_fill_manual(values = plot_colors, name = NULL,labels = c(expression(italic(A.~aceti)),expression(italic(A.~pasteurianus)),expression(italic(A.~pomorum)~DmCS),expression(italic(A.~tropicalis)~DmCS),expression(italic(A.~tropicalis)~NBRC),expression(italic(L.~brevis)~DmCS),expression(italic(L.~brevis)~gravesensis),expression(italic(L.~plantarum)~DmCS),expression(italic(L.~plantarum)~WCFS1),expression(italic(L.~fructivorans)~DmCS),expression(italic(L.~fructivorans)~KCTC),expression(italic(E.~coli)),expression(italic(P.~putida)),expression(italic(B.~subtilis)),'5-species','Axenic')) + 
  scale_linetype_manual(values = plot_lines,name = NULL,labels = c(expression(italic(A.~aceti)),expression(italic(A.~pasteurianus)),expression(italic(A.~pomorum)~DmCS),expression(italic(A.~tropicalis)~DmCS),expression(italic(A.~tropicalis)~NBRC),expression(italic(L.~brevis)~DmCS),expression(italic(L.~brevis)~gravesensis),expression(italic(L.~plantarum)~DmCS),expression(italic(L.~plantarum)~WCFS1),expression(italic(L.~fructivorans)~DmCS),expression(italic(L.~fructivorans)~KCTC),expression(italic(E.~coli)),expression(italic(P.~putida)),expression(italic(B.~subtilis)),'5-species','Axenic')) + 
  ylab("female fecundity per day") +
  xlab("time (d)") +
  guides(fill = guide_legend(ncol = 3, byrow=T), color = guide_legend(ncol = 3, byrow=T), linetype = guide_legend(ncol = 3, byrow=T)) +
  theme_cowplot() + 
  theme(legend.position = "none", legend.text.align = 0,legend.key.width = unit(.5,"in"), legend.text = element_text(size = 18), axis.text = element_text(size = 10),axis.title = element_text(size = 12),legend.justification = "center")
F


## this code is derived from loops used for supplemental figures, so in some cases it looks like it has more variables than it needs

#print 1 plot with  the same survival in each line. each line is the age class that was permuted for survival
q = 1
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
rangevals <- c(1)
m <- rangevals
n = 1
eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric())
eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
mod = n
age_class = m

subfitness <- fitness %>% filter(vial %in% subset_vec)

for(i in names(table(list(subfitness$trt)))) {
  submat <- fitness %>% 
    filter(trt==i) %>% 
    filter(date3 == min(date3)) %>% 
    mutate(date3 = 0) %>%
    mutate(fecundity = 0) %>%
    rbind(fitness %>% 
            filter(trt==i) %>% 
            arrange(date3)
    ) %>%
    rbind(fitness %>%
            filter(trt == i) %>%
            filter(date3 == max(date3)) %>% 
            mutate(date3 = max(date3)+1) %>%
            mutate(fecundity = 0)
    )
  submat$s = 1
  submat$f = 0
  
  for (j in 1:(dim(submat)[1]-1)) {
    submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
    submat$f[j] <- submat$fecundity[j] #* submat$s[j]
  }
  
  leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
  if (age_class>1) {
    for (k in 1:(age_class-1)) {
      leslie_mat <- rbind(leslie_mat,
                          matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
      )
    } 
  }
  
  for (k in age_class:min(age_class,(dim(submat)[1]-1))) {
    leslie_mat <- rbind(leslie_mat, 
                        matrix(c(rep(0,k-1),submat$s[k]*mod,rep(0,(dim(submat)[1]-k))), nrow=1)
    )
  }
  if(age_class < (dim(submat)[1]-1)) {
    for(k in (age_class+1):(dim(submat)[1]-1)) {
      leslie_mat <- rbind(leslie_mat,
                          matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
      )
    }
  }
  
  eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
  
}

eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(fitness = log(fitness)) %>% mutate(mod_death = (m)))

e3 <- eigens
eigens <- eigens %>% filter(trt%in%subset_vec) %>% droplevels()
eigens$trt <- factor(eigens$trt, levels = subset_vec)

efg <- dunn.test(eigens$fitness,eigens$trt, method = "bh", table = F)
ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05, print.comp = F) %>% 
  arrange(factor(Group, levels = c("35u","3p","4d","11c","28n","1a","71k","2b","75o","3c","7j","15h","16i","13g","Gno","Ax")))

eigens2 <- eigens %>% group_by(trt, mod_death) %>% summarize(meanfit = mean(fitness), sem = sd(fitness)/sqrt(sum(fitness>-1000))) %>% ungroup()
eigens2$trt <- factor(eigens2$trt, levels = c("35u","30p","4d","11c","28n","1a","71k","2b","75o","3c","70j","15h","16i","13g","Gno","Ax"))
eigens2 <- eigens2 %>% arrange(factor(trt, levels = c("35u","30p","4d","11c","28n","1a","71k","2b","75o","3c","70j","15h","16i","13g","Gno","Ax")))
eigens2$trt <- plyr::revalue(eigens2$trt, c("35u" = "aace","30p" = "apan","4d" = "apoc","11c"="atrc", "28n"="atrn","1a" = "lbrc","71k" = "lbga","2b" = "lplc", "75o" = "lplw","3c" = "lfrc","70j" = "lfrk","15h" = "ecok", "16i" = "pput","13g" = "bsub","Ax" = "Ax","Gno" = "5sp"))

plot_colors <- c("red","red","red","red","red","blue","blue","blue","blue","blue","blue","green","green","purple","magenta","black")
plot_lines <- c("solid","longdash","dotted","dotdash","dashed","solid","longdash","dotted","dotdash","dashed","twodash","solid","longdash","solid","solid","solid")

ggplot(eigens2, aes(x = trt, y=meanfit, ymax = meanfit+sem, ymin = meanfit-sem, col = trt, fill = trt)) + 
  geom_col(size=2.5, col=plot_colors, fill = "white", linetype = plot_lines) + 
  geom_errorbar(width = 0.5, linetype = "solid") + 
  theme_cowplot() +
  annotate(geom="text", x=c(1:length(table(list(eigens2$trt)))), y= c(eigens2$meanfit + eigens2$sem+median(eigens2$sem)), label = ghi$Letter, color = plot_colors, size=6, parse=T, hjust = 0, angle=90) + 
  ylab("mean fitness") +
  xlab("treatment") +
  coord_cartesian(ylim = c(0,2.7)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        text = element_text(size = 12)
  )
Run Figure 1D statistics
q = 1
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
rangevals <- c(1)
m <- rangevals
n = 1
eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric())
eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
mod = n
age_class = m

subfitness <- fitness %>% filter(vial %in% subset_vec)

for(i in names(table(list(subfitness$trt)))) {
  submat <- fitness %>% 
    filter(trt==i) %>% 
    filter(date3 == min(date3)) %>% 
    mutate(date3 = 0) %>%
    mutate(fecundity = 0) %>%
    rbind(fitness %>% 
            filter(trt==i) %>% 
            arrange(date3)
    ) %>%
    rbind(fitness %>%
            filter(trt == i) %>%
            filter(date3 == max(date3)) %>% 
            mutate(date3 = max(date3)+1) %>%
            mutate(fecundity = 0)
    )
  submat$s = 1
  submat$f = 0
  
  for (j in 1:(dim(submat)[1]-1)) {
    submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
    submat$f[j] <- submat$fecundity[j] #* submat$s[j]
  }
  
  leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
  if (age_class>1) {
    for (k in 1:(age_class-1)) {
      leslie_mat <- rbind(leslie_mat,
                          matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
      )
    } 
  }
  
  for (k in age_class:min(age_class,(dim(submat)[1]-1))) {
    leslie_mat <- rbind(leslie_mat, 
                        matrix(c(rep(0,k-1),submat$s[k]*mod,rep(0,(dim(submat)[1]-k))), nrow=1)
    )
  }
  if(age_class < (dim(submat)[1]-1)) {
    for(k in (age_class+1):(dim(submat)[1]-1)) {
      leslie_mat <- rbind(leslie_mat,
                          matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
      )
    }
  }
  
  eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
  
}

eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(fitness = log(fitness)) %>% mutate(mod_death = (m)))

e3 <- eigens
eigens <- eigens %>% filter(trt%in%subset_vec) %>% droplevels()
eigens$trt <- factor(eigens$trt, levels = subset_vec)

efg <- dunn.test(eigens$fitness,eigens$trt, method = "bh", table = F)
cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05, print.comp = F) %>%
  arrange(factor(Group, levels = c("35u","3p","4d","11c","28n","1a","71k","2b","75o","3c","7j","15h","16i","13g","Gno","Ax"))) %>%
  inner_join(lookup_code, by = c("Group" = "code3")) %>% select(fullname, code1, Letter)

#Figure 1 legend
## take the proper subset of data
atrc <- fecundity %>% 
  filter(vial%in%c("35u","30p","4d","11c","28n","2b","71k","1a","75o","3c","70j","15h","16i","13g","Gno","Ax")) %>% 
  filter(!is.na(fecundity)) %>% 
  filter(fecundity!="NaN") %>% 
  droplevels()

## adjust levels to maintain plot order
atrc$vial <- factor(atrc$vial, levels = c("35u","30p","4d","11c","28n","1a","71k","2b","75o","3c","70j","15h","16i","13g","Gno","Ax"))

## build vectors for plot shading
plot_colors <- c("red","red","red","red","red","blue","blue","blue","blue","blue","blue","green","green","purple","magenta","black")
plot_lines <- c("solid","longdash","dotted","dotdash","dashed","solid","longdash","dotted","dotdash","dashed","twodash","solid","longdash","solid","solid","solid")

plot(g_legend(ggplot(atrc, aes(x = date3, y=fecundity, group = vial, col = vial, fill = vial, linetype = vial)) +
                geom_smooth(aes(y = fecundity, linetype = vial), alpha = 0.07) +
                scale_color_manual(values = plot_colors, name = NULL, labels = c(expression(italic(A.~aceti)),expression(italic(A.~pasteurianus)),expression(italic(A.~pomorum)~DmCS),expression(italic(A.~tropicalis)~DmCS),expression(italic(A.~tropicalis)~NBRC),expression(italic(L.~brevis)~DmCS),expression(italic(L.~brevis)~gravesensis),expression(italic(L.~plantarum)~DmCS),expression(italic(L.~plantarum)~WCFS1),expression(italic(L.~fructivorans)~DmCS),expression(italic(L.~fructivorans)~KCTC),expression(italic(E.~coli)),expression(italic(P.~putida)),expression(italic(B.~subtilis)),'5-species','Axenic')) +
                scale_fill_manual(values = plot_colors, name = NULL,labels = c(expression(italic(A.~aceti)),expression(italic(A.~pasteurianus)),expression(italic(A.~pomorum)~DmCS),expression(italic(A.~tropicalis)~DmCS),expression(italic(A.~tropicalis)~NBRC),expression(italic(L.~brevis)~DmCS),expression(italic(L.~brevis)~gravesensis),expression(italic(L.~plantarum)~DmCS),expression(italic(L.~plantarum)~WCFS1),expression(italic(L.~fructivorans)~DmCS),expression(italic(L.~fructivorans)~KCTC),expression(italic(E.~coli)),expression(italic(P.~putida)),expression(italic(B.~subtilis)),'5-species','Axenic')) +
                scale_linetype_manual(values = plot_lines,name = NULL,labels = c(expression(italic(A.~aceti)),expression(italic(A.~pasteurianus)),expression(italic(A.~pomorum)~DmCS),expression(italic(A.~tropicalis)~DmCS),expression(italic(A.~tropicalis)~NBRC),expression(italic(L.~brevis)~DmCS),expression(italic(L.~brevis)~gravesensis),expression(italic(L.~plantarum)~DmCS),expression(italic(L.~plantarum)~WCFS1),expression(italic(L.~fructivorans)~DmCS),expression(italic(L.~fructivorans)~KCTC),expression(italic(E.~coli)),expression(italic(P.~putida)),expression(italic(B.~subtilis)),'5-species','Axenic')) +
                # coord_cartesian(ylim = c(0,0.75)) +
                ylab("female fecundity per day") +
                xlab("days") +
                guides(fill = guide_legend(ncol = 1, byrow=T), color = guide_legend(ncol = 1, byrow=T), linetype = guide_legend(ncol = 1, byrow=T)) +
                theme_cowplot() +
                theme(legend.position = "right", legend.text.align = 0, legend.text = element_text(size = 10), axis.text = element_text(size = 12),axis.title = element_text(size = 14),legend.justification = "left", legend.key.width = unit(x = .4, units = "in"),legend.background = element_rect(fill = "white"), strip.background = element_rect(fill = "white"))
)
)
Figure 2
Figure 2AC
p <- list(vec = list(vec_bac = c("28n","30p","11c","Gno","4d","35u","Ax","2b","13g","15h","75o","3c","1a","70j","71k"),
                     vec_atmut = c("SDR","PVEC","GDH","GNDH"),
                     vec_ecmut = c("7636","10862")
),
group = list("group_bac","group_atmut","group_ecmut")
)

q = 1
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric(), permuted_val = character())

rangevals <- c(1, 0.00001)

for (n in list(c(1,1))) {
  for (m in rangevals) {
    eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
    distance = n
    mod = m
    subfitness <- fitness %>% filter(vial %in% subset_vec)
    i = "11c-11"
    for(i in names(table(list(subfitness$trt)))) {
      
      ## collect the relevant data for each Leslie matrix 
      submat <- 
        # make first row
        fitness %>%  filter(trt==i) %>% filter(date3 == min(date3)) %>% mutate(date3 = 0) %>% mutate(fecundity = 0) %>% 
        # the main matrix
        rbind(fitness %>% filter(trt==i) %>% arrange(date3)) %>%
        # make last row
        rbind(fitness %>% filter(trt == i) %>% filter(date3 == max(date3)) %>%  mutate(date3 = max(date3)+1) %>% mutate(fecundity = 0))
      
      ## set the s and f coloumns to default values
      submat$s = 1
      submat$f = 0
      
      ## calculate s and f
      for (j in 1:(dim(submat)[1]-1)) {
        submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
        submat$f[j] <- submat$fecundity[j] #* submat$s[j]
      }
      
      ## begin the leslie matrix
      leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
      
      ## build in the first unmodified cells 
      if (n[1]>1) {
        for (k in 1:(n[1]-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        } 
      }
      
      ## build in cells that get modified by variable 'mod' from 'rangevals'
      for (k in n[1]:min(n[2],(dim(submat)[1]-1))) {
        leslie_mat <- rbind(leslie_mat, 
                            matrix(c(rep(0,k-1),submat$s[k]*mod,rep(0,(dim(submat)[1]-k))), nrow=1)
        )
      }
      
      ## build in the unmodified cells that occur linearly after the modified ones
      if(n[2] < (dim(submat)[1]-1)) {
        for(k in (n[2]+1):(dim(submat)[1]-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        }
      }
      
      ## the matrix
      leslie_mat
      
      ## calculate the first eigenvalue in the eigenvector of the Leslie matrix
      eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
      
    }
    
    ## add that row with the eigenvalue to a growing data frame
    eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(mod_death = (m), permuted_val = paste0(n[1],"_",n[2])))
  }
}

for(k in names(table(list(eigens$mod_death)))) {
  
  ## get rid of unnecessary rows  
  eigensA <- eigens %>% filter(trt%in%subset_vec) %>% filter(permuted_val == "1_1") %>%droplevels()
  
  ## convert trt to a factor
  eigensA$trt <- factor(eigensA$trt, levels = subset_vec)
  
  eigens2a <- eigensA %>% filter(mod_death == 1) %>% dplyr::rename(unpermuted_fitness = fitness) %>% droplevels()
  eigens2b <- eigensA %>% filter(mod_death == k) %>% dplyr::rename(permuted_fitness = fitness) %>% droplevels()
  
  eigens3 <- eigens2a %>% inner_join(eigens2b, by="vial")
  
  cortest <- cor.test(eigens3$unpermuted_fitness, eigens3$permuted_fitness, method = "spearman", exact=F)
  print(paste0("Correlation test between observed fitness and when fitness is calculated from ",1/as.numeric(as.character(k)),"-fold reduction in fly survival values. N = ",dim(eigens3 %>% filter(unpermuted_fitness > 0 & permuted_fitness>0))[1]))
  print(cortest)
  
  colors <- data.frame(table(list(eigens3$trt.x))) %>% mutate(col = as.factor(c("red","red","red","black","purple","red","red","green",rep("blue",4),"cyan",rep("blue",2))))
  eigens4 <- eigens3 %>% inner_join(colors, by=c("trt.x"="Var1"))
  plotp = round(cortest$p.value,4)
  plotrho = round(cortest$estimate,4)
  plotlabel = paste0("p = ",plotp,"; rho = ",plotrho)
  plot_colors <- eigens4$col
  print(ggplot(eigens4, aes(unpermuted_fitness, permuted_fitness, color = col, fill=col)) + 
          geom_point() +
          scale_color_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          scale_fill_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          theme_cowplot() +
          theme(legend.position = "bottom") +
          xlab(paste0("observed fitness"))+
          ylab(paste0("survival-permuted fitness"))+
          annotate(geom="label", x= min(eigens4$unpermuted_fitness), y = min(eigens4$permuted_fitness), label = deparse(plotlabel),size=6, parse=T, hjust = 0, fill = NA,label.size=0) + 
          stat_smooth(method = "lm", aes(unpermuted_fitness, permuted_fitness), inherit.aes = F)
        
  )
  
}
Figure 2BD
p <- list(vec = list(vec_bac = c("28n","30p","11c","Gno","4d","35u","Ax","2b","13g","15h","75o","3c","1a","70j","71k"),
                     vec_atmut = c("SDR","PVEC","GDH","GNDH"),
                     vec_ecmut = c("7636","10862")
),
group = list("group_bac","group_atmut","group_ecmut")
)

q = 1
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric(), permuted_val = character())

rangevals <- c(1,0.00001)

for (n in list(c(1,1))) {
  for (m in rangevals) {
    eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
    distance = n
    mod = m
    subfitness <- fitness %>% filter(vial %in% subset_vec)
    for(i in names(table(list(subfitness$trt)))) {
      
      ## collect the relevant data for each Leslie matrix 
      submat <- 
        # make first row
        fitness %>%  filter(trt==i) %>% filter(date3 == min(date3)) %>% mutate(date3 = 0) %>% mutate(fecundity = 0) %>% 
        # the main matrix
        rbind(fitness %>% filter(trt==i) %>% arrange(date3)) %>%
        # make last row
        rbind(fitness %>% filter(trt == i) %>% filter(date3 == max(date3)) %>%  mutate(date3 = max(date3)+1) %>% mutate(fecundity = 0))
      
      ## set the s and f coloumns to default values
      submat$s = 1
      submat$f = 0
      
      ## calculate s and f
      for (j in 1:(dim(submat)[1]-1)) {
        submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
        submat$f[j] <- submat$fecundity[j] #* submat$s[j]
      }
      
      ## begin the leslie matrix
      leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
      leslie_mat[(n[1]+1):(n[2]+1)] <- (leslie_mat[(n[1]+1):(n[2]+1)] * mod)
      
      ## build in the diagonal
      for (k in 1:(dim(submat)[1]-1)) {
        leslie_mat <- rbind(leslie_mat, 
                            matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
        )
      }
      
      ## the matrix
      
      
      ## calculate the first eigenvalue in the eigenvector of the Leslie matrix
      eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
      
    }
    
    ## add that row with the eigenvalue to a growing data frame
    eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(mod_death = (m), permuted_val = paste0(n[1],"_",n[2])))
  }
}

for(k in names(table(list(eigens$mod_death)))) {
  
  ## get rid of unnecessary rows  
  eigensA <- eigens %>% filter(trt%in%subset_vec) %>% filter(permuted_val == "1_1") %>%droplevels()
  
  ## convert trt to a factor
  eigensA$trt <- factor(eigensA$trt, levels = subset_vec)
  
  eigens2a <- eigensA %>% filter(mod_death == 1) %>% rename(unpermuted_fitness = fitness) %>% droplevels()
  eigens2b <- eigensA %>% filter(mod_death == k) %>% rename(permuted_fitness = fitness) %>% droplevels()
  
  table(eigens2a$vial)
  eigens3 <- eigens2a %>% inner_join(eigens2b, by="vial")
  eigens3$unpermuted_fitness
  
  cortest <- cor.test(eigens3$unpermuted_fitness, eigens3$permuted_fitness, method = "spearman", exact=F)
  print(paste0("Correlation test between observed fitness and when fitness is calculated from ",1/as.numeric(as.character(k)),"-fold reduction in fecundity values. N = ",dim(eigens3 %>% filter(unpermuted_fitness > 0 & permuted_fitness>0))[1]))
  print(cortest)
  
  colors <- data.frame(table(list(eigens3$trt.x))) %>% mutate(col = as.factor(c("red","red","red","black","purple","red","red","green",rep("blue",4),"cyan",rep("blue",2))))
  eigens4 <- eigens3 %>% inner_join(colors, by=c("trt.x"="Var1"))
  
  plotp = round(cortest$p.value,4)
  plotrho = round(cortest$estimate,4)
  plotlabel = paste0("p = ",plotp,"; rho = ",plotrho)
  ## build the color palette for the plot
  plot_colors <- eigens4$col
  
  print(ggplot(eigens4, aes(unpermuted_fitness, permuted_fitness, color = col, fill=col)) + 
          geom_point() +
          scale_color_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          scale_fill_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          theme_cowplot() +
          theme(legend.position = "bottom") +
          xlab(paste0("observed fitness"))+
          ylab(paste0("fecundity-permuted fitness"))+
          annotate(geom="label", x= min(eigens4$unpermuted_fitness), y = min(eigens4$permuted_fitness), label = deparse(plotlabel),size=6, parse=T, hjust = 0, fill = NA,label.size=0) + 
          stat_smooth(method = "lm", aes(unpermuted_fitness, permuted_fitness), inherit.aes = F)
  )
  
}
Figure 3
Figure 3A Acetobacter
## filter out the right data
atrc <- fecundity%>% filter(vial%in%c("GDH","GNDH","PVEC","SDR"))

## set the plot order
atrc$vial <- factor(atrc$vial, levels = c("PVEC","SDR","GDH","GNDH"))

## set colors
plot_colors <- c("red","red","red","red")
plot_lines <- c("solid","dashed","dotted","dotdash")

## make plot
ggplot(atrc, aes(x = date3, y=fecundity, group = vial, col = vial, fill = vial, linetype = vial)) + 
  stat_smooth(method = NULL, level= 0.95, alpha = 0.07) +
  scale_color_manual(values = plot_colors, 
                     labels = c("Empty vector","S-oxidoreductase","Glucose dehydrogenase","Gluconate dehydrogenase"), 
                     name = "") + 
  scale_fill_manual(values = plot_colors, 
                    labels = c("Empty vector","S-oxidoreductase","Glucose dehydrogenase","Gluconate dehydrogenase"), 
                    name = "")+
  scale_linetype_manual(values = plot_lines, 
                        labels = c("Empty vector","S-oxidoreductase","Glucose dehydrogenase","Gluconate dehydrogenase"), 
                        name = "") + 
  ylab("female fecundity per day") +
  xlab("days") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 3, byrow=T),
         color = guide_legend(ncol = 3, byrow=T),
         linetype = guide_legend(ncol = 3, byrow=T)
  ) +
  theme_cowplot() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        legend.key.width = unit(.5,"in"), 
        legend.justification = "center",
        text = element_text(size = 12)
  )
Figure 3A E. coli
## filter out the right data
atrc <- fecundity%>% filter(vial%in%c("7636", "10862"))

## define plot order
atrc$vial <- factor(atrc$vial, levels = c("7636", "10862"))

## set colors
plot_colors <- c("green","green")
plot_lines <- c("solid","dashed")

## plot
ggplot(atrc, aes(x = date3, y=fecundity, group = vial, col = vial, fill = vial, linetype = vial)) + 
  stat_smooth(method = NULL, level= 0.95, alpha = 0.07) +
  scale_color_manual(values = plot_colors, 
                     labels = c("WT",expression(metH)), 
                     name = "") + 
  scale_fill_manual(values = plot_colors, 
                    labels = c("WT",expression(metH)), 
                    name = "")+
  scale_linetype_manual(values = plot_lines, 
                        labels = c("WT",expression(metH)), 
                        name = "") + 
  ylab("female fecundity per day") +
  xlab("days") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 3, byrow=T),
         color = guide_legend(ncol = 3, byrow=T),
         linetype = guide_legend(ncol = 3, byrow=T)
  ) +
  theme_cowplot() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        legend.key.width = unit(.5,"in"), 
        legend.justification = "center",
        text = element_text(size = 12)
  )
Figure S1
## this uses a script in the middle
## modifies list p to exclude Pseudomonas putida
p <- list(vec = list(vec_bac = c("28n","30p","11c","Gno","4d","35u","Ax","2b","13g","15h","75o","3c","1a","70j","71k"),
                     vec_atmut = c("SDR","PVEC","GDH","GNDH"),
                     vec_ecmut = c("7636","10862")
),
group = list("group_bac","group_atmut","group_ecmut")
)

q = 1

eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric())

rangevals <- c(1:7)
n = 0.5
subset_vec = unlist(unname(p$vec[q]))
for (n in c(.5,.1,.05,.01,.005,.001,.0001,.00001)) {
  
  ## calculate the eigenvalues  
  eigen_grow <- calculate_eigen(m = rangevals, 
                                mod = n, 
                                age_class = rangevals, 
                                fitness = fitness, 
                                subset_vec = subset_vec,
                                xgroup = unlist(p$group[q])
  ) 
  
  eigens <- eigen_grow %>% 
    filter(trt%in%unlist(unname(p$vec[q]))) %>% 
    droplevels() 
  
  ## set the factor levels and colors and colors
  eigens$trt <- factor(eigens$trt, levels = subset_vec)
  eigens[1,]
  eigens2 <- eigens %>% group_by(trt, age_class_mod) %>% summarize(meanfit = mean(fitness))
  plot_colors <- colorRampPalette(c("red", "black"))
  
  
  i=1
  table(eigens$age_class_mod)
  for(i in rangevals) {
    cat('\n')
    print("######################################################")
    print(paste0("These are the statistics for age class ",i," and survival divided by ",1/n, sep = ""))
    
    
    eigen_stats <- eigens %>% filter(age_class_mod == (i)) %>% droplevels()
    
    kwoutput <- try(kruskal.test(eigen_stats$fitness~ eigen_stats$trt), silent = T)
    captured_output <- capture.output(efg <- dunn.test(eigen_stats$fitness,eigen_stats$trt, method = "bh", table = F, kw = T))
    if(kwoutput == "Error in kruskal.test.default(numeric(0), structure(integer(0), .Label = character(0), class = \"factor\")) : \n  all observations are in the same group\n") {
      print(paste0("chi-sq statistic not significant"))
    } else {
      print(paste0("chi-sq statistic: ",round(kwoutput$statistic,2)))
      print(paste0("df: ",length(table(eigen_stats$trt))-1))
      print(paste0("p: ",signif(kwoutput$p.value,3)))      
    }      
    #print(i)
    if(sum(efg$P.adjusted<0.05)>0) {
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05, print.comp = F) %>% 
        arrange(factor(Group, levels = subset_vec))
      try(assign(x = paste0("es",i), value = as.vector(ghi$Letter)),T)
    } else {
      assign(x = paste0("es",i),value = rep("a",length(subset_vec)))
    }
    assign(x = paste0("estats",i), 
           value = eigen_stats %>% 
             group_by(trt) %>% 
             summarize(mean_fitness = mean(fitness)) %>% 
             arrange(factor(trt, levels = subset_vec)) %>%
             dplyr::select(mean_fitness) %>%
             unlist() %>% 
             unname()
    )
  }
  
  #    eigens2$trt <- plyr::revalue(eigens2$trt, c("28n" = "atrn","30p" = "apan", "11c"="atrc","Gno" = "5sp","4d" = "apoc", "35u" = "aace","Ax" = "Ax","2b" = "lplc","13g" = "bsub","15h" = "ecok","75o" = "lplw","3c" = "lfrc", "1a" = "lbrc","70j" = "lfrk", "71k" = "lbga"))
  
  eigens2
  
  print(ggplot(eigens2, aes(x = trt, y=meanfit,group = as.factor(age_class_mod), col = as.factor(age_class_mod), fill = as.factor(age_class_mod))) + 
          geom_step(aes(x = trt, y = meanfit)) + #, position = position_dodge(width = widthj), size = sizej) +
          scale_color_manual(values = plot_colors(length(rangevals)),
                             name = "") +
          coord_cartesian(ylim = c(-2,3)) +
          ylab("log fitness") +
          xlab("bacterial treatment") +
          theme(legend.position = "bottom") +
          guides(color = guide_legend(ncol = 5, byrow=T, reverse = T)) +
          theme_cowplot() + 
          theme(legend.position = "bottom", 
                legend.text.align = 0,
                legend.key.width = unit(.5,"in"), 
                legend.text = element_text(size = 18), 
                axis.text = element_text(size = 18),
                axis.title = element_text(size = 24),
                legend.justification = "center", 
                text = element_text(size =14),
                axis.text.x = element_text(angle = 90, hjust = 1)
          ) + 
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[7]))+0.03, label = get(paste0("es",rangevals[7])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[6]))+0.03, label = get(paste0("es",rangevals[6])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[5]))+0.03, label = get(paste0("es",rangevals[5])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[4]))+0.03, label = get(paste0("es",rangevals[4])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[3]))+0.03, label = get(paste0("es",rangevals[3])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[2]))+0.03, label = get(paste0("es",rangevals[2])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[1]))+0.03, label = get(paste0("es",rangevals[1])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)
  )
  
}
p <- list(vec = list(vec_bac = c("28n","30p","11c","Gno","4d","35u","Ax","2b","13g","15h","75o","3c","1a","70j","71k"),
                     vec_atmut = c("SDR","PVEC","GDH","GNDH"),
                     vec_ecmut = c("7636","10862")
),
group = list("group_bac","group_atmut","group_ecmut")
)

q = 1

#print 1 plot with  the same survival in each line. each line is the age class that was permuted for survival
xgroup = unlist(p$group[q])
subsetp_vec = unlist(unname(p$vec[q]))
rangevals = 1
rangevals <- c(1:7)
n= 0.00001
for (n in c(.5,.1,.05,.01,.005,.001,.0001,.00001)) {
  eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric())
  m = 1
  for (m in rangevals) {
    eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
    mod = n
    age_class = m
    subfitness <- fitness %>% filter(vial %in% subset_vec)
    for(i in names(table(list(subfitness$trt)))) {
      submat <- fitness %>% 
        filter(trt==i) %>% 
        filter(date3 == min(date3)) %>% 
        mutate(date3 = 0) %>%
        mutate(fecundity = 0) %>%
        rbind(fitness %>% 
                filter(trt==i) %>% 
                arrange(date3)
        ) %>%
        rbind(fitness %>%
                filter(trt == i) %>%
                filter(date3 == max(date3)) %>% 
                mutate(date3 = max(date3)+1) %>%
                mutate(fecundity = 0)
        )
      submat$s = 1
      submat$f = 0
      
      for (j in 1:(dim(submat)[1]-1)) {
        submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
        submat$f[j] <- submat$fecundity[j] #* submat$s[j]
      }
      submat
      leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
      if (age_class>1) {
        for (k in 1:(age_class-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        } 
      }
      for (k in age_class:min(age_class,(dim(submat)[1]-1))) {
        leslie_mat <- rbind(leslie_mat, 
                            matrix(c(rep(0,k-1),submat$s[k]*mod,rep(0,(dim(submat)[1]-k))), nrow=1)
        )
      }
      if(age_class < (dim(submat)[1]-1)) {
        for(k in (age_class+1):(dim(submat)[1]-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        }
      }
      
      eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
      
    }
    
    eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(fitness = log(fitness)) %>% mutate(mod_death = (m)))
  }
  
  eigens <- eigens %>% filter(trt%in%subset_vec) %>% droplevels()
  eigens$trt <- factor(eigens$trt, levels = subset_vec)
  eigens2 <- eigens %>% group_by(trt, mod_death) %>% summarize(meanfit = mean(fitness))
  plot_colors <- colorRampPalette(c("red", "black"))
  
  table(list(eigens$trt))
  
  table(list(eigens$mod_death))
  for(i in rangevals) {
    cat('\n')
    print("######################################################")
    print(paste0("These are the statistics for age class ",i," and survival divided by ",1/n, sep = ""))
    
    
    eigen_stats <- eigens %>% filter(mod_death == (i)) %>% droplevels()
    #print(str(eigen_stats))
    #print(i)
    kwoutput <- try(kruskal.test(eigen_stats$fitness~ eigen_stats$trt), silent = T)
    captured_output <- capture.output(efg <- dunn.test(eigen_stats$fitness,eigen_stats$trt, method = "bh", table = F, kw = T))
    if(kwoutput == "Error in kruskal.test.default(numeric(0), structure(integer(0), .Label = character(0), class = \"factor\")) : \n  all observations are in the same group\n") {
      print(paste0("chi-sq statistic not significant"))
    } else {
      print(paste0("chi-sq statistic: ",round(kwoutput$statistic,2)))
      print(paste0("df: ",length(table(eigen_stats$trt))-1))
      print(paste0("p: ",signif(kwoutput$p.value,3)))      
    }      
    #print(i)
    if(sum(efg$P.adjusted<0.05)>0) {
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05, print.comp = F) %>% 
        arrange(factor(Group, levels = subset_vec))
      try(assign(x = paste0("es",i), value = as.vector(ghi$Letter)),T)
    } else {
      assign(x = paste0("es",i),value = rep("a",length(subset_vec)))
    }
    assign(x = paste0("estats",i), 
           value = eigen_stats %>% 
             group_by(trt) %>% 
             summarize(mean_fitness = mean(fitness)) %>% 
             arrange(factor(trt, levels = subset_vec)) %>%
             dplyr::select(mean_fitness) %>%
             unlist() %>% 
             unname()
    )
  }
  
  #    eigens2$trt <- plyr::revalue(eigens2$trt, c("28n" = "atrn","30p" = "apan", "11c"="atrc","Gno" = "5sp","4d" = "apoc", "35u" = "aace","Ax" = "Ax","2b" = "lplc","13g" = "bsub","15h" = "ecok","75o" = "lplw","3c" = "lfrc", "1a" = "lbrc","70j" = "lfrk", "71k" = "lbga"))
  
  print(ggplot(eigens2, aes(x = trt, y=meanfit,group = as.factor(mod_death), col = as.factor(mod_death), fill = as.factor(mod_death))) + 
          geom_step(aes(x = trt, y = meanfit)) + #, position = position_dodge(width = widthj), size = sizej) +
          scale_color_manual(values = plot_colors(length(rangevals)),
                             name = "") +
          coord_cartesian(ylim = c(-2,3)) +
          ylab("fitness") +
          xlab("bacterial treatment") +
          theme(legend.position = "bottom") +
          guides(color = guide_legend(ncol = 5, byrow=T, reverse = T)) +
          theme_cowplot() + 
          theme(legend.position = "bottom", 
                legend.text.align = 0,
                legend.key.width = unit(.5,"in"), 
                legend.text = element_text(size = 18), 
                axis.text = element_text(size = 18),
                axis.title = element_text(size = 24),
                legend.justification = "center", 
                text = element_text(size =14),
                axis.text.x = element_text(angle = 90, hjust = 1)
          ) + 
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[7]))+0.03, label = get(paste0("es",rangevals[7])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[6]))+0.03, label = get(paste0("es",rangevals[6])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[5]))+0.03, label = get(paste0("es",rangevals[5])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[4]))+0.03, label = get(paste0("es",rangevals[4])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[3]))+0.03, label = get(paste0("es",rangevals[3])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[2]))+0.03, label = get(paste0("es",rangevals[2])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[1]))+0.03, label = get(paste0("es",rangevals[1])),size=6, parse=T, hjust = 0, fill = NA,label.size=0)
  )
  # dev.off()
}
p <- list(vec = list(vec_bac = c("28n","30p","11c","Gno","4d","35u","Ax","2b","13g","15h","75o","3c","1a","70j","71k"),
                     vec_atmut = c("SDR","PVEC","GDH","GNDH"),
                     vec_ecmut = c("7636","10862")
),
group = list("group_bac","group_atmut","group_ecmut")
)

q = 1

#print 1 plot with  the same survival in each line. each line is the age class that was permuted for survival
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
rangevals = 1
rangevals <- c(1:7)
n= 0.00001
for (n in c(.5,.1,.05,.01,.005,.001,.0001,.00001)) {
  eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric())
  m = 1
  for (m in rangevals) {
    eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
    mod = n
    age_class = m
    subfitness <- fitness %>% filter(vial %in% subset_vec)
    for(i in names(table(list(subfitness$trt)))) {
      submat <- fitness %>% 
        filter(trt==i) %>% 
        filter(date3 == min(date3)) %>% 
        mutate(date3 = 0) %>%
        mutate(fecundity = 0) %>%
        rbind(fitness %>% 
                filter(trt==i) %>% 
                arrange(date3)
        ) %>%
        rbind(fitness %>%
                filter(trt == i) %>%
                filter(date3 == max(date3)) %>% 
                mutate(date3 = max(date3)+1) %>%
                mutate(fecundity = 0)
        )
      submat$s = 1
      submat$f = 0
      
      for (j in 1:(dim(submat)[1]-1)) {
        submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
        submat$f[j] <- submat$fecundity[j] #* submat$s[j]
      }
      submat
      leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
      if (age_class>1) {
        for (k in 1:(age_class-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        } 
      }
      for (k in age_class:min(age_class,(dim(submat)[1]-1))) {
        leslie_mat <- rbind(leslie_mat, 
                            matrix(c(rep(0,k-1),submat$s[k]*mod,rep(0,(dim(submat)[1]-k))), nrow=1)
        )
      }
      if(age_class < (dim(submat)[1]-1)) {
        for(k in (age_class+1):(dim(submat)[1]-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        }
      }
      
      eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
      
    }
    
    eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(fitness = log(fitness)) %>% mutate(mod_death = (m)))
  }
  
  eigens <- eigens %>% filter(trt%in%subset_vec) %>% droplevels()
  eigens$trt <- factor(eigens$trt, levels = subset_vec)
  eigens2 <- eigens %>% group_by(trt, mod_death) %>% summarize(meanfit = mean(fitness))
  plot_colors <- colorRampPalette(c("red", "black"))
  
  table(list(eigens$trt))
  
  table(list(eigens$mod_death))
  for(i in rangevals) {
    cat('\n')
    print("######################################################")
    print(paste0("These are the statistics for age class ",i," and survival divided by ",1/n, sep = ""))
    
    
    eigen_stats <- eigens %>% filter(mod_death == (i)) %>% droplevels()
    #print(str(eigen_stats))
    #print(i)
    kwoutput <- try(kruskal.test(eigen_stats$fitness~ eigen_stats$trt), silent = T)
    captured_output <- capture.output(efg <- dunn.test(eigen_stats$fitness,eigen_stats$trt, method = "bh", table = F, kw = T))
    if(kwoutput == "Error in kruskal.test.default(numeric(0), structure(integer(0), .Label = character(0), class = \"factor\")) : \n  all observations are in the same group\n") {
      print(paste0("chi-sq statistic not significant"))
    } else {
      print(paste0("chi-sq statistic: ",round(kwoutput$statistic,2)))
      print(paste0("df: ",length(table(eigen_stats$trt))-1))
      print(paste0("p: ",signif(kwoutput$p.value,3)))      
    }      
    #print(i)
    if(sum(efg$P.adjusted<0.05)>0) {
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05, print.comp = F) %>% 
        arrange(factor(Group, levels = subset_vec))
      try(assign(x = paste0("es",i), value = as.vector(ghi$Letter)),T)
    } else {
      assign(x = paste0("es",i),value = rep("a",length(subset_vec)))
    }
    assign(x = paste0("estats",i), 
           value = eigen_stats %>% 
             group_by(trt) %>% 
             summarize(mean_fitness = mean(fitness)) %>% 
             arrange(factor(trt, levels = subset_vec)) %>%
             dplyr::select(mean_fitness) %>%
             unlist() %>% 
             unname()
    )
  }
  
  eigens2$trt <- plyr::revalue(eigens2$trt, c("28n" = "atrn","30p" = "apan", "11c"="atrc","Gno" = "5sp","4d" = "apoc", "35u" = "aace","Ax" = "Ax","2b" = "lplc","13g" = "bsub","15h" = "ecok","75o" = "lplw","3c" = "lfrc", "1a" = "lbrc","70j" = "lfrk", "71k" = "lbga"))
  
  print(ggplot(eigens2, aes(x = trt, y=meanfit,group = as.factor(mod_death), col = as.factor(mod_death), fill = as.factor(mod_death))) + 
          geom_step(aes(x = trt, y = meanfit)) + #, position = position_dodge(width = widthj), size = sizej) +
          scale_color_manual(values = plot_colors(length(rangevals)),
                             name = "Age class") +
          coord_cartesian(ylim = c(-2,3)) +
          ylab("log fitness") +
          xlab("bacterial treatment") +
          theme(legend.position = "bottom") +
          guides(color = guide_legend(ncol = 5, byrow=T, reverse = T)) +
          theme_cowplot() + 
          theme(legend.position = "bottom", 
                legend.text.align = 0,
                legend.key.width = unit(.5,"in"), 
                legend.justification = "center", 
                text = element_text(size =12),
                axis.text.x = element_text(angle = 90, hjust = 1)
          ) + 
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[7]))+0.03, label = get(paste0("es",rangevals[7])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[6]))+0.03, label = get(paste0("es",rangevals[6])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[5]))+0.03, label = get(paste0("es",rangevals[5])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[4]))+0.03, label = get(paste0("es",rangevals[4])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[3]))+0.03, label = get(paste0("es",rangevals[3])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[2]))+0.03, label = get(paste0("es",rangevals[2])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[1]))+0.03, label = get(paste0("es",rangevals[1])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)
  )
  # dev.off()
}
Figure S2
p <- list(vec = list(vec_bac = c("28n","30p","11c","Gno","4d","35u","Ax","2b","13g","15h","75o","3c","1a","70j","71k"),
                     vec_atmut = c("SDR","PVEC","GDH","GNDH"),
                     vec_ecmut = c("7636","10862")
),
group = list("group_bac","group_atmut","group_ecmut")
)

q = 1
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric(), permuted_val = character())

rangevals <- c(1, 0.5,0.1,0.05, 0.01, 0.005,0.001,0.0001,0.00001)

for (n in list(c(1,1),c(2,2))) {
  for (m in rangevals) {
    eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
    distance = n
    mod = m
    subfitness <- fitness %>% filter(vial %in% subset_vec)
    i = "11c-11"
    for(i in names(table(list(subfitness$trt)))) {
      
      ## collect the relevant data for each Leslie matrix 
      submat <- 
        # make first row
        fitness %>%  filter(trt==i) %>% filter(date3 == min(date3)) %>% mutate(date3 = 0) %>% mutate(fecundity = 0) %>% 
        # the main matrix
        rbind(fitness %>% filter(trt==i) %>% arrange(date3)) %>%
        # make last row
        rbind(fitness %>% filter(trt == i) %>% filter(date3 == max(date3)) %>%  mutate(date3 = max(date3)+1) %>% mutate(fecundity = 0))
      
      ## set the s and f coloumns to default values
      submat$s = 1
      submat$f = 0
      
      ## calculate s and f
      for (j in 1:(dim(submat)[1]-1)) {
        submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
        submat$f[j] <- submat$fecundity[j] #* submat$s[j]
      }
      
      ## begin the leslie matrix
      leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
      
      ## build in the first unmodified cells 
      if (n[1]>1) {
        for (k in 1:(n[1]-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        } 
      }
      
      ## build in cells that get modified by variable 'mod' from 'rangevals'
      for (k in n[1]:min(n[2],(dim(submat)[1]-1))) {
        leslie_mat <- rbind(leslie_mat, 
                            matrix(c(rep(0,k-1),submat$s[k]*mod,rep(0,(dim(submat)[1]-k))), nrow=1)
        )
      }
      
      ## build in the unmodified cells that occur linearly after the modified ones
      if(n[2] < (dim(submat)[1]-1)) {
        for(k in (n[2]+1):(dim(submat)[1]-1)) {
          leslie_mat <- rbind(leslie_mat,
                              matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
          )
        }
      }
      
      ## the matrix
      leslie_mat
      
      ## calculate the first eigenvalue in the eigenvector of the Leslie matrix
      eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
      
    }
    
    ## add that row with the eigenvalue to a growing data frame
    eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(mod_death = (m), permuted_val = paste0(n[1],"_",n[2])))
  }
}

for(k in names(table(list(eigens$mod_death)))) {
  
  ## get rid of unnecessary rows  
  eigensA <- eigens %>% filter(trt%in%subset_vec) %>% filter(permuted_val == "1_1") %>%droplevels()
  
  ## convert trt to a factor
  eigensA$trt <- factor(eigensA$trt, levels = subset_vec)
  
  eigens2a <- eigensA %>% filter(mod_death == 1) %>% dplyr::rename(unpermuted_fitness = fitness) %>% droplevels()
  eigens2b <- eigensA %>% filter(mod_death == k) %>% dplyr::rename(permuted_fitness = fitness) %>% droplevels()
  
  eigens3 <- eigens2a %>% inner_join(eigens2b, by="vial")
  
  cortest <- cor.test(eigens3$unpermuted_fitness, eigens3$permuted_fitness, method = "spearman", exact=F)
  print(paste0("Correlation test between observed fitness and when fitness is calculated from ",1/as.numeric(as.character(k)),"-fold reduction in fly survival values. N = ",dim(eigens3 %>% filter(unpermuted_fitness > 0 & permuted_fitness>0))[1]))
  print(cortest)
  
  colors <- data.frame(table(list(eigens3$trt.x))) %>% mutate(col = as.factor(c("red","red","red","black","purple","red","red","green",rep("blue",4),"cyan",rep("blue",2))))
  eigens4 <- eigens3 %>% inner_join(colors, by=c("trt.x"="Var1"))
  
  plotp = round(cortest$p.value,4)
  plotrho = round(cortest$estimate,4)
  plotlabel = paste0("p = ",plotp,"; rho = ",plotrho)
  plot_colors <- eigens4$col
  print(ggplot(eigens4, aes(unpermuted_fitness, permuted_fitness, color = col, fill=col)) + 
          geom_point() +
          scale_color_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          scale_fill_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          theme_cowplot() +
          theme(legend.position = "bottom") +
          xlab(paste0("observed fitness"))+
          ylab(paste0("survival-permuted fitness"))+
          annotate(geom="label", x= min(eigens4$unpermuted_fitness), y = min(eigens4$permuted_fitness), label = deparse(plotlabel),size=6, parse=T, hjust = 0, fill = NA,label.size=0) + 
          stat_smooth(method = "lm", aes(unpermuted_fitness, permuted_fitness), inherit.aes = F)
        
  )
  
}
Figure S3
p <- list(vec = list(vec_bac = c("28n","30p","11c","Gno","4d","35u","Ax","2b","13g","15h","75o","3c","1a","70j","71k"),
                     vec_atmut = c("SDR","PVEC","GDH","GNDH"),
                     vec_ecmut = c("7636","10862")
),
group = list("group_bac","group_atmut","group_ecmut")
)

q = 1
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric(), permuted_val = character())

rangevals <- c(1, .5,.1,.05, .01, 0.005,.001,.0005,0.0001,0.00001)

#  rangevals <- c(1,0.00001)

## make up the two comparisons

for (n in list(c(1,1),c(2,7))) {
  for (m in rangevals) {
    eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
    distance = n
    mod = m
    subfitness <- fitness %>% filter(vial %in% subset_vec)
    for(i in names(table(list(subfitness$trt)))) {
      
      ## collect the relevant data for each Leslie matrix 
      submat <- 
        # make first row
        fitness %>%  filter(trt==i) %>% filter(date3 == min(date3)) %>% mutate(date3 = 0) %>% mutate(fecundity = 0) %>% 
        # the main matrix
        rbind(fitness %>% filter(trt==i) %>% arrange(date3)) %>%
        # make last row
        rbind(fitness %>% filter(trt == i) %>% filter(date3 == max(date3)) %>%  mutate(date3 = max(date3)+1) %>% mutate(fecundity = 0))
      
      ## set the s and f coloumns to default values
      submat$s = 1
      submat$f = 0
      
      ## calculate s and f
      for (j in 1:(dim(submat)[1]-1)) {
        submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
        submat$f[j] <- submat$fecundity[j] #* submat$s[j]
      }
      
      ## begin the leslie matrix
      leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
      leslie_mat[(n[1]+1):(n[2]+1)] <- (leslie_mat[(n[1]+1):(n[2]+1)] * mod)
      
      ## build in the diagonal
      for (k in 1:(dim(submat)[1]-1)) {
        leslie_mat <- rbind(leslie_mat, 
                            matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
        )
      }
      
      ## the matrix
      
      
      ## calculate the first eigenvalue in the eigenvector of the Leslie matrix
      eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
      
    }
    
    ## add that row with the eigenvalue to a growing data frame
    eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(mod_death = (m), permuted_val = paste0(n[1],"_",n[2])))
  }
}

for(k in names(table(list(eigens$mod_death)))) {
  
  ## get rid of unnecessary rows  
  eigensA <- eigens %>% filter(trt%in%subset_vec) %>% filter(permuted_val == "1_1") %>%droplevels()
  
  ## convert trt to a factor
  eigensA$trt <- factor(eigensA$trt, levels = subset_vec)
  
  eigens2a <- eigensA %>% filter(mod_death == 1) %>% rename(unpermuted_fitness = fitness) %>% droplevels()
  eigens2b <- eigensA %>% filter(mod_death == k) %>% rename(permuted_fitness = fitness) %>% droplevels()
  
  table(eigens2a$vial)
  eigens3 <- eigens2a %>% inner_join(eigens2b, by="vial")
  eigens3$unpermuted_fitness
  
  cortest <- cor.test(eigens3$unpermuted_fitness, eigens3$permuted_fitness, method = "spearman", exact=F)
  print(paste0("Correlation test between observed fitness and when fitness is calculated from ",1/as.numeric(as.character(k)),"-fold reduction in fecundity values. N = ",dim(eigens3 %>% filter(unpermuted_fitness > 0 & permuted_fitness>0))[1]))
  print(cortest)
  
  colors <- data.frame(table(list(eigens3$trt.x))) %>% mutate(col = as.factor(c("red","red","red","black","purple","red","red","green",rep("blue",4),"cyan",rep("blue",2))))
  eigens4 <- eigens3 %>% inner_join(colors, by=c("trt.x"="Var1"))
  
  plotp = round(cortest$p.value,4)
  plotrho = round(cortest$estimate,4)
  plotlabel = paste0("p = ",plotp,"; rho = ",plotrho)
  ## build the color palette for the plot
  plot_colors <- eigens4$col
  
  print(ggplot(eigens4, aes(unpermuted_fitness, permuted_fitness, color = col, fill=col)) + 
          geom_point() +
          scale_color_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          scale_fill_manual(values = c("black","blue","magenta","cyan","purple","red1"), name = NULL, labels = c("Axenic","LAB","B. subtilis","Gamma","5-sp","AAB")) + 
          theme_cowplot() +
          theme(legend.position = "bottom") +
          xlab(paste0("observed fitness"))+
          ylab(paste0("fecundity-permuted fitness"))+
          annotate(geom="label", x= min(eigens4$unpermuted_fitness), y = min(eigens4$permuted_fitness), label = deparse(plotlabel),size=6, parse=T, hjust = 0, fill = NA,label.size=0) + 
          stat_smooth(method = "lm", aes(unpermuted_fitness, permuted_fitness), inherit.aes = F)
  )
  
}
Figure S4
q = 1

#print 1 plot with  the same survival in each line. each line is the age class that was permuted for survival
xgroup = unlist(p$group[q])
subset_vec = unlist(unname(p$vec[q]))
rangevals <- c(1:7)

for (n in c(1,.5,.1,.05,.01,.005,.001,.0001,.00001)) {
  eigens <- data.frame(vial = factor(), trt = character(), exp = factor(), fitness = numeric(), mod_death = numeric())
  m = 1
  for (m in rangevals) {
    eigenvalue_df <- data.frame(vial = character(), trt = character(), exp = factor(), fitness = numeric())
    mod = n
    age_class = m
    subfitness <- fitness %>% filter(vial %in% subset_vec)
    for(i in names(table(list(subfitness$trt)))) {
      submat <- fitness %>% 
        filter(trt==i) %>% 
        filter(date3 == min(date3)) %>% 
        mutate(date3 = 0) %>%
        mutate(fecundity = 0) %>%
        rbind(fitness %>% 
                filter(trt==i) %>% 
                arrange(date3)
        ) %>%
        rbind(fitness %>%
                filter(trt == i) %>%
                filter(date3 == max(date3)) %>% 
                mutate(date3 = max(date3)+1) %>%
                mutate(fecundity = 0)
        )
      submat$s = 1
      submat$f = 0
      
      
      for (j in 1:(dim(submat)[1]-1)) {
        submat$s[j] <- (submat$cum_female[j+1] / submat$cum_female[j])
        submat$f[j] <- submat$fecundity[j] #* submat$s[j]
      }
      
      leslie_mat <- matrix(data = c(unlist(submat$f[1:dim(submat)[1]])), nrow= 1)
      leslie_mat[(m+1)] <- (leslie_mat[(m+1)] * mod)
      
      ## build in the diagonal
      for (k in 1:(dim(submat)[1]-1)) {
        leslie_mat <- rbind(leslie_mat, 
                            matrix(c(rep(0,k-1),submat$s[k],rep(0,(dim(submat)[1]-k))), nrow=1)
        )
      }
      
      eigenvalue_df <- eigenvalue_df %>% rbind(data.frame(vial = as.character(i), trt = as.character(submat$vial[1]), exp = as.factor(as.character(submat$date2[1])), fitness = as.numeric(eigen(leslie_mat)$values[1])))
      
    }
    
    eigens <- eigens %>% rbind(eigenvalue_df %>% mutate(fitness = log(fitness)) %>% mutate(mod_death = (m)))
  }
  
  eigens <- eigens %>% filter(trt%in%subset_vec) %>% droplevels()
  eigens$trt <- factor(eigens$trt, levels = subset_vec)
  eigens2 <- eigens %>% group_by(trt, mod_death) %>% summarize(meanfit = mean(fitness))
  plot_colors <- colorRampPalette(c("red", "black"))
  
  
  for(i in rangevals) {
    cat('\n')
    print("######################################################")
    print(paste0("These are the statistics for age class ",i," and survival divided by ",1/n, sep = ""))
    
    
    eigen_stats <- eigens %>% filter(mod_death == (i)) %>% droplevels()
    #print(str(eigen_stats))
    #print(i)
    kwoutput <- try(kruskal.test(eigen_stats$fitness~ eigen_stats$trt), silent = T)
    captured_output <- capture.output(efg <- dunn.test(eigen_stats$fitness,eigen_stats$trt, method = "bh", table = F, kw = T))
    if(kwoutput == "Error in kruskal.test.default(numeric(0), structure(integer(0), .Label = character(0), class = \"factor\")) : \n  all observations are in the same group\n") {
      print(paste0("chi-sq statistic not significant"))
    } else {
      print(paste0("chi-sq statistic: ",round(kwoutput$statistic,2)))
      print(paste0("df: ",length(table(eigen_stats$trt))-1))
      print(paste0("p: ",signif(kwoutput$p.value,3)))      
    }      
    #print(i)
    if(sum(efg$P.adjusted<0.05)>0) {
      ghi <- cldList(comparison = efg$comparisons, p.value = efg$P.adjusted, threshold = 0.05, print.comp = F) %>% 
        arrange(factor(Group, levels = subset_vec))
      try(assign(x = paste0("es",i), value = as.vector(ghi$Letter)),T)
    } else {
      assign(x = paste0("es",i),value = rep("a",length(subset_vec)))
    }
    assign(x = paste0("estats",i), 
           value = eigen_stats %>% 
             group_by(trt) %>% 
             summarize(mean_fitness = mean(fitness)) %>% 
             arrange(factor(trt, levels = subset_vec)) %>%
             dplyr::select(mean_fitness) %>%
             unlist() %>% 
             unname()
    )
  } 
  
  eigens2$trt <- plyr::revalue(eigens2$trt, c("28n" = "atrn","30p" = "apan", "11c"="atrc","Gno" = "5sp","4d" = "apoc", "35u" = "aace","Ax" = "Ax","2b" = "lplc","13g" = "bsub","15h" = "ecok","75o" = "lplw","3c" = "lfrc", "1a" = "lbrc","70j" = "lfrk", "71k" = "lbga"))
  print(n)
  print(m)
  #  jpeg(h=7, w = 15, filename = paste0("fitness-permute-fecundity-group_",xgroup,"_",n,".jpg"), units = "in",res = 300)
  print(ggplot(eigens2, aes(x = trt, y=meanfit,group = as.factor(mod_death), col = as.factor(mod_death), fill = as.factor(mod_death))) + 
          geom_step(aes(x = trt, y = meanfit)) + #, position = position_dodge(width = widthj), size = sizej) +
          scale_color_manual(values = plot_colors(length(rangevals)),
                             name = "Age class") +
          coord_cartesian(ylim = c(-2,3)) +
          ylab("log fitness") +
          xlab("bacterial treatment") +
          theme(legend.position = "bottom") +
          guides(color = guide_legend(ncol = 5, byrow=T, reverse = T)) +
          theme_cowplot() + 
          theme(legend.position = "bottom", 
                legend.text.align = 0,
                legend.key.width = unit(.5,"in"), 
                text = element_text(size = 12),
                axis.text.x = element_text(angle = 90),
                legend.justification = "center") + 
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[7]))+0.03, label = get(paste0("es",rangevals[7])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[6]))+0.03, label = get(paste0("es",rangevals[6])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[5]))+0.03, label = get(paste0("es",rangevals[5])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[4]))+0.03, label = get(paste0("es",rangevals[4])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[3]))+0.03, label = get(paste0("es",rangevals[3])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[2]))+0.03, label = get(paste0("es",rangevals[2])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)+
          annotate(geom="label", x= c(1.05:(length(subset_vec)+0.05)), y= get(paste0("estats",rangevals[1]))+0.03, label = get(paste0("es",rangevals[1])),size=3, parse=T, hjust = 0, fill = NA,label.size=0)
  )
  #    dev.off()
}
Figure S5
Values for flies monoassociated with Acetobacter species

Values for flies monoassociated with E. coli


