#!/usr/bin/env Rscript

"""
:Author: Ji Research Group/Stanford Genome Technology Center
:Contact: sgreer2@stanford.edu
:Creation date: 24.11.2016
:Description: 
This script plots the location of reads with HMW barcodes (using output of: count_bcs_in_windows.py)
This script requires:
- all of the R packages listed (imported) below

Revisions:
None to date
CURRENT VERSION: 1.0
"""

library(grid)
library(reshape)
library(ggplot2)
library(gridExtra)
library(lattice)

setwd(".")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (SV name)", call.=FALSE)
}

bkpt_name<-args[1]

file_1 <- paste0(bkpt_name,"_1.bc_windows.txt")
file_2 <- paste0(bkpt_name,"_2.bc_windows.txt")
file_hap <- paste0(bkpt_name,"_hap_bcs.txt")

# read in data frames
df_1 <- read.table(file_1,header=T)
df_2 <- read.table(file_2,header=T)

# read in haplotyped bc data frame
hap_bcs <- read.table(file_hap,header=T)

# obtain names of barcode columns
bc_cols<-tail(colnames(df_1),-3)

# make a vector of the first position
first_pos <- vector("numeric")
counter=3
for(i in bc_cols){
  counter = counter+1
  result<-min((df_1$window_start)[df_1[[paste(i)]]>0])
  first_pos[counter]=result
}

#################### PLOT FIRST BREAKPOINT #####################

# append vector of first positions to data frame
df_1<-as.data.frame(df_1)
first_pos_df<-as.data.frame(t(first_pos))
names(first_pos_df)<-names(df_1)
df_bind1 <- rbind(df_1,first_pos_df)


# sort data frame based on first position
df_sort1 <- df_bind1[,order(first_pos)]
df_sort1 <- head(df_sort1,-1)

# obtain new (sorted) list of df columns
bc_cols_st<-head(colnames(df_sort1),-3)

# loop through data frame -- if count is anything greater than 0 then replace value with desired y-value (+1 for each new barcode)
counter=0
for(i in bc_cols_st){
  counter = counter + 1
  df_sort1[[paste(i)]][df_sort1[[paste(i)]]>0] <- counter 
}

# melt the df in preparation for plotting -- retain only those rows with value >0 
df_melt1 <- melt(df_sort1,measure.vars = head(colnames(df_sort1),-3))
df_melt1 <-subset(df_melt1, value>0)


# add the haplotype for each barcode
df_melt1["bcs_split"] <- lapply(df_melt1["variable"], function(x) gsub("\\..*","",x) )
m1 <- merge(df_melt1, hap_bcs, by.x = "bcs_split", by.y = "bcs", all.x=TRUE)
m1$haps[is.na(m1$haps)] <- "Unassigned"

# determine breakpoint (so can plot dashed line)
brkpt1 <- min(df_1$window_start) + ((max(df_1$window_end) - min(df_1$window_start))/2)

label_calc <-min(m1$window_start)
label_low<-round(label_calc/100000)*100000
label_high<-label_low + 400000

label_y<-ceiling(max(m1$value)/20)*20

myColors <- c("1"="tomato", "2"="steelblue", "Unassigned"="grey")

theme_set(theme_bw(35))
p1 <- ggplot(m1,aes(x=(window_start/1000000),y=value, colour=haps)) +
  geom_point(size=0.6) + 
  scale_colour_manual(name = "haps", values = myColors) +
  xlab("chr10 (Mb)") +
  ylab("SV-specific barcode") +
  theme(plot.margin= unit(c(0.5, 0.1, 0.5, 0.1), "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_blank(), legend.position="top", legend.key = element_blank(), legend.text=element_text(size=28), legend.key.width=unit(2,"line"),legend.key.size = unit(10,"line"),
        panel.border = element_rect(colour = "black", fill=NA, size=3),
        axis.text=element_text(size=50),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.title.x=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_vline(xintercept = brkpt1/1000000, linetype="longdash") +
  scale_x_continuous(limits = c(min(df_1$window_start)/1000000, max(df_1$window_end)/1000000), breaks=seq(label_low/1000000,label_high/1000000,0.2)) +
  scale_y_continuous(breaks=seq(0,label_y,20))


#################### PLOT SECOND BREAKPOINT #####################

# append vector of first positions to data frame
df_2<-as.data.frame(df_2)
first_pos_df<-as.data.frame(t(first_pos))
names(first_pos_df)<-names(df_2)
df_bind2 <- rbind(df_2,first_pos_df)

# sort data frame based on first position
df_sort2 <- df_bind2[,order(first_pos)]
df_sort2 <- head(df_sort2,-1)

# obtain new (sorted) list of df columns
bc_cols_st<-head(colnames(df_sort2),-3)

# loop through data frame -- if count is anything greater than 0 then replace value with desired y-value (+1 for each new barcode)
counter=0
for(i in bc_cols_st){
  counter = counter + 1
  df_sort2[[paste(i)]][df_sort2[[paste(i)]]>0] <- counter 
}

# melt the df in preparation for plotting -- retain only those rows with value >0 
df_melt2 <- melt(df_sort2,measure.vars = head(colnames(df_sort2),-3))
df_melt2 <-subset(df_melt2, value>0)

# add the haplotype for each barcode
df_melt2["bcs_split"] <- lapply(df_melt2["variable"], function(x) gsub("\\..*","",x) )
m2 <- merge(df_melt2, hap_bcs, by.x = "bcs_split", by.y = "bcs", all.x=TRUE)
m2$haps[is.na(m2$haps)] <- "Unassigned"

# determine breakpoint (so can plot dashed line)
brkpt2 <- min(df_2$window_start) + ((max(df_2$window_end) - min(df_2$window_start))/2)

label_calc2<-min(m2$window_start)
label_low2<-round(label_calc2/100000)*100000
label_high2<-label_low2 + 400000

theme_set(theme_bw(35))
p2 <- ggplot(m2,aes(x=(window_start/1000000),y=value, colour=haps)) +
  geom_point(size=0.6) + 
  scale_colour_manual(name = "haps", values = myColors) +
  xlab("chr10 (Mb)") +
  ylab("Barcode") +
  theme(plot.margin= unit(c(0.5, 0.5, 0.5, 0.1), "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=3),
        legend.position="none",
        axis.text=element_text(size=18),
        axis.text.x = element_text(size=30),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank())+
  geom_vline(xintercept = brkpt2/1000000, linetype="longdash") +
  scale_x_continuous(limits = c(min(df_2$window_start)/1000000, max(df_2$window_end)/1000000),breaks=seq(label_low2/1000000,label_high2/1000000,0.2)) +
  scale_y_continuous(breaks=seq(0,label_y,20))

label = textGrob("chr10 (Mb)",hjust=0.5, gp = gpar(fontsize = 40))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p1)

full_plot<-grid.arrange(mylegend, arrangeGrob(p1 + theme(legend.position="none"),p2,ncol=2,nrow=1,widths=c(60, 52), heights=c(1)), label, nrow=3, heights=c(1,10,1))

ggsave(paste0(bkpt_name,"_bkpt_region.pdf"), plot=full_plot, width=15, height=6, units="in", dpi=600)

dev.off()
