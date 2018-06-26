#this script imports the large CSV and makes graphs.
setwd("C://Users//insan_000//Documents//summer 2018//simulation data//bigg data")
library(tidyverse)

dfrid = read_csv("all_averages_for_ID_sim_1000.csv")
dfrid$Kappa = as.factor(dfrid$Kappa)
dfrid$`GC%` = as.factor(dfrid$`GC%`)

dfrsub0 = subset(dfrid, dfrid$`GC%` == 0)
dfrsub1 = subset(dfrid, dfrid$`GC%` == 0.1)
dfrsub2 = subset(dfrid, dfrid$`GC%` == 0.2)
dfrsub3 = subset(dfrid, dfrid$`GC%` == 0.3)
dfrsub4 = subset(dfrid, dfrid$`GC%` == 0.4)
dfrsub5 = subset(dfrid, dfrid$`GC%` == 0.5)
dfrsub6 = subset(dfrid, dfrid$`GC%` == 0.6)
dfrsub7 = subset(dfrid, dfrid$`GC%` == 0.7)
dfrsub8 = subset(dfrid, dfrid$`GC%` == 0.8)
dfrsub9 = subset(dfrid, dfrid$`GC%` == 0.9)
dfrsub10 = subset(dfrid, dfrid$`GC%` == 1)

ggplot(dfrsub0,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 0")
ggsave("cms_vs_ID_gc0.png",device = "png")

ggplot(dfrsub1,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 10")
ggsave("cms_vs_ID_gc10.png",device = "png")

ggplot(dfrsub2,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 20")
ggsave("cms_vs_ID_gc20.png",device = "png")

ggplot(dfrsub3,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 30")
ggsave("cms_vs_ID_gc30.png",device = "png")

ggplot(dfrsub4,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 40")
ggsave("cms_vs_ID_gc40.png",device = "png")

ggplot(dfrsub5,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 50")
ggsave("cms_vs_ID_gc50.png",device = "png")

ggplot(dfrsub6,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 60")
ggsave("cms_vs_ID_gc60.png",device = "png")

ggplot(dfrsub7,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 70")
ggsave("cms_vs_ID_gc70.png",device = "png")

ggplot(dfrsub8,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 80")
ggsave("cms_vs_ID_gc80.png",device = "png")

ggplot(dfrsub9,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 90")
ggsave("cms_vs_ID_gc90.png",device = "png")

ggplot(dfrsub10,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for GC% = 100")
ggsave("cms_vs_ID_gc100.png",device = "png")

ggplot(dfrid,aes(x=`Average ID%`,y=`Average CMs`,color=Kappa)) + geom_smooth() + theme_classic() + ggtitle("CMs vs ID% for variant GC%") + facet_wrap(~`GC%`)
ggsave("cms_vs_ID_facet.png",device = "png")