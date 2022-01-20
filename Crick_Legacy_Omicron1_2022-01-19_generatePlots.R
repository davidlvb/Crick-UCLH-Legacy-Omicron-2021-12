# 
# Crick/UCLH Legacy Cohort Neutralisation Data Analysis
# Third analysis
# Against SARS-CoV-2 Variants of Concern Alpha, Delta, and Omicron
#
# 23 December 2021
#
# Edward J Carr & David LV Bauer
# Cell Biology of Infection Laboratory & RNA Virus Replication Laboratory
# The Francis Crick Institute
#
# This work is licensed under a CC BY 4.0 license and is free to use with attribution
# 

# IMPORTANT NOTES:
# Ages have been rounded *DOWN* to nearest multiple of 5
# to preserve anonymity; e.g.  44 -> 40   47 -> 45, etc. and
# elig_study_id (participant ID) is randomised and does not match IDs 
# in previous data sets to prevent identification.

#############################
#        SET UP DATA        #
#############################


### Remove all stored data from environment, reset plot, set working dir ###
rm(list = ls())
dev.off()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Load required packages ### 
library(tidyverse)
library(lubridate)
library(magrittr)
library(ggpubr)
library(patchwork)
library(cowplot)
# library(furniture)
library(gtsummary)
library(rstatix)
library(ragg)
library(systemfonts)
library(khroma)
library(infer)
library(ggvenn)

### Set up plotting parameters ###
register_variant(
  name = "Helvetica Neue Thin", 
  family = "Helvetica Neue", 
  weight = "ultralight")

vacc_types <- c("BNT162b2", "AZD1222")
centre_names <- c("crick", "uclh")

### Load Legacy Study Data table ###
load("LegacyOmicron_19jan2022_PUBLIC.Rdata")


#############################
#      DEFINE COHORTS       #
#############################

nAb_cohorted <- legacyOm

### Two-Dose Cohort ###
nAb_dose2_4wk <- nAb_cohorted %>% 
  filter(visit_date > date_dose_2+13 & visit_date < date_dose_2+43) %>% group_by(elig_study_id) %>% 
  slice_head(n=1) %>% 
  ungroup() %>%
  mutate(cohort = "2-6")
nAb_dose2_14wk <- nAb_cohorted %>% 
  filter(visit_date > date_dose_2+83 & visit_date < date_dose_2+113) %>% group_by(elig_study_id) %>% 
  slice_head(n=1) %>% 
  ungroup() %>%
  mutate(cohort = "12-16")
nAb_dose2 <- rbind(nAb_dose2_4wk, nAb_dose2_14wk) %>%
  filter( is.na(date_dose_3) |  visit_date < date_dose_3+3) %>%   # Make sure they have not received a third dose
  filter( symptomatic_afterDose2BeforeVisit == "No" ) %>% # Has not had a breakthru infection yet
  filter( posTest_afterDose2BeforeVisit == "No") %>%
  mutate(cohort=factor(cohort, levels=c("2-6", "12-16"))) %>% # Set as a factor for correct plotting order
  mutate(days_sinceEvent = days_sinceDose2)
rm(nAb_dose2_4wk, nAb_dose2_14wk)
### Two-Dose cohorts split 2-6wk/12-16w only ###
nAb_dose2_4wk <- nAb_dose2 %>% filter(cohort == "2-6")
nAb_dose2_14wk <-nAb_dose2 %>% filter(cohort == "12-16")


### Breakthrough Cohort (Dose 2 only) ###
nAb_bThru <- nAb_cohorted %>%
  filter(visit_date > date_dose_2+13) %>%   # Visit after dose 2 + 2w
  filter((visit_date < date_dose_3 + 3) | is.na(date_dose_3)) %>% #Visit before 3rd dose, or no 3rd dose yet
  filter(posTest_afterDose2BeforeVisit == "Yes" | symptomatic_afterDose2BeforeVisit == "Yes") %>%
  filter((visit_date > posTest_earliest+6 & visit_date < posTest_earliest+50) | # Visit 1-7 weeks after a posTest
           (visit_date > posTest_latest+6 & visit_date < posTest_latest+50)  | # Visit 1-7 weeks after a posTest
           (visit_date > Sx_start_date+6 & visit_date < Sx_start_date+50) ) %>% # Visit 1-7 weeks after symptoms 
  group_by(elig_study_id) %>% arrange(visit_date, .by_group = TRUE) %>% slice_head() %>% # Make sure only 1 visit/person
  ungroup() %>%
  mutate(cohort = "Breakthru") %>%
  mutate(days_sinceEvent = pmin(abs(days_sincePosTest_earliest), abs(days_sincePosTest_latest), na.rm = TRUE))  # Here, this is defined as positive test (note symp. onset is earlier!)

# Same sampling post-breakthru as dose2/3 (same overall result, but fewer datapoints, 
# incl. a few just at edges of cutoff) included here for comparison to illustrate.
nAb_bThru_1442 <- nAb_cohorted %>%
  filter(visit_date > date_dose_2+13) %>%   # Visit after dose 2 + 2w
  filter((visit_date < date_dose_3 + 3) | is.na(date_dose_3)) %>% #Visit before 3rd dose, or no 3rd dose yet
  filter(posTest_afterDose2BeforeVisit == "Yes" | symptomatic_afterDose2BeforeVisit == "Yes") %>%
  filter((visit_date > posTest_earliest+13 & visit_date < posTest_earliest+43) | # Visit 2-6 weeks after a posTest
           (visit_date > posTest_latest+13 & visit_date < posTest_latest+43)  | # Visit 2-6 weeks after a posTest
           (visit_date > Sx_start_date+13 & visit_date < Sx_start_date+43) ) %>% # Visit 2-6 weeks after symptoms 
  group_by(elig_study_id) %>% arrange(visit_date, .by_group = TRUE) %>% slice_head() %>% # Make sure only 1 visit/person
  ungroup() %>%
  mutate(cohort = "BThru+(2-6wk)") %>%
  mutate(days_sinceEvent = min(abs(days_sincePosTest_earliest), abs(days_sincePosTest_latest)))


### Three-Dose Cohort ###
nAb_dose3_before <- nAb_cohorted %>% 
  filter(visit_date > date_dose_3-28 & visit_date < date_dose_3+3) %>% group_by(elig_study_id) %>% 
  slice_head(n=1) %>% 
  ungroup() %>%
  mutate(cohort = "PRE-Boost")
nAb_dose3_after <- nAb_cohorted %>% 
  filter(visit_date > date_dose_3+13 & visit_date < date_dose_3+43) %>% group_by(elig_study_id) %>% 
  slice_head(n=1) %>% 
  ungroup() %>%
  mutate(cohort = "POST-Boost")
nAb_dose3_all <- rbind(nAb_dose3_before, nAb_dose3_after) %>%
  mutate(days_sinceEvent = days_sinceDose3) %>% 
  filter(dose_2 == "BNT162b2") # Not enough AZD1222 recipients yet!

# Exclude those with breakthrough infections prior to dose 3
nAb_dose3 <- nAb_dose3_all %>% filter(posTest_afterDose2BeforeVisit == "No" & symptomatic_afterDose2BeforeVisit == "No") %>%
  mutate(cohort = factor(cohort, levels=c("PRE-Boost", "POST-Boost")))

# There are some recent sera with only Omicon titres available.
# Create a separate cohort with no missing values, to check we have not confounded results in fold-change somehow
nAb_dose3_noNA <- nAb_dose3 %>% filter( !is.na(ic50_Alpha) & !is.na(ic50_Delta) & !is.na(ic50_Omicron) ) %>% mutate(cohort = paste(cohort, "_noNA", sep="")) 
rm(nAb_dose3_before, nAb_dose3_after)
nAb_dose3_before <- nAb_dose3 %>% filter(cohort == "PRE-Boost")
nAb_dose3_after <- nAb_dose3 %>% filter(cohort == "POST-Boost")


### Age-Matched TWO and THREE dose sub-cohort ###
# Third 'booster' dose recipients are older than general 2-dose cohort, 
# so proper comparison between doses should be age-matched
# (in this case, 45+ years old)
nAb_dose23_older <- rbind(nAb_dose2 %>% filter(cohort == "2-6") %>% filter(dose_2 == "BNT162b2"), nAb_dose3 %>% filter(cohort == "POST-Boost")) %>%
  filter(dose_2 == "BNT162b2") %>%
  filter(age > 44)
nAb_dose23_older_crick <- nAb_dose23_older %>% filter(centre == "crick")

nAb_dose23_older$cohort <- droplevels(recode(nAb_dose23_older$cohort, "2-6" = "(45+yrs) Dose 2", "POST-Boost"="(45+yrs) Dose 3"))
nAb_dose23_older_crick$cohort <- droplevels(recode(nAb_dose23_older_crick$cohort, "2-6" = "Crick::(45+yrs) Dose 2", "POST-Boost"="Crick::(45+yrs) Dose 3"))



#####################################################
### Generate Venn Diagram to illustrate sampling  ###
#####################################################

# Set up venn diagram colours
highcontrast <- colour("contrast")
vennColours <- c("#000000", as.vector(highcontrast(3)))

venn_all <- list( `2 doses\n(2 to 6wk)` = nAb_dose2 %>% filter(cohort == "2-6") %>% pull(elig_study_id),
                  `2 doses\n(12-16wk)` = nAb_dose2 %>% filter(cohort == "12-16") %>% pull(elig_study_id),
                  `PRE-Boost\n(-4 to 0wk)` = nAb_dose3 %>% filter(cohort == "PRE-Boost") %>% pull(elig_study_id),
                  `POST-Boost\n(2 to 6wk)` = nAb_dose3 %>% filter(cohort == "POST-Boost") %>% pull(elig_study_id) )

vacc <- "BNT162b2"
venn_Pf <- list( `2 doses\n(2 to 6wk)` = nAb_dose2 %>% filter(dose_2 == vacc & cohort == "2-6") %>% pull(elig_study_id),
                  `2 doses\n(12 to 16wk)` = nAb_dose2 %>% filter(dose_2 == vacc & cohort == "12-16") %>% pull(elig_study_id),
                  `PRE-Boost\n(-4 to 0wk)` = nAb_dose3 %>% filter(dose_2 == vacc & cohort == "PRE-Boost") %>% pull(elig_study_id),
                  `POST-Boost\n(2 to 6wk)` = nAb_dose3 %>% filter(dose_2 == vacc & cohort == "POST-Boost") %>% pull(elig_study_id) )

vacc <- "BNT162b2"
venn_Pf_dose23_older <- list( `2 doses\n(2 to 6wk)` = nAb_dose23_older %>% filter(dose_2 == vacc & cohort == "(45+yrs) Dose 2") %>% pull(elig_study_id),
                 `2 doses\n(12 to 16wk)` = vector(),
                 `PRE-Boost\n(-4 to 0wk)` = vector(),
                 `POST-Boost\n(2 to 6wk)` = nAb_dose23_older %>% filter(dose_2 == vacc & cohort == "(45+yrs) Dose 3") %>% pull(elig_study_id) )


vacc <- "AZD1222"
venn_AZ <- list( `2 doses\n(2 to 6wk)` = nAb_dose2 %>% filter(dose_2 == vacc & cohort == "2-6") %>% pull(elig_study_id),
                 `2 doses\n(12-16wk)` = nAb_dose2 %>% filter(dose_2 == vacc & cohort == "12-16") %>% pull(elig_study_id),
                 `PRE-Boost\n(-4 to 0wk)` = nAb_dose3 %>% filter(dose_2 == vacc & cohort == "PRE-Boost") %>% pull(elig_study_id),
                 `POST-Boost\n(2 to 6wk)` = nAb_dose3 %>% filter(dose_2 == vacc & cohort == "POST-Boost") %>% pull(elig_study_id) )


ggvenn(venn_all, digits=0, set_name_size = NA, fill_color = vennColours) + ggtitle("Whole cohort")
ggvenn(venn_Pf, digits=0, set_name_size = 4, fill_color = vennColours) + ggtitle("BNT162b2 cohort")
ggvenn(venn_Pf_dose23_older, digits=0, set_name_size = 4, fill_color = c(vennColours[1], NA, NA, vennColours[4] )) + ggtitle("Age-matched cohort: 2 dose & 3 dose ")
ggvenn(venn_AZ, digits=0, set_name_size = 4, fill_color = c(vennColours[1:2], NA, NA)) + ggtitle("AZD1222 cohort")


vennPlot <- ggvenn(venn_Pf, digits=0, set_name_size = 0, text_size = 3, stroke_size=0.27, fill_color = vennColours)
ggsave("Venn_Pfizer.pdf", vennPlot, width = 7, height=7, units="cm")

vennPlot <- ggvenn(venn_Pf_dose23_older, digits=0, set_name_size = 0, text_size = 3, stroke_size=0.27, fill_color = c(vennColours[1], NA, NA, vennColours[4] )) 
ggsave("Venn_Pfizer_older.pdf", vennPlot, width = 7, height=7, units="cm")

vennPlot <- ggvenn(venn_AZ, digits=0, set_name_size = 0, text_size = 3, stroke_size=0.27, fill_color = c(vennColours[1:2], NA, NA))  
ggsave("Venn_AZ.pdf", vennPlot, width = 7, height=7, units="cm")

######################################################
### Summarise participant & sample numbers overall ###
######################################################

# Consider all sera
legacyOm <- rbind(nAb_dose2, nAb_dose3, nAb_bThru)

# Count number of samples
legacyOm %>% group_by(dose_2) %>% nrow()

# Count number of unique participants
legacyOm_participants <- legacyOm %>% group_by(elig_study_id) %>% slice_head(n=1) %>% ungroup()
legacyOm_participants %>% group_by(dose_2) %>% nrow()


# Consider only vaccine (no breakthru counts)
legacyOm <- rbind(nAb_dose2, nAb_dose3)

# Count number of samples
legacyOm %>% group_by(dose_2) %>% furniture::table1(age)

# Count number of unique participants
legacyOm_participants <- legacyOm %>% group_by(elig_study_id) %>% slice_head(n=1) %>% ungroup()
legacyOm_participants %>% group_by(dose_2) %>% furniture::table1(age)

# Count number of breakthru samples/participants
nAb_bThru %>% group_by(dose_2) %>% furniture::table1(age)

# Count number of breakthru participants already captured in above (ie gave prior sample)
`%notin%` <- Negate(`%in%`)
sum(nAb_bThru$elig_study_id %in% legacyOm_participants$elig_study_id)
sum(nAb_bThru$elig_study_id %notin% legacyOm_participants$elig_study_id)

###################################################
### Join together all cohorts for summary table ###
###################################################

legacyOm <- rbind(nAb_dose2, nAb_bThru, nAb_dose3)

# Preview within R using furniture package
# NB NA's in days_sinceDose3 prevent its use here in the preview
legacyOm %>% group_by(cohort, dose_2) %>% furniture::table1(centre, age, sex, symptomatic_beforeVisit,days_sinceDose2,
                                                            Quant_ic50_Alpha, Quant_ic50_Delta, Quant_ic50_Omicron,
                                                            ic50_Alpha, ic50_Delta, ic50_Omicron,
                                                            Fold_ic50_AlphaOmicron,Fold_ic50_DeltaOmicron,
                                                            rounding_perc=0, NAkeep=TRUE, test = FALSE)

# Can also generate using table1 package
# NB, conflicts with furniture package, you must restart R / fully detach before
# using the below
# legacyOm %>% 
#   table1(~centre+age+sex+
#            symptomatic_beforeJoined+symptomatic_beforeVisit+
#            days_sinceDose2+days_sinceDose3+
#            Quant_ic50_Alpha+Quant_ic50_Delta+Quant_ic50_Omicron+
#            ic50_Alpha+ic50_Delta+ic50_Omicron+
#            Fold_ic50_AlphaDelta+Fold_ic50_AlphaOmicron+Fold_ic50_DeltaOmicron  | cohort-dose_2, data = .,
#          ## This removes the default summary column:
#          overall = FALSE,
#          render.continuous=c(.="Mean (SD)", "Median [IQR]"="Median [Q1-Q3]", "Geo. mean (Geo. CV%)"="GMEAN (GCV%)" ))%>% write_file(file="table1.html")

# Use gtsummary to generate HTML table for export.
# Prevent use of comma to separate thousands
# list("style_number-arg:big.mark" = "\U00B7") %>%  # Use middot as thousands separator
list("style_number-arg:big.mark" = "") %>%
  set_gtsummary_theme()

legacyOm %>%   mutate(cohort_dose = paste(dose_2, cohort, sep="::")) %>% 
  select(cohort_dose, centre, sex, age,symptomatic_beforeVisit,
         days_sinceDose2,days_sinceDose3,
         Quant_ic50_Alpha, Quant_ic50_Delta, Quant_ic50_Omicron,
         ic50_Alpha, ic50_Delta, ic50_Omicron,
         Fold_ic50_AlphaOmicron,Fold_ic50_DeltaOmicron) %>%
  tbl_summary(by=cohort_dose,
              statistic = all_continuous() ~ "{median} [{p25}-{p75}]",
              digits = list( Fold_ic50_AlphaOmicron ~ 1 ,
                             Fold_ic50_DeltaOmicron ~ 1)
  )


# We need to also generate a summary table for stratification by symptoms
nAb_dose2_4wk %>% group_by(symptomatic_beforeVisit, dose_2) %>% furniture::table1(centre, age, sex, symptomatic_beforeVisit,days_sinceDose2,
                                                            Quant_ic50_Alpha, Quant_ic50_Delta, Quant_ic50_Omicron,
                                                            ic50_Alpha, ic50_Delta, ic50_Omicron,
                                                            Fold_ic50_AlphaOmicron,Fold_ic50_DeltaOmicron,
                                                            rounding_perc=0, NAkeep=TRUE, test = FALSE)

nAb_dose2_4wk %>% mutate(symp_dose = paste(dose_2, symptomatic_beforeVisit, sep="::")) %>% 
  select(symp_dose,centre, sex, age,symptomatic_beforeVisit,
         days_sinceDose2,days_sinceDose3,
         Quant_ic50_Alpha, Quant_ic50_Delta, Quant_ic50_Omicron,
         ic50_Alpha, ic50_Delta, ic50_Omicron,
         Fold_ic50_AlphaOmicron,Fold_ic50_DeltaOmicron) %>%
  tbl_summary(by=symp_dose,
              statistic = all_continuous() ~ "{median} [{p25}-{p75}]",
              digits = list( Fold_ic50_AlphaOmicron ~ 1 ,
                             Fold_ic50_DeltaOmicron ~ 1)
  )

# And also need to also generate a summary table for age-matched BNT162b2 dose 2/3 cohort

nAb_dose23_older %>% group_by(cohort) %>% furniture::table1(centre, age, sex, symptomatic_beforeVisit,days_sinceDose2,
                                                            Quant_ic50_Alpha, Quant_ic50_Delta, Quant_ic50_Omicron,
                                                            ic50_Alpha, ic50_Delta, ic50_Omicron,
                                                            Fold_ic50_AlphaOmicron,Fold_ic50_DeltaOmicron,
                                                            rounding_perc=0, NAkeep=TRUE, test = FALSE)

nAb_dose23_older %>%
  select(cohort, centre, sex, age,symptomatic_beforeVisit,
         days_sinceDose2,days_sinceDose3,
         Quant_ic50_Alpha, Quant_ic50_Delta, Quant_ic50_Omicron,
         ic50_Alpha, ic50_Delta, ic50_Omicron,
         Fold_ic50_AlphaOmicron,Fold_ic50_DeltaOmicron) %>%
  tbl_summary(by=cohort,
              statistic = all_continuous() ~ "{median} [{p25}-{p75}]",
              digits = list( Fold_ic50_AlphaOmicron ~ 1 ,
                             Fold_ic50_DeltaOmicron ~ 1)
  )



#################################
###       BEGIN PLOTS           #
#################################

#### #### #### ##
#### PANEL B ####
#### #### #### ##


source("9-dlvb-plotFNs_Omicron.R")

relevantData_nAb_dose2 <-  nAb_dose2 %>% pivot_longer(cols = starts_with("ic50_"), 
                                                      names_to = "strain", 
                                                      values_to = "ic50") %>%
  mutate(strain = str_remove(strain, "ic50_")) %>%
  filter(strain %in% c("Alpha", "Delta", "Omicron")) %>%
  mutate(strain = factor(strain, levels = strainOrder))

outplot2dose <-ggplot(relevantData_nAb_dose2,aes(x=cohort, y=ic50, color=strain) ) + 
  facet_wrap(~dose_2+strain, ncol =6) +
  labs(title="Neutralisation following 2nd Dose", x="Weeks post-second dose")
panel.b <- ggbauer_violin(outplot2dose) + theme(axis.title.x = element_text(size=12))


panel.b.barplot <- barplot_bauer_prop(relevantData_nAb_dose2) + 
  facet_wrap(~dose_2+cohort+strain, ncol = 3)

panel.b
panel.b.barplot

panel.b.barplot.nolabs <- panel.b.barplot + theme(strip.background = element_blank(), strip.text.x = element_blank()) + ylab(element_blank())
panel.b.nolabs <- panel.b + theme(strip.background = element_blank(), strip.text.x = element_blank(), title = element_blank())

row1 <- plot_grid(panel.b.nolabs, NULL, panel.b.barplot.nolabs, nrow=1, rel_widths = c(3,0.1,1.5) )
row1

ggsave("row1.png", plot = row1, height=90, width=300, units="mm")

#### Panel B stats ###
### Calculate boostrapped stats data re. fold-changes between variants, median and 95%CI
FCdata <- relevantData_nAb_dose2 %>% filter(dose_2 == "BNT162b2", cohort == "2-6")

inferFC_AO <- FCdata %>%
  infer::specify(response = Fold_ic50_AlphaOmicron) %>%  
  infer::generate(reps = 5000, type = "bootstrap") %>%
  infer::calculate(stat = "median")

inferFC_DO <- FCdata %>%
  infer::specify(response = Fold_ic50_DeltaOmicron) %>%  
  infer::generate(reps = 5000, type = "bootstrap") %>%
  infer::calculate(stat = "median")

median(FCdata$Fold_ic50_AlphaOmicron)
median(inferFC_AO$stat)
inferFC_AO %>% infer::get_confidence_interval(level = 0.95)

median(FCdata$Fold_ic50_DeltaOmicron)
median(inferFC_DO$stat)
inferFC_DO %>% infer::get_confidence_interval(level = 0.95)




#### #### #### ##
#### PANEL C ####
#### #### #### ##
# 2dose, 2-6wk, Stratified by symptoms

relevantData_nAb_dose2_bySymp <- relevantData_nAb_dose2 %>%  filter(cohort=="2-6") %>% filter(!is.na(priorSxAtFirstVisit))
relevantData_nAb_dose2_bySymp$priorSxAtFirstVisit <- recode(relevantData_nAb_dose2_bySymp$priorSxAtFirstVisit, Y="Yes", N="No")

panel.c <- relevantData_nAb_dose2_bySymp %>%
  filter(cohort=="2-6") %>%
  filter(!is.na(priorSxAtFirstVisit)) %>%
  ggplot(aes(x=priorSxAtFirstVisit, y=ic50, color=strain #, group=elig_study_id
  )) + 
  facet_wrap(~dose_2+strain, ncol=6) +
  labs(title="2-dose Vaccine Recipients",
       subtitle="Neutralisation titres at 4 weeks after 2nd dose\nstratified by prior COVID symptoms",
       x= "Prior COVID Symptoms")
panel.c <- ggbauer_violin(panel.c) + 
  theme(axis.title.x = element_text(size=12))
panel.c

# Chi-squared for AZD1222
AZOmi <- relevantData_nAb_dose2_bySymp %>% filter(dose_2 == "AZD1222" & strain == "Omicron") %>% select(elig_study_id, priorSxAtFirstVisit, Quant_ic50_Omicron)
M <- rbind( summary(AZOmi %>% filter(priorSxAtFirstVisit == "Yes") %>% pull(Quant_ic50_Omicron)),
       summary(AZOmi %>% filter(priorSxAtFirstVisit == "No") %>% pull(Quant_ic50_Omicron)) ) %>% as_tibble() %>% select(1:2)
M
stats::chisq.test(M) %>% tidy() 

# Wilcoxon UNPAIRED for BNT162b2
PfOmi <- relevantData_nAb_dose2_bySymp %>% filter(dose_2 == "BNT162b2" & strain == "Omicron") %>% select(elig_study_id, priorSxAtFirstVisit, strain, ic50)
PfOmi_SxYES <- PfOmi %>% filter(priorSxAtFirstVisit == "Yes")
PfOmi_SxNO <- PfOmi %>% filter(priorSxAtFirstVisit == "No")
rbind( wilcox.test(PfOmi_SxYES$ic50, PfOmi_SxNO$ic50, paired = FALSE) %>% tidy(),
       ks.test(PfOmi_SxYES$ic50, PfOmi_SxNO$ic50, paired = FALSE) %>% tidy()  )



#### #### #### ##
#### PANEL D ####
#### #### #### ##
# Breakthru

relevantData_nAb_bThru <-  nAb_bThru %>% pivot_longer(cols = starts_with("ic50_"), 
                                                      names_to = "strain", 
                                                      values_to = "ic50") %>%
  mutate(strain = str_remove(strain, "ic50_")) %>%
  filter(strain %in% c("Alpha", "Delta", "Omicron")) %>%
  mutate(strain = factor(strain, levels = strainOrder))

outplotBT <-ggplot(relevantData_nAb_bThru,aes(x=strain, y=ic50, color=strain) ) + 
  facet_wrap(~dose_2, ncol =6) +
  labs(title="Breakthough after 2 vaccine doses", x="Variant")
panel.d <- ggbauer_violin(outplotBT) + theme(axis.title.x = element_text(size=12))
panel.d


# We report median and IQR of ALL breakthrus together in main text 
summary(nAb_bThru$ic50_Omicron)



#### #### #### ##
#### PANEL E ####
#### #### #### ##
# Booster

relevantDataBOOST <-  nAb_dose3 %>% 
  pivot_longer(cols = starts_with("ic50_"), names_to = "strain", values_to = "ic50") %>%
  mutate(strain = str_remove(strain, "ic50_")) %>%
  filter(strain %in% c("Alpha", "Delta", "Omicron")) %>%
  mutate(strain = factor(strain, levels = strainOrder)) %>%
  filter(dose_2 == "BNT162b2")


outplotBOOST <- ggplot(relevantDataBOOST,aes(x=cohort, y=ic50, color=strain,)) + 
  facet_wrap(~strain,ncol=6) +
  labs(title="3-dose Vaccine Recipients",
       subtitle="Neutralisation titres before (-28d to +2d) and after\n (+14d to +42d) third dose of BNT162b2")

panel.e <- ggbauer_violin(outplotBOOST)
panel.e

#### Panel E stats ###
### Calculate boostrapped stats data re. fold-changes between variants, median and 95%CI
FCdata <- relevantDataBOOST %>% filter(dose_2 == "BNT162b2", cohort == "POST-Boost")

inferFC_AO <- FCdata %>%
  infer::specify(response = Fold_ic50_AlphaOmicron) %>%  
  infer::generate(reps = 5000, type = "bootstrap") %>%
  infer::calculate(stat = "median")

inferFC_DO <- FCdata %>%
  infer::specify(response = Fold_ic50_DeltaOmicron) %>%  
  infer::generate(reps = 5000, type = "bootstrap") %>%
  infer::calculate(stat = "median")

median(FCdata$Fold_ic50_AlphaOmicron, na.rm = TRUE)
median(inferFC_AO$stat)
inferFC_AO %>% infer::get_confidence_interval(level = 0.95)

median(FCdata$Fold_ic50_DeltaOmicron, na.rm = TRUE)
median(inferFC_DO$stat)
inferFC_DO %>% infer::get_confidence_interval(level = 0.95)


#### #### #### #####
#### OUTPUT C-E ####
#### #### #### #####

panel.c.nolabs <- panel.c + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  title = element_blank()
)
panel.d.nolabs <- panel.d + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  title = element_blank()
)

panel.e.nolabs <- panel.e + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  title = element_blank()
)


row2 <- plot_grid(panel.c.nolabs, panel.d.nolabs, panel.e.nolabs, nrow=1, rel_widths = c(2,1.0,1.5) , align = "hv")
row2

ggsave("row2.png", plot = row2, height=90, width=300, units="mm")



#### #### #### ##
#### PANEL SUPP ####
#### #### #### ##
# d2/3 age-matched


relevantData_older <-  nAb_dose23_older %>% 
  pivot_longer(cols = starts_with("ic50_"), names_to = "strain", values_to = "ic50") %>%
  mutate(strain = str_remove(strain, "ic50_")) %>%
  filter(strain %in% c("Alpha", "Delta", "Omicron")) %>%
  mutate(strain = factor(strain, levels = strainOrder)) %>%
  filter(dose_2 == "BNT162b2")

outplot_older <- ggplot(relevantData_older,aes(x=cohort, y=ic50, color=strain,)) + 
  facet_wrap(~strain,ncol=6) +
  labs(title="NAbT following second and third doses",
       subtitle="Age-matched cohort (45 years and over")
panel.s2 <- ggbauer_violin(outplot_older)
panel.s2

panel.s2.nolabs <- panel.s2 + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  title = element_blank()
)

row3 <- plot_grid(NA, NA, panel.s2.nolabs, nrow=1, rel_widths = c(2,1.0,1.5) , align = "hv")
row3

ggsave("row3.png", plot = row3, height=90, width=300, units="mm")



# Stats comparing VOCS
FCdata <- relevantData_older %>% filter(dose_2 == "BNT162b2", cohort == "(45+yrs) Dose 3")

inferFC_AO <- FCdata %>%
  infer::specify(response = Fold_ic50_AlphaOmicron) %>%  
  infer::generate(reps = 5000, type = "bootstrap") %>%
  infer::calculate(stat = "median")

inferFC_DO <- FCdata %>%
  infer::specify(response = Fold_ic50_DeltaOmicron) %>%  
  infer::generate(reps = 5000, type = "bootstrap") %>%
  infer::calculate(stat = "median")

summary(inferFC_AO)
summary(inferFC_DO)

median(FCdata$Fold_ic50_AlphaOmicron, na.rm = TRUE)
median(FCdata$Fold_ic50_DeltaOmicron, na.rm = TRUE)

inferFC_AO %>% infer::get_confidence_interval(level = 0.95)
inferFC_DO %>% infer::get_confidence_interval(level = 0.95)

# With well-matched data, it is reasonable to now compare across
# use boot package to use bootstratp stats to estimate
# the fold-boost in median NAbTs after 3rd dose


library(boot)
currStrain <- "Alpha"
testTable <- relevantData_older %>% ungroup %>% select(sample_barcode, strain, cohort, ic50) %>% filter(strain == currStrain)
CIbootstrap_unpaired_d23(testTable, currStrain)
currStrain <- "Delta"
testTable <- relevantData_older %>% ungroup %>% select(sample_barcode, strain, cohort, ic50) %>% filter(strain == currStrain)
CIbootstrap_unpaired_d23(testTable, currStrain)
currStrain <- "Omicron"
testTable <- relevantData_older %>% ungroup %>% select(sample_barcode, strain, cohort, ic50) %>% filter(strain == currStrain)
CIbootstrap_unpaired_d23(testTable, currStrain)

