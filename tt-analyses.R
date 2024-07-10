library(nhanesA)
library(tidyverse)
library(survey)


#https://wwwn.cdc.gov/nchs/data/Tutorials/Code/DB303_R.r for how to estimate data for the U.S. population given the complex sampling scheme of NHANES

#Per NHANES, if doing age group analyses for stable estimates in age groups for CDC recommends bracketing by: 5 years and under, 6-11 years, 12-19 years, 20-39 years, 40-59 years, 60 years and over -- we do not do age group analyses so not a concern


test_b <- nhanes("SSCHL_B") #2001 and 2002 #ng/mL measurement
test_c <- nhanes("SSCHL_C") #2003 and 2004 #ng/mL measurement 

test_g <- nhanes("TST_G") #2011 and 2012 #ng/dL measurement
test_h <- nhanes("TST_H") #2013 and 2014 #ng/dL measurement
test_i <- nhanes("TST_I") #2015 and 2016 #ng/dL measurement

test_b <- test_b %>% select(1:2) %>% as_tibble()
test_c <- test_c %>% select(1:2) %>% as_tibble()
test_g <- test_g %>% select(1:2) %>% as_tibble()
test_h <- test_h %>% select(1:2) %>% as_tibble()
test_i <- test_i %>% select(1:2) %>% as_tibble()

hsq_b <- nhanes("HSQ_B") %>% select("SEQN", "HSD010") %>% as_tibble()
hsq_c <- nhanes("HSQ_C") %>% select("SEQN","HSD010") %>% as_tibble()
hsq_g <- nhanes("HSQ_G") %>% select("SEQN","HSD010") %>% as_tibble()
hsq_h <- nhanes("HSQ_H") %>% select("SEQN","HSD010") %>% as_tibble()
hsq_i <- nhanes("HSQ_I") %>% select("SEQN","HSD010") %>% as_tibble()


demo_b <- nhanes("DEMO_B") %>% select(SEQN, RIAGENDR, RIDAGEYR, SDMVPSU, SDMVSTRA, WTMEC2YR)

demo_c <- nhanes("DEMO_C") %>% select(SEQN, RIAGENDR, RIDAGEYR, SDMVPSU, SDMVSTRA, WTMEC2YR)

demo_g <- nhanes("DEMO_G") %>% select(SEQN, RIAGENDR, RIDAGEYR, SDMVPSU, SDMVSTRA, WTMEC2YR)

demo_h <- nhanes("DEMO_H") %>% select(SEQN, RIAGENDR, RIDAGEYR, SDMVPSU, SDMVSTRA, WTMEC2YR)

demo_i <- nhanes("DEMO_I") %>% select(SEQN, RIAGENDR, RIDAGEYR, SDMVPSU, SDMVSTRA, WTMEC2YR)

combine_data <- function(test, demo, hsq){
	
	df <- left_join(demo, test, by=c("SEQN"="SEQN")) 
	df <- left_join(df, hsq, by=c("SEQN"="SEQN"))

	colnames(df) <- c("SEQN", "Gender", "Age", "SDMVPSU", "SDMVSTRA", "WEIGHT", "TT", "HSTATUS")

	#converting ng/ml to ng/dl
	if (mean(df$TT, na.rm = TRUE) < 100){
		df$TT <- round(df$TT * 100, 2)
	}

	df <- df %>% mutate(LowT = ifelse(TT < 300, 1, 0))

	#only adults
	df <- df %>% mutate(age_indicator = (Age >= 18))

	df %>% as_tibble()
}


getSummary <- function(varformula, byformula, design){
  # Get mean, stderr, and unweighted sample size
  c <- svyby(varformula, byformula, design, unwtd.count ) 
  p <- svyby(varformula, byformula, design, svymean ) 
  outSum <- left_join(select(c,-se), p) 
  outSum
}


getQuantiles <- function(varformula, byformula, design){
  m <- svyby(varformula, byformula, design, svyquantile, quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95))
  c <- svyby(varformula, byformula, design, unwtd.count )
  m <- m %>% as_tibble() %>% select(-starts_with("se"))  %>% mutate(n = c$counts)


  return(m)
}


get_estimate <- function(test, demo, hsq){

	#get males, and self-reported excellent or very good health status adults
	data <- combine_data(test, demo, hsq) %>% mutate(inAnalysis = (Gender == "Male" & HSTATUS %in% c("Excellent", "Excellent,", "Very good,") & !is.na(TT) & age_indicator))


	NHANES_data <- svydesign(data = data, id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WEIGHT, nest = TRUE)
	NHANES <- subset(NHANES_data, inAnalysis)

	tt <- getQuantiles(~TT, ~Gender, NHANES)
	lowt <- getSummary(~LowT, ~Gender, NHANES)

	list(nhanes = (data %>% filter(inAnalysis)), tt = tt, lowt = lowt)

}


data_b <- get_estimate(test_b, demo_b, hsq_b)

data_c <- get_estimate(test_c, demo_c, hsq_c)

data_g <- get_estimate(test_g, demo_g, hsq_g) 

data_h <- get_estimate(test_h, demo_h, hsq_h)

data_i <- get_estimate(test_i, demo_i, hsq_i) 


collect_raw_data <- function(obj, cycle_years){
	obj$nhanes %>% mutate(cycle = cycle_years)
}

collate_data <- function(estimate_obj, cycle_years, type = "TT"){

	if (type == "TT"){
		estimate_obj$tt %>% mutate(type = "TT", cycle = cycle_years) %>% magrittr::set_colnames(c("gender", "whisker_low", "box_low", "median", "box_high", "whisker_high" ,"n", "type", "cycle"))
	} else{
		estimate_obj$lowt %>% mutate(type = "LowT_Rate", cycle = cycle_years) %>% magrittr::set_colnames(c("gender", "n", "value", "se", "type", "cycle")) %>% mutate(value = round(value*100, 2)) %>% mutate(se = round(se *100, 2))
	}

}



data_TT <- rbind(collate_data(data_b, "2001-2002", "TT"),
collate_data(data_c, "2003-2004", "TT"),
collate_data(data_g, "2011-2012", "TT"),
collate_data(data_h, "2013-2014", "TT"),
collate_data(data_i, "2015-2016", "TT"))

data_lowT <- rbind(collate_data(data_b, "2001-2002", "lowT"),
collate_data(data_c, "2003-2004", "lowT"),
collate_data(data_g, "2011-2012", "lowT"),
collate_data(data_h, "2013-2014", "lowT"),
collate_data(data_i, "2015-2016", "lowT"))

#median TT levels
print(data_TT)

#low T percentages in each cycle
print(data_lowT)

sizing_theme <- theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size =12),axis.title=element_text(size=14), legend.text=element_text(size=10), legend.title=element_text(size=14), plot.title=element_text(size=18, hjust=0.5))
panel_theme <- theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank())

figa <- data_TT %>% mutate(Assay = ifelse(cycle %in% c("2001-2002", "2003-2004"), "Roche Elecsys", "LC-MS/MS")) %>% mutate(Assay = factor(Assay, levels = c("Roche Elecsys", "LC-MS/MS"))) %>% mutate(cycle = factor(cycle, levels = c("2001-2002", "2003-2004", "2011-2012", "2013-2014", "2015-2016"))) %>% ggplot(., aes(x=cycle, color = Assay, fill = Assay)) + geom_boxplot(aes(ymin = whisker_low, ymax = whisker_high, middle = median, lower = box_low, upper = box_high), stat = "identity", alpha = 0.5) + theme_minimal() + ylab("Total Testosterone (ng/dl)") + xlab("Year") + panel_theme + sizing_theme + scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, 1000)) + scale_color_manual(values = c("Roche Elecsys" = "#374E5599", "LC-MS/MS" = "#Df8F4499")) + scale_fill_manual(values = c("Roche Elecsys" = "#374E5599", "LC-MS/MS" = "#Df8F4499")) + theme(legend.position = "none") 



x <- data_lowT %>% mutate(Assay = ifelse(cycle %in% c("2001-2002", "2003-2004"), "Roche Elecsys", "LC-MS/MS")) %>% mutate(Assay = factor(Assay, levels = c("Roche Elecsys", "LC-MS/MS")))

x <- x %>% mutate(cycle = factor(cycle, levels = c("2001-2002", "2003-2004", "2011-2012", "2013-2014", "2015-2016")))


figb <- x %>% filter(type == "LowT_Rate")  %>% ggplot(aes(x=factor(cycle), y=value, color = Assay, fill = Assay)) + geom_col(alpha = 0.6, position = position_dodge()) + theme_minimal() +  ylab("TT < 300 ng/dl (%)") + xlab("Year") + panel_theme + sizing_theme + scale_y_continuous(breaks = scales::pretty_breaks(n = 8), limits = c(0, 25)) + scale_color_manual(values = c("Roche Elecsys" = "#374E5599", "LC-MS/MS" = "#Df8F4499")) + scale_fill_manual(values = c("Roche Elecsys" = "#374E5599", "LC-MS/MS" = "#Df8F4499"))

cowplot::plot_grid(figa, figb, labels = c("A", "B"), label_size = 16, nrow = 1, rel_widths = c(1, 1.15)) %>% ggsave(filename = "TT-fig.pdf", plot = ., height = 5.4, width = 12, units = "in", device = cairo_pdf)
