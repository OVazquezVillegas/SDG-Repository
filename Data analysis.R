# Data analysis

# Set working directory
setwd("C:/Users/6674828/OneDrive - Universiteit Utrecht/SDG/SDG-Repository")


# Load required packages
library(data.table)
library(stringr)
library(dplyr)

# Load file
NL = read.table("NL.csv", header = T, sep = ",")
string = NL$C1



# Clean data
# Eindhoven University of Technology (28874)
EI_p =  "TECH UNIV EINDHOVEN | EINDHOVEN UNIV TECHNOL"
EI = as.data.frame(str_extract(string,))
sum(is.na(EI$`str_extract(string, "TECH UNIV EINDHOVEN | EINDHOVEN UNIV TECHNOL")`)) # 953365
EI = cbind(EI, NL)
EI = EI[complete.cases(EI$`str_extract(string, "TECH UNIV EINDHOVEN | EINDHOVEN UNIV TECHNOL")`),]


# Wageningen University & Research
WA_p = c("UNIV WAGENINGEN | WAGENINGEN UNIV")
WA = as.data.frame(str_extract(string, ))
sum(is.na(WA$`str_extract(string, "TECH UNIV WANDHOVEN | WANDHOVEN UNIV TECHNOL")`)) # 953365
WA = cbind(WA, NL)
WA = WA[complete.cases(WA$`str_extract(string, "TECH UNIV WANDHOVEN | WANDHOVEN UNIV TECHNOL")`),]

