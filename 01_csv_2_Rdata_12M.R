#############
# file IO
#############

lumen_file_1 = "./data_12M/12M_Lumen_01.csv"
lumen_file_2 = "./data_12M/12M_Lumen_02.csv"
lumen_file_3 = "./data_12M/12M_Lumen_03.csv"
lumen_file_4 = "./data_12M/12M_Lumen_04.csv"

mucosa_file_1 = "./data_12M/12M_Mucosa_01.csv"
mucosa_file_2 = "./data_12M/12M_Mucosa_02.csv"
mucosa_file_3 = "./data_12M/12M_Mucosa_03.csv"
mucosa_file_4 = "./data_12M/12M_Mucosa_04.csv"

blood_file = "./data_12M/12M_Blood.csv"
rectum_file = "./data_12M/12M_Rectum.csv"
feces_file = "./data_12M/12M_Feces.csv"

lumen_T1 <- read.csv(lumen_file_1)
save(lumen_T1, file="./Rdata_12M/12M_lumen_T1.Rdata")
lumen_T2 <- read.csv(lumen_file_2)
save(lumen_T2, file="./Rdata_12M/12M_lumen_T2.Rdata")
lumen_T3 <- read.csv(lumen_file_3)
save(lumen_T3, file="./Rdata_12M/12M_lumen_T3.Rdata")
lumen_T4 <- read.csv(lumen_file_4)
save(lumen_T4, file="./Rdata_12M/12M_lumen_T4.Rdata")

mucosa_T1 <- read.csv(mucosa_file_1)
save(mucosa_T1, file="./Rdata_12M/12M_mucosa_T1.Rdata")
mucosa_T2 <- read.csv(mucosa_file_2)
save(mucosa_T2, file="./Rdata_12M/12M_mucosa_T2.Rdata")
mucosa_T3 <- read.csv(mucosa_file_3)
save(mucosa_T3, file="./Rdata_12M/12M_mucosa_T3.Rdata")
mucosa_T4 <- read.csv(mucosa_file_4)
save(mucosa_T4, file="./Rdata_12M/12M_mucosa_T4.Rdata")

blood <- read.csv(blood_file)
save(blood, file="./Rdata_12M/12M_blood.Rdata")

rectum <- read.csv(rectum_file)
save(rectum, file="./Rdata_12M/12M_rectum.Rdata")

feces <- read.csv(feces_file)
save(feces, file="./Rdata_12M/12M_feces.Rdata")
