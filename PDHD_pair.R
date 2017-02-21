setwd('/home/Yulong/RESEARCH/WangTao/Rawdata/')

##################### deal with data ###################
PD <- read.csv('PD.csv')
PD <- PD[, -2]
PD[, 1] <- paste0('Name', 1:nrow(PD))
colnames(PD) <- c('PatientName', 'Time', 'HydroxyVitaminD',
                  'BlootProtein', 'AlkalinePhosphatase', 'BloodUreaNitrogen',
                  'Creatinine', 'BloodUricAcid', 'Calcium', 'Phosphorum',
                  'iPHT', 'CReactiveProtein', 'BNP')

HD <- read.csv('HD.csv')
HD <- HD[, -2]
HD[, 1] <- paste0('Name', 1:nrow(HD))
colnames(HD) <- c('PatientName', 'Time', 'HydroxyVitaminD',
                  'BlootProtein', 'AlkalinePhosphatase', 'BloodUreaNitrogen',
                  'Creatinine', 'BloodUricAcid', 'Calcium', 'Phosphorum',
                  'iPHT', 'CReactiveProtein', 'BNP')

save(PD, HD, file = 'PDHDprocess.RData')
########################################################

#######################
