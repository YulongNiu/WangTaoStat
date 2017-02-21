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

######################VD comparison#####################
library('ggplot2')

load('PDHDprocess.RData')

vd <- data.frame(VD = c(PD[, 'HydroxyVitaminD'],
                        HD[, 'HydroxyVitaminD']),
                 Group = rep(c('PD', 'HD'), c(nrow(PD), nrow(HD))))

ggplot(vd, aes(Group, VD)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = '血透和腹透的25羟维生素D对比',
       x = '透析方法',
       y = '25羟维生素D水平') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(family='GB1', 'Vd_average.pdf')

t.test(PD[, 'HydroxyVitaminD'], HD[, 'HydroxyVitaminD'])
wilcox.test(PD[, 'HydroxyVitaminD'], HD[, 'HydroxyVitaminD'])
########################################################


#########################VD proportion test##################
library('ggplot2')
library('devtools')

source_url('https://raw.githubusercontent.com/YulongNiu/FunFunc/master/test_two_proportion.R')

load('PDHDprocess.RData')

PDvd <- PD[, 'HydroxyVitaminD']
PDn <- length(PDvd)

HDvd <- HD[, 'HydroxyVitaminD']
HDn <- length(HDvd)

## lack <15nm/ml
## decrease 15ng/ml ~ 30ng/ml
## normal >30ng/ml
lack <- c(sum(PDvd < 15),
          sum(HDvd < 15))

decrease <- c(sum(PDvd >= 15 & PDvd <= 30),
              sum(HDvd >= 15 & HDvd <= 30))

normal <- c(sum(PDvd > 30),
            sum(HDvd > 30))

patients <- c(PDn, HDn)

## R
prop.test(lack, patients)
prop.test(decrease, patients)
prop.test(normal, patients)

## pooled
TestPropor(lack[1], patients[1], lack[2], patients[2])
TestPropor(decrease[1], patients[1], decrease[2], patients[2])
TestPropor(normal[1], patients[1], normal[2], patients[2])

## binomial
TestPropor(lack[1], patients[1], lack[2], patients[2], method = 'binomial')
TestPropor(decrease[1], patients[1], decrease[2], patients[2], method = 'binomial')
TestPropor(normal[1], patients[1], normal[2], patients[2], method = 'binomial')

VDvalue <- c(PDvd, HDvd)
VDlevel <- rep('decrease', PDn + HDn)
VDlevel[VDvalue < 15] <- 'lack'
VDlevel[VDvalue > 30] <- 'normal'

VDLevelMat <- data.frame(VD = VDvalue,
                         Group = rep(c('PD', 'HD'), c(PDn, HDn)),
                         Level = VDlevel)
#############################################################
