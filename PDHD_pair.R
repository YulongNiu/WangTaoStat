setwd('/home/Yulong/RESEARCH/WangTao/Rawdata/')

##################### deal with data ###################
PD <- read.csv('PD.csv')
PD <- PD[, -2]
PD[, 1] <- paste0('Name', 1:nrow(PD))
colnames(PD) <- c('PatientName', 'age', 'Time', 'HydroxyVitaminD',
                  'BloodProtein', 'AlkalinePhosphatase', 'BloodUreaNitrogen',
                  'Creatinine', 'BloodUricAcid', 'Calcium', 'Phosphorum',
                  'iPHT', 'CReactiveProtein', 'BNP')

HD <- read.csv('HD.csv')
HD <- HD[, -2]
HD[, 1] <- paste0('Name', 1:nrow(HD))
colnames(HD) <- c('PatientName', 'age', 'Time', 'HydroxyVitaminD',
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
                 Group = rep(c('腹透', '血透'), c(nrow(PD), nrow(HD))))

ggplot(vd, aes(Group, VD)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = '腹透和血透的25羟维生素D对比',
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
VDlevel <- rep('降低', PDn + HDn)
VDlevel[VDvalue < 15] <- '缺乏'
VDlevel[VDvalue > 30] <- '正常'

VDLevelMat <- data.frame(VD = VDvalue,
                         Group = rep(c('腹透', '血透'), c(PDn, HDn)),
                         Level = factor(VDlevel))

ggplot(VDLevelMat, aes(Group, VD)) +
  geom_boxplot(aes(colour = Level))+
  labs(title = '血透和腹透的25羟维生素D分组对比',
       x = '透析方法',
       y = '25羟维生素D水平') +
  scale_colour_discrete(name = '25羟维生素D水平') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(family='GB1', 'Vd_group_box.pdf')

VDPerMat <- data.frame(Level = rep(c('缺乏', '降低', '正常'), each = 2),
                       Percentage = c(lack/patients, decrease/patients, normal/patients),
                       Group = rep(c('腹透', '血透'), 3))

ggplot(VDPerMat, aes(x = Level, y = Percentage)) +
  geom_bar(aes(fill = Group), position = 'dodge', stat = 'identity') +
  labs(title = '血透和腹透的25羟维生素D比例对比',
       x = '25羟维生素D水平') +
  scale_y_continuous('比例', labels = scales::percent) +
  scale_fill_discrete(name = '透析方法') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(family='GB1', 'Vd_group_bar.pdf')
#############################################################

############################lm anlaysis######################
library('ggplot2')
library('reshape2')
library('gridExtra')

load('PDHDprocess.RData')

colName <- c('25羟维生素D', '白蛋白', '碱性磷酸酶', '肌酐', '钙', '磷', 'iPHT', 'C反应蛋白', 'BNP')

PDa <- PD[, c(-1, -2, -3, -7, -9)]
PDaObj <- lapply(2:length(colName), function(x) {
  ggplot(PDa, aes_string(colnames(PDa)[1], colnames(PDa)[x])) +
    geom_point() +
    geom_smooth(span = 0.3) +
    labs(x = colName[1],
         y = colName[x]) +
    theme(plot.title = element_text(hjust = 0.5))
    ## geom_smooth(method = 'lm', se = FALSE)
})

PDaObjPlot <- marrangeGrob(grobs = PDaObj, nrow = 4, ncol = 2)
ggsave(family='GB1', filename = 'PD_reg.pdf', plot = PDaObjPlot, width = 12, height = 10)

HDa <- HD[, c(-1, -2, -3, -7, -9)]
HDaObj <- lapply(2:length(colName), function(x) {
  ggplot(HDa, aes_string(colnames(HDa)[1], colnames(HDa)[x])) +
    geom_point() +
    geom_smooth(span = 0.3) +
    labs(x = colName[1],
         y = colName[x]) +
    theme(plot.title = element_text(hjust = 0.5))
    ## geom_smooth(method = 'lm', se = FALSE)
})

HDaObjPlot <- marrangeGrob(grobs = HDaObj, nrow = 4, ncol = 2)
ggsave(family='GB1', filename = 'HD_reg.pdf', plot = HDaObjPlot, width = 12, height = 10)

summary(lm(PDa[, 1] ~ PDa[, 2]))
#############################################################
