library(ggplot2)
library(tidyr)
###plotting the time difference
#graph with number of nodes at bottom
#3 lines:
#	time of MRPC
#	time of MRPC with FDR correction
#	time of MRdualPC
time_MRPC <- c(3982, 1049, 3454, 6366, 8351, 4971, 5455, 642, 2434, 5911, 7651, 5722, 7588, 5215, 919, 7089, 1961, 2112, 6058, 1969, 1944, 1886, 818, 9362, 468, 3428, 2354, 257, 4178, 2438, 2455, 243, 868, 2824, 7641, 1653, 3483, 3189, 767, 4703, 5058, 1141, 2507, 3801, 1308, 661, 2159, 3551, 2740, 2631, 5340, 4899, 1018, 2393, 4443, 6159, 1504, 787, 4264, 1634, 1967, 3086, 977, 3248, 1894, 843, 1510, 5002, 2843, 1928, 982, 2946, 5297, 3070, 9.7, 8167, 7624, 1829, 5452, 5573, 5468, 4000, 2704, 806, 1332, 1973, 2159, 3259, 4948, 1044, 3446.46, 5373, 4060, 367, 1502, 490, 4302, 763, 800, 1293, 536, 5404, 6037, 1169, 107, 650, 986, 3460, 4641)
time_dual <- c(43.93, 8.54, 23.90, 54.69, 109.47, 68.85, 52.45, 7.66, 21.59, 42.52, 144.98, 42.29, 30.65, 28.62, 5.27, 45.73, 16.20, 24.13, 82.51, 19.20, 20.87, 17.59, 4.69, 83.92, 6.94, 20.90, 23.61, 7.48, 30.22, 39.04, 17.21, 5.66, 10.16, 30.88, 37.06, 13.78, 35.59, 16.74, 12.56, 19.08, 21.23, 12.51, 21.34, 27.31, 14.26, 9.60, 24.56, 21.65, 17.01, 23.62, 53.54, 33.83, 10.13, 29.81, 46.30, 64.66, 13.55, 10.72, 59.4, 19.60, 20.22, 27.47, 13.38, 30.11, 19.39, 6.81, 10.54, 71.88, 21.51, 15.42, 10.19, 24.67, 41.98, 20.71, 5.78, 83.53, 97.11, 22.81, 38.77, 155.74, 38.70, 60.47, 21.91, 13.76, 20.17, 13.09, 30.08, 29.60, 112.31, 5.92, 29.95, 59.18, 18.64, 6.69, 17.46, 7.03, 23.30, 12.02, 9.32, 8.67, 5.56, 38.99, 36.78, 10.26, 14.62, 6.14, 15.55, 20.06, 34.19)
size <- c(767, 499, 689, 864, 944, 888, 880, 512, 710, 863, 919, 839, 819, 811, 432, 908, 575, 690, 989, 630, 616, 578, 443, 958, 390, 721, 676, 455, 866, 840, 700, 400, 541, 738, 887, 742, 899, 674, 564, 801, 793, 633, 701, 867, 604, 556, 801, 822, 661, 747, 985, 786, 517, 646, 835, 927, 602, 521, 823, 675, 573, 757, 556, 707, 592, 470, 551, 906, 716, 628, 525, 667, 790, 694, 139, 903, 974, 603, 868, 884, 813, 794, 707, 549, 667, 628, 684, 684, 986, 506, 812, 807, 729, 489, 655, 526, 792, 568, 546, 511, 480, 769, 943, 533, 308, 484, 646, 778, 833)
#run MRPC on larger than 1000 to see which work with MRPC (and to get time for dualMRPC that will run)


times_fast <- (time_MRPC/time_dual)
data <- data.frame(size, time_MRPC, time_dual, times_fast)

#use quadratic term
theme_set(theme_bw())
data %>%
  ggplot( aes(x=size, color = grp)) +
    geom_point(aes(y = log(time_MRPC)), color = "#ffa600") +
    geom_point(aes(y = log(time_dual)), color = "#003f5c") +
    geom_smooth(aes(y = log(time_MRPC), color = "#ffa600"), method = "lm", se = F, formula = y ~ x + I(x^2)) +  #this line is included for the legend
    geom_smooth(aes(y = log(time_MRPC)), method = "lm", se = T, color = "#ffa600", formula = y ~ x + I(x^2)) +
    geom_smooth(aes(y = log(time_dual), color = "#003f5c"), method = "lm", se = F, formula = y ~ x + I(x^2)) +
    geom_smooth(aes(y = log(time_dual)), method = "lm", color = "#003f5c", se = T, formula = y ~ x + I(x^2)) +
    labs(x = "Nodes in module", y = "log(seconds)", color = "Legend") +
    scale_colour_manual(name="Method", values=c("MRPC" = "#ffa600", "MRdualPC" = "#003f5c")) +
    ggtitle("Computation time as a function of number of nodes") +
    theme(text = element_text(size=15))

data %>%
  ggplot( aes(x=size)) +
    geom_point(aes(y = log(times_fast)), color = "#003f5c") +
    geom_smooth(aes(y = log(times_fast)), method = "lm", se = T, color = "#003f5c", formula = y ~ x + I(x^2)) +
    xlab("Nodes in module") + ylab("log(times faster)") +
    ggtitle("Log(Times faster) MRdualPC is than MRPC for varying module size")


time_extra_dual <- c(115.13, 43.98, 73.99, 55.47, 309.11, 109.70, 107.97, 87.42, 299.18, 63.96, 342.62, 175.26, 1695.36, 1111.58, 290.59, 123.32, 4586.36, 605.49, 85.66, 179.14, 539.10, 390.68, 929.65, 192.61, 95.87, 82.32, 507.81, 1113.38, 554.66, 986.58, 77.08, 842.59, 148.92, 170.74, 778.85, 4445.37, 51.01, 269.84)
extra_size <- c(993, 972, 924, 973, 1288, 1129, 1072, 1011, 1190, 1010, 1230, 1305, 1727, 1537, 1121, 1032, 1967, 1229, 1041, 1257, 1464, 1334, 1594, 1166, 1061, 1011, 1492, 1610, 1256, 1413, 1081, 1494, 1047, 1161, 1488, 1924, 1006, 1028)
time_extra_MRPC <- c(time_MRPC, rep(NA, length(time_extra_dual)))
time_extra_dual <- c(time_dual, time_extra_dual)
extra_size <- c(size, extra_size)

data1 <- data.frame(extra_size, time_extra_MRPC, time_extra_dual)
data1 %>%
  ggplot( aes(x=extra_size, color = grp)) +
    geom_point(aes(y = log(time_extra_MRPC)), color = "#ffa600") +
    geom_point(aes(y = log(time_extra_dual)), color = "#003f5c") +
    geom_smooth(aes(y = log(time_extra_MRPC), color = "#ffa600"), method = "lm", se = F, formula = y ~ x + I(x^2)) +  #this line is included for the legend
    geom_smooth(aes(y = log(time_extra_MRPC)), method = "lm", se = T, color = "#ffa600", formula = y ~ x + I(x^2)) +
    geom_smooth(aes(y = log(time_extra_dual), color = "#003f5c"), method = "lm", se = F, formula = y ~ x + I(x^2)) +
    geom_smooth(aes(y = log(time_extra_dual)), method = "lm", color = "#003f5c", se = T, formula = y ~ x + I(x^2)) +
    labs(x = "Nodes in module", y = "log(seconds)", color = "Legend") +
    scale_colour_manual(name="Method", values=c("MRPC" = "#ffa600", "MRdualPC" = "#003f5c")) +
    ggtitle("Speed of MRPC vs MRdualPC")
