#########
final_points7 <- generate_new_runs(c(ems_wave7, ems_wave6, ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 1000, targets, verbose=TRUE)

final_results7 <- setNames(data.frame(t(apply(final_points7, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


final_wave7 <- cbind(final_points7, final_results7)

write_xlsx(final_wave7,"C:\\Users\\LaraGosce\\Documents\\ASCENT\\Ethiopia\\Ethiopia code\\results\\calibration results\\final_wave7.xlsx")

all_points7 <- list(final_wave7)
simulator_plot(all_points7, targets)

#########
final_r7 <- data.frame()
row_to_keep <- 1:1000
for (i in 1:1000) {
  if (final_wave7[i,21]>=targets$Dtb1[1]&&final_wave7[i,21]<=targets$Dtb1[2]&&
      final_wave7[i,22]>=targets$Dtb2[1]&&final_wave7[i,22]<=targets$Dtb2[2]&&
      final_wave7[i,23]>=targets$Dtb3[1]&&final_wave7[i,23]<=targets$Dtb3[2]&&
      final_wave7[i,24]>=targets$Dtb4[1]&&final_wave7[i,24]<=targets$Dtb4[2]&&
      final_wave7[i,25]>=targets$Dtb5[1]&&final_wave7[i,25]<=targets$Dtb5[2]&&
      final_wave7[i,26]>=targets$Dtb6[1]&&final_wave7[i,26]<=targets$Dtb6[2]&&
      final_wave7[i,27]>=targets$Dtb7[1]&&final_wave7[i,27]<=targets$Dtb7[2]&&
      final_wave7[i,28]>=targets$Dtb8[1]&&final_wave7[i,28]<=targets$Dtb8[2]&&
      final_wave7[i,29]>=targets$Inc1[1]&&final_wave7[i,29]<=targets$Inc1[2]&&
      final_wave7[i,30]>=targets$Inc2[1]&&final_wave7[i,30]<=targets$Inc2[2]&&
      final_wave7[i,31]>=targets$Inc3[1]&&final_wave7[i,31]<=targets$Inc3[2]&&
      final_wave7[i,32]>=targets$Inc4[1]&&final_wave7[i,32]<=targets$Inc4[2]&&
      final_wave7[i,33]>=targets$Inc5[1]&&final_wave7[i,33]<=targets$Inc5[2]&&
      final_wave7[i,34]>=targets$Inc6[1]&&final_wave7[i,34]<=targets$Inc6[2]&&
      final_wave7[i,35]>=targets$Inc7[1]&&final_wave7[i,35]<=targets$Inc7[2]&&
      final_wave7[i,36]>=targets$Inc8[1]&&final_wave7[i,36]<=targets$Inc8[2])
  {
    #final_r[i,] <- final_wave[i,]#rbind(final_r, final_wave[i,]) #final_wave[i,]#
    # new_element <- final_wave[i,]
    # final_r[[length(final_r) + 1]] <- new_element
    row_to_keep[i] = TRUE
  } else {
    row_to_keep[i] = FALSE
  }
}
final_r7 <- final_wave7
keep <- which(row_to_keep==TRUE)
final_r72 <- final_r7[keep,]
write_xlsx(final_r72,"C:\\Users\\LaraGosce\\Documents\\ASCENT\\Ethiopia\\Ethiopia code\\results\\calibration results\\final_results7.xlsx")

all_points7 <- list(final_r72)
simulator_plot(all_points7, targets)

all_points7 <- list(wave1, wave2, wave3, wave4, wave5, final_r72) #wave0,  wave6, wave7, wave8, wave9, wave10, wave11, wave12, wave13, 
simulator_plot(all_points7, targets)