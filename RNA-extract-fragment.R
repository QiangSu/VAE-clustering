library("stringr")
m=list()
s <- c("CGACCACAACATCCAGTACCAGTTCCGCACAGAGACAAATGGAGGACAGGCTGTGATCCAAAATCCCTTCAGCAATGGTGGCAGTCCGGCGGCCGAGG")

k <- str_length(s)

for (i in 1:k-49) {
  m[i] <- substr(s, i, i+49)
}

write.csv(m,"J://shenzhen Uui//project//circular RNA//code//USF2_NM_207291.3-50mer-withoutgap-0925.csv",row.names = FALSE)

