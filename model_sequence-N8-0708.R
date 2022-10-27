library(stringr)

m1=list()
m2=list()
m3=list()
m4=list()
m5=list()
m6=list()
m7=list()
m8=list()

k <- c("A","T","C","G")
for (i in 1:4) {
  a <- c("GGTCGCTGCGNGTCNCACNAGCNCCGNGCGNGCCNGCANGCCTGCTCTGC")
  str_sub(a,11,11) <- k[i]
  m1[[i]] <- a
}
#print(m1)
for (i in 1:4) {
  str_sub(m1,15,15) <- k[i]
  m2[[i]] <- m1
}
#print(m2)

l <- unlist(m2)


for (i in 1:4) {
  str_sub(l,19,19) <- k[i]
  m3[[i]] <- l
}
#print(m3)

l <- unlist(m3)


for (i in 1:4) {
  str_sub(l,23,23) <- k[i]
  m4[[i]] <- l
}
#print(m4)

l <- unlist(m4)
for (i in 1:4) {
  str_sub(l,27,27) <- k[i]
  m5[[i]] <- l
}
#print(m5)
l <- unlist(m5)
for (i in 1:4) {
  str_sub(l,31,31) <- k[i]
  m6[[i]] <- l
}
#print(m6)

l <- unlist(m6)
for (i in 1:4) {
  str_sub(l,35,35) <- k[i]
  m7[[i]] <- l
}
#print(m7)

l <- unlist(m7)
for (i in 1:4) {
  str_sub(l,39,39) <- k[i]
  m8[[i]] <- l
}
#print(m8)


write.csv(m8,"J://shenzhen Uui//project//circular RNA//code//N8-GC1-sequence-0708.csv",row.names = FALSE)














