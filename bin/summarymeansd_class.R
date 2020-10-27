library(Rmisc)

# Import data 
read.csv("soilvariables.csv")
soilvariables<-read.csv("soilvariables.csv", row.names = 1)

###summary(data)

pH <- summarySE(soilvariables, measurevar= "pH", groupvars=c("Site"), na.rm = TRUE)
Pdis <-summarySE(soilvariables, measurevar= "Pdis", groupvars=c("Site"), na.rm = TRUE)
Ca <-summarySE(data, measurevar= "Ca", groupvars=c("Site"), na.rm = TRUE)
Mg <-summarySE(data, measurevar= "Mg", groupvars=c("Site"), na.rm = TRUE)
K <-summarySE(data, measurevar= "K", groupvars=c("Site"), na.rm = TRUE)
Na <-summarySE(data, measurevar= "Na", groupvars=c("Site"), na.rm = TRUE)
H <-summarySE(data, measurevar= "H", groupvars=c("Site"), na.rm = TRUE)
Al <-summarySE(data, measurevar= "Al", groupvars=c("Site"), na.rm = TRUE)
SoilM <-summarySE(data, measurevar= "SoilM", groupvars=c("Site"), na.rm = TRUE)

soilvar= (data.frame(pH,Pdis,Ca,Mg,K,Na,H,Al,SoilM))
