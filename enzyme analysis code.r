#Read in datafiles.  Datafile must be set up as specified for this code.
raw.dat=read.delim("DataTxtFiles/hydroenzymes.txt", header=FALSE)
weight.dat=read.delim("DataTxtFiles/sampleweights.txt")

#This code rearranges the data so that all data for each plate is on one row.
n.row.raw=nrow(raw.dat)
n.plate=n.row.raw/9
plate.dat=matrix(nrow=n.plate, ncol=108)
#The iterator i is set by the for loop command below, and keeps track of the row in the raw data.
#This iterator keeps track of what plate we are on.
k=1
#This iterator keeps track of what row we are on within a plate.
j=1
#These iterators keep track of the correct columns within a plate row to fill in.
m=1
n=12
for(i in 1:n.row.raw){
	plate.dat[k,m:n]=as.matrix(raw.dat[i,1:12])
	m=ifelse(j==9,1,n+1)
	n=ifelse(j==9,12,n+12)
	k=ifelse(j==9,k+1,k)
	j=ifelse(j==9,1,j+1)
}

plate.dat=as.data.frame(plate.dat)

# This code provides a summary of each column in each plate.
# There are 12 columns in the plate, each containing a different type of assay (blanks, sample+MUB, etc.).
# Look carefully at the values here, because it is very important to look 
# for outliers in the data that enlarge a column standard deviation.
# If you see any suspicious columns, go back to the original data and look at the well values for that column.
plate.summary=matrix(nrow=n.plate,ncol=36)
p=0
#The iterator here is going through sample plate columns (i.e., blanks, then buffer+MUB, etc.)
for(i in 0:11){
	vals=matrix(as.numeric(as.matrix(plate.dat[,c(13+i,25+i,37+i,49+i,61+i,73+i,85+i,97+i)])),nrow=n.plate)
	plate.summary[,i*3+1]=apply(vals,1,mean)
	plate.summary[,i*3+2]=apply(vals,1,sd)
	plate.summary[,i*3+3]=100*(plate.summary[,i*3+2]/plate.summary[,i*3+1])
}
plate.summary=data.frame(plate.summary)
colnames(plate.summary)=c("Blank_Ave","Blank_SD","Blank_CV","RefMUB_Ave","RefMUB_SD","RefMUB_CV","Neg_Ave","Neg_SD","Neg_CV","S1Buf_Ave","S1Buf_SD","S1Buf_CV","S1MUB_Ave","S1MUB_SD","S1MUB_CV","S1Sub_Ave","S1Sub_SD","S1Sub_CV","S2Buf_Ave","S2Buf_SD","S2Buf_CV","S2MUB_Ave","S2MUB_SD","S2MUB_CV","S2Sub_Ave","S2Sub_SD","S2Sub_CV","S3Buf_Ave","S3Buf_SD","S3Buf_CV","S3MUB_Ave","S3MUB_SD","S3MUB_CV","S3Sub_Ave","S3Sub_SD","S3Sub_CV")
round(plate.summary,digits=1)
plate.summary


#########################################
# Constants that depend on your laboratory procedures.  
# Make sure these match what you did in the lab!!!!!!!!!!!!!!

enz.const=matrix(nrow=6,ncol=4)
colnames(enz.const)=c("Enz","hrs","vol.pip","nmol")
# Enzymes tested.  Names must match what is in datafile.  Case-sensitive.
enz.const[,1]=c("nag","bglu","bxyl","bcell","phos","lap")
# Time of incubtation for each enzyme.  Time is shown in hours.
enz.const[,2]=c(2, 0.75, 2, 2,0.75,2)
# Volume of sample used in each well, in mL.
enz.const[,3]=c(0.2, 0.2, 0.2, 0.2,0.2,0.2)
# Nanomoles of substrate and MUB pipetted into assay wells and controls.
enz.const[,4]=c(0.5, 0.5, 0.5, 0.5,0.5,0.5)
edit(enz.const)

# Nanomoles of substrate and MUB pipetted into assay wells and controls.
nmol=0.5

#########################################


################# Calculations
# Calculate average of Blank wells in each plate.  
# These are the wells with only buffer (column 1 on plate).
blanks=matrix(as.numeric(as.matrix(plate.dat[,c(13,25,37,49,61,73,85,97)])),nrow=n.plate)
blank.ave=apply(blanks,1,mean)

# Calculate fluorescence that would result from complete hydrolysis of the substrate.  
# These are the wells with buffer+MUB (column 2 on plate).
MUB.only.raw=matrix(as.numeric(as.matrix(plate.dat[,c(14,26,38,50,62,74,86,98)])),nrow=n.plate)
MUB.only.raw.ave=apply(MUB.only.raw,1,mean)
# Subtract the blanks from unhydrolyzed substrate fluorescence.
MUB.only.ave=MUB.only.raw.ave-blank.ave

# Calculate fluorescence due to unhydrolyzed substrate.  
# These are the wells with buffer+substrate (column 3 on plate).
# Blank is subtracted below.
substrate.only.raw=matrix(as.numeric(as.matrix(plate.dat[,c(15,27,39,51,63,75,87,99)])),nrow=n.plate)
substrate.only.raw.ave=apply(substrate.only.raw,1,mean)
# Subtract the blanks from unhydrolyzed substrate fluorescence.
substrate.only.ave=substrate.only.raw.ave-blank.ave

# Calculate fluorescence due to sample alone (autofluorescence).  
# These are the wells with sample+buffer (columns 5, 8, 11 on plate).
# Blank is subtracted below.
sample.1.buf.raw=matrix(as.numeric(as.matrix(plate.dat[,c(16,28,40,52,64,76,88,100)])),nrow=n.plate)
sample.1.buf.raw.ave=apply(sample.1.buf.raw,1,mean)
sample.2.buf.raw=matrix(as.numeric(as.matrix(plate.dat[,c(19,31,43,55,67,79,91,103)])),nrow=n.plate)
sample.2.buf.raw.ave=apply(sample.2.buf.raw,1,mean)
sample.3.buf.raw=matrix(as.numeric(as.matrix(plate.dat[,c(22,34,46,58,70,82,94,106)])),nrow=n.plate)
sample.3.buf.raw.ave=apply(sample.3.buf.raw,1,mean)
# Subtract the blanks.
sample.1.buf.ave=sample.1.buf.raw.ave-blank.ave
sample.2.buf.ave=sample.2.buf.raw.ave-blank.ave
sample.3.buf.ave=sample.3.buf.raw.ave-blank.ave

# Calculate fluorescence of MUB with quenching due to sample.    
# These are the wells with sample+MUB (column 4, 7, 10 on plate).
# We substract off both the blank and the sample autofluorescence below.
sample.1.que.raw=matrix(as.numeric(as.matrix(plate.dat[,c(17,29,41,53,65,77,89,101)])),nrow=n.plate)
sample.1.que.raw.ave=apply(sample.1.que.raw,1,mean)
sample.2.que.raw=matrix(as.numeric(as.matrix(plate.dat[,c(20,32,44,56,68,80,92,104)])),nrow=n.plate)
sample.2.que.raw.ave=apply(sample.2.que.raw,1,mean)
sample.3.que.raw=matrix(as.numeric(as.matrix(plate.dat[,c(23,35,47,59,71,83,95,107)])),nrow=n.plate)
sample.3.que.raw.ave=apply(sample.3.que.raw,1,mean)
# Subtract the blanks and autofluorescence.
sample.1.que.ave=sample.1.que.raw.ave-blank.ave-sample.1.buf.ave
sample.2.que.ave=sample.2.que.raw.ave-blank.ave-sample.2.buf.ave
sample.3.que.ave=sample.3.que.raw.ave-blank.ave-sample.3.buf.ave

# Calculate due to hydrolysis of substrate by sample.    
# These are the wells with sample+substrate (column 6, 9, 12 on plate).
# We substract off the blank, the sample autofluorescence, and the fluorescence due to unhydrolyzed substrate below.
sample.1.act.raw=matrix(as.numeric(as.matrix(plate.dat[,c(18,30,42,54,66,78,90,102)])),nrow=n.plate)
sample.1.act.raw.ave=apply(sample.1.act.raw,1,mean)
sample.2.act.raw=matrix(as.numeric(as.matrix(plate.dat[,c(21,33,45,57,69,81,93,105)])),nrow=n.plate)
sample.2.act.raw.ave=apply(sample.2.act.raw,1,mean)
sample.3.act.raw=matrix(as.numeric(as.matrix(plate.dat[,c(24,36,48,60,72,84,96,108)])),nrow=n.plate)
sample.3.act.raw.ave=apply(sample.3.act.raw,1,mean)
# Subtract the blanks, autofluorescence, and unhydrolyzed substrate. Second, subtract substrate (negative) control.
###COLLEEN CHANGED THIS v###
sample.1.act.ave=sample.1.act.raw.ave-blank.ave-sample.1.buf.ave
sample.2.act.ave=sample.2.act.raw.ave-blank.ave-sample.2.buf.ave
sample.3.act.ave=sample.3.act.raw.ave-blank.ave-sample.3.buf.ave

# Calculate quantity of fluorescence per nmol MUB.
emission.coef=MUB.only.ave/nmol

# Calculate quench coefficient for each sample.
# This is the sample+MUB divided by MUB alone, but with sample+MUB corrected for autofluorescence as above.
sample.1.que=sample.1.que.ave/MUB.only.ave
sample.2.que=sample.2.que.ave/MUB.only.ave
sample.3.que=sample.3.que.ave/MUB.only.ave

# Calculate nmol MUB hydrolyzed per well.
###COLLEEN CHANGED THIS v###
sample.1.nmol.well=((sample.1.act.ave/sample.1.que)-substrate.only.ave)/emission.coef
sample.2.nmol.well=((sample.2.act.ave/sample.2.que)-substrate.only.ave)/emission.coef
sample.3.nmol.well=((sample.3.act.ave/sample.3.que)-substrate.only.ave)/emission.coef

# Calculate nmol/hr/mL assay.
sample.nmol.mL=matrix(nrow=n.plate*3, ncol=3)
for(i in 1:n.plate){
	sample.nmol.mL[(3*(i-1)+1):(3*(i-1)+3),1]=enz.const[plate.dat[i,1]==enz.const[,1],1]
	sample.nmol.mL[3*(i-1)+1,2]=as.character(plate.dat[i,2])
	sample.nmol.mL[3*(i-1)+2,2]=as.character(plate.dat[i,3])
	sample.nmol.mL[3*(i-1)+3,2]=as.character(plate.dat[i,4])
	hrs.i=as.numeric(enz.const[plate.dat[i,1]==enz.const[,1],2])
	vol.pip.i=as.numeric(enz.const[plate.dat[i,1]==enz.const[,1],3])
	sample.nmol.mL[3*(i-1)+1,3]=sample.1.nmol.well[i]/(hrs.i*vol.pip.i)
	sample.nmol.mL[3*(i-1)+2,3]=sample.2.nmol.well[i]/(hrs.i*vol.pip.i)
	sample.nmol.mL[3*(i-1)+3,3]=sample.3.nmol.well[i]/(hrs.i*vol.pip.i)
}
colnames(sample.nmol.mL)=c("Enzyme","Sample","nmol_per_ml")
sample.nmol.mL

# Calculate nmol/hr/g sample (in column called "activity").
sample.nmol.g=cbind(sample.nmol.mL,weight.dat[match(sample.nmol.mL[,2],weight.dat$SampleID),])
sample.nmol.g=na.omit(sample.nmol.g)
sample.nmol.g$activity=as.numeric(as.matrix(sample.nmol.g$nmol_per_ml))*sample.nmol.g$volume/as.numeric(as.character(sample.nmol.g$enz.wt))
sample.nmol.g

write.csv(sample.nmol.g,"DataTxtFiles/enzymeactivity.csv",row.names=FALSE)
