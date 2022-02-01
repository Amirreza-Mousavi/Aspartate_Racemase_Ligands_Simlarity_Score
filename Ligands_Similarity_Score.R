###Initialization
library(Rcpi)
Jaccardindex=function(A,B){
  int=table(A&B)[2]
  return(as.numeric(int/(length(A)+length(B)-int)))
  
}

###Load molecules
x=readMolFromSmi("Smile Files.smi")

###MACCS166 calculation
y=extractDrugMACCSComplete(x)

###Compute similarity matrix
m=matrix(0,131,131)

for(i in 1:131){
  for(j in 1:131){
   m[i,j]=Jaccardindex(y[i,],y[j,])
  }}

###Print the similarity matrix
print(m)

###Cluster the ligands based on their similarity
###1/m is a measure of dissimilarity(to be updated later)
pdf("Result.pdf",width = 16,height = 9,paper = "special")
plot(hclust(as.dist(1/m)))
dev.off()
