### function called `S.obs` to calculate observed richness

S.obs <- function(x = ""){
  rowSums(x > 0) * 1
}

### function to calculate Good's Coverage: How well did you sample your site?

C <- function(x = ""){
  1 - (rowSums(x==1) / rowSums(x))
}

### Richness estimators

# function to calculate **Chao1**

S.chao1 <- function(x = ""){
  S.obs(x) + (sum(x==1)^2) / (2 * sum(x==2))
}

# function to calculate **Chao2** 

S.chao2 <- function(site = "", SbyS = ""){
  SbyS = as.data.frame(SbyS)
  x = SbyS[site,]
  SbyS.pa <- (SbyS > 0) * 1 #convert SbyS to presence/absence
  Q1 = sum(colSums(SbyS.pa) ==1) #spp. observed once
  Q2 = sum(colSums(SbyS.pa) ==2) #spp. observed twice
  S.chao2 = S.obs(x) + (Q1^2)/(2*Q2)
  return(S.chao2)
}

# function to calculate **ACE**

S.ace <- function(x = "", thresh=10){
  x <- x[x>0]                           #excludes zero-abundance taxa
  S.abund <- length(which(x > thresh)) #richness of abundant taxa
  S.rare <- length(which(x <= thresh)) #richness of rare taxa
  singlt <- length(which(x == 1)) #number of singletons
  N.rare <- sum(x[which(x <= thresh)]) #abundance of rare taxa
  C.ace <- 1 - (singlt/N.rare) #coverage (proportion non-singleton rare taxa)
  i <- c(1:thresh)            #threshold abundance range
  count <- function(i,y){     #counter to go through i range
    length(y[y==i])
  }
  a.1 <- sapply(i, count, x) #number of individuals in richness i richness classes
  f.1 <- (i * (i - 1)) * a.1 #k(k-1)kf sensu Gotelli
  G.ace <- (S.rare/C.ace)*(sum(f.1)/(N.rare*(N.rare-1)))
  S.ace <- S.abund + (S.rare/C.ace) + (singlt/C.ace) * max(G.ace,0)
  return(S.ace)
}

## SPECIES EVENNESS
# how abundance varies among species, that is, **species evenness**.

### Visualizing evenness: the rank abundance curve (RAC)

RAC <- function(x="")+{
  x.ab = x[x > 0]
  x.ab.ranked = x.ab[order(x.ab, decreasing=TRUE)]
  as.data.frame(lapply(x.ab.ranked, unlist))
  return(x.ab.ranked)
}

# function to quantify unevenness using Simpson's evenness ($E_{1/D}$) 

SimpE <- function(x = ""){
  S <- S.obs(x)
  x = as.data.frame(x)
  D <- diversity(x, "inv")
  E <- (D)/S
  return(E)
}

# function to quantify unevenness and Smith and Wilson's evenness index ($E_{var}$).


### Smith and Wilson's evenness index ($E_{var}$)

Evar <- function(x){
  x <- as.vector(x[x > 0])
  1 - (2/pi) * atan(var(log(x)))
}

### Shannon's diversity (a.k.a., Shannon's entropy)
ShanH <- function(x=""){
  H = 0
  for (n_i in x){
    if(n_i > 0){
      p = n_i / sum(x)
      H = H - p*log(p)
    }
  }
  return(H)
}
                                        
### Simpson's diversity (or dominance)

SimpD <- function(x=""){
  D = 0
  N = sum(x)
  for (n_i in x){
    D = D + (n_i^2)/(N^2)
  }
  return(D)
}             


### HILL NUMBERS 
# can calculate Hill numbers for $q$ exponent 0, 1 and 2 to interpret the effect of rare species in your community based on the response of diversity to increasing exponent $q$

# richness
# q.zero <- specnumber(x) 

# exponential shannons entropy
# q.one <- diversity(x, "shannon") 

# inverse Simpsons diversity
# q.two <- diversity(x, "invsimpson") 


profile <- function(C){
  cbind(seq(0, 7, by=0.11),
        unlist(lapply(seq(0, 7, by=0.11), function(q) sum(apply(C, 1, function(x) (x/sum(x))^q))^(1/(1-q)))))
}


# add.spec.scores function

`add.spec.scores.class` <-
  function(ordi,comm,method="cor.scores",multi=1,Rscale=F,scaling="1") {
    ordiscores <- scores(ordi,display="sites")
    n <- ncol(comm)
    p <- ncol(ordiscores)
    specscores <- array(NA,dim=c(n,p))
    rownames(specscores) <- colnames(comm)
    colnames(specscores) <- colnames(ordiscores)
    if (method == "cor.scores") {
      for (i in 1:n) {
        for (j in 1:p) {specscores[i,j] <- cor(comm[,i],ordiscores[,j],method="pearson")}
      }
    }
    if (method == "wa.scores") {specscores <- wascores(ordiscores,comm)}
    if (method == "pcoa.scores") {
      rownames(ordiscores) <- rownames(comm)
      eigenv <- ordi$eig
      accounted <- sum(eigenv)
      tot <- 2*(accounted/ordi$GOF[2])-(accounted/ordi$GOF[1])
      eigen.var <- eigenv/(nrow(comm)-1)
      neg <- length(eigenv[eigenv<0])
      pos <- length(eigenv[eigenv>0])
      tot <- tot/(nrow(comm)-1)
      eigen.percen <- 100*eigen.var/tot
      eigen.cumpercen <- cumsum(eigen.percen)
      constant <- ((nrow(comm)-1)*tot)^0.25
      ordiscores <- ordiscores * (nrow(comm)-1)^-0.5 * tot^-0.5 * constant
      p1 <- min(p, pos)
      for (i in 1:n) {
        for (j in 1:p1) {
          specscores[i,j] <- cor(comm[,i],ordiscores[,j])*sd(comm[,i])/sd(ordiscores[,j])
          if(is.na(specscores[i,j])) {specscores[i,j]<-0}
        }
      }
      if (Rscale==T && scaling=="2") {
        percen <- eigen.var/tot
        percen <- percen^0.5
        ordiscores <- sweep(ordiscores,2,percen,"/")   
        specscores <- sweep(specscores,2,percen,"*")
      }
      if (Rscale==F) {
        specscores <- specscores / constant
        ordiscores <- ordi$points
      }        
      ordi$points <- ordiscores
      ordi$eig <- eigen.var
      ordi$eig.percen <- eigen.percen
      ordi$eig.cumpercen <- eigen.cumpercen
      ordi$eigen.total <- tot
      ordi$R.constant <- constant
      ordi$Rscale <- Rscale
      ordi$scaling <- scaling
    }
    specscores <- specscores * multi    
    ordi$cproj <- specscores
    return(ordi)
  }
