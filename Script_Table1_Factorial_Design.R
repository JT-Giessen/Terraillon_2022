################################################################################
# Factorial-Design
# Table: Correlations 
################################################################################


{ # Startup and function definitions
    library ("SelectionTools")
    library(sqldf)
    st.set.num.threads(48)
    st.set.info.level(-2)
    st.input.dir  <- "input"
    st.output.dir <- "output"
  # Load Data to default data set
    st.read.marker.data ("rgs_marker_2.mpo",format="m")
    st.read.map         ("rgs-gerste.map",skip=1, format="mcp")
    y <- st.read.performance.data ("LP18-02-yld-a.dta")
  # Estmate effects
    st.copy.marker.data ( "C0", source.data.set="default" )
    gs.esteff.rr ( method="BLUP", data.set="C0")
  # Load parental lines to a simulation population
    pNames <- c("101", "146", "D33", "D37", "OTT", "ANT", "MER", "ETI",
           "JEN", "QUA" )
    pCode <- c("RGS251", "RGS252", "RGS253", "RGS254", "RGS255",
           "RGS256", "RGS257", "RGS258", "RGS259", "RGS260")
    st.copy.marker.data("p","C0")
    st.restrict.marker.data( ind.list=pCode, data="p")
   # Prepare simulation fuctions
    reset.all()
    st.set.sim.gp("p")
    st.set.sim.mp("p")
    p.names <- st.marker.data.statistics(data="p")$individual.list$Name
    for (n in p.names) {
        st.copy.marker.data(n,"p")
        st.restrict.marker.data( ind.list=n, data=n)
        st.set.sim.pp(n, n)
    }
factorial.cross <- function ( lines.per.family = 10)
  {   
    families <- NULL
    remove.population("DHlines")
    for (i in c(1,2,3,4,6)) for (j in c(5,7,8,9,10)) {
      cross ("F1",pCode[i],pCode[j],1)
      family.name <- sprintf("%s%s",pCode[i],pCode[j])
      families <- c(families,family.name)
      dh (family.name,"F1",lines.per.family)
      append.population("DHlines",family.name)
    }
  }
diallel.cross <- function (lines.per.family = 25)
  {  
    families <- NULL
    remove.population("DHlines")
    for (i in 6:9) for (j in ((i+1):10) ) {
      cross ("F1",pCode[i],pCode[j],1)
      family.name <- sprintf("%s%s",pCode[i],pCode[j])
      families <- c(families,family.name)
      dh (family.name,"F1",lines.per.family)
      append.population("DHlines",family.name)
    }
  }    
random.cross <- function ( lines.per.family = 10)
  {   
    families <- NULL
    remove.population("DHlines")
    lst <- NULL
    for (i in 1:9) for (j in ((i+1):10) ) lst <- rbind(lst,c(i,j))
    idx <- lst[sample(1:nrow(lst),25),]
    for (k in 1:nrow(idx)) {
        i <- idx[k,1]
        j <- idx[k,2]
        cross ("F1",pCode[i],pCode[j],1)
        family.name <- sprintf("%s%s",pCode[i],pCode[j])
        families <- c(families,family.name)
        dh (family.name,"F1",lines.per.family)
        append.population("DHlines",family.name)
    }
  }    
genotypic.values <- function()
  { # Genotypic values of the new DH lines / Training set T0
    st.get.simpop("DHlines","DHlines")
    h <- gs.predict.genotypes(training.set="C0",prediction.set="DHlines")
    write.table (h,file="input/DHlines.yld", sep=";", quote=FALSE, row.names=FALSE)
    st.read.performance.data ("DHlines.yld",data="DHlines")
    g <- h$yhat; names(g) <- h$i
    g
  }
adj.entry.means <- function(                # Simulates a field trial
                            g,              # Vector of genotypic values 
                            n.fams,         # Number of families
                            fam.size,       # Size of families 
                            var.fam.env,    # Family x env variance
                            var.gen.env,    # Line x env variance
                            var.err,        # Error variance
                            n.env,          # No of envrionments
                            n.rep)          # No pf replications
  {
    n.g    <- length(g)
    p <- g +
        rnorm (n=n.fams, mean=0, sd=sqrt(var.fam.env/n.env) ) %x% rep(1,fam.size) +
        rnorm (n=n.g, mean=0, sd=sqrt(var.err/(n.env*n.rep))) +
        rnorm (n=n.g, mean=0, sd=sqrt(var.gen.env/n.env))
    names(p) <- names(g)
    return(p)
  }
gs.plot.corr <- function (c1, title = "",sub="") 
  {
    par(mar = c(4, 4, 3, 1))
    l <- c(min(c1), max(c1))
    l <- c(65,110)
    p1 <- l[1] + 0.65 * (l[2] - l[1])
    p2 <- l[1] + 0.1 * (l[2] - l[1])
    p3 <- l[1] + 0.05 * (l[2] - l[1])
    p4 <- l[1] + 0.9 * (l[2] - l[1])
    plot(c1, ylim = l, xlim = l, col = "blue",main = title,cex.main=1)
    text(p1, p2, sprintf("r = %4.2f\nrho = %4.2f", cor(c1)[2], 
        cor(c1, method = "spearman")[2]), pos = 4)
    text(p3, p4, sprintf("%s", sub), pos = 4)
    lines(l, l)
    lines(c(mean(c1[, 1]), mean(c1[, 1])), c(min(c1[, 2]), max(c1[,2])))
    lines(c(min(c1[, 1]), max(c1[, 1])), c(mean(c1[, 2]), mean(c1[,2])))
  }
  predicted.values <- function(p)
  {
  # Save phenotypic values
    p.tab <- cbind(rownames(p),p)
    write.table (p.tab,file="input/TSPS.yld", sep=";", quote=FALSE, row.names=FALSE)
  # Estimate effects 
    st.copy.marker.data("TS","DHlines")
    st.read.performance.data ("TSPS.yld",data="TS")
    st.restrict.marker.data (ExHet.MIN=0.1, data.set="TS")
    gs.esteff.rr ( method="BLUP", data.set="TS")
  # Predict the genotypic values of the line in the training set
    j <- gs.predict.genotypes(training.set="TS", prediction.set="TS")
    h <- j$yhat; names(h) <- j$i
    h  
  }
  eval.correlations <- function (m.s,
                               n.c,
                               f.s,
                               v.ce,
                               v.le,
                               v.er,
                               n.L,
                               n.R  ){
         if ("f" == m.s)  factorial.cross( f.s )
    else if ("d" == m.s)  diallel.cross( f.s )
    else if ("r" == m.s)  random.cross( f.s )
    else { print("MÄÄH MÄÄH MÄÄH"); return() }
    g <- genotypic.values()
    p <- adj.entry.means( g,
                          n.c, f.s, v.ce, v.le, v.er,
                          n.L, n.R )
    h <- predicted.values(p)
    data.frame (
          ms = m.s , 
          nc = n.c , 
          fs = f.s , 
          vce= v.ce,
          vle= v.le,
          ver= v.er,
          nL = n.L , 
          nR = n.R ,
          vg  = var(g),
          cpg = cor(p,g),
          cph = cor(p,h),
          chg = cor(h,g))
  }
  do.it <- function(s=""){
    for (r in 1:n.rep) for (i in 1:sc) 
    {
      cat(sprintf("%s rep: %03i sc: %02i\r",s,r,i))
      cor.results <- eval.correlations(
        m.s  [i],
        n.c  [i],
        f.s  [i],
        v.ce [i],
        v.le [i],
        v.er [i],
        n.L  [i],    
        n.R  [i]
      )
    conn <- dbConnect(RSQLite::SQLite(), dbfile)
    dbWriteTable(conn, "results", cor.results, append=TRUE)
    dbDisconnect(conn)
   }
  }
}

################################################################################
# Carry out calculations
################################################################################

dbfile <- "results/Factorial_runs.sqldb"

{
    m.s  <- c ( rep("f",10), rep("f",10) )
    n.c  <- c ( rep(25, 10), rep(25, 10) )
    f.s  <- c ( rep(10, 10), rep( 6, 10) )
    v.ce <- c(1,1) %x% c ( 10 , 18, 02, 10,  20,  10 , 18, 02, 10,  20)
    v.le <- c(1,1) %x% c ( 10 , 02, 18, 10,  20,  10 , 02, 18, 10,  20)
    v.er <- c(1,1) %x% c ( 30 , 30, 30, 60, 120,  30 , 30, 30, 60, 120)
    n.L  <- c(1,1) %x% c (  1 ,  1,  1,  1,   1,   3 ,  3,  3,  3,   3)
    n.R  <- c(1,1) %x% c (  1 ,  1,  1,  1,   1,   2 ,  2,  2,  2,   2)
    n.rep <- 2  	#For test run
    #n.rep <- 1000  	#For actual run
    sc    <- length(m.s)
    do.it()
}

conn <- dbConnect(RSQLite::SQLite(), dbfile)
Q <- dbGetQuery(conn, "SELECT * FROM results")
dbDisconnect(conn)
#
res.tab <-
sqldf ("SELECT ms, nc, fs, nL , nR, vce, vle, ver , 
               AVG(vg) , AVG(cpg), AVG(cph), AVG(chg), COUNT(*)
        FROM Q
        GROUP BY nc, fs, nL , nR, vce, vle, ver")
#
names(res.tab) <- c( "m.s","n.c","f.s", "nL", "nR", "v.ce", "v.le", "v.er",
                     "v.g", "c.pg", "c.ph", "c.hg", "n.sim")
res.tab[,2:ncol(res.tab)] <- round(res.tab[,2:ncol(res.tab)],2)
srt <- c(12, 11, 14, 13, 15, 17, 16, 19, 18, 20, 2, 1, 4, 3, 5, 7, 6, 9, 8, 10)
format(res.tab,digits=2)[srt,]



