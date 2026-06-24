


# install packages if not present in R (R version >= 3.4.1)

if(!("gtools" %in% .packages(all.available=TRUE))) install.packages("gtools")
if(!("corpcor" %in% .packages(all.available=TRUE))) install.packages("corpcor")
if(!("qgraph" %in% .packages(all.available=TRUE))) install.packages("qgraph")
if(!("mgm" %in% .packages(all.available=TRUE))) install.packages("mgm")
if(!("igraph" %in% .packages(all.available=TRUE))) install.packages("igraph")


if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
	BiocManager::install("RBGL")
	}
if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
	BiocManager::install("graph")
	}

require(gtools)
require(corpcor)
require(qgraph)
require(mgm)
require(igraph)
require(RBGL)
require(graph)

if(!("pcalg" %in% .packages(all.available=TRUE))) install.packages("pcalg")
require(pcalg)


sin.ag <-
function (pcor, n, plot = TRUE, alpha = 0.1, beta = 0.5) 
{
    p <- dim(pcor)[2]
	pval.pcor <- function(pcor, n){
		pcor.unic <- pcor[lower.tri(pcor, diag=FALSE)]
		z <- 0.5*log((1+pcor.unic)/(1-pcor.unic))
		pval <- 2*(1-pnorm(sqrt(n-3)*abs(z)))
		return(p.adjust(pval, method="holm"))
	}   
	sin.lt <- pval.pcor(pcor, n)
	sin.amat <- matrix(0,ncol=p,nrow=p)
	sin.amat[lower.tri(sin.amat)] <- sin.lt
	sin.amat <- sin.amat + t(sin.amat)
	sin.p <- sin.amat
	sin.amat[sin.p < 0.1] <- 1
	sin.amat[sin.p >= 0.1] <- 0
	diag(sin.amat) <- 0
	if (plot) {
    	connect <- combinations(p, 2)
    	lc <- dim(connect)[1]
    	make.name <- function(a) paste(a[1],paste("-",a[2],sep=""),sep="")
    	leg <- apply(connect,1,make.name)
        plot(sin.lt, pch = 16, bty = "n", axes = FALSE, xlab = "edge", 
            ylab = "P-value")
        axis(1, at = 1:lc, labels = leg, las = 2)
        axis(2, at = c(0.2, 0.4, 0.6, 0.8, 1), labels = TRUE)
        lines(c(1, lc), c(alpha, alpha), col = "gray")
        lines(c(1, lc), c(beta, beta), col = "gray")
        name <- deparse(substitute(data))
        title(main = paste(attr(data, "cond")), font.main = 1)
        text(1,0.15,"0.1",col="gray")
        text(1,0.55,"0.5",col="gray")        
    }
    names(sin.lt) <- leg
    invisible(list(pval=sin.lt,amat=sin.amat))
}

plot.fci.igraph <- function(fci.out){
	amat <- fci.out@amat
	amat.dir <- amat
	amat.dir[amat.dir==1 | amat.dir==3] <- 0
	amat.dir[amat.dir==2 & t(amat.dir)==3] <- 1
	amat.dbl <- amat
	amat.dbl[amat.dbl==2 | amat.dbl==3] <- 0
	amat.tot <- amat.dir+amat.dbl
	amat.tot[amat.tot>1] <- 1
	amat.col <- amat.tot
	amat.col[amat.col==1 & t(amat.col)==0] <- 3
	graph <- graph_from_adjacency_matrix(amat.tot,weighted=TRUE)
	graph.col <- graph_from_adjacency_matrix(amat.col,weighted=TRUE)
	E(graph)$color <- ifelse(E(graph.col)$weight==3, "blue", "red") 
	l <- layout_with_fr(graph)
	l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
	plot(graph,edge.arrow.size=0.4,vertex.size=20,vertex.frame.color="#ffffff",
			vertex.color="gray80",vertex.label.color="blue",edge.curved=FALSE,rescale=FALSE,layout=l*1)
}

