#https://github.com/enricogiudice/dualPC/blob/main/dualPC.R


# This function inverts a matrix and puts the inverse in correlation form
psolve <- function(M) {
  M1 <- solve(M, tol = 5e-26) # should use pseudoinverse
  scale <- 1/sqrt(diag(M1))
  t(M1*scale)*scale 
}


# this function generates all subsets
combinations <- function (n, r, v = 1:n) {
  v0 <- vector(mode(v), 0)
  if (r == 0) 
    v0
  else if (r == 1) 
    matrix(v, n, 1)
  else if (r == n) 
    matrix(v, 1, n)
  else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])), 
             Recall(n - 1, r, v[-1]))
}

# this function just generates the next subset
nextSubS <- function(SubS, max) {
  k <- length(SubS)
  if (SubS[k] == max) { # need to shift
    rmost_vec <- which(SubS != 1:k + max - k)
    if (length(rmost_vec) > 0) {
      change_point <- max(which(SubS != 1:k + max - k))
      SubS[change_point:k] <- SubS[change_point] + 1:(k - change_point + 1)
    } else {
      warning("no more sets!")
    }
  } else {
    SubS[k] <- SubS[k] + 1
  }
  SubS
}

# This function performs the Fisher z-test for conditional independence
#info here https://online.stat.psu.edu/stat505/lesson/6/6.3
Fztest <- function(x, Nmk, T_star) {
  T <- abs(sqrt(Nmk - 3) * 0.5 * pcalg:::log.q1pm(x))
  #if (is.na(T)){ #added in this code so that if is NA doesn't give error and instead adds 0 (as is done with zStat)
  #  T <- 0
  #}
  T[is.na(T)] = 0 #replaced above code with this as T can be more than single item
  (T > T_star)
  #(2 * pnorm(abs(T), lower.tail = FALSE)) < 0.05 #alternative to finding sig. pvalues
}

# This function finds the nodes connected to x
Find_Nbhd <- function(G, x, y) {
  if (G[x, y] == 1) { # if there is an edge to test
    edge_vec <- G[x, ]
    edge_vec[c(x, y)] <- 0 # don't want these
    which(edge_vec == 1)
  } else { # otherwise
    c()
  }
}

# Compute the conditional correlation coefficient
Compute_rho <- function(M) {
  chol_mat <- chol(M[-c(1:2), -c(1:2), drop = FALSE]) # Cholesky decomposition
  x_mat <- backsolve(chol_mat, M[-c(1:2), 1:2, drop = FALSE], transpose = TRUE)
  local_mat <- M[1:2, 1:2] - t(x_mat) %*% x_mat
  rho <- local_mat[1, 2]/sqrt(local_mat[1, 1]*local_mat[2, 2])
  if(is.na(rho)){
    rho <- 0
  }
  rho
}



dual_pc <- function(dat, cor_mat, N, alpha, ord_ind = TRUE, skeleton = FALSE, pattern_graph = FALSE, max_ord = NULL) {
  n <- ncol(cor_mat) # number of variables
  if (is.null(max_ord)) { # the maximum subset size to test
    max_ord <- n
  }
  if (length(alpha) == 1) {
    T_c <- abs(qnorm(alpha/2)) # p-value cut-off in z-space for normal tests
    T_p <- T_c # p-value cut-off in z-space for dual tests
  } else { # in case we want different thresholds for normal and dual tests
    T_c <- abs(qnorm(alpha[1]/2))
    T_p <- abs(qnorm(alpha[2]/2))
  }
  c_mat <- cor_mat # local correlation matrix 
  ord <- 0 # size of conditioning sets, 0th level is correlation/precision matrix
  # T statistics on correlation space
  # keep edges which are significant, ie delete those which aren't
  #Gc_0 <- Fztest(cor_mat*upper.tri(c_mat), N, T_c)
  #instead of using fishers z tranformation use cor.test
  x <- psych::corr.test(suffStat_C1$C, adjust = "none") #corr.test runs cor.test on matrix
  Gc_0 <- x$p <= 0.05
  
  # precision matrix
  #p_mat <- psolve(c_mat)
  #this returns errors - try using partial correlation
  ##GeneNet 
  p_mat <- pcor(dat)
  p_mat <- p_mat$estimate
  #p_mat <- ggm.estimate.pcor(dat)
  pcor = sm2vec(p_mat)
  indexes = sm.index(p_mat)
  colnames(indexes) = c("node1", "node2")
  w = cbind(pcor, indexes)
  #cat("Estimate (local) false discovery rates (partial correlations):\n")  
  fdr.out = fdrtool(w[, 1], statistic = "correlation", cutoff.method = "locfdr")
  pval = fdr.out$pval
  qval = fdr.out$qval
  prob = 1 - fdr.out$lfdr
  result = cbind(w, pval, qval, prob)
  sort.idx = order(-abs(result[, 1]))
  result = as.data.frame(result[sort.idx, ])
  num.significant.1 <- sum(result$pval <= 0.05)
  sig_node <- result[1:num.significant.1,]
  temp <- matrix(ncol=nrow(p_mat),nrow=nrow(p_mat))
  temp[is.na(temp)] <- FALSE
  for(i in 1:nrow(sig_node)){
    info <- sig_node[i,]
    pos1 <- info$node1
    pos2 <- info$node2
    temp[pos1, pos2] <- TRUE
  }
  colnames(temp) <- colnames(p_mat)
  rownames(temp) <- rownames(p_mat)
  Gp_0 <- temp  

  # keep track of sepsets of precision matrix
  pres_sepsets <- Gp_0 | t(Gp_0)
  # combined 
  G_0 <- (Gc_0 & Gp_0)*upper.tri(Gc_0)
  # we keep track of the skeleton in a symmetric matrix
  G_cur <- 1*(G_0 | t(G_0)) # current graph
  if (ord_ind) { # to track the edges we need to delete
    del_mat <- matrix(1, n, n)
  }
  seq_p <- seq_len(n)
  sepsets <- lapply(seq_p, function(.) vector("list", n))
  done_flag <- FALSE # track when we have performed all tests needed
  while (ord < max_ord && done_flag == FALSE) { # keep track of order
    done_flag <- TRUE
    ord <- ord + 1
    edges <- which(G_cur == 1, arr.ind = TRUE)
    for (ii in 1:nrow(edges)) {
      x <- edges[ii, 1]
      y <- edges[ii, 2]
      S <- Find_Nbhd(G_cur, x, y) #taking edges that are also sig. for x that aren't y
      nbhd_size <- length(S)
      if (ord <= nbhd_size) { # only need dual tests up to when ord is half nbhd_size
        # we could track which larger tests have already been performed and skip those
        #c_mat <- cor_mat[c(x, y, S), c(x, y, S)] # local correlation matrix
        #p_mat <- psolve(c_mat) # local precision matrix
        if(((length(S)+3) >= N) == TRUE){ #this is +3 as n-N-3 for degrees of freedom (ensures doesn't hit 0 or below so can sqrt())
          #p_mat <- ggm.estimate.pcor(dat[, c(x, y, S)], verbose = FALSE)
          #test.edges <- network.test.edges(p_mat, verbose = FALSE) #calculates edges that are sig. higher than 0
          p_mat <- pcor(dat[, c(x, y, S)])
          p_mat <- p_mat$estimate
          pcor = sm2vec(p_mat)
          indexes = sm.index(p_mat)
          colnames(indexes) = c("node1", "node2")
          w = cbind(pcor, indexes)
          cat("Estimate (local) false discovery rates (partial correlations):\n")
          fdr.out = fdrtool(w[, 1], statistic = "correlation", cutoff.method = "locfdr")
          pval = fdr.out$pval
          qval = fdr.out$qval
          prob = 1 - fdr.out$lfdr
          result = cbind(w, pval, qval, prob)
          sort.idx = order(-abs(result[, 1]))
          test.edges = as.data.frame(result[sort.idx, ])
          #in test.edges, the posistion of nodes are how given in p_mat
          #ie. node x is in pos 1, node y is pos 2 and S are following this
          test.edges1 <- test.edges[test.edges$node1 == 1, ]
          test.edges1 <- test.edges1[test.edges1$node2 == 2, ] 
          #test to see if they are swapped (will not be, but used as a failsafe for some buggy behaviour)
          test.edges2 <- test.edges[test.edges$node1 == 2, ] 
          test.edges2 <- test.edges2[test.edges2$node2 == 1, ]
          test.results <- rbind(test.edges1, test.edges2)
          test_flag <- test.results$pval <= 0.05
          if(nrow(test.results) == 0){
            test_flag <- FALSE
          }
        } else {
          c_mat <- pcor(dat[, c(x, y, S)])
          c_mat <- c_mat$estimate
          #c_mat <- ggm.estimate.pcor(dat[, c(x, y, S)], verbose = FALSE)
          colnames(c_mat) <- rownames(c_mat) <- colnames(dat[, c(x, y, S)])
          p_mat <- c_mat*(-1)
          diag(p_mat) <- 1
          c_mat <- cor_mat[c(x, y, S), c(x, y, S)] # local correlation matrix
          test_flag <- Fztest(p_mat[1, 2], N-nbhd_size, T_p)
        }
        if (ord == nbhd_size) {
          n_subsets <- 0 # nothing else to test
        } else {
          # only need subsets up to half the size - the dual part takes care of the others!
          subset_size <- min(ord, nbhd_size - ord)
          n_subsets <- choose(nbhd_size, subset_size)
          if (ord < nbhd_size - 1) { # we need another round
            done_flag <- FALSE
          }
        }
        jj <- 0
        while (jj < n_subsets && test_flag == TRUE) {
          jj <- jj + 1 # loop
          if (jj == 1) {
            subset <- 1:subset_size
          } else {
            subset <- nextSubS(subset, nbhd_size)
          }
          cond_set <- subset + 2 # which rows/columns to condition on
          if (ord <= nbhd_size/2) {
            # normal test
            rho <- Compute_rho(c_mat[c(1, 2, cond_set), c(1, 2, cond_set)])
            test_flag <- Fztest(rho, N-ord, T_c)
          } else {
            # normal test, but use dual space to compute more efficiently
            rho <- Compute_rho(p_mat[c(1, 2, cond_set), c(1, 2, cond_set)])
            test_flag <- Fztest(rho, N-ord, T_c)
          }
          if (test_flag == FALSE) {  # record separating subset
            sepsets[[x]][[y]] <- sepsets[[y]][[x]] <- S[subset]
          }
          if (ord < nbhd_size/2 && test_flag == TRUE) {
            # dual test as well
            rho <- Compute_rho(p_mat[c(1, 2, cond_set), c(1, 2, cond_set)])
            test_flag <- Fztest(rho, N-nbhd_size+ord, T_p)
          }
          if (test_flag == FALSE) {  # record separating subset
            sepsets[[x]][[y]] <- sepsets[[y]][[x]] <- S[-subset]
          }
        }
        if (test_flag == FALSE) { # a test was rejected
          if (ord_ind) { # we only delete at the end of each main loop
            del_mat[x, y] <- 0
            del_mat[y, x] <- 0
          } else {
            G_cur[x, y] <- 0 # delete edges
            G_cur[y, x] <- 0
          }
          if (jj == 0) { # record separating subsets
            sepsets[[x]][[y]] <- sepsets[[y]][[x]] <- S
          }
        }
      }
    }
    if (ord_ind) { # we delete edges now instead
      G_cur <- 1*(G_cur & del_mat)
    }
  }
  if (skeleton) {
    new("MRPCclass", graph = as(G_cur, "graphNEL"), sepset = sepsets)
  } else {
    orient_vstructures(G_cur, sepsets, pres_sepsets, pattern_graph)
  }
  
}


EdgeOrientation1 <- function (gInput, GV, suffStat, FDR, alpha, indepTest,
                             FDRcontrol, tau = 0.5, lambda = 0.25,
                             verbose = FALSE) {
  
  g <- as(gInput@graph, "matrix") # g ia an adjacency from undirected graph (skleton)
  g1 <- g
  p <- nrow(g)
  # tarmat (adjacency matrix for directed graph) is updated in every step, 
  # and contains the output of final topology
  tarmat <- matrix(0, nrow(g),ncol(g)) #same row and column from g
  rownames(tarmat) <- rownames(g)     #same row names from g
  colnames(tarmat) <- colnames(g)     #same column names from g
  # extract all edges
  g[lower.tri(g)] <- 0  # Use only upper triangular because g is symmetric matrix
  edges <- which (g==1, arr.ind = TRUE)
  
  #Step-1 start
  if (GV>0) {
    # identify edges involving Vs
    edgesWithBothVs <- edges[which (edges[,1] <= GV & edges[,2] <= GV), ]
    edgesWithVFirst <- edges[which (edges[,1] <= GV & edges[,2] > GV), ]
    edgesWithVSecond <- edges[which (edges[,1] > GV & edges[,2] <= GV), ]
    
    # assign 1s to corresponding edges in tarmat
    # edges between two Vs are undirected (or bidirected)
    if (length (edgesWithBothVs) > 2) {
      tarmat[edgesWithBothVs] <- 1
      tarmat[edgesWithBothVs[,2:1]] <- 1
    } else {
      tarmat[edgesWithBothVs[1], edgesWithBothVs[2]] <- 1
      tarmat[edgesWithBothVs[2], edgesWithBothVs[1]] <- 1
    }
    # Edges involving one V go from V to the other node
    # Need to distinguish whether edgesWithVFirst or edgesWithVSecond is 
    # a matrix or a vector.
    if (length (edgesWithVFirst) > 2) {
      tarmat[edgesWithVFirst] <- 1      
    } else {
      tarmat[edgesWithVFirst[1], edgesWithVFirst[2]] <- 1
    }
    
    if (length(edgesWithVSecond) > 2) {
      tarmat[edgesWithVSecond[,2:1]] <- 1
    } else {
      tarmat[edgesWithVSecond[2], edgesWithVSecond[1]] <- 1
    }
  }
  #Step-1 end
  #Step-2 start
  #Start to orient v-structures
  ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
  V <- colnames(g)

  for (i in seq_len(nrow(ind))) {
    x <- ind[i, 1]
    y <- ind[i, 2]
    allZ <- setdiff(which(g1[y, ] == 1), x)
    for (z in allZ) {
      # Triplet x-y-z is directed x-->y<--z if x and z conditionally dependent given y
      if ((g1[x, z] == 0 & g1[x, y] == 1 & g1[y, z] == 1)  & !(tarmat[y, x] ==1) & !(tarmat[z, y] ==1) & !(tarmat[y, z] ==1) & 
          !(y %in% gInput@sepset[[x]][[z]] || y %in% gInput@sepset[[z]][[x]])) 
      {       
        if(indepTest=="gaussCItest") #if indepTest=gaussCItest
        {
          pval <- gaussCItest(x, z, y, suffStat)
        }
        if(indepTest=="disCItest") #if indepTest=disCItest
        {
          pval <- disCItest(x, z, y, suffStat) #additional
        }
        #pval=disCItest(x, z, y, suffStat) #additional conditional test
        Alpha <- alpha
        if (pval <= Alpha) {  #Reject H0 (H0:nodes are independent)
          tarmat[x, y] <- tarmat[z, y] <- 1 #directed x-->y<--z
        }       
      }
    }
  }
  #Step-2 end to orient v-structures
  
  #Step-3
  #Orient remaining edges, weather gMR is applicable or not based on the following process:
  #Pullout all edges from g1 matrix and define 1st and 2nd node (goal to make a triplet for the remaining edges)
  #extract all edges with involving those 1st and 2nd nodes (from row and column of g)
  #ignore the nodes that already have direction
  #form a triplet
  #determine the direction using test results from part 1 based on gMR
  #Repeat until all undirected edges to directed
  # m <- m
  # Alpha <- Alpha
  # R <- R
  #start 
  #when data contain genetic variants
  if (any(tarmat == 1) & GV>0) #if at least one edge directed already made so far
  {
    WW1 <- unique(which(tarmat==1,arr.ind = T)[,2]) #pullout the all canidate genes for v-structure
    #if(WW1<GV){
    #WW1=WW1[-c(1:GV)] #ignor canidate genes for GV
    WW1 <- setdiff(WW1,1:GV)
    #}
    if(length(WW1)!=0)
    {
      #WW1=setdiff(WW1,1:GV)
      for (v1 in 1:length(WW1))
      {
        WW2 <- unique(which(tarmat[,WW1[v1]]==1,arr.ind = T)) #edge between v-structure  
        WW3 <- unique(which(g1[,WW1[v1]]==1,arr.ind = T)) #edges between canidate genes and others  
        WW3 <- setdiff(WW3,WW2)   #ignore if already direction
        if(length(WW3)!=0)
        {
          for (v2 in 1:length(WW3))
          {
            x <- WW2[1]
            y <- WW1[v1]
            z <- WW3[v2]
            
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]]))
              {
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest for continuous data
                {
                  pval <- gaussCItest(x, z, y, suffStat) #additional pval
                }
                if(indepTest=="disCItest") #if indepTest=disCItest for discrete data 
                {
                  pval <- disCItest(x, z, y, suffStat) #additional pval
                } 
                Alpha <- alpha
                if (pval <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
            
          }
          
        }
      }
      #if(length(WW3)!=0)
      #{
      #Orient remaining edges,
      T1 <- which(tarmat==1,arr.ind = T)
      G1 <- which(g==1,arr.ind = T)
      D1 <- setdiff(G1,T1)
      #
      if(length(D1)!=0)
      {
        for (d1 in 1:length(D1))
        {
          Rem1_row <- which(g[,D1[d1]]==1,arr.ind = T)
          Rem1_col <- which(g[D1[d1],]==1,arr.ind = T)
          Rem1 <- c(Rem1_row,Rem1_col)
          Rem2 <- which(tarmat[,Rem1]==1,arr.ind = T)
          if(length(Rem1)!=0 & length(Rem2)!=0)
          {
            for (d11 in 1:length(Rem1))
            {
              x <- Rem2[d11]
              y <- Rem1[d11]
              z <- D1[d1]
            }
            
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]]))
                
              {
                
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest
                {
                  pval <- gaussCItest(x, z, y, suffStat)
                }
                if(indepTest=="disCItest") #if indepTest=gaussCItest
                {
                  pval <- disCItest(x, z, y, suffStat) #additional
                }                
                Alpha <- alpha
                if (pval <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  R[m] <- 1
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
          }
        }
      }
      #Check the reamning edge orientation
      ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        
        if(tarmat[x, y]==0 & tarmat[y, x]==0) #bidirected if still no edge in tarmat
        {
          tarmat[x, y] <- 1
          tarmat[y, x] <- 1
        }
      }
    }
    
  }#end when data contain genetic variants
  #Start to orient remaining edges, when data not necessary to contain genetic variants.
  #First identify all the canidate genes for v-structure
  #Then identify the edges (undirected) with canidate genes 
  #ignore the nodes that already have direction
  #Make a triplet
  #determine the direction using test results from part 1 based on gMR
  #Repeat until all undirected edges to directed
  
  #Start when data not necessary to contain genetic variants.
  if (any(tarmat == 1) & GV==0) #if at least one edge directed already made so far
  {
    WW1 <- unique(which(tarmat == 1,arr.ind = T)[,2]) #pullout the all canidate genes for v-structure
    #WW1=WW1[-c(1:GV)]
    if(length(WW1)!=0)
    {
      #WW1=setdiff(WW1,1:GV)
      for (v1 in 1:length(WW1))
      {
        WW2 <- unique(which(tarmat[,WW1[v1]]==1,arr.ind = T)) #edge between v-structure  
        WW3 <- unique(which(g1[,WW1[v1]]==1,arr.ind = T)) #edges between canidate genes and others  
        WW3 <- setdiff(WW3,WW2)   #ignore if already direction
        if(length(WW3)!=0)
        {
          for (v2 in 1:length(WW3))
          {
            x <- WW2[1]
            y <- WW1[v1]
            z <- WW3[v2]
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]])) 
              {       
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest
                {
                  pval <- gaussCItest(x, z, y, suffStat)
                }
                if(indepTest=="disCItest") #if indepTest=disCItest
                {
                  pval <- disCItest(x, z, y, suffStat) #additional
                }
                Alpha <- alpha
                #Alpha=SeqFDR(m,FDR,a=2,R) #Alpha valued from sequential FDR test
                if (pval <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
            
          }
          
        }
      }
      #if(length(WW3)!=0)
      #{
      #Orient remaining edges,
      T1 <- which(tarmat==1,arr.ind = T)
      G1 <- which(g==1,arr.ind = T)
      D1 <- setdiff(G1,T1)
      #
      if(length(D1)!=0)
      {
        for (d1 in 1:length(D1))
        {
          Rem1 <- which(g[,D1[d1]]==1,arr.ind = T)
          Rem2 <- which(tarmat[,Rem1]==1,arr.ind = T)
          if(length(Rem1)!=0 & length(Rem2)!=0)
          {
            for (d11 in 1:length(Rem1))
            {
              x <- Rem2[d11]
              y <- Rem1[d11]
              z <- D1[d1]
            }
            
            if ((!is.na(x) & !is.na(y) & !is.na(z)) & (x!=0 & y!=0 & z!=0) & (x!=y & y!=z & z!=x))
            {
              #Case-2: If,y and z are adjacent, x and z conditionally independent given y,
              #then the edge direction will be y-->z
              if (g1[y, z] == 1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 & ((y %in% gInput@sepset[[x]][[z]]) || (y %in% gInput@sepset[[z]][[x]])))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 0
              }
              
              #Case-3: If,y and z are adjacent, x and z conditionally dependent given y,
              #then the edge direction will be z-->y.
              if (g1[y, z] == 1 & g1[x, z] != 1 & tarmat[x,y]==1 & tarmat[y, x]!=1 & tarmat[y, z]!=1 & tarmat[z, y]!=1 &!(y %in% gInput@sepset[[x]][[z]]) & !(y %in% gInput@sepset[[z]][[x]]))
                
              {
                if(indepTest=="gaussCItest") #if indepTest=gaussCItest
                {
                  pval <- gaussCItest(x, z, y, suffStat)
                }
                if(indepTest=="disCItest") #if indepTest=gaussCItest
                {
                  pval <- disCItest(x, z, y, suffStat) #additional
                }                
                Alpha <- alpha
                if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
                  tarmat[z, y] <- 1 #directed z-->y
                } 
                else {
                  tarmat[y, z] <- 1
                }
              }
              #Case-4:If, y & z have relation and x & y conditionally dependent given z and x & z conditionally independent given y,
              #then edge direction will be z<-->y
              if (g1[y, z] == 1 & g1[x, z] == 1 & g1[x, y] == 1 & !(z %in% gInput@sepset[[x]][[y]]) &!(y %in% gInput@sepset[[x]][[z]]))
              {
                tarmat[y, z]  <- 1
                tarmat[z, y]  <- 1
              }
            }
          }
        }
      }
      #Check the reaming edges
      ind <- which(g1 == 1, arr.ind = TRUE)  #Pullout the all relation in adjacency matrix from undirected graph
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        if(tarmat[x, y]==0 & tarmat[y, x]==0) #bidirected if still no edge in tarmat
        {
          tarmat[x, y] <- 1
          tarmat[y, x] <- 1
        }
      }
    }
    
  }#end when data not necessary to contain genetic variants.
  
  #Step 4
  #Produce as a same skeleton if no edges involves with genetic variants and no v-structures
  if(all(tarmat==0))
  {
    tarmat <- g1
  }
  #Step 5
  #If found the direction with only genetic variants, no v-structures and 
  #all the other nodes still undirected, then remaining edges will be bidirected 
  if(GV>0 & (any(tarmat[1:GV,]==1) || any(tarmat[,1:GV]==1)) & all(tarmat[-c(1:GV),-c(1:GV)]==0))
  {
    tarmat1 <- g1
    tarmat1[1:GV,] <- tarmat[1:GV,]
    tarmat1[,1:GV] <- tarmat[,1:GV]
    tarmat <- tarmat1
  }
  
  gInput@graph <- as(tarmat, "graphNEL")
  gInput
  
}

ModiSkeleton1 <- function (data, suffStat, FDR, alpha, filtered, indepTest = c("gaussCItest", "disCItest","citest"), labels, p, method = c("stable", "original", "stable.fast"), m.max = Inf, fixedGaps = NULL,
                          fixedEdges = NULL, NAdelete = TRUE, FDRcontrol = c("LOND", "ADDIS"), 
                          tau = 0.5, lambda = 0.25, verbose = FALSE) {
  cl <- match.call()
  if (!missing(p))
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  }
  else {
    stopifnot(is.character(labels))
    if (missing(p))
      p <- length(labels)
    else if (p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)))
    stop("fixedGaps must be symmetric")
  else G <- !fixedGaps
  diag(G) <- FALSE
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)))
    stop("fixedEdges must be symmetric")
  {
    sepset <- lapply(seq_p, function(.) vector("list", p))
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    m <- 0L    #Current test number
    Alpha <- 0L
    K <- 0
    
    # Initialize vectors/scalars that are specific to the addis function.
    R <- pval <- kappai <- Ci <- Si <- Ci_plus <- numeric(dim(data)[2]^2)
    gammai <- kappai_star <- alphai <- numeric(dim(data)[2]^2)
    
    # Create objects for the numerator and exponent of the gamma series. This
    # is a p-series whose infinite sum is 1.
    normalizer <- 0.4374901658
    exponent <- 1.6
    
    # Initialize the gammai_sum to be zero. This will be updated after the first
    # rejection occurs.
    gammai_sum <- 0
    
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L

    if (!missing(filtered)){
      G <- filtered
    }

    n.edgetests <- numeric(1)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      #ind <- which(upper.tri(G), arr.ind = TRUE) # if not using own input use this instead?
      ind <- which(G, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      # Order refers to the number of nodes being conditioned on.
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,
            "\n", sep = "")
      if (method == "stable") {
        G.l <- split(G, gl(p, p))
      }
      
      # Loop through all possible edges.
      for (i in 1:remEdges) {
        
        # Print every 100th index of the edges considered and the number of all
        # possible edges. If a test is performed, the details of this test will
        # also be printed later (starting at line 161).
        if (verbose && (verbose >= 2 || i%%100 == 0))
          cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable")
            
            nbrsBool1 <- G.l[[x]]
          nbrsBool2 <- G.l[[y]] 
          # G.l[[x]]
          #else G[, x]
          nbrsBool1[y] <- FALSE
          nbrs_x <- seq_p[nbrsBool1]
          #G.l[[y]]
          #else G[, y]
          nbrsBool2[x] <- FALSE
          nbrs_y <- seq_p[nbrsBool2]
          
          nbrs <- unique(union(nbrs_x,nbrs_y))
          
          #G.l[[x]]
          #else G[, x]
          #nbrsBool[y] <- FALSE
          #nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              
              m <- m+1  #Total number of the test
              
              # Increase the length of the R, pval, and other ADDIS vectors if
              # their length is greater than or equal to the current iteration.
              if (m >= length(R)) {
                
                R <- c(R, numeric(dim(data)[2]^2))
                pval <- c(pval, numeric(dim(data)[2]^2))
                kappai <- c(kappai, numeric(dim(data)[2]^2))
                kappai_star <- c(kappai_star, numeric(dim(data)[2]^2))
                Ci <- c(Ci, numeric(dim(data)[2]^2))
                Si <- c(Si, numeric(dim(data)[2]^2))
                Ci_plus <- c(Ci_plus, numeric(dim(data)[2]^2))
                gammai <- c(gammai, numeric(dim(data)[2]^2))
                alphai <- c(alphai, numeric(dim(data)[2]^2))
                
              }
              
              ##Start to calculate P-value using ci.test and gaussCItest
              
              if(indepTest == "citest") #if indepTest=ci.test
              {   
                x <- data[,ind[i,1]]
                y <- data[,ind[i,2]]
                z <- data[,nbrs[S]]
                if (length(S)==0) {
                  P <- ci.test(x, y)
                  pval[m] <- P$p.value
                } else {
                  P <- ci.test(x, y, z)
                  pval[m] <- P$p.value  #P-Value
                }
                x <- ind[i, 1]
                y <- ind[i, 2]
              }
              if(indepTest=="gaussCItest") #if indepTest=gaussCItest
              {
                pval[m] <- gaussCItest(x, y, nbrs[S], suffStat)
              }
              if(indepTest=="disCItest") #if indepTest=disCItest
              {
                pval[m] <- disCItest(x, y, nbrs[S], suffStat)
              }
              
              #End to calculate P-value using ci.test and gaussCItest
              
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S],"\n")
              if (is.na(pval[m]))
                pval[m] <- as.numeric(NAdelete)
              if (pMax[x, y] < pval[m])
                pMax[x, y] <- pval[m]
              if (verbose)
                cat("Test number =", m, "\n")
              if (verbose)
                cat("pval =", pval[m], "\n")
              
              if (FDRcontrol == 'LOND') { #if want to control sequential FDR 
                
                # Calculate alpha using the LOND method.
                alphai[m] <- SeqFDR(m,FDR,a=2,R)
                
                Alpha <- alphai[m]
                
              } else if (FDRcontrol == 'ADDIS') {
                
                # Calculate alpha using the ADDIS algorithm.
                
                # Initialize all the vectors and values for the first iteration.
                if (m == 1) {
                  
                  # Calculate w0 from tau, lambda, and alpha. This value is used
                  # in calculating alpha_t at each iteration.
                  w0 <- tau * lambda * FDR/2
                  
                  # Calculate the sum of the candidate p-values.
                  Ci_sum <- 0
                  
                  # Calculate the sum of the selected tests for each p-value tested.
                  Si_sum <- 0
                  
                  # The total number of rejections so far.
                  K <- 0
                  
                  # Calculate the first element in the gamma sequence.
                  gammai[1] <- normalizer / 1^exponent
                  
                  # Calculate alphai for the first test.
                  alphai[1] <- w0 * gammai[1]
                  
                  # Update the Alpha value with alphai[1]
                  Alpha <- alphai[1]
                  
                  # Determine if the first test should be rejected
                  R[1] <- pval[1] <= alphai[1]
                  
                } else {
                  
                  # Run ADDIS on the current iteration and update all the vectors
                  # and other values.
                  run_addis <- addis(alpha = FDR,
                                     tau = tau,
                                     lambda = lambda,
                                     iter = m,
                                     w0 = w0,
                                     pval = pval,
                                     alphai = alphai,
                                     gammai = gammai,
                                     kappai = kappai,
                                     kappai_star = kappai_star,
                                     K = K,
                                     Ci = Ci,
                                     Si = Si,
                                     Ri = R,
                                     Ci_plus = Ci_plus,
                                     Ci_sum = Ci_sum,
                                     Si_sum = Si_sum,
                                     gammai_sum = gammai_sum,
                                     normalizer = normalizer,
                                     exponent = exponent)
                  
                  # Update all values and vectors output from the addis function.
                  alphai[[m]] <- run_addis[[1]]
                  gammai[[m]] <- run_addis[[2]]
                  K <- run_addis[[3]]
                  R[[m]] <- run_addis[[5]]
                  Si[[m - 1]] <- run_addis[[6]]
                  Ci[[m - 1]] <- run_addis[[7]]
                  Ci_sum <- run_addis[[9]]
                  Si_sum <- run_addis[[10]]
                  gammai_sum <- run_addis[[12]]
                  
                  # Only update the kappai and Ci_plus vectors if K is greater
                  # than one. If K is zero then the first element in kappai will
                  # remain zero until the first rejection.
                  if (K != 0) {
                    
                    kappai[[K]] <- run_addis[[4]]
                    kappai_star[[K]] <- run_addis[[11]]
                    Ci_plus[1:K] <- run_addis[[8]]
                    
                  }
                  
                  # Update the Alpha value with the current alphai value.
                  Alpha <- alphai[[m]]
                  
                }
                
                
              } else {
                Alpha <- alpha #if want to use fixed significance level 
              }
              if (verbose)
                cat("Alpha value =", Alpha, "\n")
              
              if (pval[m] <= Alpha) {  #Reject H0 (H0:nodes are independent)
                R[m] <- 1
                if (verbose)
                  cat("Since pval<Alpha,test is rejected: Nodes are dependent", "\n")
              } else {
                R[m] <- 0  #Accept H0
                if (verbose)
                  cat("Since pval>Alpha,test is accepted:Nodes are independent", "\n")
              }
              if (pval[m] >= Alpha) {
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                
                break
              } else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
              
            }
            
          }
        }
        
      }
      
      ord <- ord + 1L
    }
    for (i in 1:(p - 1)) {
      for (j in 2:p) {
        pMax[i, j] <- pMax[j, i] <- max(pMax[i,j], pMax[j, i])
      }
    }
    }
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  } else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  
  #temp<-new("pcAlgo",graph = Gobject,call = cl, n = integer(0),
  # max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
  # sepset = sepset,pMax = pMax, zMin = matrix(NA, 1, 1))
  
  new("MRPCclass",
      graph = Gobject,
      call = cl,
      n = integer(0),
      max.ord = as.integer(ord - 1),
      n.edgetests = n.edgetests,
      sepset = sepset,
      pMax = pMax,
      zMin = matrix(NA, 1, 1),
      test = m,
      alpha = Alpha,
      R = R,
      K = K,
      pval = pval,
      normalizer = normalizer,
      exponent = exponent,
      alphai = alphai,
      kappai = kappai,
      kappai_star = kappai_star,
      Ci = Ci,
      Si = Si,
      Ci_plus = Ci_plus,
      gammai = gammai,
      gammai_sum = gammai_sum)
  
  #return(list(obj=temp,test=m,alpha=Alpha,R=R))
  #else
  #{
  #return(list(graph = Gobject,call = cl, n = integer(0),
  # max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
  # sepset = sepset,pMax = pMax, zMin = matrix(NA, 1, 1),test=m,alpha=Alpha,R=R))
  
  #}
}


find_filtered <- function (data, alpha, cor_mat) {
  #x <- psych::corr.test(data, adjust = "none") #corr.test runs cor.test on matrix
  #Gc_0 <- x$p <= alpha
  #cor_mat = cor(data)
  #add in ability to find cor_mat in here if not give
  #add in ability to include N with the function if wanted
  N = nrow(data)
  T_c <- abs(qnorm(alpha/2))
  Gc_0 <- Fztest(cor_mat*upper.tri(cor_mat), N, T_c)
  ##GeneNet 
  p_mat <- pcor(data)
  p_mat <- p_mat$estimate
  pcor = sm2vec(p_mat)
  indexes = sm.index(p_mat)
  colnames(indexes) = c("node1", "node2")
  w = cbind(pcor, indexes)
  temp <- w[, 1]
  temp[temp < -0.999999] <- temp[temp < -0.999999] + 0.0000001 #if any pcor are -1, increases to slightly higher to avoid error in fdrtool
  temp[temp > 0.999999] <- temp[temp > 0.999999] - 0.0000001 #if any pcor are 1, reduces to slightly lower to avoid error in fdrtool
  w[, 1] <- temp
  fdr.out = fdrtool(w[, 1], statistic = "correlation", cutoff.method = "locfdr")
  #fdr.out = fdrtool(w[, 1], statistic = "correlation", cutoff.method = "fndr")
  pval = fdr.out$pval
  qval = fdr.out$qval
  prob = 1 - fdr.out$lfdr
  result = cbind(w, pval, qval, prob)
  sort.idx = order(-abs(result[, 1]))
  result = as.data.frame(result[sort.idx, ])
  num.significant.1 <- sum(result$pval <= alpha)
  sig_node <- result[1:num.significant.1,]
  temp <- matrix(ncol=nrow(p_mat),nrow=nrow(p_mat))
  temp[is.na(temp)] <- FALSE
  for(i in 1:nrow(sig_node)){
    info <- sig_node[i,]
    pos1 <- info$node1
    pos2 <- info$node2
    temp[pos1, pos2] <- TRUE
  }
  colnames(temp) <- colnames(p_mat)
  rownames(temp) <- rownames(p_mat)
  Gp_0 <- temp
  # keep track of sepsets of precision matrix
  pres_sepsets <- Gp_0 | t(Gp_0)
  # combined 
  G_0 <- (Gc_0 & Gp_0)*upper.tri(Gc_0)
  # we keep track of the skeleton in a symmetric matrix
  G_cur <- 1*(G_0 | t(G_0)) # current graph
  G <- sapply(as.data.frame(G_cur), as.logical)
}

