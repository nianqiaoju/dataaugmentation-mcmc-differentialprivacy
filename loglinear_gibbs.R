rm(list = ls());
iexp <- 1
noiselevel <- 1;
niter <- 10000; ## number of gibbs iterations
invisible(eval(parse(text=commandArgs(TRUE))));
set.seed(iexp);
epsvec <- rev(c(0.3,0.1,1,3,10));
epsilon <- epsvec[noiselevel];
## read non-private data n 
load(file = "loglinear_datageneratingparameters.RData");
load(file = "loglinear_notdpdata_N100.RData");
library("extraDistr");

## add noise
mijk <- list();
for(k in 1 : K){
	mijk[[k]] <- nijk[[k]] + VGAM::rlaplace(n = Jk[k] * I, location = 0, scale = 2 * K / epsilon);
}

## some helper functions for gibbs sampler-------
### summarize group identity into n
get_n <- function(g){
  return(as.vector(table(factor(g, levels = 1 : I))));
}

### summarize individual survey answers into nijk
get_nijk <- function(g,ans){
  nijk <- list();
  for(k in 1 : K){
    nijk[[k]] <- as.matrix(table(factor(ans[,k], levels = 1 : Jk[k]), factor(g, levels = 1 : I)));
  }
  return(nijk);
}

### sum squared error from nijk to mijk
get_sse <- function(nijk){
  dis <- 0;
  for(k in 1 : K){
    dis <- dis + sum((nijk[[k]] - mijk[[k]])^2);
  }
  return(dis);
}

### initialize the gibbs sampler to have p, pijk, g, ans, n, and nijk
rgibbs_init <- function(){
  p <- rgamma(n = I, shape = alp, rate = 1);
  p <- p / sum(p);
  pijk <- list();
  for(k in 1: K){
    pijk_ <- matrix(0, nrow = Jk[k], ncol= I);
    for(i in 1 : I){
      pijk_[,i] <- rgamma(n = Jk[k], shape = alpha[[k]][,i], rate = 1);
      pijk_[,i] <- pijk_[,i] / sum(pijk_[,i]);
    }
    pijk[[k]] <- pijk_;
  }
  state <- list();
  state$p <- p;
  state$pijk <- pijk;
  ## sample the individual group identities
  g <- sample.int(I, N, TRUE, p);
  ans <- matrix(nrow = N, ncol = K);
  for(k in 1 : K){
    for(id in 1 : N){
      ans[id,k] <- extraDistr::rcat(1, pijk[[k]][,g[id]]);
    }
  }
  state$g <-  g;
  state$ans <- ans;
  state$n <- get_n(g);
  state$nijk <- get_nijk(g,ans);
  return(state);
}

gibbs_iter <- function(state){
  ## Sample p given n from dirichlet
  state$p <- extraDistr::rdirichlet(1, as.vector(state$n) + alp); 
  for(k in 1 : K){
    for(i in 1 : I){
      state$pijk[[k]][,i] <- extraDistr::rdirichlet(1, state$nijk[[k]][,i] + alpha[[k]][,i]); 
    }
  }
  ## independent MH update for each g[id] and ans[id,]
  acnt <- 0; 
  for(id in 1 : N){## for each person
    lratio <- 0;
    newgid <- extraDistr::rcat(1,state$p);
    newans <- rep(NA, K);
    for(k in 1 : K){## for each feature
      newans[k] <- rcat(1, state$pijk[[k]][,newgid]); ## sample new answer xk
      ## update ratio of eta's, in log scale
      lratio <- lratio + extraDistr::dlaplace(state$nijk[[k]][newans[k], newgid] - mijk[[k]][state$ans[newans[k]], newgid] + 1, 0, 2 * K / epsilon, TRUE); 
      lratio <- lratio + extraDistr::dlaplace(state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]] - 1, 0, 2 * K / epsilon, TRUE);
      lratio <- lratio - extraDistr::dlaplace(state$nijk[[k]][newans[k], newgid] - mijk[[k]][state$ans[newans[k]], newgid], 0, 2 * K / epsilon, TRUE);
      lratio <- lratio - extraDistr::dlaplace(state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]] , 0, 2 * K / epsilon, TRUE);

    }
    if(log(runif(1)) < lratio){## if accept then book keeping
      for(k in 1 : K){
        state$nijk[[k]][state$ans[id,k],state$g[id]] <- state$nijk[[k]][state$ans[id,k],state$g[id]] - 1;
        state$nijk[[k]][newans[k], newgid] <- state$nijk[[k]][newans[k], newgid] + 1;
      }
      state$n[state$g[id]] <- state$n[state$g[id]] - 1;
      state$n[newgid] <- state$n[newgid] + 1;
      state$g[id] <- newgid;
      state$ans[id,] <- newans;
      acnt <- acnt + 1; ## increment acceptance count
    }
    rm(newgid); rm(newans);lratio <- 0;
  }
  state$accept <- acnt / N;
  state$dis <- get_sse(state$nijk);
  return(state);
}



## initialize gibbs sampler
state <- rgibbs_init();
state
p_gibbs <- matrix(NA, ncol = I, nrow = niter);
acc_gibbs <- rep(NA, niter);
dis_gibbs <- rep(NA, niter);
## run gibbs sampler
for(iter in 1 : niter){
  state <- gibbs_iter(state);
  p_gibbs[iter, ] <- state$p;
  acc_gibbs[iter] <- state$accept;
  dis_gibbs[iter] <- state$dis;
  if(iter %% 100 == 0){
    cat(paste(iter, "out of", niter, "iterations done.\n"));
    cat(state$p, "\n");
  }
} 
## save posterior samples into csv file
dfsample <- data.frame(p_gibbs);
names(dfsample) <- paste("p", c(1:5), sep = "");
dfsample$iter <- 1 : niter;
dfsample$eps <- epsilon;
dfsample$iexp <- iexp;
head(dfsample)

## save acceptance rate into csv file
dfacc <- data.frame(iexp = iexp, eps = epsilon, acc = acc_gibbs);
head(dfacc);