# Function to generate mean reliability diagrams + simulation example
options(scipen = 9999) # to avoid scientific notation

meanreldiag = function(fcast,y,lim = NULL,
                       resampling = FALSE,n_resamples = 999,digits = 3,region_level = 0.9,
                       inset_hist = FALSE,main = "Mean Reliability",xlab = NULL,ylab = NULL){
  n = 1
  if(is.null(lim)){
    lim = range(sapply(n,function(n) range(fcast$m[,n])))
    lim[1] = lim[1]
    lim[2] = lim[2]
    adj_summand = c(-1,1)*max(abs(lim))*0.2
  }
  else adj_summand = 0
  if(is.null(xlab)) xlab = expression(m[1](F))
  if(is.null(ylab)) ylab = expression(m[group("",list(1,rc),"")](F))
  
  plot(NULL,xlim = lim+adj_summand,ylim = lim+adj_summand,main = main,
       xlab = xlab,ylab = ylab)
  
  m_rc = isoreg(fcast$m[,n],y)$yf # correspond to ORDERED forecast values!
  
  score = function(x,y) mean((x-y)^2)
  s = score(fcast$m[,n],y)
  s_rc = score(m_rc,y[order(fcast$m[,n])])
  s_mg = score(mean(y),y)
  mcb = s - s_rc
  dsc = s_mg - s_rc
  unc = s_mg
  
  if(resampling){
    set.seed(99)
    
    n_samples = n_resamples + 1 # total number of samples including observed sample
    low = floor(n_samples * (1-region_level)/2)
    up = n_samples - low
    pval_digits = ceiling(log(n_samples,10))
    
    resamples = sapply(1:n_resamples,function(i) fcast$sample())
    
    m = fcast$m[,n]
    
    m_rc_resamples = apply(resamples, 2, function(y) isoreg(m,y)$yf)
    m_rc_resamples_sorted = apply(cbind(m_rc,m_rc_resamples),1,sort) # includes observed values!
    
    # color_step_poly(lim,sort(m),sort(m),m_rc_resamples_sorted[up,],m_rc_resamples_sorted[low,])
    polygon(c(lim[1],sort(m),lim[2],rev(sort(m)),lim[1]),
            c(lim[1],m_rc_resamples_sorted[up,],lim[2],rev(m_rc_resamples_sorted[low,]),lim[1]),
            border = NA,col = "lightblue1")
    points(sort(m),m_rc_resamples_sorted[low,],type = "l",lty = 1,col = "lightblue2")
    points(sort(m),m_rc_resamples_sorted[up,],type = "l",lty = 1,col = "lightblue2")
    box()
    
    mcb_resamples = sapply(1:n_resamples,function(i) score(m,resamples[,i]) - score(m_rc_resamples[,i],resamples[order(m),i]))
    mcb_bounds = sort(c(mcb,mcb_resamples))[c(low,up)]
    
    rank_obs = tail(rank(c(mcb_resamples,mcb)),1)
    pval = 1 - (rank_obs - 1)/(n_resamples + 1)
    
    text(x = (lim + adj_summand)[1],y = (lim + adj_summand)[2],
         paste0(c("MCB ","DSC ","UNC "),
                bquote(.(format(round(c(mcb,dsc,unc),digits = digits)),nsmall = digits)),
                c(paste0(" [","p = ", bquote(.(format(round(pval,digits = pval_digits),nsmall = pval_digits))),"]"),"",""),
                collapse = "\n"),
         adj = c(0,1))
  }
  else text(x = (lim + adj_summand)[1],y = (lim + adj_summand)[2],
            paste0(c("MCB ","DSC ","UNC "),
                   bquote(.(format(round(c(mcb,dsc,unc),digits = digits)),nsmall = digits)),
                   collapse = "\n"),
            adj = c(0,1))
  
  abline(a = 0,b = 1,col = "grey",lty = 2)
  points(sort(fcast$m[,n]),m_rc,type = "l")
  
  if(inset_hist){
    par.old = par(fig = c(0.7,0.9,0.25,0.45),mar = c(1,0,0,0),new = TRUE)
    plot(hist(fcast$m[,n],main = "",yaxt = "n",xlab = "",ylab = ""),add = TRUE)
    par(par.old)
  }
  
  return(list(MDU = c(mcb,dsc,unc),if(resampling) list(pval,mcb_bounds) else NULL))
}

# Simulation example

setup_normal = function(n,tau0,eta0,seed = 99){
  set.seed(seed)
  mu = rnorm(n)
  tau = sample(tau0*c(-1,1),n,replace = TRUE)
  eta = sample(eta0*c(-1,1),n,replace = TRUE)
  y = rnorm(n) + mu
  
  F_perf = function(x) pnorm(x,mu)
  F_unf = function(x) 0.5*(pnorm(x,mu) + pnorm(x,mu+tau))
  F_lop = function(x) ifelse(x <= mu,
                             (1-eta)*pnorm(x,mu),
                             (1+eta)*pnorm(x,mu)-eta)
  
  m_perf = cbind(mu, mu^2 + 1, mu^3 + 3*mu)
  m_unf = cbind(mu + 0.5*tau,
                mu^2 + tau*mu + 0.5*tau^2 + 1,
                mu^3 + 1.5*tau*mu^2 + 3*(0.5*tau^2 + 1)*mu + 0.5*(tau^3 + 3*tau))
  phi0 = dnorm(0)
  m_lop = cbind(mu + 2*eta*phi0,
                mu^2 + 1 + 4*eta*phi0*mu,
                mu^3 + 3*mu + 2*eta*phi0*(3*mu^2 + 2))
  
  sample_perf = function(cases = 1:n) rnorm(length(cases)) + mu[cases]
  sample_unf = function(cases = 1:n){
    mix = sample(c(0,1),length(cases),replace = TRUE)
    return(rnorm(length(cases)) + mu[cases] + mix*tau[cases])
  }
  sample_lop = function(cases = 1:n){
    mix = sample(c(-1,1),length(cases),replace = TRUE,prob = c(1-eta0,1+eta0))*sign(eta[cases])
    return(mu[cases] + mix*abs(rnorm(length(cases))))
  }
  
  return(list(mu = mu, tau = tau, eta = eta, y = y,
              perf = list(F = F_perf, m = m_perf, sample = sample_perf),
              unf = list(F = F_unf, m = m_unf, sample = sample_unf),
              lop = list(F = F_lop, m = m_lop, sample = sample_lop)))
}

sim = setup_normal(400,1.5,0.7)

# Mean reliability diagrams for perfect, unfocused and lopsided forecast
meanreldiag(sim$perf,sim$y,resampling = TRUE,inset_hist = TRUE)
meanreldiag(sim$unf,sim$y,resampling = TRUE,inset_hist = TRUE)
meanreldiag(sim$lop,sim$y,resampling = TRUE,inset_hist = TRUE)
