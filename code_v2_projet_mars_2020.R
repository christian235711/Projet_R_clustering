
library(readr)
query <- read_csv("~/Downloads/query_6.csv")

data1 = query[,3:2]
names(data1)= c("Longitude", "Latitude")
plot(data1, xlab = "Longitude", ylab="Latitude")


###################################################################################################
###################################################################################################
###################################################################################################
## KNNC 

library(mclust)
library(prabclus) # 3 data

nnc = NNclean(data= data1, k = 5 , distances = NULL, edge.correct = FALSE, wrap = 0.1, convergence = 0.001, plot=F, quiet=TRUE)
plot(data1, col= 1+nnc$z )

longueur_1 = length(nnc$z[nnc$z==1])
longueur_1

data_apres_nnc =data.frame(data1, nnc$z)
names(data_apres_nnc)= c("xx", "yy","cc")
head(data_apres_nnc)

data_sans_bruit= subset(data_apres_nnc, cc==1, select = -cc) # dim 2
plot(data_sans_bruit, xlab="Longitude", ylab = "Latitude", col="blue" )

data_blackground_noise =subset(data_apres_nnc, cc==0, select = -cc) # dim 2
plot(data_blackground_noise, xlab="Longitude", ylab = "Latitude")

##########################################################
##########################################################
##########################################################
## MCLUST

mod1 <- Mclust(data_sans_bruit )   
summary(mod1)
plot(mod1, what="BIC")
plot(mod1,  what= "classification", ylab = "Latitude", xlab="Longitude")


##########################################################
##########################################################
##########################################################
## HPCC

library(princurve)

#######
#There are 5 functions here: hpcc,clust.var.spline, plot.pclust,penalty,vdist

hpcc<-function(x,clust.init,alpha=.4,plot.iter=F,plot.result=T,
               df=5,force=T,forcenum=1) {
  
  # This function returns the final classification and the list of
  # total variance at each merge (which aids determination of the
  # correct number of clusters without use of a penalty function).
  
  #x is the data in two columns (each row is a data point)
  #clust.init is the initial clustering; it has the same length as the
  #   number of rows in x, with each entry giving the number of the cluster
  #   to which that point initially belongs.  
  #alpha control the relative weight of variance along the curve and
  #   variance about the curve
  #plot.iter determines whether to plot intermediate results at each
  #   iteration
  #plot.result determines whether to produce a plot of the final result
  #df specifies the number of degrees of freedom to use in each feature
  #   cluster 
  #force=T causes clustering to continue until the number of clusters
  #   is equal to forcenum, instead of stopping when the total variance
  #   increases
  
  # plot.pclust is called to plot the result if plot.result=T.
  # functions clust.var.spline, vdist, and penalty are also called.	
  # trace/debugging info can be printed by uncommenting the cat
  #   statements
  
  
  # x is an n x 2 matrix of points
  n<-length(x[,1])
  p<-length(x[1,])  #p must be 2
  
  #uses clust.init as initial clustering.
  class1<-clust.init
  
  if(plot.result==T) {plot.pclust(x,class1,dg=df )}
  
  #initialize variables.
  lastmerge<-c(0,100)
  clvars<-0     
  mergevars<-0  
  totalvar.track<-rep(0,50)
  totalvar.counter<-0
  done<-F
  
  while(done==F) {	
    
    classes<-unique(class1)
    clnum<-length(classes)
    
    lastm<-min(lastmerge)
    lastn<-max(lastmerge)
    oldclvars<-clvars
    clvars<-rep(0,clnum)
    
    # Find variance of each cluster
    for(i in 1:clnum)  {
      #cat("Clvars for cluster ",i,":   ")
      if (i<lastm) {clvars[i]<-oldclvars[i]} #for optimizing
      else {if(i>lastn) {clvars[i]<-oldclvars[i+1]}
        
        else {
          
          #singletons are labeled 0 by mclust.  We arbitrarily give them
          #a large variance to force their inclusion in other clusters.
          
          if(classes[i]==0) {clvars[i]<-100}
          else {
            clvars[i]<-clust.var.spline(x[class1==classes[i],],plot.me=plot.iter,
                                        dg=df,alpha)}
        }}
      #cat(clvars[i],"\n") 
    }
    
    #calculate variance for each possible merge
    #alternate approaches: 
    #    1. calculate merge variances only until one is found that will
    #       reduce total variance (could be adapted to simulated annealing
    #       for large numbers of clusters) 
    #    2. could randomly choose merges to evaluate
    #    3. only look at merges of clusters which are "close"
    
    oldmergevars<-mergevars
    
    mergevars<-matrix(0,clnum,clnum)
    for(i in 1:clnum) {
      for(j in i:clnum) { 
        #cat("mergevars of ",i,j," : ")
        
        if(i==j){mergevars[i,j]<-clvars[i]}
        else {
          
          if(i<lastm&&j<lastm) {mergevars[i,j]<-oldmergevars[i,j]}
          else { if(i>lastn&&j>lastn) {mergevars[i,j]<-oldmergevars[i+1,j+1]}
            else { 
              mergevars[i,j]<- clust.var.spline(rbind(x[class1==classes[i],],
                                                      x[class1==classes[j],]),plot.me=plot.iter,
                                                dg=df,alpha) 
            }
          }
        }
        
        #cat(mergevars[i,j],"\n")                            
        
      }}
    
    #choose merge to minimize total variance and update class1
    #if no merge will reduce variance and force=F, we are done
    
    nomergevar<-sum(clvars)+penalty(clnum)
    currentmin<-rep(0,3)  #variance, clust num 1, clust num 2
    currentmin[1]<-1000   # this is just for initialization
    
    totalvars<-mergevars  #initialize totalvars 
    
    for(i in 1:clnum) {
      for(j in i:clnum) {
        
        if (j==i) {totalvars[i,j]<-nomergevar}
        else {
          
          totalvars[i,j]<-mergevars[i,j]+sum(clvars)-clvars[i]-clvars[j]+penalty(clnum-1)
          if(totalvars[i,j]<currentmin[1]) {currentmin[1]<-totalvars[i,j]
          currentmin[2]<-i
          currentmin[3]<-j }
        }}}
    totalvar.counter<-totalvar.counter+1
    totalvar.track[totalvar.counter]<-nomergevar
    
    if(currentmin[1]>nomergevar) {
      if(force==F||clnum==forcenum) {  cat("Done!\n")
        done<-T}
      else { 
        #cat("Forcing a merge.\n")
        #cat("If we stop merging now, the total variance is: ",nomergevar,"\n")
        #cat("Merging ",currentmin[2],currentmin[3],
        #         " will yield total variance of ",currentmin[1],"\n")
        lastmerge<-c(currentmin[2],currentmin[3])
        class1[class1==classes[currentmin[2]]]<-classes[currentmin[3]]
        if(plot.result==T) {plot.pclust(x,class1, dg=df)
        }}}
    
    
    else {
      if(clnum==forcenum) {cat(" Princlust Done!/n")
        done<-T }
      else{
        #cat("If we stop merging now, the total variance is: ",nomergevar,"\n")
        #cat("Merging ",currentmin[2],currentmin[3]," will yield total
        #     variance of ",currentmin[1],"\n")
        lastmerge<-c(currentmin[2],currentmin[3])
        class1[class1==classes[currentmin[2]]]<-classes[currentmin[3]]
        if(plot.result==T) {plot.pclust(x,class1, dg=df)}}}
  }
  
  if(plot.result==T) {plot.pclust(x,class1, dg=df) }
  
  classes<-class1
  
  l1= list(classes,totalvar.track,x)
  names(l1)= c("classes","totalvar.track","x")
  
  return (l1)
}

#------------------------------------------------------------------

clust.var.spline<-function(x,plot.me=F,dg=5,alpha=.4) {
  #This is the function which calls principal.curve
  #
  #If plot.me is T, then the iterative plots of the principal curve
  #will be displayed as they are being fit.
  
  #if there are fewer than 7 points in a cluster, we fit a principal component
  #line instead of a principal curve.  
  
  
  #fewer than seven points
  if(length(x[,1])<7) {
    
    temp<-prcomp(cbind(x[,1],x[,2]))
    temp.slope<-temp$rot[2,1]/temp$rot[1,1]
    temp.int<-mean(x[,2])-(temp.slope)*(mean(x[,1]))
    temp.coef<-c(temp.int,temp.slope)
    a<-projpoints(x,temp.coef)
    
    if(plot.me) {
      plot(x[,1],x[,2])
      abline(temp.coef, col="blue")# AJOUTEE
      points(a,pch=2)
    }
    
    vabout<-vdist(x,a)
    nn<-length(x[,1])
    epsilon<-a[1:(nn-1),]-a[2:nn,]
    mu<-apply(epsilon,2,mean)
    mu<-matrix(mu,byrow=T,ncol=2,nrow=nn-1)
    valong<-(.5)*vdist(epsilon,mu)
  }
  
  #more than seven points
  else {
    #temp<-principal.curve(x,plot.true=plot.me,df=dg,maxit=3) ## on a changé principal.curve par principal_curve
    temp<-principal_curve(x, plot_iterations = plot.me,df =dg, maxit = 5, smoother="smooth_spline") # argument inutilisé plot.true=FALSE
    vabout<-temp$dist
    nn<-length(temp$s[,1])
    
    # For closed principal curves, the nn index in epsilon should wrap
    # around, so length(epsilon) = nn.  For open curves, we have nn-1 segments.
    # We must put the points in the s matrix in the correct order before
    # calculating variance along the curve.  The ordering is provided by
    # the $tag value from principal.curve()
    
    
    # news<-temp$s[temp$tag,] # changé
    news<-temp$s[temp$ord,]
    
    epsilon<-news[1:(nn-1),]-news[2:nn,]
    mu<-apply(epsilon,2,mean)
    mu<-matrix(mu,byrow=T,ncol=2,nrow=nn-1)
    valong<-(.5)*vdist(epsilon,mu)
  }
  total<-alpha*vabout+valong
  #total<-vabout+alpha*valong
  total}

#------------------------------------------------------------------

plot.pclust<-function(x,classif, dg) {
  
  #x is the data that was used by hpcc (an n x 2 matrix)
  #classif is the classification (output of hpcc)
  #a graphics device must already be open
  
  classes<-unique(classif)
  clnum<-length(classes)
  plot(x,type="n",ylab = "Latitude", xlab="Longitude") # image blanche car type="n"
  
  for(i in 1:clnum) {
    points(x[classif==classes[i],],pch=i ,col= i+1 ) # pch=as.character(i)
    lines(principal_curve(x[(classif==classes[i]),],maxit=5, df=dg, smoother="smooth_spline" ),col="magenta", lwd=2) # AJOUTEE
  } 
  fnord<-23
  fnord}

#------------------------------------------------------------------
penalty<-function(k) {
  #A penalty can be specified here.  k is the number of clusters.
  final<-0
  final
}

#------------------------------------------------------------------

vdist<-function(x,y)  {
  #x, y are n x p matrices.  Each row is a point.
  total<-0
  n<-length(x[,1])
  p<-length(x[1,])
  for(i in 1:n) {
    for(j in 1:p) {total<-total+((x[i,j]-y[i,j])^2) }   
  }
  total}

#------------------------------------------------------------------

# HPCC APPLIQUEE POUR K=8 ET K=3

m = hpcc(x=as.matrix(data_sans_bruit) ,clust.init=mod1$classification, alpha= .1, plot.iter = F, plot.result = T,df=8, force=T,forcenum=3)

m$classes
table(m$classes)

##########################################################
##########################################################
##########################################################

# CEM-PCC




#This file has two functions: cem and calc.likelihood

cem<-function(x, df=5, clust.init = rep(0,(length(x[,1]))), trace=F,
              max.iter=5,plot.me=F,trim=0,mix=T,var.bound=-1)  {
  
  #x is the data in two columns (each row is a data point)
  #clust.init is the initial clustering; it has the same length as the
  #   number of rows in x, with each entry giving the number of the cluster
  #   to which that point initially belongs.  The feature clusters are 
  #   numbered by consecutive integers starting with 1; 0 indicates the
  #   noise cluster. 
  #trace=T causes debugging information to be printed
  #max.iter is the number of EM iterations to do
  #plot.me=T causes results to be plotted
  #trim can be used to compute a robust estimate of the variance about the
  #   curve for any feature cluster.  The trim amount should be between
  #   0 and 0.4.
  # mix determines whether to use the likelihood or the mixture
  #   likelihood (recommended) for each point when forming the new
  #   clustering at each iteration.  The default value T uses the
  #   mixture likelihood.
  #var.bound is a lower bound on the estimate of the variance about the
  #   curve for any feature cluster.  If var.bound < 0, it is computed
  #   automatically based on the range of the data.  Setting var.bound
  #   equal to zero removes any constraint from the variance estimate.
  
  
  n<-length(x[,1])
  area<-(max(x[,1])-min(x[,1]))*(max(x[,2])-min(x[,2]))
  done<-F
  iterations<-0
  
  # initialize the clustering and the variables which will keep track of
  # results at each iteration
  clust<-clust.init
  clust.log<-rep(0,max.iter)
  overall.loglikelihood<-rep(0,max.iter)
  all.curvevars<-matrix(rep(0,max.iter*length(unique(clust))),
                        ncol=length(unique(clust)))
  all.varalong<-matrix(rep(0,max.iter*length(unique(clust))),
                       ncol=length(unique(clust)))
  all.meanalong<-matrix(rep(0,max.iter*length(unique(clust))),
                        ncol=length(unique(clust)))
  all.curvelengths<-matrix(rep(0,max.iter*length(unique(clust))),
                           ncol=length(unique(clust)))
  all.mixprop<-matrix(rep(0,max.iter*length(unique(clust))),
                      ncol=length(unique(clust)))
  
  #compute the variance bound, if needed
  if (var.bound<0) { 
    var.bound<- (max(c(range(x[,1])[2]-range(x[,1])[1],
                       range(x[,2])[2]-range(x[,2])[1]))/4000)^2 }
  
  
  
  ### MAIN PROGRAM LOOP ###	
  while (done==F) {
    
    #initialize loop variables
    last.clust<-clust
    classes<-unique(clust)
    num.of.classes<-length(classes)
    
    curvelengths<-rep(0,num.of.classes)
    curvevars<-rep(0,num.of.classes)
    varalong<-rep(0,num.of.classes)
    meanalong<-rep(0,num.of.classes)
    sqdists<-matrix(rep(0,n*num.of.classes),ncol=num.of.classes)
    distalong<-matrix(rep(NA,n*num.of.classes),ncol=num.of.classes)
    likelihoods<-matrix(rep(0,n*num.of.classes),ncol=num.of.classes)
    likelihoods.mix<-matrix(rep(0,n*num.of.classes),ncol=num.of.classes)
    like.along<-matrix(rep(0,n*num.of.classes),ncol=num.of.classes)
    like.about<-matrix(rep(0,n*num.of.classes),ncol=num.of.classes)
    
    if (trace) {cat("Clustering\n",clust,"\n")}
    if (plot.me) {
      plot(x,xlab="X",ylab="Y",main="CEM Clustering",pch= 1, col=1 )
    }
    
    #------------------------------------------------------
    #fit curve and get parameter estimates for each cluster
    
    for ( j in 1:num.of.classes) {
      
      if (classes[j]>0) {  #this is a feature cluster
        filter<-(clust==classes[j])
        n.curve <- length(x[filter, 1])
        if (n.curve!=sum(filter)) {stop("Filter error.")}
        curve <- principal_curve(x[filter,],maxit=5, df=df,  smoother= "smooth_spline") ###### changé
        if (plot.me) {
          plot(x,xlab="X",ylab="Y", col=0) # AJOUTEE 
          points(x[clust==classes[j],],pch=classes[j], col= j+1 ) # car ici col = classes[j] ne contient pas 0
          points(x[clust!=classes[j],], col =0) # AJOUTEE 
          lines(curve, col = "magenta")  }
        curvelengths[j] <- (curve$lambda[curve$ord])[n.curve] # on a changé tag par ord
        
        if (trace) {
          cat("\n","j= ",j," \n")
          cat("classes[j]= ",classes[j]," \n")
          cat("n.curve= ",n.curve," \n")
          cat("curvelengths[j]= ",curvelengths[j]," \n")
          cat("length(curve$lambda) ",length(curve$lambda)," \n")
          cat("length(curve$tag)",length(curve$ord)," \n") # on a changé tag par ord
          cat("curve$lambda[curve$tag]",curve$lambda[curve$ord]," \n") # on a changé tag par ord
        }
        
        #naive variance estimate: curvevars[j] <- (curve$dist)/(n.curve - 1)
        #or robust estimate from trimmed variance
        temp.dists<-x[filter,]
        temp.dists<-((temp.dists[,1]-curve$s[,1])^2)+
          ((temp.dists[,2]-curve$s[,2])^2)
        temp.dists<-temp.dists[order(temp.dists)]
        
        if (trim>0) {
          trimlow<-ceiling((trim/2)*length(temp.dists))
          trimhigh<-floor((1-(trim/2))*length(temp.dists))
          temp.dists<-temp.dists[c(trimlow:trimhigh)]
        }	
        
        #check variance bound
        curvevars[j] <-sum(temp.dists)/(length(temp.dists)-1)
        if (curvevars[j]<var.bound) { curvevars[j]<-var.bound }
        
        for (i in 1:n) {
          sqdists[i,j]<-project_to_curve(cbind(x[i, 1], x[i, 2]), 
                                         s = curve$s,stretch=0)$dist
        } # on a supprimé tag = curve$tag
        
        distincurve<-rep(NA,(n.curve-1))
        
        for (i in 2:n.curve) {
          distincurve[(i-1)]<-(curve$lambda[curve$ord])[i] -
            (curve$lambda[curve$ord])[i-1]
        } # on a changé tag par ord
        #ad hoc method to deal with distincurve at endpoints
        distincurve<-c(distincurve,distincurve[(n.curve-1)])
        
        meanalong[j]<-mean(distincurve)
        varalong[j]<-var(distincurve)
        
        for (i in 1:n) {
          #find order on curve
          temp.filter<-filter
          filter.ind<-0
          if (temp.filter[i]) {filter.ind<-1}
          temp.filter[i]<-F
          this.proj<-NA
          this.proj<-project_to_curve(rbind(x[i,],x[temp.filter,]),s=curve$s,
                                      stretch=0) # on a supprimé tag=curve$tag
          #find index of this point. It is point i, but it will be the 
          #first point in this.proj 
          this.n<-sum(temp.filter)+1
          this.point<-(c(1:this.n))[this.proj$ord==1] # on a changé tag par ord
          #if (plot.me) {segments(x[i,1],x[i,2],this.proj$s[1,1], 
          #                      this.proj$s[1,2])} ###### ON EFFACE POUR ELIMINER LES LIGNES DE LA PROJECTIONS
          if (this.point==1)  {
            #now we know that point i projects to the first endpoint
            #we want to extend the line from the last few projection points
            #to get distalong and sqdists.  This would be hard and
            #computationally expensive, so we approximate:
            distalong[i,j]<-sqrt(sqdists[i,j])
            #this will err on the side of making the curve too
            #conservative at the endpoints
            if (trace) {cat("Point ",i," first in cluster ",j,
                            " distalong: ",distalong[i,j],"\n")}
          }
          
          else { if (this.point==this.n) {
            #now we know that point i projects to the last endpoint
            #we want to extend the line from the last two projection points
            #to get distalong and sqdists
            
            distalong[i,j]<-sqrt(sqdists[i,j])
            if (trace) {cat("Point ",i," last in cluster ",j,
                            " distalong: ",distalong[i,j],"\n")}
          }
            
            else { #now we know that point i projects inside curve j
              
              #put the lambdas in order along the curve
              current.lams<-this.proj$lambda[this.proj$ord] # on a changé tag par ord
              dist1<-current.lams[this.point+1]-current.lams[this.point]
              dist2<-current.lams[this.point]-current.lams[this.point-1]
              distalong[i,j]<-dist1
              
              if (trace) {
                cat("Point ",i," is ",this.point," in cluster ",j,
                    " distalong: ",distalong[i,j],"\n")}
            }#end of else
          }#end of else
        }} #end of feature cluster routine
      
      
      if (classes[j]==0) { #now deal with noise cluster @#######@@@@@@@@ changé avec else if (classes[j]==0)
        curvelengths[j]<-0
        curvevars[j]<-0
        for (i in 1:n) {
          sqdists[i,j]<-0
        }
        if (plot.me) { points(x[clust==classes[j],],pch=classes[j] ) }
      } #end of noise cluster routine
      
    } #end of for loop
    
    if (plot.me) { # AJOUTEE
      plot(x,ylab = "Latitude", xlab="Longitude",pch= 1, col=0 )
      for ( j in 1:num.of.classes){
        if(classes[j]>0){
          points(x[clust==classes[j],],pch=classes[j], col= j+1 )
          lines(principal_curve(x[clust==classes[j],],maxit=5, df=df, smoother="smooth_spline"), col = "magenta", lwd=2) 
        }
        else{ points(x[clust==classes[j],], col= 1 ) }
      }
      
    }
    
    
    #---------------------------------------------------------
    #calculate likelihood of each point being in each cluster
    
    for (i in 1:n) {
      for (j in 1:num.of.classes) {
        if(classes[j]>0) {
          
          ### uniform*normal method:
          like.about[i,j]<-(1/(sqrt(2*pi*curvevars[j])))*
            (exp((-1/2)*(sqdists[i,j]/curvevars[j])))
          like.along[i,j]<-(1/curvelengths[j])
          mix.prop<-(sum(clust==classes[j]))/n
          likelihoods.mix[i,j]<- mix.prop*like.about[i,j]*like.along[i,j]
          likelihoods[i,j]<- like.about[i,j]*like.along[i,j]
        }
        
        if (classes[j]==0) { mix.prop<-(sum(clust==classes[j]))/n
        likelihoods[i,j]<-1/area
        likelihoods.mix[i,j]<-mix.prop/area
        like.along[i,j]<-0
        like.about[i,j]<-0
        
        }	                   
      }
    }
    #---------------------------------------------------------
    #save all parameter estimates
    
    all.curvevars[iterations+1,]<- curvevars # on a jaouté t
    all.varalong[iterations+1,]<-varalong
    all.meanalong[iterations+1,]<-meanalong
    all.curvelengths[iterations+1,]<-curvelengths
    mix.prop<-rep(NA,num.of.classes)
    for (j in 1:num.of.classes) {mix.prop[j]<-(sum(clust==classes[j]))/n}
    all.mixprop[iterations+1,]<-mix.prop
    
    
    #---------------------------------------------------------
    # make new clustering
    
    if (trace) {
      for (j in 1:num.of.classes) {
        cat("Likelihoods for cluster ",j,"\n",round(likelihoods[,j],dig=4),"\n\n")
      }
      cat("\n\nVar along: ",varalong,"\n")
      cat("\nVar about: ",curvevars,"\n")
      cat("\nMean along: ",meanalong,"\n")
    }
    
    if (mix) {
      for (i in 1:n) {
        largest<-max(likelihoods.mix[i,])
        for (j in 1:num.of.classes) { 
          if (likelihoods.mix[i,j]==largest) {clust[i]<-classes[j]}
        }
      }}
    
    if (!mix) {
      for (i in 1:n) {
        largest<-max(likelihoods[i,])
        for (j in 1:num.of.classes) { 
          if (likelihoods[i,j]==largest) {clust[i]<-classes[j]}
        }
      }}
    
    #---------------------------------------------------------
    # test for convergence of clustering or max number of iterations
    
    iterations<-iterations+1
    sameness.counter<-0
    for(i in 1:n) {if (last.clust[i]==clust[i]) 
    {sameness.counter<-sameness.counter+1}}
    clust.log[iterations]<-sameness.counter
    if(sameness.counter==n) {done<-T}
    if(iterations>=max.iter) {done<-T}
    temp<-0
    for (i in 1:n) {
      temp<-temp+log(sum(likelihoods.mix[i,]))
    }
    overall.loglikelihood[iterations]<-temp
    if
    (overall.loglikelihood[iterations]>=max(overall.loglikelihood[1:iterations])) 
    { clust.best<-clust }
    
  } # end of while loop -- END OF MAIN PROGRAM LOOP
  
  
  
  #---------------------------------------------------------
  l2=  list(varalong,meanalong,curvelengths,curvevars,likelihoods,likelihoods.mix,
            clust,iterations,clust.log,x,df,var.bound,like.along,like.about,
            distalong,
            sqdists,clust.best,overall.loglikelihood,
            all.curvevars,all.varalong,all.meanalong,all.curvelengths,all.mixprop) 
  names(l2)=c("varalong","meanalong","curvelengths","curvevars","likelihoods","likelihoods.mix","clust",
              "iterations","clust.log","x","df","var.bound", "like.along","like.about","distalong",
              "sqdists","clust.best","overall.loglikelihood",
              "all.curvevars","all.varalong","all.meanalong","all.curvelengths","all.mixprop")
  
  return (l2)
  
}

#---------------------------------------------------------
#---------------------------------------------------------

calc.likelihood <- function(cem.obj, mix.best = T)
{
  #cem.obj is the output from cem().  it includes:
  #x is the data
  #df is the degrees of freedom
  #clust is the classification (0 indicates the noise cluster)
  #likelihoods is a matrix of likelihoods (one column per cluster)
  #if trim was used in cem, then the likelihood calculated here will
  #reflect it
  
  #mix.best is used to indicate if the best overall loglikelihood (over
  #all iterations of CEM that were run) should be used
  
  n <- length(cem.obj$x[, 1])
  if (mix.best == F) {
    num.clusters <- length(unique(cem.obj$clust))
  }
  if (mix.best == T) {
    num.clusters <- length(unique(cem.obj$clust.best))
  }
  cluster.sizes <- rep(0, num.clusters)
  score <- 0
  if (mix.best == F) {
    for (i in 1:num.clusters) {
      cluster.sizes[i] <- sum(cem.obj$clust == unique(cem.obj$clust)[i])
    }
    
    #now the entries of cluster.sizes correspond to each column of the
    #likelihood matrix (regardless of ordering)
    for (i in 1:n) {
      this.lik <- 0
      for (j in 1:num.clusters) {
        this.lik <- this.lik + ((cluster.sizes[j] / n) * cem.obj$likelihoods[i, j])
      }
      score <- score + log(this.lik)
    }
  }
  if (mix.best == T) {
    score <- max(cem.obj$overall.loglikelihood[cem.obj$overall.loglikelihood != 0])
  }
  #assume n-1 curve clusters, 1 noise cluster.
  #estimate curves, mixture proportions, noise
  num.parameters <- ((num.clusters - 1) * (cem.obj$df + 2)) + (num.clusters - 1) + 1
  bic <- (2 * score) - (num.parameters * log(n))
  bic2 <- (2 * score) - (2 * num.parameters * log(n))
  #2= num of dimensions (bic2 is experimental)
  
  l3 = list(score, bic, bic2)
  names(l3)= c("score", "bic", "bic2")
  
  return (l3)
}

#------------------------------------------------------------------

#### CEM-PCC APPLIQUE SUR 3 CLUSTERS AVEC 8 DEGRES DE LIBERTE

data_features_points_etiquette = data.frame(data_sans_bruit, m[1] ) 
names(data_features_points_etiquette) = c("xx", "yy", "cc")
data_features_points_etiquette

data_blackground_noise_etiquette = data.frame(data_blackground_noise, rep(0, length(data_blackground_noise[,1] ))  )
names(data_blackground_noise_etiquette) = c("xx", "yy", "cc")
data_blackground_noise_etiquette


data_final_pour_cem_pcc =  rbind(data_blackground_noise_etiquette,data_features_points_etiquette)

algo_cem_pcc =  cem( x= as.matrix(data_final_pour_cem_pcc[c(1,2)]) ,clust.init= as.matrix(data_final_pour_cem_pcc[3] ),df =8 ,plot.me = T, mix=T )
algo_cem_pcc$overall.loglikelihood


vraisem = calc.likelihood(algo_cem_pcc, mix.best=T)
vraisem$bic


##########################################################
##########################################################
#### BOUCLE POUR AFFICHER LES VALEURS DU BIC 

degre_libre =20
nb_cluster = 9

matrice_bic = matrix( rep(0,degre_libre*nb_cluster), ncol =nb_cluster)
k=0

for( i in 2: degre_libre ){
  
  
  for( j in 1:nb_cluster ){
    
    hpcc_boucle = hpcc(x=as.matrix(data_sans_bruit) ,clust.init=mod1$classification, alpha= .1, plot.iter = F, df=i, plot.result = F, force=T, forcenum=j) # pas besoin de dégrés de liberté
    
    data_features_points_etiquette_dd = data.frame(data_sans_bruit, hpcc_boucle[1] ) 
    names(data_features_points_etiquette_dd) = c("xx", "yy", "cc")
    data_features_points_etiquette_dd
    
    data_blackground_noise_etiquette = data.frame(data_blackground_noise, rep(0, length(data_blackground_noise[,1] ))  )
    names(data_blackground_noise_etiquette) = c("xx", "yy", "cc")
    data_blackground_noise_etiquette
    
    data_final_pour_cem_pcc_dd =  rbind(data_blackground_noise_etiquette, data_features_points_etiquette_dd)
    
    algo_cem_pcc =  cem( x= as.matrix(data_final_pour_cem_pcc_dd[c(1,2)]) ,clust.init= as.matrix(data_final_pour_cem_pcc_dd[3] ),df=i, plot.me = F, mix=T )
    
    matrice_bic[i,j]= calc.likelihood(algo_cem_pcc, mix.best=T)$bic
    
    print(calc.likelihood(algo_cem_pcc, mix.best=T)$bic)
    k=k+1
    print(k)
  }
  
}



###########
# POUR LE CLUSTER 0

algo_cem_pcc =  cem( x= as.matrix(data_final_pour_cem_pcc[c(1,2)]) ,clust.init= rep(0, 1034 ),df =20 ,plot.me = T, mix=T )
algo_cem_pcc$overall.loglikelihood

vraisem = calc.likelihood(algo_cem_pcc, mix.best=T)
vraisem$bic




