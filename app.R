#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            TrapSim Tool - Multi Method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Code developed by Andrew Gormley, Manaaki Whenua - Landcare Research
#Developed for the NSC Eco-economics project

#Includes Audrey's grid based approach as well - Audrey says it works with cells of 500m and even 200m...(?)
#28/8 - some testing would suggest that for nightly checking, setting the grid square dimension at 4*sigma is best when compared to IBM


library("shiny")
library("shinythemes")
library("leaflet")
library("maptools")
library("spatstat") #For the owin command...and for runifpoint
library("RColorBrewer")
library("leaflet")
library("rgdal")
library("rgeos")
library("proxy")
library("sf")
library("raster")
library("DT")

# setwd("C:\\Users\\gormleya\\OneDrive - MWLR\\Documents\\CAEM\\IslandConservation\\TrapSimFeasibility\\Shiny")

# def.shp<-"Robinson_Coati"  #The default shape.
def.shp<-"WhakatipuMahia"
# shp.zones<-"Robinson_Coati_Workzones"
ras.2<-raster("habitat_specific.asc")
# ras.z2<-raster("huntingeffort.asc")

# cols.tst<-brewer.pal(6,"Blues")
cols.vec<-brewer.pal(8,"Paired")
cols.eff<-brewer.pal(8,"Reds")
proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
#4326 is the EPSG for WGS84...

#~~~~~~~~~~~~~~~ A bunch o' functions ~~~~~~~~~~~~~~~~~~~~~~~~
#Define the resampling function
resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 

trap.cost.func<-function(a,b,c,d,e){
  #e.g. trap.cost.func(a=number of checks, b=n.traps, c=input$traps.per.day, d=input$day.rate,e=cost.per.trap)
  fixed.cost<-b*e
  
  # a<-2
  # b<-100.5
  # c<-50
  # d<-100
  # ceiling(2*b/c)/2*a*d
  
  labor.cost<-ceiling(2*b/c)/2*a*d #- scales up to a half day, mostly works 
  trap.cost<-as.integer(labor.cost+fixed.cost)
  # trap.cost<-as.integer((((b*checks)/c)*d)+(b*e))
  return(trap.cost)
}

hunt.cost.func<-function(a,b){
  day.rate<-a
  days.zone<-b
  hunt.cost<-day.rate*sum(days.zone)
  return(hunt.cost)
}

get.alpha.beta<-function(g0.mean, g0.sd){
  g0.var<-g0.sd^2
  paren<-    (g0.var + g0.mean^2 - g0.mean)
  alpha<- -g0.mean*paren/g0.var
  beta<- paren*(g0.mean-1)/g0.var
  return(list(alpha=alpha, beta=beta))
}

get.loc.shape<-function(sigma.mean, sigma.sd){
  m<-sigma.mean
  s<-sigma.sd
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  return(list(location=location, shape=shape))
}

#Function to make animal locations from a raster of relative abundance/habitat.
get.pest.locs<-function(ras, n.poss, shp){
  
  gridpolygon<-rasterToPolygons(ras)
  grids<-length(gridpolygon@data[,1])
  pest.per.grid<-round(gridpolygon@data[,1]*n.poss/grids,0)   #We want to sample more than we need becase we are going to remove some that are outside the shapefile
  coords<-vector("list", grids) #Set up to save the coordinates - cant work out how to get spsample to sample across all grids
  for (i in 1: grids){
    if(pest.per.grid[i]>0){
      coords[[i]]<-spsample(gridpolygon@polygons[[i]], n =pest.per.grid[i] , type = "random")@coords
    }
  }
  coords<-as.data.frame(do.call("rbind", coords)) #Put them altogether from the list
  coords<-coords[inside.owin(coords[,1], coords[,2], shp),]  #Remove the ones from outside the shapefile...
  coords<-coords[sample(1:dim(coords)[1], size=n.poss, replace=FALSE),]  #Then sample to get the desired actual sample size
  
  return(coords)
  
}


#  Make the trap locations from a spacing, buffer and shapfile
make.trap.locs<-function(x.space,y.space,buff,shp){
  
  b.box<-bbox(shp)
  
  traps.x<-seq(from=b.box[1,1], by=x.space, to=b.box[1,2])
  traps.y<-seq(from=b.box[2,1], by=y.space, to=b.box[2,2])
  
  traps<-as.data.frame((expand.grid(traps.x, traps.y)))
  colnames(traps)<-c("X","Y")
  shp.buff<-gBuffer(shp,width=-buff)
  #Remove traps that are outside the window,,,
  traps<-traps[inside.owin(traps[,1], traps[,2], shp.buff),]
  
  return(traps)
}


# Function to split the input boxes for scenarios
split.str<-function(a){
  b<-as.numeric((strsplit(as.character(a),split= "/"))[[1]])
  return(b)
}
#~~~~~~~~~~~~~~~~~  End functions ~~~~~~~~~~~~~~~~~~~~~~


server<-function(input, output, session) {
  
  
  #~~~~~~ Set up the scenarios
  scenParam <- reactiveVal()
  
  #Delete all scenarios 
  observeEvent(input$deleteAllRows,{
    t<-scenParam()
    t <- t[-(1:dim(t)[1]),]
    scenParam(t)
  })
  
  #Delete selected scenarios
  observeEvent(input$deleteRows,{
    t<-scenParam()
    if (!is.null(input$tableDT_rows_selected)) {
      t <- t[-as.numeric(input$tableDT_rows_selected),]
      rownames(t)<-1:dim(t)[1]
    }
    scenParam(t)
  })
  
  observeEvent(input$update,{
    # validate(
    #   need(input$trap_methods==TRUE | input$hunt_methods==TRUE,"You need at least trapping or hunting")
    # )
    x.space.a = NA
    y.space.a = NA
    trap.start.a = NA
    trap.nights.a = NA
    check.interval.a = NA
    g.mean.a=NA
    g.zero.a=NA    #This is ther proportion that is untrappable, not g0....
    x.space.b = NA
    y.space.b = NA
    trap.start.b = NA
    trap.nights.b = NA
    check.interval.b = NA
    g.mean.b=NA
    g.zero.b=NA
    
    if(input$trap_methods==1){    
      x.space.a = input$traps.x.space.a
      y.space.a = input$traps.y.space.a
      trap.start.a = input$trap.start.a
      trap.nights.a = input$trap.nights.a
      check.interval.a = input$n.check.a
      g.mean.a=input$g0.mean.a
      g.zero.a=input$g0.zero.a
      if(input$show_trap_b==1){
        x.space.b = input$traps.x.space.b
        y.space.b = input$traps.y.space.b
        trap.start.b = input$trap.start.b
        trap.nights.b = input$trap.nights.b
        check.interval.b = input$n.check.b
        g.mean.b=input$g0.mean.b
        g.zero.b=input$g0.zero.b
      }
    }
    
    
    
    # #Hunting    
    # hunt.start.a = NA 
    # hunt.days.a = NA
    # hunt.eff.a = NA
    # hunt.start.b = NA 
    # hunt.days.b = NA
    # hunt.eff.b = NA
    # if(input$hunt_methods==1){    
    #   hunt.start.a = input$hunt.start.a
    #   hunt.days.a = input$hunt.days.a
    #   hunt.eff.a = input$effort.a 
    #   
    #   if(input$show_hunt_b==1){
    #     # if(is.na(hunt.start.b)==FALSE){
    #     hunt.start.b = input$hunt.start.b
    #     hunt.days.b = input$hunt.days.b
    #     hunt.eff.b = input$effort.b 
    #   }
    # }
    
    # 
    
    
    bait.start.a = NA
    bait.nights.a = NA
    bait.check.a = NA
    bait.x.space.a = NA
    bait.y.space.a = NA
    bait.g.mean.a = NA
    bait.g.zero.a = NA
    if(input$bait_methods==1){
      bait.start.a = input$bait.start.a
      bait.nights.a = input$bait.nights.a
      bait.check.a = input$bait.check.a
      bait.x.space.a = input$bait.x.space.a
      bait.y.space.a = input$bait.y.space.a
      bait.g.mean.a = input$bait.g0.mean.a
      bait.g.zero.a = input$bait.g.zero.a
    }
    
    to_add <- data.frame(
      x.space.a = x.space.a,
      y.space.a = y.space.a,
      trap.start.a = trap.start.a,
      trap.nights.a = trap.nights.a,
      check.interval.a = check.interval.a,
      g.mean.a=g.mean.a,
      g.zero.a=g.zero.a,
      x.space.b = x.space.b,
      y.space.b = y.space.b,
      trap.start.b = trap.start.b,
      trap.nights.b = trap.nights.b,
      check.interval.b = check.interval.b,
      g.mean.b=g.mean.b,
      g.zero.b=g.zero.b,
      bait.start.a = bait.start.a,
      bait.nights.a = bait.nights.a,
      bait.check.a = bait.check.a,
      bait.x.space.a = bait.x.space.a,
      bait.y.space.a = bait.y.space.a,
      bait.g.mean.a = bait.g.mean.a,
      bait.g.zero.a = bait.g.zero.a
      # hunt.start.a = hunt.start.a,
      # hunt.days.a = hunt.days.a,
      # hunt.eff.a = hunt.eff.a,
      # hunt.start.b = hunt.start.b,
      # hunt.days.b = hunt.days.b,
      # hunt.eff.b = hunt.eff.b
      
    )
    newScenParam <- rbind(scenParam(),to_add) # adding new data
    scenParam(newScenParam) # updating data
    
    #Test for duplicated
    t<-scenParam()
    t<-t[!duplicated(t), ]        #Dont add duplicates
    t<-t[!rowSums(is.na(t))==21,]  #Dont add blank scenrios
    scenParam(t)
    
    # scenParam<-rbind(scenParam,to_add)
    return(list(scenParam=scenParam))
  })
  
  
  #Read in the shapefile  
  mydata.shp<-reactive({
    shp<-readOGR("Shapefiles",def.shp)
    
    if(input$area_type=="Map"){
      if(is.null(input$shp.file)==FALSE){
        
        myshape<-input$shp.file
        dir<-dirname(myshape[1,4])
        for ( i in 1:nrow(myshape)) {
          file.rename(myshape[i,4], paste0(dir,"/",myshape[i,1]))
        }
        getshp <- list.files(dir, pattern="*.shp", full.names=TRUE)
        
        shp<-readShapePoly(getshp,proj4string=CRS(proj4string))
        shp<-gBuffer(shp, width=1) #Can fix some orphaned holes issues
      }
    }
    proj4string<-crs(shp)  #get the proj string from the shapefile itself
    shp<-gBuffer(shp, width=1) #Can fix some orphaned holes issues
    ha<-sapply(slot(shp, "polygons"), slot, "area")/10000  
    return(list(shp=shp, p4s=proj4string, ha=ha))
  })
  
  
  # #Two hunting methods. Calculate the probability of kill based on k/rho and effort.
  # #! Currently not used.
  # mydata.hunt.prob<-reactive({
  #   effort.zone<-c(input$effort.a, input$effort.b)
  #   hunt.rho<-c(input$hunt.rho.a, input$hunt.rho.b)
  #   hunt.k<-c(input$hunt.k.a, input$hunt.k.b)
  #   theta.hat<-1-exp(-((hunt.rho*log(effort.zone))^hunt.k))
  #   p.hunt.day<-theta.hat
  #   
  #   return(list(p.hunt.day=p.hunt.day))
  #   
  # })
  
  
  
  # Make the traps and animals on the interactive map. Not linked to the actual simulation
  mydata.map<-reactive({
    shp<-mydata.shp()$shp    
    traps.x.space<-as.numeric(input$traps.x.space.i)
    traps.y.space<-as.numeric(input$traps.y.space.i)
    buff<-as.numeric(input$traps.buff.i)
    traps<-make.trap.locs(traps.x.space, traps.y.space, buff, shp)
    
    #Temporary pest animal coordinates...
    n.poss<-input$numb.poss#.i
    #Temporary commented out.    
    # if(is.null(input$ras.1)==TRUE){
    if((input$ras.1)=="Random"){
      #1. Random locations.
      n.poss.tmp<-(runifpoint(n.poss,shp))
      animals.xy.ini<-as.data.frame(n.poss.tmp)
    }else{
      #2. Grid specific densities
      # infile.ras <-input$ras.1
      # ras.2<-raster(infile.ras$datapath)
      #Call the function
      animals.xy.ini<-get.pest.locs(ras.2, n.poss, shp)
    }
    colnames(animals.xy.ini)<-c("X","Y")
    
    return(list(traps=traps, animals.xy.ini=animals.xy.ini))
  })
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~  Now simulate the actual trapping...~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #    When the simulations get run, automatically switch to the results tab. Noice!
  observeEvent(input$act.btn.trapsim,{
    updateTabsetPanel(session, "inTabset",selected = "4. Results")
  })
  
  datab<-eventReactive(input$act.btn.trapsim,{
    buffer<-100  #A single buffer
    #Some validation stuff - need to fill this out more completely.
    validate(
      # need(input$traps.buff.a != "", "Please enter a value for the trap buffer"),
      need(input$traps.x.space.a != "", "Please enter a value for the X trap spacing"),
      need(input$traps.y.space.a != "", "Please enter a value for the Y trap spacing"),
      need(input$n.nights != "", "Please enter a value for the number of Nights"),
      need(input$n.check.a != "", "Please enter a value for the Check interval"),
      need(input$p.bycatch.a != "", "Please enter a value for the Daily bycatch rate"),
      need(input$max.catch.a != "", "Please enter a value for the prob of Max catch")
    )
    

    #Get the parameter values for the scenarios

    params<-as.data.frame(scenParam())
    #Set up some places to store results 
    n.scen<-dim(params)[1]    
    pop.size.list<-vector("list",n.scen)
    pop.zone.list<-vector("list",n.scen)
    trap.catch.list<-vector("list",n.scen)
    hunt.catch.list<-vector("list",n.scen)
    bait.catch.list<-vector("list",n.scen)
    # pop.size.zone.mat[[ii]][,t+1]
    
    withProgress(message="Running simulation ",value=0,{
      
      for(kk in 1:n.scen){  #For each scenario...
        incProgress(kk/n.scen, detail = paste("Doing scenario ", kk," of", n.scen))
        
        #Pass the parameters from params to a parameter name.
        # #Trapping period
        trap.start.a<-params$trap.start.a[kk]
        trap.nights.a<-params$trap.nights.a[kk]
        n.check.a<-params$check.interval.a[kk]
        x.space.a<-params$x.space.a[kk]
        y.space.a<-params$y.space.a[kk]
        buffer.a<-buffer
        g0.mean.a<-params$g.mean.a[kk]
        g.zero.a<-params$g.zero.a[kk]
        
        trap.start.b<-params$trap.start.b[kk]
        trap.nights.b<-params$trap.nights.b[kk]
        n.check.b<-params$check.interval.b[kk]
        x.space.b<-params$x.space.b[kk]
        y.space.b<-params$y.space.b[kk]
        buffer.b<-buffer
        g0.mean.b<-params$g.mean.b[kk]
        g.zero.b<-params$g.zero.b[kk]
        
        
        bait.start.a<-params$bait.start.a[kk]
        bait.nights.a<-params$bait.nights.a[kk]
        bait.check.a<-params$bait.check.a[kk]
        bait.x.space.a<-params$bait.x.space.a[kk]
        bait.y.space.a<-params$bait.y.space.a[kk]
        bait.buff.a<-buffer
        bait.g0.mean.a<-params$bait.g.mean.a[kk]
        bait.g.zero.a<- params$bait.g.zero.a[kk]
        
        
        # When we have traps//baits etc, then calculate the interval and the checking interval and bycatch/failure
        if(is.na(trap.start.a)==FALSE){
          trap.period.a<-seq(from=trap.start.a, to=(trap.start.a+trap.nights.a-1), by=1)
          # #This sets the trap checking interval. i.e. traps are cleared and reset on these nights only...
          check.vec.a<-seq(from=trap.start.a, to=(trap.start.a+trap.nights.a), by=n.check.a)
          p.bycatch.a<-input$p.bycatch.a
          
          # trap.start.a<-1
          # trap.nights.a<-10
          # n.check.a<-11
          # (trap.period.a<-seq(from=trap.start.a, to=(trap.start.a+trap.nights.a-1), by=1))
          # (check.vec.a<-seq(from=trap.start.a, to=(trap.start.a+trap.nights.a), by=n.check.a))
          
        }
        # if(input$show_trap_b==1){
        if(is.na(trap.start.b)==FALSE){
          trap.period.b<-seq(from=trap.start.b, to=(trap.start.b+trap.nights.b-1), by=1)
          check.vec.b<-seq(from=trap.start.b, to=(trap.start.b+trap.nights.b), by=n.check.b)
          p.bycatch.b<-input$p.bycatch.b
        }
        
        if(is.na(bait.start.a)==FALSE){
          bait.period.a<-seq(from=bait.start.a, to=(bait.start.a+bait.nights.a-1), by=1)
          # #This sets the trap checking interval. i.e. traps are cleared and reset on these nights only...
          bait.check.vec.a<-seq(from=bait.start.a, to=(bait.start.a+bait.nights.a), by=bait.check.a)
          p.failure.a<-0#input$p.bycatch.a
        }
        
        
        #How long to run the simulation for.
        n.nights<-input$n.nights

        #The g0 uncertainty, and sigma values for traps and bait - not in params - but should be, 
        g0.sd.a<-input$g0.sd.a
        g0.sd.b<-input$g0.sd.b
        bait.g0.sd.a<-input$bait.g0.sd.a
        
        sigma.mean<-input$sigma.mean 
        sigma.sd<-input$sigma.sd
        rmax.poss<-input$rmax.poss
        K.poss<-10 #input$K.poss
        #When the reproductive period starts and how long it lasts (i.e. spread out over...)
        rep.start<-input$rep.start
        rep.nights<-input$rep.nights
        rep.tmp<-seq(from=rep.start, by=1, to=(rep.nights+rep.start-1))
        aa<-ceiling(n.nights/365) #basically the number of years
        bb<-((1:aa)-1)*365
        repro.interval<-rep(rep.tmp,aa)+rep(bb, each=rep.nights)
        rep.start.vec<-(bb+rep.start)
        
        ha<-mydata.shp()$ha
        shp<-mydata.shp()$shp
        # shp.2<-mydata.zone()$shp.2

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #Make the trap & bait station locations 
        if(is.na(trap.start.a)==FALSE){
          traps.a<-make.trap.locs(x.space.a, y.space.a,buffer.a,shp)
          coordinates(traps.a) <- c( "X", "Y" )
          proj4string(traps.a) <- CRS(proj4string)
          traps.xy.a<-as.data.frame(traps.a)
          n.traps.a<-dim(traps.xy.a)[1]
        }
        
        # if(input$show_trap_b==1){
        if(is.na(trap.start.b)==FALSE){
          traps.b<-make.trap.locs(x.space.b, y.space.b, buffer.b,shp)
          coordinates(traps.b) <- c( "X", "Y" )
          proj4string(traps.b) <- CRS(proj4string)
          traps.xy.b<-as.data.frame(traps.b)
          n.traps.b<-dim(traps.xy.b)[1]
        }
        
        #Make the bait station locations
        if(is.na(bait.start.a)==FALSE){
          baits.a<-make.trap.locs(bait.x.space.a, bait.y.space.a, bait.buff.a, shp)
          coordinates(baits.a) <- c( "X", "Y" )
          proj4string(baits.a) <- CRS(proj4string)
          baits.xy.a<-as.data.frame(baits.a)
          n.baits.a<-dim(baits.xy.a)[1]
        }
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        #Carrying capacity for the area
        K.tot<-K.poss*ha
        
        n.poss<-input$numb.poss
        
        max.catch.a<-input$max.catch.a
        max.catch.b<-input$max.catch.b
        
        #~~~~Calculate the costs~~~~
        trap.cost.sim<-0
        bait.cost.sim<-0
        if(is.na(trap.start.a)==FALSE){
          checks<-ceiling(input$trap.nights.a/input$n.check.a)+1
          # trap.cost.sim<-trap.cost.func(a=check.vec.a, b=n.traps.a, c=input$traps.per.day.a, d=input$day.rate.a, e=input$cost.per.trap.a)
          trap.cost.sim<-trap.cost.func(a=checks, b=n.traps.a, c=input$traps.per.day.a, d=input$day.rate.a, e=input$cost.per.trap.a)
        }
        # if(input$show_trap_b==1){
        if(is.na(trap.start.b)==FALSE){
          checks<-ceiling(input$trap.nights.b/input$n.check.b)+1
          # trap.cost.sim<-trap.cost.sim+trap.cost.func(a=check.vec.b, b=n.traps.b, c=input$traps.per.day.b, d=input$day.rate.b, e=input$cost.per.trap.b)
          trap.cost.sim<-trap.cost.sim+trap.cost.func(a=checks, b=n.traps.b, c=input$traps.per.day.b, d=input$day.rate.b, e=input$cost.per.trap.b)
        }
        
        if(is.na(bait.start.a)==FALSE){
          checks<-ceiling(input$bait.nights.a/input$bait.check.a)+1
          # bait.cost.sim<-trap.cost.func(a=bait.check.vec.a, b=n.baits.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a)
          bait.cost.sim<-trap.cost.func(a=checks, b=n.baits.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a)        
          }
        
        
        
        if (input$sim_type=='grid'){
          cell.width<-input$cell.width
          cell.area.m2<-cell.width^2
          
          #Make the grid...
          shp<-st_as_sf(shp)
          shp_grid<-shp%>%st_make_grid(cellsize=cell.width, what="polygons")%>%st_intersection(shp)
          grid.traps.master<-lengths(st_intersects(shp_grid, st_as_sf(traps)))#*max.catch  #Number in each cell
          
        }
        
        n_its<-input$n.its  #Number of iterations
        pop.size.mat<-matrix(NA,nrow=n_its,ncol=n.nights+1)
        trap.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of trapping - for each iteration - how many that night/day
        bait.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of baiting - for each iteration - how many that night/day
        hunt.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of hunting
        pop.size.zone.vec<-vector("list",n_its)
        
        
        for(ii in 1:n_its){
          pop.size.mat[ii,1]<-n.poss
          
          #~~~~~~~~~Make some animals~~~~~~~~~~
          # if(is.null(input$ras.1)==TRUE){
          if((input$ras.1)=="Random"){
            #1. Random locations.
            n.poss.tmp<-(runifpoint(n.poss,shp))
            animals.xy<-as.data.frame(n.poss.tmp)
          }else{
            #2. Grid specific densities
            animals.xy<-get.pest.locs(ras.2, n.poss, shp)
          }

          colnames(animals.xy)<-c("X","Y")
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
          # n.poss.tmp<-(runifpoint(n.poss,shp))
          # animals.xy<-as.data.frame(n.poss.tmp)
          # colnames(animals.xy)<-c("X","Y")
          n.animals<-dim(animals.xy)[1]
          animals.xy$Dead<-0
          coordinates(animals.xy) <- ~ X+Y
          proj4string(animals.xy) <- proj4string(shp)
          
          #g0 for traps
          if(is.na(trap.start.a)==FALSE){
          alpbet.a<-get.alpha.beta(g0.mean.a, g0.sd.a)
          animals.xy$g0.a<-rbeta(n.animals, alpbet.a$alpha, alpbet.a$beta)
          animals.xy$g0.a[animals.xy$g0.a<0]<-0 #Probably not needed
          animals.xy$g0.a[sample(x=n.animals,size=round(n.animals*g.zero.a),replace=F)]<-0  #Make some of the g0 values zero - untrappable animals...
          }else{
            animals.xy$g0.b<-0
          }
          #g0 for trap type 2
          if(is.na(trap.start.b)==FALSE){
            alpbet.b<-get.alpha.beta(g0.mean.b, g0.sd.b)
            animals.xy$g0.b<-rbeta(n.animals, alpbet.b$alpha, alpbet.b$beta)
            animals.xy$g0.b[animals.xy$g0.b<0]<-0 #Probably not needed
            # if(is.na(g.zero.b)==FALSE){
              animals.xy$g0.b[sample(x=n.animals,size=round(n.animals*g.zero.b),replace=F)]<-0
            # }
          }else{
            animals.xy$g0.b<-0
          }
          
          #g0 for trap type 2
          if(is.na(bait.start.a)==FALSE){
            alpbet.bait<-get.alpha.beta(bait.g0.mean.a, bait.g0.sd.a)
            animals.xy$g0.bait<-rbeta(n.animals, alpbet.bait$alpha, alpbet.bait$beta)
            animals.xy$g0.bait[animals.xy$g0.bait<0]<-0 #Probably not needed
            animals.xy$g0.bait[sample(x=n.animals,size=round(n.animals*bait.g.zero.a),replace=F)]<-0
            
            # if(is.na(g.zero.b)==FALSE){
            #   animals.xy$g0.b[sample(x=n.animals,size=round(n.animals*g.zero.b),replace=F)]<-0
            # }
          }else{
            animals.xy$g0.bait<-0
          }
          # bait.g0.mean.a
          
          
          locshp<-get.loc.shape(sigma.mean, sigma.sd)
          animals.xy$Sigma<-rlnorm(n.animals ,meanlog=locshp$location, sdlog=locshp$shape)
          
          #The first one initialises for a grid based simulation - i.e. Audrey's model.
          #The seond is the trap+animal pairwise 
          if (input$sim_type=='grid'){
            animals.SP<-animals.xy
            coordinates(animals.SP) <- c( "X", "Y" )
            proj4string(animals.SP) <- CRS(proj4string)
            which.grid<-st_within(st_as_sf(animals.SP), shp_grid)
            animals.xy$CellIndex<-as.data.frame(which.grid)$col.id
            animals.xy$PreProb<-(2*pi*animals.xy$g0.a*animals.xy$Sigma^2)/cell.area.m2   #The pre probability...
            grid.traps<-grid.traps.master
            grid.traps.capacity<-grid.traps.master*max.catch.a
          }else{
            if(is.na(trap.start.a)==FALSE){
              dist2.xy.a<-matrix(NA,n.traps.a,n.animals)
              prob.xy.a<-matrix(0,n.traps.a,n.animals)
              dist.xy.a<-dist(as.data.frame(traps.xy.a), as.data.frame(animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
              prob.xy.a<-exp(-(dist.xy.a^2)/(2*animals.xy$Sigma^2))*animals.xy$g0.a #Use the g0u for sampling...
              rm(dist.xy.a)
            }
            # if(input$show_trap_b==1){
            if(is.na(trap.start.b)==FALSE){
              dist2.xy.b<-matrix(NA,n.traps.b,n.animals)
              prob.xy.b<-matrix(0,n.traps.b,n.animals)
              dist.xy.b<-dist(as.data.frame(traps.xy.b), as.data.frame(animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
              prob.xy.b<-exp(-(dist.xy.b^2)/(2*animals.xy$Sigma^2))*animals.xy$g0.b #Use the g0u for sampling...
              rm(dist.xy.b)
            }
            #Bait stations
            if(is.na(bait.start.a)==FALSE){
              dist2.xy.c<-matrix(NA,n.baits.a,n.animals)
              prob.xy.c<-matrix(0,n.baits.a,n.animals)
              dist.xy.c<-dist(as.data.frame(baits.xy.a), as.data.frame(animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
              prob.xy.c<-exp(-(dist.xy.c^2)/(2*animals.xy$Sigma^2))*animals.xy$g0.bait #Use the g0u for sampling...
              rm(dist.xy.c)
            }
          }
          
          
          #~~~Now run the simulation of trapping the animals...
          # withProgress(message="Running Simulation",value=0,{
          trap.catch.vec<-rep(0, n.nights)
          if(is.na(trap.start.a)==FALSE){
            trap.catch.a<-matrix(0,n.traps.a,n.nights)  			#Trap.catch stores the captures in each trap each night
            trap.remain.a<-rep(max.catch.a, n.traps.a)
          }
          # if(input$show_trap_b==1){
          if(is.na(trap.start.b)==FALSE){
            trap.catch.b<-matrix(0,n.traps.b,n.nights)  			#Trap.catch stores the captures in each trap each night
            trap.remain.b<-rep(max.catch.b, n.traps.b)
          }
          # trap.remain<-rep(T,n.traps)						        #Record if the trap remains in operation...TRUE, FALSE
          if(is.na(bait.start.a)==FALSE){
            bait.catch.a<-matrix(0,n.baits.a,n.nights)  			#Trap.catch stores the captures in each trap each night
            # bait.remain.a<-rep(max.catch.b, n.traps.b)
          }
          
          
          
          pop.size.zone.vec[[ii]]<-matrix(NA,nrow=4,ncol=n.nights+1)
          # pop.size.zone.vec[[ii]][,1]<-table(over(animals.xy[animals.xy$Dead==0,], shp.2))
          
          
          for (t in 1:n.nights){								#For each night
            not.caught<-(1:n.animals)[animals.xy$Dead==0]			#Animals not already caught
            
            if(is.na(trap.start.a)==FALSE){
              if(t%in%trap.period.a==TRUE){
                if(t%in%check.vec.a==TRUE){#If it is a trap clearance day...then reset the traps to T *before* trappig starts!
                  # trap.remain<-rep(T,n.traps)		
                  trap.remain.a<-rep(max.catch.a,n.traps.a)
                  if (input$sim_type=='grid'){
                    grid.traps<-grid.traps.master
                  }
                  
                }
                #Turn off some of the traps according to the random probability
                # trap.remain[rbinom(n.traps,1, p.bycatch)==1]<-FALSE #Nedd to only turn off those that are on...Might be okay...
                trap.remain.a<-trap.remain.a-rbinom(n.traps.a, trap.remain.a, p.bycatch.a) #This modifcation deals with multiple capture traps 
                trap.remain.a[trap.remain.a<0]<-0
                
                if(sum(not.caught)>0){
                  for (j in not.caught){ 							#For each animal not already caught
                    if (input$sim_type=='grid'){
                      #Based on Audrey's model...
                      pcap<-1-exp(-(animals.xy$PreProb[j]*grid.traps[animals.xy$CellIndex[j]]))
                      if(rbinom(1,1,prob=pcap)==1){ #If animal gets caught
                        animals.xy$Dead[j]<-1
                        grid.traps[animals.xy$CellIndex[j]]<-grid.traps[animals.xy$CellIndex[j]]-1
                        trap.catch.vec[t]<-trap.catch.vec[t]+1
                      }
                    }else{
                      #New code based on Multinomial (from Dean...
                      prob.tmp<-prob.xy.a[,j]*(trap.remain.a>0)				#Adjust the probabilities with trap.remain so previous traps cannot catch anything
                      cumulative.capture.prob<-1-prod(1-prob.tmp) #Cumulative probability of possum i getting captured
                      
                      if(rbinom(1,1,prob=cumulative.capture.prob)==1){#If the animal is going to get caught...
                        trap.id<-match(1,rmultinom(1,1,prob.tmp))#Which was the successful trap...
                        trap.catch.a[trap.id,t]<-trap.catch.a[trap.id,t]+1
                        # trap.animals[j]<-1					#Set trapped animal to 1
                        animals.xy$Dead[j]<-1
                        # animals.xy$Day2[j]<-t
                        # trap.remain[trap.id]<-(trap.catch[trap.id,t]<max.catch)			#Calculate whether the trap is full or not!!
                        trap.remain.a[trap.id]<-trap.remain.a[trap.id]-1
                      }
                    }
                  }
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Module for immigration to go in here...
                }
              } #End of if t %in% trap.period
            }
            
            
            # if(input$show_trap_b==1){
            if(is.na(trap.start.b)==FALSE){
              if(t%in%trap.period.b==TRUE){
                if(t%in%check.vec.b==TRUE){#If it is a trap clearance day...then reset the traps to T *before* trappig starts!
                  # trap.remain<-rep(T,n.traps)		
                  trap.remain.b<-rep(max.catch.b,n.traps.b)
                }
                #Turn off some of the traps according to the random probability
                # trap.remain[rbinom(n.traps,1, p.bycatch)==1]<-FALSE #Nedd to only turn off those that are on...Might be okay...
                trap.remain.b<-trap.remain.b-rbinom(n.traps.b, trap.remain.b, p.bycatch.b) #This modifcation deals with multiple capture traps 
                trap.remain.b[trap.remain.b<0]<-0
                
                if(sum(not.caught)>0){
                  for (j in not.caught){ 							#For each animal not already caught
                    
                    #New code based on Multinomial (from Dean...
                    prob.tmp<-prob.xy.b[,j]*(trap.remain.b>0)				#Adjust the probabilities with trap.remain so previous traps cannot catch anything
                    cumulative.capture.prob<-1-prod(1-prob.tmp) #Cumulative probability of possum i getting captured
                    
                    if(rbinom(1,1,prob=cumulative.capture.prob)==1){#If the animal is going to get caught...
                      trap.id<-match(1,rmultinom(1,1,prob.tmp))#Which was the successful trap...
                      trap.catch.b[trap.id,t]<-trap.catch.b[trap.id,t]+1
                      # trap.animals[j]<-1					#Set trapped animal to 1
                      animals.xy$Dead[j]<-1
                      # animals.xy$Day2[j]<-t
                      # trap.remain[trap.id]<-(trap.catch[trap.id,t]<max.catch)			#Calculate whether the trap is full or not!!
                      trap.remain.b[trap.id]<-trap.remain.b[trap.id]-1
                    }
                    # }
                  }
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Module for immigration to go in here...
                }
              } #End of if t %in% trap.period
            }
            
            #Now do some baiting....
            
            if(is.na(bait.start.a)==FALSE){
              if(t%in%bait.period.a==TRUE){

                if(sum(not.caught)>0){
                  for (j in not.caught){ 							#For each animal not already caught
                    
                    #New code based on Multinomial (from Dean...
                    prob.tmp<-prob.xy.c[,j]#*(trap.remain.b>0)				#Adjust the probabilities with trap.remain so previous traps cannot catch anything
                    cumulative.capture.prob<-1-prod(1-prob.tmp) #Cumulative probability of possum i getting captured
                    
                    if(rbinom(1,1,prob=cumulative.capture.prob)==1){#If the animal is going to get caught...
                      trap.id<-match(1,rmultinom(1,1,prob.tmp))#Which was the successful trap...
                      bait.catch.a[trap.id,t]<-bait.catch.a[trap.id,t]+1
                      animals.xy$Dead[j]<-1
                    }
                  }
                }
              } #End of if t %in% trap.period
            }
            
            
            
            
            #Hunting
            
            
            
            hunt.vec<-c(0,0)
            # zones<-c("A","B","C","D")
            
            # zz<-ifelse(input$show_hunt_b==TRUE,2,1)
            
            
            
            # if(is.na(hunt.start.a)==FALSE){
            #   zz<-ifelse(is.na(hunt.start.b)==TRUE,1,2)
            #   for (z in 1:zz){
            #     if(days.zone[z]>0){
            #       if(t%in% hunt.periods[[z]] ){
            #         # idx.animal.zone<-inside.owin(animals.xy[,1], animals.xy[,2], shp.2[shp.2$WorkZone=='C',])
            #         # idx.animal.zone<-(inside.owin(as.data.frame(animals.xy)[,1], as.data.frame(animals.xy)[,2], shp.2[shp.2$WorkZone==zones[z],])&animals.xy$Dead==0)
            #         # N.tmp<-sum(idx.animal.zone)
            #         idx.animal<-(animals.xy$Dead==0)
            #         N.tmp<-sum(idx.animal)
            #         
            #         if(N.tmp>0){
            #           n.kill<-rbinom(n=1, p=p.hunt.day[z], size=N.tmp)
            #           hunt.vec[z]<-n.kill
            #           if(n.kill>0){
            #             idx.dead<-sample(which(idx.animal,TRUE), size=n.kill, replace=FALSE) 
            #             animals.xy$Dead[idx.dead]<-1
            #           }
            #         }
            #       }
            #     }
            #   }
            #   hunt.catch.mat[ii,t]<-sum(hunt.vec)
            # }else{
              hunt.catch.mat[ii,t]<-0
            # }
            
            #~~~~~~~~~~~~~~~~~~~~~Reproduction~~~~~~~~~~~~~~~~~~~~~~~~~
            #Get the reproductive vector when we hit the start of the interval...
            #If t is one of the start days of the breeding season...
            if(t%in%rep.start.vec){
              #Discrete version of rmax
              rd<-exp(rmax.poss)-1
              N0<-sum(animals.xy$Dead==0) #Current population
              K<-K.tot
              new.animals<-rd*N0*((K-N0)/K)#Number of new animals
              mu.animals<-new.animals/rep.nights #new animals per night
              if(mu.animals>0){ #Draw animals when mu is positive
                new.animal.vec<-(rpois(rep.nights, lambda=mu.animals))
                new.animal.vec[(N0+cumsum(new.animal.vec))>(K*1.1)]<-0
              }else{
                new.animal.vec<-rep(0,rep.nights) #new animals is 0 when mu is negative
              }
              yr<-match(t, rep.start.vec) #What year are we in...?
            }
            
            
            if(t%in%repro.interval==TRUE){
              #Get the index from the number of new animals
              t.idx<-t-(365*(yr-1)) #TO adjust for multi years
              N.new<-new.animal.vec[match(t.idx,repro.interval)]
              if(N.new>0){  
                #Make the new locations randomly...
                
                # if(is.null(input$ras.1)==TRUE){
                if((input$ras.1)=="Random"){
                  
                  #1. Random locations.
                  
                  new.animals.xy<-as.data.frame(runifpoint(N.new,shp))
                }else{
                  #2. Grid specific densities
                  # infile.ras <-input$ras.1
                  # ras.2<-raster(infile.ras$datapath)
                  #Call the function
                  new.animals.xy<-get.pest.locs(ras.2, N.new, shp)
                }
                
                # new.animals.xy<-as.data.frame(runifpoint(N.new,shp))
                colnames(new.animals.xy)<-c("X","Y")
                new.animals.SP<-new.animals.xy  #Why dis?
                coordinates(new.animals.SP) <- c( "X", "Y" )
                proj4string(new.animals.SP) <- CRS(proj4string)
                # new.animals.xy$g0<-g0.mean
                # new.animals.xy$Sigma<-sigma.mean
                new.animals.xy$Dead<-0
                new.animals.xy$g0.a<-rbeta(N.new, alpbet.a$alpha, alpbet.a$beta)
                new.animals.xy$g0.a[new.animals.xy$g0.a<0]<-0.000001
                new.animals.xy$g0.a[sample(x=N.new,size=round(N.new*g.zero.a),replace=F)]<-0
                
                if(is.na(trap.start.b)==FALSE){
                  new.animals.xy$g0.b<-rbeta(N.new, alpbet.b$alpha, alpbet.b$beta)
                  new.animals.xy$g0.b[new.animals.xy$g0.b<0]<-0.000001
                  if(is.na(g.zero.b)==FALSE){
                    new.animals.xy$g0.b[sample(x=N.new,size=round(N.new*g.zero.b),replace=F)]<-0
                  }
                }else{
                  new.animals.xy$g0.b<-0
                }
                
                if(is.na(bait.start.a)==FALSE){
                  # alpbet.bait<-get.alpha.beta(bait.g0.mean.a, bait.g0.sd.a)
                  new.animals.xy$g0.bait<-rbeta(N.new, alpbet.bait$alpha, alpbet.bait$beta)
                  new.animals.xy$g0.bait[new.animals.xy$g0.bait<0]<-0 #Probably not needed
                }else{
                  new.animals.xy$g0.bait<-0
                }
                
                
                new.animals.xy$Sigma<-rlnorm(N.new ,meanlog=locshp$location, sdlog=locshp$shape)        
                
                
                
                
                
                if (input$sim_type=='grid'){
                  new.which.grid<-st_within(st_as_sf(new.animals.SP), shp_grid)
                  new.animals.xy$CellIndex<-as.data.frame(new.which.grid)$col.id
                  new.animals.xy$PreProb<-(2*pi*new.animals.xy$g0.a*new.animals.xy$Sigma^2)/cell.area.m2   #The pre probability...
                }else{
                  # #Calculate the probabilities for these new ones...
                  if(is.na(trap.start.a)==FALSE){
                    dist.xy.new<-matrix(NA,n.traps.a,N.new)
                    prob.xy.new<-matrix(0,n.traps.a,N.new)
                    dist.xy.new<-dist(as.data.frame(traps.xy.a), as.data.frame(new.animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
                    prob.xy.new<-exp(-(dist.xy.new^2)/(2*new.animals.xy$Sigma^2))*new.animals.xy$g0.a #Use the g0u for sampling...
                    rm(dist.xy.new)
                    prob.xy.a<-cbind(prob.xy.a, prob.xy.new)
                  }
                  
                  if(is.na(trap.start.b)==FALSE){
                    dist2.xy.new.b<-matrix(NA,n.traps.a,N.new)
                    prob.xy.new.b<-matrix(0,n.traps.b,N.new)
                    dist.xy.new.b<-dist(as.data.frame(traps.xy.b), as.data.frame(new.animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
                    prob.xy.new.b<-exp(-(dist.xy.new.b^2)/(2*new.animals.xy$Sigma^2))*new.animals.xy$g0.b #Use the g0u for sampling...
                    rm(dist.xy.new.b)
                    prob.xy.b<-cbind(prob.xy.b, prob.xy.new.b)
                  }
                  
                  if(is.na(bait.start.a)==FALSE){
                    dist2.xy.new.c<-matrix(NA,n.baits.a,N.new)
                    prob.xy.new.c<-matrix(0,n.baits.a,N.new)
                    dist.xy.new.c<-dist(as.data.frame(baits.xy.a), as.data.frame(new.animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
                    prob.xy.new.c<-exp(-(dist.xy.new.c^2)/(2*new.animals.xy$Sigma^2))*new.animals.xy$g0.bait #Use the g0u for sampling...
                    rm(dist.xy.new.c)
                    prob.xy.c<-cbind(prob.xy.c, prob.xy.new.c)
                  }
                  
                }
                
                coordinates(new.animals.xy) <- ~ X+Y
                proj4string(new.animals.xy) <- proj4string(shp)
                animals.xy<-rbind(animals.xy,new.animals.xy)
                
                n.animals<-dim(animals.xy)[1]
              }}
            
            #Update the age-class...
            pop.size.mat[ii,t+1]<-sum(animals.xy$Dead==0)
            #This calculates the number of alive animals in each zone.  
            # pop.size.zone.vec[[ii]][,t+1]<-table(over(animals.xy[animals.xy$Dead==0,], shp.2))
            
            
            
          }	 #End of the night
          if (input$sim_type=='grid'){
            trap.catch.mat[ii,]<-(trap.catch.vec)
          }else{
            if(is.na(trap.start.a)==FALSE){
              trap.catch.mat[ii,]<-colSums(trap.catch.a)
            }else{
              trap.catch.mat[ii,]<-rep(0, n.nights)
            }
            # if(input$show_trap_b==1){
            if(is.na(trap.start.b)==FALSE){
              trap.catch.mat[ii,]<-colSums(trap.catch.a)+colSums(trap.catch.b)
            }
            if(is.na(bait.start.a)==FALSE){
              bait.catch.mat[ii,]<-colSums(bait.catch.a)
            }else{
              bait.catch.mat[ii,]<-rep(0, n.nights)
            }
          }
        }#End of iteration ii
        
        
        # pop.zone.list[[kk]]<-apply(simplify2array(pop.size.zone.vec),c(1,2), mean)
        params$TrapCost[kk]<-trap.cost.sim
        params$BaitCost[kk]<-bait.cost.sim
        # params$HuntCost[kk]<-hunt.cost.sim
        # params$TotalCost[kk]<-hunt.cost.sim + trap.cost.sim
        params$TotalCost[kk]<-trap.cost.sim+bait.cost.sim
        
        params$MeanPopSize[kk]<-round(mean(pop.size.mat[,n.nights+1]),2)
        
        pop.size.list[[kk]]<-pop.size.mat
        hunt.catch.list[[kk]]<-hunt.catch.mat
        trap.catch.list[[kk]]<-trap.catch.mat
        bait.catch.list[[kk]]<-bait.catch.mat
        # cat(file=stderr(), "drawing histogram with", kk, "bins", "\n")
      } #End kk    
      
      
      #Add a comment 
    })  #End of progress
    
    
    
    params<-params[,c(25,24,22,23,1:21)]
    
    return(list(trap.catch.mat=trap.catch.mat, bait.catch.mat=bait.catch.mat, pop.size.mat=pop.size.mat, animals.xy=animals.xy, hunt.catch.mat=hunt.catch.mat, params=params, pop.size.list=pop.size.list, trap.catch.list=trap.catch.list, bait.catch.list=bait.catch.list))#, pop.zone.list=pop.zone.list))#, animals.done.xy=animals.xy))    
  }
  )
  
  
  
  output$mymap<-renderLeaflet({
    shp<-mydata.shp()$shp
    # shp<-mydata.zone()$shp.2
    
    shp.proj<-spTransform(shp,CRS("+proj=longlat +datum=WGS84"))
    # ha<-mydata()$ha
    
    m<-leaflet() %>%
      addTiles(group="Default")%>%
      addProviderTiles("Esri.WorldTopoMap", group = "Topo")%>%
      addProviderTiles("Esri.WorldImagery", group = "Aerial")%>%
      addPolygons(data=shp.proj, weight=2, fillColor="grey40", color="black", fillOpacity=0.2)
    m
    
  })
  
  
  
  observe({
    traps<-mydata.map()$traps
    animals.xy<-mydata.map()$animals.xy.ini
    proj4string<-mydata.shp()$p4s
    
    traps.proj<-proj4::project(traps, proj=proj4string, inverse=T) 
    tmp<-(proj4::project(animals.xy[,1:2], proj=proj4string, inverse=T))
    animals.xy$Lat<-tmp$y
    animals.xy$Lon<-tmp$x
    map<-leafletProxy("mymap")
    map%>%clearMarkers()
    
    map%>%addCircleMarkers(lng=traps.proj$x,lat=traps.proj$y, radius=3, color="black", weight=1, fill=TRUE, fillColor="red", fillOpacity=1, stroke=TRUE, group="Traps")
    map%>%addCircleMarkers(lng=animals.xy$Lon,lat=animals.xy$Lat, radius=5, color="black", weight=1, fill=TRUE, fillColor=cols.vec[2], fillOpacity=1, stroke=TRUE, group="Animals")
    
    map%>%addLayersControl(
      baseGroups = c("Default","Topo","Aerial"),
      overlayGroups=c("Traps","Animals"),
      options = layersControlOptions(collapsed = FALSE)
    )
    
  })
  
  
  
  
  output$plot1<-renderPlot({
    
    res<-datab()$params
    par(mar=c(4,4,1,1), mgp=c(2.5,1,0), tcl=-0.25)
    plot(res$MeanPopSize, res$TotalCost, xlab="Final population", ylab="Total cost")
    
    
  })
  
  
  
  output$plot2<-renderPlot({
    trap.catch.list<-datab()$trap.catch.list
    trap.catch.mat<-trap.catch.list[[as.numeric(input$result_scenario)]]
    bait.catch.list<-datab()$bait.catch.list
    bait.catch.mat<-bait.catch.list[[as.numeric(input$result_scenario)]]
    
    nights.vec<-1:input$n.nights
    ymax.t<-max(cumsum(colMeans(trap.catch.mat)))
    ymax.b<-max(cumsum(colMeans(bait.catch.mat)))
    ymax<-max(ymax.b, ymax.t)
    
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,ymax), type='n', xlab="Nights", ylab="Cumulative kills", las=1)
    for(i in 1:input$n.its){
      lines(nights.vec,cumsum(trap.catch.mat[i,]), col=cols.vec[1])
      lines(nights.vec,cumsum(bait.catch.mat[i,]), col=cols.vec[3])
    }
    lines(nights.vec,cumsum(colMeans(trap.catch.mat)), col=cols.vec[2])
    points(nights.vec,cumsum(colMeans(trap.catch.mat)), bg=cols.vec[2], pch=21)
    
    lines(nights.vec,cumsum(colMeans(bait.catch.mat)), col=cols.vec[4])
    points(nights.vec,cumsum(colMeans(bait.catch.mat)), bg=cols.vec[4], pch=22)
    
    legend("topleft", legend=c("Traps","Bait station"), pch=c(21,22), pt.bg=cols.vec[c(2,4)], bty="n")
    mtext("Cumulative captures",3, cex=1.5, line=1)
    # mtext("Trapped per night",3, cex=1.5, line=1)
  })
  
  
  # output$plot2.hunt<-renderPlot({
  #   
  #   hunt.catch.list<-datab()$hunt.catch.list
  #   trap.catch.mat<-hunt.catch.list[[as.numeric(input$result_scenario)]]
  #   
  #   nights.vec<-1:input$n.nights
  #   ymax<-max(cumsum(colMeans(trap.catch.mat)))
  #   
  #   par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
  #   plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,ymax), type='n', xlab="Nights", ylab="Cumulative captures", las=1)
  #   for(i in 1:input$n.its){
  #     lines(nights.vec,cumsum(trap.catch.mat[i,]), col="grey")
  #   }
  #   lines(nights.vec,cumsum(colMeans(trap.catch.mat)), col="black")
  #   points(nights.vec,cumsum(colMeans(trap.catch.mat)), bg="black", pch=21)
  #   mtext("Cumulative captures: Hunting",3, cex=1.5, line=1)
  #   # mtext("Trapped per night",3, cex=1.5, line=1)
  # })
  
  
  output$plot4<-renderPlot({
    
    pop.size.list<-datab()$pop.size.list
    pop.size.mat<-pop.size.list[[as.numeric(input$result_scenario)]]
    # pop.size.mat<-datab()$pop.size.mat
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,max(pop.size.mat)), type='n', xlab="Nights", ylab="Population Size", las=1)
    grid()
    for(i in 1:input$n.its){
      lines(0:input$n.nights,pop.size.mat[i,], col="grey")
    }
    lines(0:input$n.nights,colMeans(pop.size.mat), col="black")
    points(0:input$n.nights,colMeans(pop.size.mat), bg="black", pch=21)
    mtext("Population size",3, cex=1.5, line=1)
    
  })
  
  observe({
    # updateSelectInput(session=session, inputId="result_scenario",choices=mydata.scen()$params$Scenario)
    updateSelectInput(session=session, inputId="result_scenario",choices=rownames(scenParam()))
    
  })
  output$plot7<-renderPlot({
    par(mar=c(4,3,1,1),tcl=-.2, mgp=c(2,0.5,0))
    locshp<-get.loc.shape(input$sigma.mean, input$sigma.sd)
    hist(rlnorm(10000 ,meanlog=locshp$location, sdlog=locshp$shape), xlab="", ylab="", main="", col="grey")
    mtext("Sigma (m)",1, cex=1, line=2)
    
  })
  
  
  output$plot6<-renderPlot({
    par(mar=c(4,3,1,1),tcl=-.2, mgp=c(2,0.5,0))
    alpbet<-get.alpha.beta(input$g0.mean.a, input$g0.sd.a)
    hist(rbeta(10000, alpbet$alpha, alpbet$beta), xlim=c(0,1), main="", col="grey", xlab="", ylab="")
    mtext("g0",1, cex=1, line=2)
    
    
  })
  
  
  
  output$plot.hunt<-renderPlot({
    eff<-seq(from=50, to=2000, by=50)
    theta.hat.a<-1-exp(-((input$hunt.rho.a*log(eff))^input$hunt.k.a))
    theta.hat.b<-1-exp(-((input$hunt.rho.b*log(eff))^input$hunt.k.b))
    plot(eff, theta.hat.a, ylim=c(0,1), ylab="Pr.Kill",las=1, xlab="Distance (m)", type='l', col="blue", main="P.kill as a function of effort/day")
    lines(eff, theta.hat.b, col="red")
    legend("topleft", legend=c("Method 1","Method 2"), col=c(4,2), lty=1)
  })
  
  
  output$text5<-renderText({
    traps<-mydata.map()$traps
    ntraps<-dim(traps)[1]
    ha.area<-mydata.shp()$ha
    
    paste(ntraps, " Traps (One per ",round(ha.area/ntraps,1)," ha)", sep="")
  })
  
  
  output$text6<-renderText({
    hr.radius<-input$sigma.mean*2.45
    # ha.area<-input$area.ha
    return(paste0("Home range size: ",round(((pi*hr.radius^2)/10000),2), " ha"))
  })
  
  
  output$text9<-renderText({
    
    sigma.mean<-input$sigma.mean
    grid.size<-sigma.mean*4
    
    return(paste0("Suggested grid size: ", round(grid.size,0)," m"))
  })
  
  
  output$text8<-renderText({
    
    g0.mean<-input$g0.mean.a
    max.sd<-sqrt((g0.mean^2 - g0.mean)*-1)
    return(paste0("Max SD: ", round(max.sd,2)))
  })
  
  
  output$text10<-renderText({
    ha<-mydata.shp()$ha
    return(paste0("Total area (ha): ", round(ha,0)))
  })

  
  
  output$text_trap_a_cost<-renderText({
    shp<-mydata.shp()$shp    
    traps.a<-make.trap.locs(input$traps.x.space.a, input$traps.y.space.a, 100, shp)
    n.traps.a<-dim(traps.a)[1]
    # check.vec.a<-seq(from=input$trap.start.a, to=(input$trap.start.a+input$trap.nights.a), by=input$n.check.a)
    checks<-ceiling(input$trap.nights.a/input$n.check.a)+1  #The number of checks - copes with check intervals that dont fit nealy into the duration
    cost<-trap.cost.func(a=checks, b=n.traps.a, c=input$traps.per.day.a, d=input$day.rate.a, e=input$cost.per.trap.a)
    
    return(paste0(n.traps.a," Traps\nCost = $", cost ))
  })

  output$text_trap_b_cost<-renderText({
    shp<-mydata.shp()$shp    
    traps.b<-make.trap.locs(input$traps.x.space.b, input$traps.y.space.b, 100, shp)
    n.traps.b<-dim(traps.b)[1]
    # check.vec.b<-seq(from=input$trap.start.b, to=(input$trap.start.b+input$trap.nights.b), by=input$n.check.b)
    checks<-ceiling(input$trap.nights.b/input$n.check.b)+1
    cost<-trap.cost.func(a=checks, b=n.traps.b, c=input$traps.per.day.b, d=input$day.rate.b, e=input$cost.per.trap.b)
    
    return(paste0(n.traps.b," Traps\nCost = $", cost ))
  })
    
  
  
  output$text_bait_a_cost<-renderText({
    shp<-mydata.shp()$shp    
    bait.a<-make.trap.locs(input$bait.x.space.a, input$bait.y.space.a, 100, shp)
    n.bait.a<-dim(bait.a)[1]
    # bait.check.vec.a<-seq(from=input$bait.start.a, to=(input$bait.start.a+input$bait.nights.a), by=input$bait.check.a)
    checks<-ceiling(input$bait.nights.a/input$bait.check.a)+1
    cost<-trap.cost.func(a=checks, b=n.bait.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a)
    
    return(paste0(n.bait.a," Bait Stations\nCost = $", cost ))
  })
  
output$text_density<-renderText({
  ha<-mydata.shp()$ha
  numb<-input$numb.poss
  
  
  return(paste0(round(numb/ha,2)," per ha"))
})
  
  
  # output$text11<-renderText({
  #   ha<-mydata()$ha
  #   hunt.cell.size<-input$hunt.cell.size
  #   return(paste0("Number of cells ", round(ha/hunt.cell.size,0)))
  # })
  
  
  # output$scenario_dropdown<-renderUI({
  #   
  #   params<-mydata.c()$params
  #   
  #   selectInput("result_scenario","Choose Scenario to plot",params$Scenario)
  #   
  # })
  output$texthunta<-renderText({
    p.hunt.day<-mydata.hunt.prob()$p.hunt.day
    return(paste0("Daily kill prob ", round(p.hunt.day[1],2)))
  })
  output$texthuntb<-renderText({
    p.hunt.day<-mydata.hunt.prob()$p.hunt.day
    return(paste0("Daily kill prob ", round(p.hunt.day[2],2)))
  })
  # 
  
  
  output$results.table<-renderDataTable({
    
    datab()$params
    
  })
  
  # output$scenarios.table<-renderTable({
  #   
  #   mydata.scen()$params
  #   
  # })
  
  output$tableDT <- DT::renderDataTable(
    scenParam()
  )
  
  
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         The user interface
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui<-fluidPage(theme=shinytheme("flatly"),
              title="Pest Control",
              tags$head(
                tags$style(HTML("
                                .shiny-output-error-validation {
                                color: blue;
                                }
                                "))),
              fluidRow(
                column(width=4,  
                       h1("Pest Control DSS"),
                       checkboxInput("showinfo","Background information and instructions"),

                       
                       # h3("Trapping Simulation Tool")
                ),
                column(8,
                       img(src="manaaki_logo.png", height = 90, align="right", hspace=20,vspace=10),
                       # img(src="ari_logo.jpg", height = 90, align="right", hspace=20,vspace=10),
                       img(src="bhnsc.png", height = 90, align="right", hspace=20,vspace=10)

                       # img(src="IC_logo.png", height = 90, align="right", hspace=20,vspace=10)
                )
              ),
              
              
              fluidRow(
                column(width=5, 
                       conditionalPanel(
                         condition = "input.showinfo == 1",
                         "This simulation app is intended to provide guidance as to the approximate amount of trapping/baiting that you may require to achieve a various levels of pest reduction.",
                         # "You can specify an area of a specified size, or upload a shapefile of the area of interest.",
                         # "When you have set up the area, trap layout and pest parameters, click the", strong("Run Trap Sim"), "button at the bottom of the page to run the trapping simulation.",
                         p(),"It is based on an earlier tool called TrapSim which was focused on simulating trapping only.",
                         "For a background report on TrapSim, ", 
                         tags$a(href="https://www.pfhb.nz/assets/Document-Library/Gormley-and-Warburton-2017-TrapSim-a-decision-support-tool.pdf", "Click here.", target="_blank"), p(),
                         
                         h3("Instructions"),                        
                         "1. Set the area and pest parameters" ,br(),
                         "2. Control Methods - set the various scenarios by selecting the control methods and the corresponding parameters. Click 'Add Scenario' to build up a list of scenarios", br(),
                         "3. Run Sceanrios - check the scenarios, enter the number of iterations for each and the simulation length. Click 'Run TrapSim'", br(),
                         "4. Results - explore the results  - graphs and tables of animals remaining, costs etc", br(),
                         
                         
                         # "For input parameters with red labels, enter values separated by a slash, e.g. 1000/500 ",
                         p(),
                         strong("Hover cursor over each of the input boxes for pop-up help.")

                       )
                )
              ),
              
              
              
              fixedRow(
                column(width=12,
                       tabsetPanel(id="inTabset",
                                   tabPanel("1. Area and Pest Parameters",
                                            column(width=4,
                                                   h4("Area"),
                                                   wellPanel(
                                                     #   div(style="display:inline-block",
                                                     #       tags$div(title="'Specify-Size': specify a total size for a random square; 'Upload Shapefile': upload a shapefile of your study area",
                                                     #                
                                                     radioButtons(inputId="area_type", label="Chose area",choices=c("Mahia Peninsula"="RC","Upload Shapefile"="Map"), selected="RC"),
                                                     #   
                                                     #   conditionalPanel(
                                                     #     condition="input.area_type=='Area'",
                                                     #     numericInput(inputId = "area.ha", label="Area (ha)", value=10000, width="120px")
                                                     #   ),
                                                     #   
                                                     conditionalPanel(
                                                       condition="input.area_type=='Map'",
                                                       
                                                       tags$div(title="Be sure to select all the components (.shp, .dbf etc)",
                                                                fileInput(inputId = "shp.file", label="Chose the VCZ Shapefile", accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj"), multiple=TRUE, width="200px")
                                                       )
                                                     ),
                                                     verbatimTextOutput("text10")
                                                   ),
                                                   h4("Pest parameters"),
                                                   wellPanel(
                                                     div(style="display:inline-block;vertical-align:top",
                                                         numericInput(inputId = "numb.poss", label="Number (total)", value=100, width="180px")
                                                     ),
                                                     div(style="display:inline-block;vertical-align:top",
                                                         # fileInput(inputId = "ras.1", label="Ascii file of relative abundance", accept=c('.asc'), multiple=FALSE, width="250px")
                                                         # radioButtons(inputId = "ras.1", label="Relative abundance", choices=c("Random","Habitat Specific"), selected="Random",width="250px")
                                                         radioButtons(inputId = "ras.1", label="Relative abundance", choices=c("Random"), selected="Random",width="250px")
                                                     ),
                                                     div(style="display:inline-block;vertical-align:bottom",
                                                     verbatimTextOutput("text_density")),
                                                     tags$div(title="Sigma x 2.45 is the radius of a circle  where an indivudual spends 95% of its time.",
                                                              h5(strong("Home range (sigma)"))),
                                                     div(style="display:inline-block",
                                                         numericInput(inputId = "sigma.mean", label="Mean", value=100, width="120px")),
                                                     div(style="display:inline-block",
                                                         numericInput(inputId="sigma.sd", label='StdDev', value=5, width="120px")),
                                                     div(style="display:inline-block;vertical-align:bottom",
                                                         checkboxInput(inputId = "show_sigma",label="Show Sigma dist.", value=FALSE)
                                                     ),
                                                     conditionalPanel(
                                                       condition="input.show_sigma==1", 
                                                       plotOutput(outputId = "plot7", width = "250px", height="200px")),
                                                     verbatimTextOutput("text6"),
                                                     
                                                     
                                                     
                                                     h5(strong("Reproductive parameters")),
                                                     div(style="display:inline-block",
                                                         tags$div(title="This is the maximum rate of increase.",
                                                                  numericInput(inputId = "rmax.poss", label="Rmax", value=0.4, width="120px")
                                                         )),
                                                     div(style="display:inline-block",
                                                         tags$div(title="This is the start day of the reproductive period (i.e. start of the birth pulse).",
                                                                  numericInput(inputId = "rep.start", label="Start day", value=110, width="120px")
                                                         )),
                                                     div(style="display:inline-block",
                                                         tags$div(title="This is the number of days the reproductive period lasts (i.e. duration of the birth pulse).",
                                                                  numericInput(inputId = "rep.nights", label="Length (days)", value=60, width="120px")
                                                         ))
                                                   )
                                            ),
                                            column(width=8,
                                                   "This map is provided to visualise the area. The parameter values specified here are ", strong("not"), " included in the simulations",
                                                   p(),
                                                   # div(style="display:inline-block;vertical-align:bottom",
                                                   #     numericInput(inputId = "numb.poss.i", label="Population size", value=10, width="180px")
                                                   # ),
                                                   div(style="display:inline-block;vertical-align:bottom",
                                                       tags$div(title="The trap spacing in the east-west direction",
                                                                numericInput(inputId = "traps.x.space.i", label="Trap space E-W (m)", value="1000",width="135px"))),
                                                   div(style="display:inline-block;vertical-align:bottom",
                                                       tags$div(title="The trap spacing in the north-south direction",
                                                                numericInput(inputId = "traps.y.space.i", label="Trap space N-S (m)", value="1000",width="135px"))),
                                                   div(style="display:inline-block;vertical-align:bottom",
                                                       tags$div(title="The buffer from the edge",
                                                                numericInput(inputId = "traps.buff.i", label="Buffer (m)", value="100", min=0, max=1000, width="110px"))),
                                                   div(style="width:300px;display:inline-block;vertical-align:bottom",
                                                       verbatimTextOutput("text5")),
                                                   # p(),
                                                   leafletOutput(outputId = "mymap", height = "1000px", width="1200px")
                                                   
                                            )
                                   ),
                                   tabPanel("2. Control Methods",

                                            # h3("Choose 'Trapping' and/or 'Hunting' and then specify inputs for up to two methods of each. Click 'Add Scenario' to add them to the 'Run Scenarios' tab "),
                                            fluidRow(
                                              column(width=2,
                                                     
                                                     wellPanel(
                                                       h5("Choose control method(s)"),
                                                       checkboxInput(inputId="trap_methods", label="Trapping", value=TRUE),
                                                       checkboxInput(inputId="bait_methods", label="Bait stations", value=TRUE)
                                                     )),
                                              column(width=1,
                                                     wellPanel(
                                                       actionButton("update", "Add Scenario"))
                                              )
                                            ),
                                            
                                            
                                            conditionalPanel(
                                              condition="input.trap_methods==1",
                                              
                                              h4("Trapping"),
                                              wellPanel(

                                                div(style="display:inline-block;vertical-align:middle",h5("Trapping Method 1")),p(),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The trap spacing in the east-west direction",
                                                             numericInput(inputId = "traps.x.space.a", label="Spacing E-W (m)", value="1000",width="135px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The trap spacing in the north-south direction",
                                                             numericInput(inputId = "traps.y.space.a", label="Spacing N-S (m)", value="1000",width="135px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Nightly probability of by-catch, false triggers etc  ",
                                                             numericInput(inputId = "p.bycatch.a", label="Daily bycatch", value=0, width="120px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Maximum catch per trap  ",
                                                             numericInput(inputId = "max.catch.a", label="Max catch", value=1, width="135px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    numericInput(inputId = "g0.mean.a", label="Trap Probability (g0) Mean", value=0.1, width="120px")),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    numericInput(inputId="g0.sd.a", label='StdDev', value=.01, width="120px")),
                                                # div(style="display:inline-block;vertical-align:bottom",
                                                #     verbatimTextOutput("text8")),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    numericInput(inputId = "g0.zero.a", label="Proportion untrappable", value=0.05, width="120px")),
                                                #Timings
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="Start night of trapping.",
                                                             numericInput(inputId = "trap.start.a", label="Start night", value="5", width="120px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="Number of nights traps are set for.",
                                                             numericInput(inputId = "trap.nights.a", label="Duration (nights)", value="20", width="120px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                             numericInput(inputId = "n.check.a", label="Check interval", value="1", width="120px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Number of traps checked per day",
                                                             numericInput(inputId = "traps.per.day.a", label="Traps checked per day", value=40, width="120px")
                                                    )),
                                                
                                                
                                                #Costs
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Labour cost - day rate ($)",
                                                             numericInput(inputId = "day.rate.a", label="Day rate ($)", value=400, width="135px")
                                                    )),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Fixed cost ($) per trap",
                                                             numericInput(inputId = "cost.per.trap.a", label="Fixed cost per trap ($)", value=20, width="120px")
                                                    )),
                                                
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    verbatimTextOutput("text_trap_a_cost")),
                                                
                                                
                                                
                                                # div(style="display:inline-block;vertical-align:bottom",
                                                #     checkboxInput(inputId = "show_g0",label="Show g0 dist.", value=FALSE)),
                                                # conditionalPanel(
                                                #   condition="input.show_g0==1", 
                                                #   plotOutput(outputId = "plot6", width = "250px", height="200px")),
                                                
                                                p(),
                                                
                                                
                                                tags$style(type="text/css", "#redtitle {color: black}"),
                                                # ),
                                                
                                                p(),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    checkboxInput(inputId = "show_trap_b",label="Second Trap Method", value=FALSE)
                                                ),
                                                conditionalPanel(
                                                  condition="input.show_trap_b==1",
                                                  # wellPanel(
                                                  div(style="display:inline-block;vertical-align:bottom",h5("Trapping Method 2")),p(),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(id="redtitle",title="The trap spacing in the east-west direction",
                                                               
                                                               numericInput(inputId = "traps.x.space.b", label="Spacing E-W (m)", value="200",width="135px"))),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(id="redtitle",title="The trap spacing in the north-south direction",
                                                               numericInput(inputId = "traps.y.space.b", label="Spacing N-S (m)", value="200",width="135px"))),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(title="Nightly probability of by-catch, false triggers etc  ",
                                                               numericInput(inputId = "p.bycatch.b", label="Daily bycatch", value=0, width="120px"))),
                                                  
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(title="Maximum catch per trap  ",
                                                               numericInput(inputId = "max.catch.b", label="Max catch", value=1, width="135px"))),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      numericInput(inputId = "g0.mean.b", label="Trap Probability (g0) Mean", value=0.2, width="120px")),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      numericInput(inputId="g0.sd.b", label='StdDev', value=.01, width="120px")),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      numericInput(inputId = "g0.zero.b", label="Proportion untrappable", value=0.05, width="120px")),

                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(id="redtitle",title="Start night of trapping.",
                                                               numericInput(inputId = "trap.start.b", label="Start night", value="50", width="120px"))),
                                                  
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(id="redtitle",title="Number of nights traps are set for.",
                                                               numericInput(inputId = "trap.nights.b", label="Duration (nights)", value="20", width="120px"))),
                                                  
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(id="redtitle",title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                               numericInput(inputId = "n.check.b", label="Check interval", value="10", width="120px"))),
                                                  
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(title="Number of traps checked per day",
                                                               numericInput(inputId = "traps.per.day.b", label="Traps checked per day", value=40, width="120px")
                                                      )),
                                                  
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(title="Labour cost - day rate ($)",
                                                               numericInput(inputId = "day.rate.b", label="Day rate ($)", value=400, width="135px")
                                                      )),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      tags$div(title="Fixed cost ($) per trap",
                                                               numericInput(inputId = "cost.per.trap.b", label="Fixed cost per trap ($)", value=20, width="120px")
                                                      )),
                                                  div(style="display:inline-block;vertical-align:bottom",
                                                      verbatimTextOutput("text_trap_b_cost"))
                                                  
                                                  
                                                ) #End of Well Panel
                                              )
                                            ),
                                            
                                            
                                            conditionalPanel(
                                            condition="input.bait_methods==1",
                                            
                                            h4("Bait Stations"),
                                            wellPanel(
                                              
                                              # div(style="display:inline-block;vertical-align:middle",h5("Baiting Method 1")),p(),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="The bait station spacing in the east-west direction",
                                                           numericInput(inputId = "bait.x.space.a", label="Spacing E-W (m)", value="500",width="135px"))),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="The bait station spacing in the north-south direction",
                                                           numericInput(inputId = "bait.y.space.a", label="Spacing N-S (m)", value="500",width="135px"))),
                                              # div(style="display:inline-block;vertical-align:bottom",
                                              #     tags$div(title="The buffer from the edge",
                                              #              numericInput(inputId = "bait.buff.a", label="Edge buffer (m)", value="100", min=0, max=1000, width="120px"))),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Nightly probability of station failure  ",
                                                           numericInput(inputId = "p.failure.a", label="Daily rate of failure", value=0, width="120px"))),
                                              
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  numericInput(inputId = "bait.g0.mean.a", label="Bait Sation Probability (g0) Mean", value=0.1, width="120px")),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  numericInput(inputId="bait.g0.sd.a", label='StdDev', value=.01, width="120px")),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  numericInput(inputId = "bait.g.zero.a", label="Proportion untrappable", value=0.05, width="120px")),

                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="Start night of baiting.",
                                                           numericInput(inputId = "bait.start.a", label="Start night", value="5", width="120px"))),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="Number of nights baits are set for.",
                                                           numericInput(inputId = "bait.nights.a", label="Duration (nights)", value="50", width="120px"))),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="The checking interval of the bait stations. For stations that are not cleared, set equal to Nights ",
                                                           numericInput(inputId = "bait.check.a", label="Check interval", value="10", width="120px"))),
                                              
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Number of stations checked per day",
                                                           numericInput(inputId = "bait.per.day.a", label="Bait stations checked per day", value=40, width="120px")
                                                  )),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Labour cost - day rate ($)",
                                                           numericInput(inputId = "bait.day.rate.a", label="Day rate ($)", value=400, width="135px")
                                                  )),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Fixed cost ($) per bait station",
                                                           numericInput(inputId = "cost.per.bait.a", label="Fixed cost per bait station ($)", value=60, width="120px")
                                                  )),
                                              
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  verbatimTextOutput("text_bait_a_cost")),
                                              # text_bait_a_cost
                                              
                                              
                                              tags$style(type="text/css", "#redtitle {color: black}")
                                              # ),
                                            )
                                            ),#End of Well Panel
                                            
                                            
                                            
                                            # ), #End of TabPanel
                                            # tabPanel("Hunting",
                                            conditionalPanel(
                                              condition="input.hunt_methods==1",
                                              
                                              h4("Hunting"),
                                              wellPanel(
                                                
                                                div(style="display:inline-block;vertical-align:middle",h5("Hunting Method 1")),p(),
                                                div(style="display:inline-block;vertical-align:middle",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "hunt.start.a", label="Start day", value="30", width="120px"))),
                                                div(style="display:inline-block;vertical-align:middle",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "hunt.days.a", label="Days hunted", value="10", width="120px"))),
                                                div(style="display:inline-block;vertical-align:middle",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "effort.a", label="Effort/day (m)", value="250", width="120px"))),
                                                div(style="display:inline-block;vertical-align:middle",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "day.rate.hunt.a", label="Day rate ($)", value=500, width="120px"))),
                                                div(style="display:inline-block;vertical-align:middle",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "hunt.k.a", label="K", value="3.2", width="100px"))),
                                                div(style="display:inline-block;vertical-align:middle",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "hunt.rho.a", label="rho", value="0.1", width="100px"))),
                                                div(style="display:inline-block;vertical-align:middle",
                                                    verbatimTextOutput("texthunta")),
                                                
                                                p(),
                                                
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    checkboxInput(inputId = "show_hunt_b",label="Second Hunting Method ?", value=FALSE)
                                                ),
                                                conditionalPanel(
                                                  condition="input.show_hunt_b==1",
                                                  
                                                  div(style="display:inline-block;vertical-align:middle",h5("Hunting Method 2")),p(),
                                                  div(style="display:inline-block;vertical-align:middle",
                                                      tags$div(title="help text ",
                                                               numericInput(inputId = "hunt.start.b", label="Start day", value="50", width="120px"))),
                                                  div(style="display:inline-block;vertical-align:middle",
                                                      tags$div(title="help text ",
                                                               numericInput(inputId = "hunt.days.b", label="Days hunted", value="20", width="120px"))),
                                                  div(style="display:inline-block;vertical-align:middle",
                                                      tags$div(title="help text ",
                                                               numericInput(inputId = "effort.b", label="Effort/day (m)", value="150", width="120px"))),
                                                  div(style="display:inline-block;vertical-align:middle",
                                                      tags$div(title="help text ",
                                                               numericInput(inputId = "day.rate.hunt.b", label="Day rate ($)", value=450, width="120px"))),
                                                  div(style="display:inline-block;vertical-align:middle",
                                                      tags$div(title="help text ",
                                                               numericInput(inputId = "hunt.k.b", label="K", value="2.2", width="100px"))),
                                                  div(style="display:inline-block;vertical-align:middle",
                                                      tags$div(title="help text ",
                                                               numericInput(inputId = "hunt.rho.b", label="rho", value="0.12", width="100px"))),
                                                  div(style="display:inline-block;vertical-align:middle",
                                                      verbatimTextOutput("texthuntb"))
                                                  
                                                )
                                              )
                                            )
                                            
                                   ),

                                   tabPanel("3. Run Scenarios",
                                            # tableOutput('scenarios.table'),
                                            h3("Use this tab to check scenarios, and delete them if needed (click on the rows to delete), then click 'Run TrapSim' to run the simulations."),
                                            fluidRow(
                                              column(width=2,
                                                     wellPanel(
                                                       
                                                       radioButtons(inputId="sim_type", label="Simulation Type",choices=c("Individual traps"="individ"), selected="individ", width="200px"),
                                                       div(style="display:inline-block;vertical-align:bottom",
                                                           numericInput(inputId = "n.nights",label="Simulation length", value=100, width="120px")),
                                                       div(style="display:inline-block;vertical-align:bottom",
                                                           numericInput(inputId = "n.its",label="Iterations", value=5, width="120px"))
                                                       
                                                     )),
                                              column(width=1,
                                                     wellPanel(
                                                       div(style="display:inline-block;vertical-align:bottom",
                                                           actionButton("act.btn.trapsim","Run TrapSim"))
                                                       
                                                     )
                                              )
                                            ),
                                            
                                            # fluidRow(
                                            actionButton("deleteRows", strong("Delete Selected Rows")),
                                            actionButton("deleteAllRows", strong("Delete All Rows")),
                                            p(),
                                            # ),
                                            
                                            DT::dataTableOutput("tableDT")
                                   ),
                                   

                                   
                                   # ),
                                   tabPanel("4. Results",
                                            fluidRow(
                                              column(width=2,
                                                     selectInput("result_scenario","Choose Scenario to plot",choices=NULL),
                                                     plotOutput(outputId = "plot1", width = "350px", height="300px")
                                                     
                                              ),
                                              column(width=3,
                                                     plotOutput(outputId = "plot2", width = "550px", height="350px")  #Cumulative Captures
                                              ),
                                              # column(width=3,
                                              #        plotOutput(outputId = "plot2.hunt", width = "550px", height="350px")
                                              # ),
                                              column(width=3,
                                                     plotOutput(outputId = "plot4", width = "550px", height="350px") #Population Size
                                              )
                                              
                                              # )
                                            ),
                                            fluidRow(
                                              # column(width=6,
                                              dataTableOutput('results.table')
                                              # uiOutput("scenario_dropdown")
                                              
                                            )
                                            
                                   )
                                   
                       ))
              ), #End of Row
              
              
              h6("v0.9: December 2021"),
              h6("email: gormleya@landcareresearch.co.nz")
              # h6("TrapSim was originally developed using funding from Centre for Invasive Species Solutions (CISS)"),
              # img(src="ciss_logo.jpg", height = 90, align="right", hspace=20,vspace=10)
              # h6("Developed using funding from Centre for Invasive Species Solutions (CISS), MBIE (New Zealand), and Island Conservation")
)


shinyApp(ui=ui, server=server)

#1.2 includes ability to read in a csv file of traps...
#1.3 - multi-capture traps...
#1.4 - costs and rearrange interface...
#1.5 - run scenarios across some parameter combinations
#1.5.2 - outputs by zones as well. 
#1.6 - removed zones - too fricking complicated, and got two methods of hunting and trapping in there. But no scenarios....
#1.7.1 - scenarios are back! - with add and delete buttons.

#Need to include bait costs...?
#Group costs by labour and fixed costs...?





