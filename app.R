#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            TrapSim Tool - Multi Method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Code developed by Andrew Gormley, Manaaki Whenua - Landcare Research
#Developed for the NSC Eco-economics project and CISS

#Includes Audrey's grid based approach as well - Audrey says it works with cells of 500m and even 200m...(?)
#28/8 - some testing would suggest that for nightly checking, setting the grid square dimension at 4*sigma is best when compared to IBM

#Doesn't currently have the grid-based part in there...


#Testing to see that things push to the correct place...
#I think it works....

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
library("shinyvalidate")
library("dplyr")
library("shinycssloaders")

# setwd("C:\\Users\\gormleya\\OneDrive - MWLR\\Documents\\CAEM\\IslandConservation\\TrapSimFeasibility\\Shiny")

# def.shp<-"Robinson_Coati"  #The default shape.
def.shp<-"WhakatipuMahia"
# def.shp1<-"AshleyForestVCZ"
# shp.zones<-"Robinson_Coati_Workzones"
# ras.2<-raster("habitat_specific.asc")

# cols.tst<-brewer.pal(6,"Blues")
cols.vec<-brewer.pal(8,"Paired")
cols.eff<-brewer.pal(8,"Reds")
proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
#4326 is the EPSG for WGS84...

#~~~~~~~~~~~~~~~ A bunch of functions ~~~~~~~~~~~~~~~~~~~~~~~~
#Define the resampling function
resamp <- function(x,...){if(length(x)==1) x else sample(x,...)} 

#~~~~~~~~~~~~~~~~COST FUNCTIONS~~~~~~~~~~~~~
#Calculate the cost of trapping
trap.cost.func<-function(a,b,c,d,e){
  #e.g. trap.cost.func(a=number of checks, b=n.traps, c=input$traps.per.day, d=input$day.rate,e=cost.per.trap)
  fixed.cost<-b*e  #n.traps x cost.per.trap
  labor.cost<-ceiling(2*b/c)/2*a*d #- scales up to a half day, mostly works 
  trap.cost<-as.integer(labor.cost+fixed.cost)
  # trap.cost<-as.integer((((b*checks)/c)*d)+(b*e))
  return(trap.cost)
}
# bait.cost.a

#Cost of bait-stations
bait.cost.func<-function(a,b,c,d,e,f){
  #e.g. trap.cost.func(a=number of checks, b=n.traps, c=input$traps.per.day, d=input$day.rate, e=cost.per.trap, fbait cost per station per day)
  fixed.cost<-b*e  #n.traps x cost.per.trap
  bait.cost<-b*a*f
  labor.cost<-ceiling(2*b/c)/2*a*d #- scales up to a half day, mostly works 
  trap.cost<-as.integer(labor.cost+fixed.cost+bait.cost)
  # trap.cost<-as.integer((((b*checks)/c)*d)+(b*e))
  return(trap.cost)
}


#Calculate the cost of hunting
hunt.cost.func<-function(a,b){
  day.rate<-a
  days.hunt<-b
  hunt.cost<-day.rate*days.hunt
  return(hunt.cost)
}

#Calculate the cost of poisoning
pois.cost.func<-function(a,b,c){
  pois.per.ha<-a  #Cost per ha
  hectares<-b   #Total hectares
  pois.prop<-c/100  #Proportion of area treated.
  pois.cost<-round(pois.per.ha*hectares*pois.prop,0)
  return(pois.cost)
}


#For a given mean/sd, calculate the alpha & beta parameters
get.alpha.beta<-function(g0.mean, g0.sd){
  g0.var<-g0.sd^2
  paren<-    (g0.var + g0.mean^2 - g0.mean)
  alpha<- -g0.mean*paren/g0.var
  beta<- paren*(g0.mean-1)/g0.var
  return(list(alpha=alpha, beta=beta))
}

#For a given mean/sd, calculate the location and shape parameters for a log-normal distribution
get.loc.shape<-function(sigma.mean, sigma.sd){
  m<-sigma.mean
  s<-sigma.sd
  location <- log(m^2 / sqrt(s^2 + m^2))
  shape <- sqrt(log(1 + (s^2 / m^2)))
  return(list(location=location, shape=shape))
}

#Function to make animal locations from a raster of relative abundance/habitat.
# Takes the raster, the number of possums and the shapefile
get.pest.locs<-function(ras, n.poss, shp){
  # 
  # gridpolygon<-rasterToPolygons(ras)
  # grids<-length(gridpolygon@data[,1])
  # pest.per.grid<-round(gridpolygon@data[,1]*n.poss/grids,0)   #We want to sample more than we need because we are going to remove some that are outside the shapefile
  # coords<-vector("list", grids) #Set up to save the coordinates - cant work out how to get spsample to sample across all grids
  # for (i in 1: grids){
  #   if(pest.per.grid[i]>0){
  #     coords[[i]]<-spsample(gridpolygon@polygons[[i]], n =pest.per.grid[i] , type = "random")@coords
  #   }
  # }
  # coords<-as.data.frame(do.call("rbind", coords)) #Put them altogether from the list
  # 
  buff<-n.poss*5
  # Get the celll resolution
  res.x<-res(ras)[1]
  res.y<-res(ras)[1]
  
  df <- as.data.frame(ras, xy = T, na.rm = T)
  # sampled_cells <- sample(1:nrow(df), size = n.poss+buff, prob = df$habitat_specific, replace = T)
  sampled_cells <- sample(1:nrow(df), size = n.poss+buff, prob = df[,3], replace = T)
  x.val<-df$x[(sampled_cells)]+res.x*(runif(n.poss+buff)-0.5)
  y.val<-df$y[(sampled_cells)]+res.y*(runif(n.poss+buff)-0.5)
  coords<-data.frame('x'=x.val, 'y'=y.val)
  coords<-coords[inside.owin(coords[,1], coords[,2], shp),]  #Remove the ones from outside the shapefile...
  # coords<-coords[sample(1:dim(coords)[1], size=n.poss, replace=FALSE),]  #Then sample to get the desired actual sample size
  
  
  if(dim(coords)[1]>0){
    coords<-coords[sample(1:dim(coords)[1], size=n.poss, replace=FALSE),]  #Then sample to get the desired actual sample size
  }else{
    showModal(modalDialog(
      title = "Error","Raster and shapefile do not overlap.",easyClose = TRUE, fade=FALSE, size="s",
      footer =  modalButton("OK")
    ))
    
    validate(need(dim(coords)[1]>0,"Raster does not overlap with the shapefile"))
    
    # coords<-""
  }
  
  
  
  
  return(coords)
  
}





#  Make the trap locations from a x/y spacing, buffer and shapefile
make.trap.locs<-function(x.space,y.space,buff,shp, ras=NULL){
  
  b.box<-bbox(shp)
  
  traps.x<-seq(from=b.box[1,1], by=x.space, to=b.box[1,2])
  traps.y<-seq(from=b.box[2,1], by=y.space, to=b.box[2,2])
  
  traps<-as.data.frame((expand.grid(traps.x, traps.y)))
  colnames(traps)<-c("X","Y")
  shp.buff<-gBuffer(shp,width=-buff)   #Buffer in from the shapefuile edge
  #Remove traps that are outside the window,,,
  traps<-traps[inside.owin(traps[,1], traps[,2], shp.buff),]
  if(is.null(ras)){
  }else{
    cell.trap<-which(values(ras)%in%1)  #Which cells are in the raster mask
    cell.idx<-cellFromXY(ras,traps) #which raster cells are the traps in...?
    traps<-traps[which(cell.idx%in%cell.trap),]
  }
  return(traps)
}



#~~~~~~~~~~~~~~~~~  End functions ~~~~~~~~~~~~~~~~~~~~~~


server<-function(input, output, session) {

  
  iv <- InputValidator$new()
  iv$add_rule("scenname", sv_required())
  iv$add_rule("numb.poss", sv_required())
  iv$enable()

  
  #~~~~~~ Set up the scenarios
  scenParam <- reactiveVal()  #Can't recall why needed..
  
  #Code to delete all scenarios when the Delete All Rows is pressed
  observeEvent(input$deleteAllRows,{
    t<-scenParam()
    t <- t[-(1:dim(t)[1]),]
    scenParam(t)
  })
  
  #Code to delete selected scenarios when Delete Selected Rows is pressed
  observeEvent(input$deleteRows,{
    t<-scenParam()
    if (!is.null(input$tableDT_rows_selected)) {
      t <- t[-as.numeric(input$tableDT_rows_selected),]
      rownames(t)<-1:dim(t)[1]
    }
    scenParam(t)
  })
  
  #Code to add a scenario when the Add Scenario button is pressed
  observeEvent(input$update,{
    shp<-mydata.shp()$shp

    if(input$scenname==""){
      showModal(modalDialog(
        title = "Error","Add a scenario name",easyClose = TRUE, fade=FALSE, size="s",
        footer =  modalButton("OK")
      ))
    }
    
    validate(need(!input$scenname=="", "Provide a scenario name"))
    scen.name=input$scenname
    #~~~ Trapping ~~~
    x.space.a = NA
    y.space.a = NA
    trap.start.a = NA
    trap.nights.a = NA
    check.interval.a = NA
    g.mean.a=NA
    g.zero.a=NA    #This is the proportion that is untrappable, not g0....
    trap.mask=NA
    
    x.space.b = NA
    y.space.b = NA
    trap.start.b = NA
    trap.nights.b = NA
    check.interval.b = NA
    g.mean.b=NA
    g.zero.b=NA
    
    
    
    #Add the trap stuff
    if(input$trap_methods==1){    
      x.space.a = input$traps.x.space.a
      y.space.a = input$traps.y.space.a
      trap.start.a = input$trap.start.a
      trap.nights.a = input$trap.nights.a
      check.interval.a = input$n.check.a
      g.mean.a=input$g0.mean.a
      g.zero.a=input$g0.zero.a
      
      if(input$trap_mask==1){ #Work out the trapped area and save the file path to the raster that will be used.
        trap.mask = mydata.trap()$myraster
      }
      
      #For a second type of trapping....Might not be needed...
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
    hunt.start.a = NA
    hunt.days.a = NA
    hunt.eff.a = NA
    hunt.rho.a = NA
    hunt.mask = NA
    if(input$hunt_methods==1){
      hunt.start.a = input$hunt.start.a
      hunt.days.a = input$hunt.days.a
      hunt.eff.a = input$hunt.effort.a
      hunt.rho.a = input$hunt.rho.a
      
      if(input$hunt_mask==1){ #Work out the hunted area and save the file path to the raster that will be used.
        hunt.mask = mydata.hunt()$myraster
      }
      
    }
    

    #~~~ Bait Stations ~~~
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

    
    #~~~ Poisoning ~~~
    pois.start.a = NA
    pois.days.a = NA
    # pois.prop.a = NA
    pois.pkill.a = NA
    pois.area.a = NA
    pois.mask = NA
    if(input$pois_methods==1){
      pois.start.a = input$pois.start.a
      pois.days.a = input$pois.days.a
      # pois.prop.a = input$pois.prop.a
      pois.pkill.a = input$pois.pkill.a
      pois.area.a = round(gArea(shp)/10000,0)
      if(input$pois_mask==1){ #Work out the poisoned area and save the file path to the raster that will be used.
        r.tmp<-mask(disaggregate(mydata.pois()$ras.pois,100), shp)
        pois.area.a<-sum(values(r.tmp), na.rm=TRUE)*prod(res(r.tmp))/10000
        pois.mask = mydata.pois()$myraster
      }
      
    }
        
    #The letter codes used to help build up the 'Methods' part of the scenario names
    scen.code = ""
    if(input$trap_methods==1){
      scen.code<-paste0(scen.code,"T") #trapping
    }
    if(input$bait_methods==1){
      scen.code<-paste0(scen.code,"B") #baiting
    }
    if(input$hunt_methods==1){
      scen.code<-paste0(scen.code,"H") #hunting
    }
    if(input$pois_methods==1){
      scen.code<-paste0(scen.code,"A") #aerial toxin
    }
    
    
    
    #Full scenario to add.
    to_add <- data.frame(
      scen.name = scen.name,
      x.space.a = x.space.a,
      y.space.a = y.space.a,
      trap.start.a = trap.start.a,
      trap.nights.a = trap.nights.a,
      check.interval.a = check.interval.a,
      g.mean.a=g.mean.a,
      g.zero.a=g.zero.a,
      trap.mask=trap.mask,
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
      bait.g.zero.a = bait.g.zero.a,
      pois.start.a = pois.start.a,
      pois.days.a = pois.days.a,
      # pois.prop.a = pois.prop.a,
      pois.area.a = pois.area.a,
      pois.pkill.a = pois.pkill.a,
      pois.mask=pois.mask,
      hunt.start.a = hunt.start.a,
      hunt.days.a = hunt.days.a,
      hunt.eff.a = hunt.eff.a,
      hunt.rho.a = hunt.rho.a,
      hunt.mask = hunt.mask,
      methods = scen.code
      # hunt.start.b = hunt.start.b,
      # hunt.days.b = hunt.days.b,
      # hunt.eff.b = hunt.eff.b
      
    )
    newScenParam <- rbind(scenParam(),to_add) # adding new data
    scenParam(newScenParam) # updating data
    
    #Test for duplicates and blank scenarios
    t<-scenParam()
    if(sum(duplicated(t[,2:33]))==1){
      t<-t[!duplicated(t[,2:33]), ]        #Dont add duplicates
      showModal(modalDialog(
        title = "Error: Duplicate scenario!","This control scenario has already been added",easyClose = TRUE, fade=FALSE, size="s",
        footer =  modalButton("Cancel", icon=icon("exclamation"))
      ))
    }else if (sum(rowSums(is.na(t))==dim(t)[2]-1)==1){
      
      t<-t[!rowSums(is.na(t))==dim(t)[2]-1,]  #Dont add blank scenarios
      showModal(modalDialog(
        title = "Error: Blank scenario!","No control methods were selected",easyClose = TRUE, fade=FALSE, size="s",
        footer =  modalButton("Cancel", icon=icon("exclamation"))
      ))
    }else{
      showModal(modalDialog(
        title = "Success",paste0(input$scenname," added"),easyClose = TRUE, fade=TRUE, size="s",
        footer =  modalButton("OK", icon=icon("smile"))
      ))
      id<-dim((t))[1]+1
      updateTextInput(session,"scenname", "Scenario name", value=sprintf("S%03d", id))
    }
    scenParam(t)
    

    
    # scenParam<-rbind(scenParam,to_add)
    return(list(scenParam=scenParam))
  })
  
  

  #Read in the trapping mask. And output the raster as well as the path to it.
  mydata.trap<-reactive({
    if(input$trap_mask==1){
      if(is.null(input$trap_asc)==FALSE){
        myraster<-input$trap_asc$datapath   #basename for filename, dirname
        # myraster<-input$hunt_asc$basename
        ras.trap<-raster(myraster)
        
      }}
    return(list(ras.trap=ras.trap, myraster=myraster))
  })
  
  #Read in the hunting mask. And output the raster as well as the path to it.
  mydata.hunt<-reactive({
    if(input$hunt_mask==1){
      if(is.null(input$hunt_asc)==FALSE){
      myraster<-input$hunt_asc$datapath   #basename for filename, dirname
      # myraster<-input$hunt_asc$basename
      ras.hunt<-raster(myraster)

    }}
    return(list(ras.hunt=ras.hunt, myraster=myraster))
  })
  
  #Read in the aerial mask
  mydata.pois<-reactive({
    if(input$pois_mask==1){
      if(is.null(input$pois_asc)==FALSE){
        myraster<-input$pois_asc$datapath
        ras.pois<-raster(myraster)
      }}
    return(list(ras.pois=ras.pois, myraster=myraster))
  })
  
  
  #Read in the habitat mask
  mydata.habitat<-reactive({
    if(input$ras_hab=='Hab'){
      if(is.null(input$habitat_asc)==FALSE){
        myraster<-input$habitat_asc$datapath
        ras.habitat<-raster(myraster)
      }}
    return(list(ras.habitat=ras.habitat, myraster=myraster))
  })
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Read in the main shapefile - default has it loaded with Mahia (area_type=RC)
  mydata.shp<-reactive({
    
    #1. The default shape - Mahia Peninsula - defined at the very top
    shp<-readOGR("Shapefiles",def.shp)
    proj4string<-crs(shp)
    if(input$area_type=="MP"){
      shp<-readOGR("Shapefiles",def.shp)
    }
    #2. 'Upload Shapefile, then read in all the components, 
    if(input$area_type=="Map"){
      
      if(is.null(input$shp.file)==FALSE){
        myshape<-input$shp.file
        
        
        if(nrow(myshape)==1){  #zipped files
          dir<- dirname(myshape$datapath[1])
          unzip(myshape$datapath, exdir = dir)
        }else{
          dir<-dirname(myshape[1,4])
          for ( i in 1:nrow(myshape)) {
            file.rename(myshape[i,4], paste0(dir,"/",myshape[i,1]))
          }}
        
        
        getshp <- list.files(dir, pattern="*.shp", full.names=TRUE)
        shp<-readOGR(getshp)
        # shp<-gBuffer(shp, width=1) #Fixes some issues with orphaned holes
      }
    }
    
    #3. An indicative area - currently not an option.
    if(input$area_type=="Area"){
      ha.tmp<-input$area.ha
      size<-sqrt(ha.tmp*10000)
      xy<-c(1730000,5620000)
      a<-cbind(c(0,size,size,0, 0)+xy[1],c(0,0,size,size,0)+xy[2])
      ID<-"a"
      shp<-SpatialPolygons(list(Polygons(list(Polygon(a)),ID)), proj4string=CRS(proj4string))
    }
    
    proj4string<-crs(shp)  #get the proj string from the shapefile itself
    shp<-gBuffer(shp, width=1) #Can fix some orphaned holes issues
    ha<-sapply(slot(shp, "polygons"), slot, "area")/10000  
    return(list(shp=shp, p4s=proj4string, ha=ha))
  })
  

  
  
  # Make the traps and animals on the interactive map. Not linked to the actual simulation
  mydata.map<-reactive({
    shp<-mydata.shp()$shp    
    traps.x.space<-as.numeric(input$traps.x.space.i)
    traps.y.space<-as.numeric(input$traps.y.space.i)
    buff<-as.numeric(input$traps.buff.i)
    traps<-make.trap.locs(traps.x.space, traps.y.space, buff, shp)
    
    #Temporary pest animal coordinates...
    n.poss<-input$numb.poss#.i
    if(is.na(n.poss)){
      n.poss<-0
    }


    #Temporary commented out.    
    # if(is.null(input$ras.1)==TRUE){
    if((input$ras_hab)=="Ran"){
      #1. Random locations.
      n.poss.tmp<-(runifpoint(n.poss,shp))
      animals.xy.ini<-as.data.frame(n.poss.tmp)
    }else{
      #2. Grid specific densities
      # infile.ras <-input$ras.1
      # ras.2<-raster(infile.ras$datapath)
      #Call the function
      validate(need(input$habitat_asc !="","Upload a habitat raster"))
      # need(input$traps.x.space.a != "", "Please enter a value for the X trap spacing"),
      ras.habitat<-raster(mydata.habitat()$myraster)
      
      # ras.2<-raster(mydata.habitat()$myraster)
      animals.xy.ini<-get.pest.locs(ras.habitat, n.poss, shp)
    }
    colnames(animals.xy.ini)<-c("X","Y")
    
    return(list(traps=traps, animals.xy.ini=animals.xy.ini))
  })
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~  Now simulate the actual trapping...~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #When the simulations get set to run, automatically switch to the results tab. Noice!
  observeEvent(input$act.btn.trapsim,{
    updateTabsetPanel(session, "inTabset",selected = "4. Results")
  })
  
  #After pushing the RunTrapSim button,...
  datab<-eventReactive(input$act.btn.trapsim,{
    buffer<-100  #A single buffer
    #Some validation stuff - need to fill this out more completely. - actually this wont work - now that we make a table of scenarios...
    validate(
      # need(input$traps.buff.a != "", "Please enter a value for the trap buffer"),
      need(input$traps.x.space.a != "", "Please enter a value for the X trap spacing"),
      need(input$traps.y.space.a != "", "Please enter a value for the Y trap spacing"),
      need(input$n.nights != "", "Please enter a value for the number of Nights"),
      need(input$n.check.a != "", "Please enter a value for the Check interval"),
      need(input$p.bycatch.a != "", "Please enter a value for the Daily bycatch rate"),
      need(input$max.catch.a != "", "Please enter a value for the prob of Max catch"),
      need(input$numb.poss != "", "Please enter a value for the number of animals ")
    )

    #Read in all the parameter values for the scenarios
    params<-as.data.frame(scenParam())
    #Set up some places to store results 
    n.scen<-dim(params)[1]    
    pop.size.list<-vector("list",n.scen)
    trap.catch.list<-vector("list",n.scen)
    hunt.catch.list<-vector("list",n.scen)
    bait.catch.list<-vector("list",n.scen)
    pois.catch.list<-vector("list",n.scen)
    # pop.size.zone.mat[[ii]][,t+1]
    
    withProgress(message="Running simulation ",value=0,{
      #Go through the scenarios one by one...
      for(kk in 1:n.scen){  #For each scenario...
        incProgress(kk/n.scen, detail = paste("Doing scenario ", kk," of", n.scen))
        
        #Pass the parameters from params to a parameter name that will be used..
        # Traps 1
        trap.start.a<-params$trap.start.a[kk]
        trap.nights.a<-params$trap.nights.a[kk]
        n.check.a<-params$check.interval.a[kk]
        x.space.a<-params$x.space.a[kk]
        y.space.a<-params$y.space.a[kk]
        buffer.a<-buffer
        g0.mean.a<-params$g.mean.a[kk]
        g.zero.a<-params$g.zero.a[kk]
        trap.mask<-params$trap.mask[kk]
        
        #Traps 2 - why do we have two trap types...?
        trap.start.b<-params$trap.start.b[kk]
        trap.nights.b<-params$trap.nights.b[kk]
        n.check.b<-params$check.interval.b[kk]
        x.space.b<-params$x.space.b[kk]
        y.space.b<-params$y.space.b[kk]
        buffer.b<-buffer
        g0.mean.b<-params$g.mean.b[kk]
        g.zero.b<-params$g.zero.b[kk]
        
        #Bait 1
        bait.start.a<-params$bait.start.a[kk]
        bait.nights.a<-params$bait.nights.a[kk]
        bait.check.a<-params$bait.check.a[kk]
        bait.x.space.a<-params$bait.x.space.a[kk]
        bait.y.space.a<-params$bait.y.space.a[kk]
        bait.buff.a<-buffer
        bait.g0.mean.a<-params$bait.g.mean.a[kk]
        bait.g.zero.a<- params$bait.g.zero.a[kk]
        
        hunt.start.a<-params$hunt.start.a[kk]
        hunt.days.a<-params$hunt.days.a[kk]
        hunt.eff.a<-params$hunt.eff.a[kk]
        hunt.rho.a<-params$hunt.rho.a[kk]
        hunt.mask<-params$hunt.mask[kk]
      
        pois.start.a<-params$pois.start.a[kk]
        pois.days.a<-params$pois.days.a[kk]
        # pois.prop.a<-params$pois.prop.a[kk]
        pois.pkill.a<-params$pois.pkill.a[kk]
        pois.mask<-params$pois.mask[kk]
        
        
        
        #How long to run the simulation for.
        n.nights<-input$n.nights

        #The g0 uncertainty, and sigma values for traps and bait. The sds are not in params - but should be!! 
        g0.sd.a<-input$g0.sd.a
        g0.sd.b<-input$g0.sd.b
        bait.g0.sd.a<-input$bait.g0.sd.a
        
        sigma.mean<-input$sigma.mean 
        sigma.sd<-input$sigma.sd
        # rmax.poss<-input$rmax.poss
        ann.growth.poss<-input$ann.growth.poss
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
          if(is.na(trap.mask)==TRUE){
            traps.a<-make.trap.locs(x.space.a, y.space.a,buffer.a,shp)
          }else{
            traps.a<-make.trap.locs(x.space.a, y.space.a,buffer.a,shp, ras=raster(trap.mask))            
          }

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
        
        #Mask stuff for pois and hunting
        if(is.na(pois.mask)==FALSE){
          ras.pois<-raster(pois.mask)
          cell.pois<-which(values(ras.pois)%in%1)  #which cells are poisoning ones...?
        }
        
        if(is.na(hunt.mask)==FALSE){
          ras.hunt<-raster(hunt.mask)
          cell.hunt<-which(values(ras.hunt)%in%1)  #which cells are hunting ones...?
        }
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
        # When we have traps//baits etc, then calculate the interval and the checking interval and bycatch/failure
        #for the traps and the bait stations
        if(is.na(trap.start.a)==FALSE){
          trap.period.a<-seq(from=trap.start.a, to=(trap.start.a+trap.nights.a-1), by=1)
          # #This sets the trap checking interval. i.e. traps are cleared and reset on these nights only...
          check.vec.a<-seq(from=trap.start.a, to=(trap.start.a+trap.nights.a), by=n.check.a)
          p.bycatch.a<-input$p.bycatch.a
          
        }
        # if(input$show_trap_b==1){
        if(is.na(trap.start.b)==FALSE){
          trap.period.b<-seq(from=trap.start.b, to=(trap.start.b+trap.nights.b-1), by=1)
          check.vec.b<-seq(from=trap.start.b, to=(trap.start.b+trap.nights.b), by=n.check.b)
          p.bycatch.b<-input$p.bycatch.b
        }
        
        if(is.na(bait.start.a)==FALSE){
          bait.period.a<-seq(from=bait.start.a, to=(bait.start.a+bait.nights.a-1), by=1)
          # #This sets the checking interval. - cleared and reset on these nights only...
          bait.check.vec.a<-seq(from=bait.start.a, to=(bait.start.a+bait.nights.a), by=bait.check.a)
          p.failure.a<-0#input$p.bycatch.a
        }
        
        if(is.na(pois.start.a)==FALSE){
          pois.period.a<-seq(from=pois.start.a, to=(pois.start.a+pois.days.a-1), by=1)
          # pois.pkill.actual<-pois.pkill.a*pois.prop.a/(10000)  #Check why this is in here...
          pois.pkill.actual<-pois.pkill.a/(100)  
          pois.daily.pkill<-1-((1-pois.pkill.actual)^(1/pois.days.a))
        }
        
        
        if(is.na(hunt.start.a)==FALSE){
          hunt.period.a<-seq(from=hunt.start.a, to=(hunt.start.a+hunt.days.a-1), by=1)
          
          shp<-mydata.shp()$shp
          #Need to work out the prob kill based on the hunting area not necessairly the entire shapefile...
          if(is.na(hunt.mask)==FALSE){
            r.tmp<-mask(disaggregate(ras.hunt,100), shp)
            ha<-sum(values(r.tmp), na.rm=TRUE)*prod(res(r.tmp))/10000
          }else{
            ha<-mydata.shp()$ha
          }
          
          Eff<-hunt.eff.a/(ha/100)
          
          # H<-log(Eff+1)
          hunt.daily.pkill<-1-exp(-(hunt.rho.a*Eff))
          
        }
        
        
        
        
        
        
                
        #Carrying capacity for the area in terms of total number of animals - should be grid based?
        K.tot<-K.poss*ha
        
        n.poss<-input$numb.poss #- should this be a parameter...? Or better to leave - maybe leave cause can re-run with same params.
        if(is.na(n.poss)){
          n.poss<-0
        }
        max.catch.a<-input$max.catch.a
        max.catch.b<-input$max.catch.b
        
        #~~~~Calculate the costs~~~~
        trap.cost.sim<-0
        bait.cost.sim<-0
        hunt.cost.sim<-0
        pois.cost.sim<-0
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
          bait.cost.sim<-bait.cost.func(a=checks, b=n.baits.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a, f=input$bait.cost.a)        
          }
        
        
        
          if(is.na(hunt.start.a)==FALSE){
            hunt.cost.sim<-hunt.cost.func(a=input$day.rate.hunt.a, b=hunt.days.a)    
          }
        
        
        if(is.na(pois.start.a)==FALSE){
          pois.ha<-ha
          # if(input$pois_mask==1){
            if(is.na(pois.mask)==FALSE){
              ras.pois<-raster(pois.mask)
              r.tmp<-mask(disaggregate(ras.pois,100), shp)
            # r.tmp<-mask(disaggregate(mydata.pois()$ras.pois,100), shp)
            pois.ha<-sum(values(r.tmp), na.rm=TRUE)*prod(res(r.tmp))/10000
          }
          
          pois.cost.sim<-pois.cost.func(a=input$pois.per.ha.a, b=pois.ha, c=100)#pois.prop.a)
        }
        
        # if (input$sim_type=='grid'){
        #   cell.width<-input$cell.width
        #   cell.area.m2<-cell.width^2
        #   
        #   #Make the grid...
        #   shp<-st_as_sf(shp)
        #   shp_grid<-shp%>%st_make_grid(cellsize=cell.width, what="polygons")%>%st_intersection(shp)
        #   grid.traps.master<-lengths(st_intersects(shp_grid, st_as_sf(traps)))#*max.catch  #Number in each cell
        #   
        # }
        
        n_its<-input$n.its  #Number of iterations
        pop.size.mat<-matrix(NA,nrow=n_its,ncol=n.nights+1)
        trap.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of trapping - for each iteration - how many that night/day
        bait.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of baiting - for each iteration - how many that night/day
        hunt.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of hunting
        pois.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of poisoning
        # pop.size.zone.vec<-vector("list",n_its)
        
        
        for(ii in 1:n_its){
          pop.size.mat[ii,1]<-n.poss
          
          #~~~~~~~~~Make some animals~~~~~~~~~~
          # if(is.null(input$ras.1)==TRUE){
          if((input$ras_hab)=="Ran"){
            #1. Random locations.
            n.poss.tmp<-(runifpoint(n.poss,shp))
            animals.xy<-as.data.frame(n.poss.tmp)
          }else{
            #2. Grid specific densities
            ras.habitat<-raster(mydata.habitat()$myraster)
            
            animals.xy<-get.pest.locs(ras.habitat, n.poss, shp)
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
          
          #Calculate the g0 values...
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
          
          #g0 for bait stations
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
          #The second is the trap+animal pairwise  - i.e. it works out the probability of capture for each device and animal, by device type.
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
          
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          #~~~Now run the simulation of trapping the animals...
          # withProgress(message="Running Simulation",value=0,{
          #Set up some storage matrices
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

          pois.catch.vec<-rep(0, n.nights)
          hunt.catch.vec<-rep(0, n.nights)
          # pop.size.zone.vec[[ii]]<-matrix(NA,nrow=4,ncol=n.nights+1)
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
            
            #Now the poison operation...
            if(is.na(pois.start.a)==FALSE){
              if(t%in%pois.period.a==TRUE){
                
                # if(input$pois_mask==1){
                  # ras.pois<-mydata.pois()$ras.pois
                  if(is.na(pois.mask)==FALSE){
                  # ras.pois<-raster(pois.mask)  #Moved up earlier
                  
                  cell.idx<-cellFromXY(ras.pois,animals.xy[,1:2])
                  # cell.idx
                  # cell.pois<-which(values(ras.pois)%in%1)  #which cells are hunting ones...?
                  # cell.hunt
                  idx.target<-which(cell.idx%in%cell.pois & animals.xy$Dead==0) #Which animals are targets...
                  n.pois<-rbinom(1,length(idx.target),prob=pois.daily.pkill)#hunt.daily.pkill) #How many will be killed
                  idx.kill<-sample(idx.target,n.pois, replace=FALSE)     #Which ones will be killed
                  
                }else{
                alive<-which(animals.xy$Dead==0)
                n.pois<-rbinom(1,length(alive),prob=pois.daily.pkill) #How many will be killed
                idx.kill<-sample(alive,n.pois, replace=FALSE)     #Which ones will be killed
                }
                animals.xy$Dead[idx.kill]<-1
                pois.catch.vec[t]<-n.pois
              }}
            
            #Hunting
            if(is.na(hunt.start.a)==FALSE){
              if(t%in%hunt.period.a==TRUE){

                # if(input$hunt_mask==1){
                  # ras.hunt<-mydata.hunt()$ras.hunt
                  if(is.na(hunt.mask)==FALSE){
                    # ras.hunt<-raster(hunt.mask)
                  
                  cell.idx<-cellFromXY(ras.hunt,animals.xy[,1:2])
                  # cell.hunt<-which(values(ras.hunt)%in%1)  #which cells are hunting ones...?
                  
                  # cell.idx and cell.hunt from earlier
                  idx.target<-which(cell.idx%in%cell.hunt & animals.xy$Dead==0) #Which animals are targets...
                  n.hunt<-rbinom(1,length(idx.target),prob=hunt.daily.pkill)#hunt.daily.pkill) #How many will be killed
                  idx.kill<-sample(idx.target,n.hunt, replace=FALSE)     #Which ones will be killed
                }else{
                
                alive<-which(animals.xy$Dead==0)
                n.hunt<-rbinom(1,length(alive),prob=hunt.daily.pkill) #How many will be killed
                idx.kill<-sample(alive,n.hunt, replace=FALSE)     #Which ones will be killed
                }

                animals.xy$Dead[idx.kill]<-1
                
                hunt.catch.vec[t]<-n.hunt

              }}

            
            #~~~~~~~~~~~~~~~~~~~~~Reproduction~~~~~~~~~~~~~~~~~~~~~~~~~
            #Get the reproductive vector when we hit the start of the interval...
            #If t is one of the start days of the breeding season...
            if(t%in%rep.start.vec){
              #Discrete version of rmax
              # rd<-exp(rmax.poss)-1
              rd<-ann.growth.poss/100
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
                if((input$ras_hab)=="Ran"){
                  
                  #1. Random locations.
                  
                  new.animals.xy<-as.data.frame(runifpoint(N.new,shp))
                }else{
                  #2. Grid specific densities
                  # infile.ras <-input$ras.1
                  # ras.2<-raster(infile.ras$datapath)
                  #Call the function
                  ras.habitat<-raster(mydata.habitat()$myraster)
                  
                  new.animals.xy<-get.pest.locs(ras.habitat, N.new, shp)
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
            if(is.na(pois.start.a)==FALSE){
              pois.catch.mat[ii,]<-pois.catch.vec #This is a single number for each night
            }else{
              pois.catch.mat[ii,]<-rep(0, n.nights)
            }
            if(is.na(hunt.start.a)==FALSE){
              hunt.catch.mat[ii,]<-hunt.catch.vec #This is a single number for each night
            }else{
              hunt.catch.mat[ii,]<-rep(0, n.nights)
            }
          }
        }#End of iteration ii
        
        
        # pop.zone.list[[kk]]<-apply(simplify2array(pop.size.zone.vec),c(1,2), mean)
        params$TrapCost[kk]<-trap.cost.sim
        params$BaitCost[kk]<-bait.cost.sim
        params$HuntCost[kk]<-hunt.cost.sim
        params$PoisCost[kk]<-pois.cost.sim
        # params$HuntCost[kk]<-hunt.cost.sim
        # params$TotalCost[kk]<-hunt.cost.sim + trap.cost.sim
        params$TotalCost[kk]<-trap.cost.sim+bait.cost.sim+hunt.cost.sim+pois.cost.sim
        
        params$MeanPopSize[kk]<-round(mean(pop.size.mat[,n.nights+1]),2)
        
        pop.size.list[[kk]]<-pop.size.mat
        hunt.catch.list[[kk]]<-hunt.catch.mat   #This will be all 0s at the moment
        trap.catch.list[[kk]]<-trap.catch.mat
        bait.catch.list[[kk]]<-bait.catch.mat
        pois.catch.list[[kk]]<-pois.catch.mat
        # cat(file=stderr(), "drawing histogram with", kk, "bins", "\n")
      } #End kk    
      
      
      #Add a comment 
    })  #End of progress
    
    
    
    # params<-params[,c(36,35,31,32,33,34,1:30)]
    # params<-params[,c(37,36,32,33,34,35,1:31)]
    # params<-params[,c(37,36,32,33,34,35,1:7,15:31)]
    # params<-params[,c(1,33,39,38,34,35,36,37,2:8,16:32)]
    params<-params[,c(1,34,40,39,35,36,37,38,2:9,17:33)]

    return(list(trap.catch.mat=trap.catch.mat, bait.catch.mat=bait.catch.mat, hunt.catch.list=hunt.catch.list, pois.catch.list=pois.catch.list, pop.size.mat=pop.size.mat, animals.xy=animals.xy, hunt.catch.mat=hunt.catch.mat, params=params, pop.size.list=pop.size.list, trap.catch.list=trap.catch.list, bait.catch.list=bait.catch.list))#, pop.zone.list=pop.zone.list))#, animals.done.xy=animals.xy))    
  }
  )
  
  
  
  output$mymap<-renderLeaflet({
    shp<-mydata.shp()$shp
    # shp<-mydata.zone()$shp.2
      
    shp.proj<-spTransform(shp,CRS("+proj=longlat +datum=WGS84"))
    # ha<-mydata()$ha
    
    m<-leaflet() %>%
      addTiles(group="Default")%>%
      # addProviderTiles("Esri.WorldTopoMap", group = "Topo")%>%
      addProviderTiles("OpenTopoMap", group = "Topo")%>%
      addProviderTiles("Esri.WorldImagery", group = "Aerial")%>%
      # addProviderTiles("Esri.WorldStreetMap", group = "Street")%>%
      addPolygons(data=shp.proj, weight=2, fillColor="grey40", color="black", fillOpacity=0.2)
    m
    
  })
  
  
  
  observe({
    # validate(need(dim(coords)[1]>0,"Raster does not overlap with the shapefile"))
    traps<-mydata.map()$traps
    animals.xy<-mydata.map()$animals.xy.ini
    proj4string<-mydata.shp()$p4s
    
    traps.proj<-proj4::project(traps, proj=proj4string, inverse=T) 
    tmp<-(proj4::project(animals.xy[,1:2], proj=proj4string, inverse=T))
    animals.xy$Lat<-tmp$y
    animals.xy$Lon<-tmp$x
    map<-leafletProxy("mymap")
    map%>%clearMarkers()
    
    map%>%addCircleMarkers(lng=traps.proj$x,lat=traps.proj$y, radius=3, color="black", weight=1, fill=TRUE, fillColor="red", fillOpacity=1, stroke=TRUE, group="Devices")
    map%>%addCircleMarkers(lng=animals.xy$Lon,lat=animals.xy$Lat, radius=5, color="black", weight=1, fill=TRUE, fillColor=cols.vec[2], fillOpacity=1, stroke=TRUE, group="Animals")
    
    map%>%addLayersControl(
      baseGroups = c("Default","Topo","Aerial"),
      overlayGroups=c("Devices","Animals"),
      options = layersControlOptions(collapsed = FALSE)
    )
    
  })
  
  
  
  
  output$plot1<-renderPlot({
    #This is the plot of cost vs number left...First plot on Results tab
    # validate(need(input$results.table_rows_selected !="","Select a results row"))
    idx<-as.numeric(input$results.table_rows_selected)
    
    res<-datab()$params
    # par(mar=c(4,4,1,1), mgp=c(2.5,1,0), tcl=-0.25)
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(res$MeanPopSize, res$TotalCost, xlab="Final population", ylab="Total cost", type='n')
    points(res$MeanPopSize, res$TotalCost, pch=21, bg=cols.vec[2], cex=1.5)    
    points(res$MeanPopSize[idx], res$TotalCost[idx], pch=21, bg=cols.vec[6], cex=2)    
    mtext("Total Cost vs Animals Remaining",3, cex=1.5, line=1)
  })
  
  
  
  output$plot2<-renderPlot({
    validate(need(input$results.table_rows_selected !="","Select a results row from the table below"))
    idx<-as.numeric(input$results.table_rows_selected)
    
    trap.catch.list<-datab()$trap.catch.list
    trap.catch.mat<-trap.catch.list[[idx]]
    bait.catch.list<-datab()$bait.catch.list
    bait.catch.mat<-bait.catch.list[[idx]]
    pois.catch.list<-datab()$pois.catch.list
    pois.catch.mat<-pois.catch.list[[idx]]
    hunt.catch.list<-datab()$hunt.catch.list
    hunt.catch.mat<-hunt.catch.list[[idx]]
    
    nights.vec<-1:input$n.nights
    ymax.t<-max(cumsum(colMeans(trap.catch.mat)))
    ymax.b<-max(cumsum(colMeans(bait.catch.mat)))
    ymax.p<-max(cumsum(colMeans(pois.catch.mat)))
    ymax.h<-max(cumsum(colMeans(hunt.catch.mat)))
    ymax<-max(ymax.b, ymax.t, ymax.p, ymax.h)
    
    par(mar=c(4,4,3,2), tcl=-.2, mgp=c(2.5,1,0))
    plot(1,1,xlim=c(0,input$n.nights), ylim=c(0,ymax), type='n', xlab="Nights", ylab="Cumulative kills", las=1)
    for(i in 1:input$n.its){
      lines(nights.vec,cumsum(trap.catch.mat[i,]), col=cols.vec[1])
      lines(nights.vec,cumsum(bait.catch.mat[i,]), col=cols.vec[3])
      lines(nights.vec,cumsum(pois.catch.mat[i,]), col=cols.vec[5])
      lines(nights.vec,cumsum(hunt.catch.mat[i,]), col=cols.vec[7])
    }
    lines(nights.vec,cumsum(colMeans(trap.catch.mat)), col=cols.vec[2])
    points(nights.vec,cumsum(colMeans(trap.catch.mat)), bg=cols.vec[2], pch=21)
    
    lines(nights.vec,cumsum(colMeans(bait.catch.mat)), col=cols.vec[4])
    points(nights.vec,cumsum(colMeans(bait.catch.mat)), bg=cols.vec[4], pch=22)

    lines(nights.vec,cumsum(colMeans(pois.catch.mat)), col=cols.vec[6])
    points(nights.vec,cumsum(colMeans(pois.catch.mat)), bg=cols.vec[6], pch=23)
    
    lines(nights.vec,cumsum(colMeans(hunt.catch.mat)), col=cols.vec[8])
    points(nights.vec,cumsum(colMeans(hunt.catch.mat)), bg=cols.vec[8], pch=24)

    legend("topleft", legend=c("Traps","Bait station","Aerial poison", "Hunting"), pch=c(21,22,23, 24), pt.bg=cols.vec[c(2,4,6,8)], bty="n")
    mtext("Cumulative kills",3, cex=1.5, line=1)
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
    validate(need(input$results.table_rows_selected !="","Select a results row"))
    
    idx<-as.numeric(input$results.table_rows_selected)
    # idx<-as.numeric(input$result_scenario)
    
    pop.size.list<-datab()$pop.size.list
    # pop.size.mat<-pop.size.list[[as.numeric(input$result_scenario)]]
    pop.size.mat<-pop.size.list[[idx]]
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
  
  #This updates the dropdown box in the results section
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
  
  
output$plot_hunt<-renderPlot({  
  
  
  hunt.eff<-seq(from=0.1,to=20, by=0.1)
  
  # ha.hunt<-mydata.shp()$hunt.area.a
  #Calculate the Effort - i.e. the distance /area to be hunted
  if(input$hunt_mask==1){
    validate(need(input$hunt_asc != "", "Upload a tif or asc file to mask the hunting area"),)
    shp<-mydata.shp()$shp
    r.tmp<-mask(disaggregate(mydata.hunt()$ras.hunt,100), shp)
    ha<-sum(values(r.tmp), na.rm=TRUE)*prod(res(r.tmp))/10000
  }else{
    ha<-mydata.shp()$ha  
  }
  
  
  
  
  ha.km2<-ha/100
  
  Eff<-hunt.eff/ha.km2
  hunt.daily.pkill<-1-exp(-(input$hunt.rho.a*Eff))
  plot(hunt.eff, hunt.daily.pkill, xlab="Km per day", ylab="Prob. of daily kill", las=1)
  abline(v=input$hunt.effort.a, col="red")
  
  
  abline(h=hunt.daily.pkill[match(input$hunt.effort.a,hunt.eff)], col="blue")
  mtext(round(ha,0),3, line=1)
  # mtext(hunt.ha,3, line=2)
})
  
  
  
  
  
  #Plots of the asc maps
output$plot.trap.asc<-renderPlot({
  
  validate(
    need(input$trap_asc != "", "upload a tif or asc file to mask the trapping area"),
  )
  
  ras.trap<-raster(mydata.trap()$myraster)
  plot(ras.trap)
  # plot(mydata.pois()$ras.pois)
  plot(mydata.shp()$shp, add=TRUE)
  
})



  output$plot.pois.asc<-renderPlot({
    
    validate(
           need(input$pois_asc != "", "upload a tif or asc file to mask the aerial area"),
    )
      
    ras.pois<-raster(mydata.pois()$myraster)
    plot(ras.pois)
    # plot(mydata.pois()$ras.pois)
    plot(mydata.shp()$shp, add=TRUE)
    
  })
  output$plot.hunt.asc<-renderPlot({
    
    validate(
      need(input$hunt_asc != "", "Upload a tif or asc file to mask the hunting area"),
    )
    ras.hunt<-raster(mydata.hunt()$myraster)
    plot(ras.hunt)
    # plot(mydata.hunt()$ras.hunt)
    plot(mydata.shp()$shp, add=TRUE)
    
  })
  
  output$plot.habitat.asc<-renderPlot({
    # validate(need(dim(coords)[1]>0,"Raster does not overlap with the shapefile"))
    validate(
      need(input$habitat_asc != "", "Upload a tif or asc file of the relative abundance"),
    )
    ras.habitat<-raster(mydata.habitat()$myraster)
    plot(ras.habitat)
    # plot(mydata.hunt()$ras.hunt)
    plot(mydata.shp()$shp, add=TRUE)
    
  })
  
  

  
  #Text on the indicative map
  output$text5<-renderText({
    traps<-mydata.map()$traps
    ntraps<-dim(traps)[1]
    ha.area<-mydata.shp()$ha
    
    paste(ntraps, " Traps (One per ",round(ha.area/ntraps,1)," hectare)", sep="")
  })
  
  
  output$text6<-renderText({
    hr.radius<-input$sigma.mean*2.45
    # ha.area<-input$area.ha
    hr.ha<-((pi*hr.radius^2)/10000)
    hr.km<-hr.ha/100
    return(paste0("Home range size: ",round(hr.ha,1), " hectares"))
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
  
  #Area in hectares on Tab 1
  output$text10<-renderText({
    ha<-mydata.shp()$ha
    return(paste0("Total area (ha): ", round(ha,0)))
  })

  
  #Cost of traps on the Control Methods tab
  output$text_trap_a_cost<-renderText({
    
    if(input$trap_mask==1){
      validate(
        need(input$trap_asc != "", "NA"),
      )
      trap.mask<-mydata.trap()$ras.trap
      traps.a<-make.trap.locs(input$traps.x.space.a, input$traps.y.space.a, 100, shp, ras=(trap.mask))
    }else{
      traps.a<-make.trap.locs(input$traps.x.space.a, input$traps.y.space.a, 100, shp)  
    }
    
    
    
    shp<-mydata.shp()$shp    
    
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
    cost<-bait.cost.func(a=checks, b=n.bait.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a, f=input$bait.cost.a)
    
    return(paste0(n.bait.a," Bait Stations\nCost = $", cost ))
  })
  
  output$text_hunt_a_cost<-renderText({
    hunt.cost<-hunt.cost.func(a=input$hunt.days.a, b=input$day.rate.hunt.a)
    return(paste0("Cost = $", hunt.cost ))
  })
  
  output$text_pois_a_cost<-renderText({
    ha<-mydata.shp()$ha  
    shp<-mydata.shp()$shp  
    pois.ha<-ha
    

    
    if(input$pois_mask==1){
      validate(
        need(input$pois_asc != "", "NA"),
      )
      r.tmp<-mask(disaggregate(mydata.pois()$ras.pois,100), shp)
      pois.ha<-sum(values(r.tmp), na.rm=TRUE)*prod(res(r.tmp))/10000
    }
    
    pois.cost<-pois.cost.func(a=input$pois.per.ha.a, b=pois.ha, c=100)#pois.prop.a)
    
    # pois.cost<-pois.cost.func(a=input$pois.per.ha.a, b=ha, c=100)#input$pois.prop.a)
    
    return(paste0("Cost = $", pois.cost ))
  })
  
#The number of animals per hectare shown in Tab 1 under Pest Parameters  
output$text_density<-renderText({
  need(input$numb.poss != "", "Please enter a value for the number of animals ")
  ha<-mydata.shp()$ha
  # numb<-200
  numb<-input$numb.poss
  return(paste0(round(numb/ha,2)," per ha"))
})
  

  #This contains the results 
  # output$results.table<-renderDataTable({
  #   datab()$params
  # })
  
  output$results.table<-renderDataTable(
    datab()$params, selection='single'
  )
  
  #This contains the scenarios on Tab 3: Run Scenarios
  output$tableDT <- DT::renderDataTable(
    scenParam()[c(1,34,2:9,17:33)]
  )
  
  
  # https://shiny.rstudio.com/articles/generating-reports.html
  
  
  output$report <- downloadHandler(
    #Name of the report.
    filename = function() {
      paste("Results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      
      data<-datab()$params
      
      write.csv(data, file)
    }
  )
  
  
  
  
  
  
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#         The user interface
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui<-fluidPage(theme=shinytheme("flatly"),
              title="EradSim",
              tags$head(
                tags$style(HTML("
                                .shiny-output-error-validation {
                                color: blue;
                                }
                                "))),
              fluidRow(
                column(width=6,  
                       h1("Eradication Feasibility"),
                       strong("How much control?"),
                       "This app is intended to provide guidance as to the approximate level of control required to achieve a desired level of pest reduction.",p(),
                       
                       # "v1"
                       # checkboxInput("showinfo","Background information and instructions"),

                       # h3("Trapping Simulation Tool")
                ),
                column(6,
                       img(src="manaaki_logo.png", height = 90, align="right", hspace=20,vspace=10),
                       img(src="ciss_logo.jpg", height = 80, align="right", hspace=20,vspace=10),
                       # img(src="bhnsc.png", height = 80, align="right", hspace=20,vspace=10)
                       img(src="ari_logo.jpg", height = 80, align="right", hspace=20,vspace=10)

                       # img(src="IC_logo.png", height = 90, align="right", hspace=20,vspace=10)
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
                                                     # radioButtons(inputId="area_type", label="Chose area",choices=c("Mahia Peninsula (Default)"="RC","Upload Shapefile"="Map", "Indicative Area"="Area"), selected="RC"),
                                                     # radioButtons(inputId="area_type", label="Chose area",choices=c("Ashley Forest"="RC","Mahia Peninsula"="MP","Upload Shapefile"="Map"), selected="RC"),
                                                     radioButtons(inputId="area_type", label="Chose area",choices=c("Mahia Peninsula"="MP","Upload Shapefile"="Map"), selected="MP"),
                                                     #   
                                                       conditionalPanel(
                                                         condition="input.area_type=='Area'",
                                                         numericInput(inputId = "area.ha", label="Area (ha)", value=10000, width="120px")
                                                       ),
                                                     #   
                                                     conditionalPanel(
                                                       condition="input.area_type=='Map'",
                                                       
                                                       tags$div(title="Be sure to select all four components (.dbf, .prj, .shp, .shx), or upload a .zip file",
                                                                fileInput(inputId = "shp.file", label="Select the shapefile (.zip or .dbf/.prj/.shp/.shx components)", accept=c('.shp','.dbf','.sbn','.sbx','.shx',".prj",".zip"), multiple=TRUE, width="300px")
                                                       )
                                                     ),
                                                     verbatimTextOutput("text10")
                                                   ),
                                                   h4("Pest parameters"),
                                                   wellPanel(
                                                     div(style="display:inline-block;vertical-align:top",title="Number of animals",
                                                         numericInput(inputId = "numb.poss", label="Number (total)", value=100, width="180px")
                                                     ),
                                                     div(style="display:inline-block;vertical-align:top",title="Pests distributed randomly or according to habitat (which requires reading in a raster of relative distribution)",
                                                         # fileInput(inputId = "ras.1", label="Ascii file of relative abundance", accept=c('.asc'), multiple=FALSE, width="250px")
                                                         radioButtons(inputId = "ras_hab", label="Pest distribution", choices=c("Random"="Ran","Habitat Specific"="Hab"), selected="Ran",width="250px")
                                                         # radioButtons(inputId = "ras.1", label="Relative abundance", choices=c("Random"), selected="Random",width="250px")
                                                     ),
                                                     

                                                     div(style="display:inline-block;vertical-align:top",
                                                         conditionalPanel(condition="input.ras_hab == 'Hab'",
                                                                          div(style="display:inline-block;vertical-align:top", 
                                                                              fileInput(inputId = "habitat_asc", label="Upload a raster of habitat (.asc or .,tif)", accept=c(".tif",".asc"), multiple=FALSE, width="200px")),
                                                                          div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.habitat.asc", width = "450px", height="350px"))
                                                                          # plot.habitat.asc
                                                                          # div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.hunt.asc", width = "450px", height="350px"))
                                                         )),
                                                     
                                                     
                                                     
                                                     div(style="display:inline-block;vertical-align:bottom",
                                                     verbatimTextOutput("text_density")),
                                                     tags$div(title="Sigma x 2.45 is the radius of a circle  where an indivudual spends 95% of its time.",
                                                              h5(strong("Home range (Sigma)"))),
                                                     div(style="display:inline-block",title="Sigma x 5 is the diameter of the home range where an indivudual spends 95% of its time.",
                                                         numericInput(inputId = "sigma.mean", label="Mean (metres)", value=100, width="120px")),
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
                                                         tags$div(title="This is the annual percentage increase of the population.",
                                                                  # numericInput(inputId = "rmax.poss", label="Rmax", value=0.4, width="120px")
                                                                  numericInput(inputId = "ann.growth.poss", label="Annual population growth (%)", value=20, width="120px", min=0, max=200)
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
                                                   "This map is provided to visualise the area. The device spacings specified here are ", strong("not"), " included in the simulations",
                                                   p(),
                                                   # div(style="display:inline-block;vertical-align:bottom",
                                                   #     numericInput(inputId = "numb.poss.i", label="Population size", value=10, width="180px")
                                                   # ),
                                                   div(style="display:inline-block;vertical-align:bottom",
                                                       tags$div(title="The trap spacing in the east-west direction",
                                                                numericInput(inputId = "traps.x.space.i", label="Spacing E-W (m)", value="1000",width="135px"))),
                                                   div(style="display:inline-block;vertical-align:bottom",
                                                       tags$div(title="The trap spacing in the north-south direction",
                                                                numericInput(inputId = "traps.y.space.i", label="Spacing N-S (m)", value="1000",width="135px"))),
                                                   div(style="display:inline-block;vertical-align:bottom",
                                                       tags$div(title="The buffer from the edge",
                                                                numericInput(inputId = "traps.buff.i", label="Edge buffer (m)", value="100", min=0, max=1000, width="110px"))),
                                                   div(style="width:300px;display:inline-block;vertical-align:bottom",
                                                       verbatimTextOutput("text5")),
                                                   # p(),
                                                   leafletOutput(outputId = "mymap", height = "1000px", width="1200px")%>% withSpinner(type=7)
                                                   
                                            )
                                   ),
                                   tabPanel("2. Control Methods",
                                            h4("Chose one or more control methods to build a control scenario. Click Add Scenario to add it to the Run Scenarios tab."),
                                            h4("Repeat this for all control scenarios you want to run"),
                                            
                                            # h3("Choose 'Trapping' and/or 'Hunting' and then specify inputs for up to two methods of each. Click 'Add Scenario' to add them to the 'Run Scenarios' tab "),
                                            fluidRow(
                                              column(width=2,
                                                     
                                                     wellPanel(
                                                       h4("Control method(s)"),
                                                       div(style="width:300px;height:20px;vertical-align:top",
                                                           title="Ground trapping including kill traps and/or capture traps",
                                                           checkboxInput(inputId="trap_methods", label="Trapping", value=FALSE)),
                                                       div(style="width:300px;height:20px;vertical-align:top",
                                                           title="Bait locations including bait stations, or point-based hand-laid baits",
                                                           checkboxInput(inputId="bait_methods", label="Bait stations", value=FALSE)),
                                                       div(style="width:300px;height:20px;vertical-align:top",
                                                           title="Hunting including ground or aerial hunting. The kill rate parameters will need to be changed!",
                                                           checkboxInput(inputId="hunt_methods", label="Hunting", value=FALSE)),
                                                       div(style="width:300px;height:20px;vertical-align:top",
                                                           title="Aerial poisoning  - main parameter is the %kill.",
                                                           checkboxInput(inputId="pois_methods", label="Aerial poisoning", value=FALSE))
                                                     )),
                                              column(width=1,
                                                     # wellPanel(
                                                       actionButton("update", "Add Scenario", icon("plus"),style="background-color:green")
                                              )
                                            ),
                                            
                                            h4("Scenario name"),
                                            textInput(inputId = "scenname",label="", value="S001"),
                                            conditionalPanel(
                                              condition="input.trap_methods==1",
                                              
                                              h4("Trapping"),
                                              wellPanel(

                                                div(style="display:inline-block;vertical-align:middle",h5("Trapping Method 1")),p(),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The trap spacing in the east-west direction",
                                                             numericInput(inputId = "traps.x.space.a", label="Spacing E-W (m)", value="500",width="135px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The trap spacing in the north-south direction",
                                                             numericInput(inputId = "traps.y.space.a", label="Spacing N-S (m)", value="500",width="135px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Nightly probability of by-catch, false triggers etc  ",
                                                             numericInput(inputId = "p.bycatch.a", label="Daily bycatch", value=0.01, min=0.01, max=0.99, step=0.01, width="120px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Maximum catch per trap  ",
                                                             numericInput(inputId = "max.catch.a", label="Max catch", value=1, width="135px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    numericInput(inputId = "g0.mean.a", label="Trap Probability (g0) Mean", value=0.1, min=0.01, max=0.99, step=0.01, width="120px")),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    numericInput(inputId="g0.sd.a", label='StdDev', value=.01, width="120px")),
                                                # div(style="display:inline-block;vertical-align:bottom",
                                                #     verbatimTextOutput("text8")),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    numericInput(inputId = "g0.zero.a", label="Proportion untrappable", value=0.05, min=0.01, max=0.99, step=0.01, width="120px")),
                                                #Timings
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="Start night of trapping.",
                                                             numericInput(inputId = "trap.start.a", label="Start night", value="5", width="120px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="Number of nights traps are set for.",
                                                             numericInput(inputId = "trap.nights.a", label="Duration (nights)", value="30", width="120px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The checking interval of the traps. For traps that are not cleared, set equal to Nights ",
                                                             numericInput(inputId = "n.check.a", label="Check interval", value="7", width="120px"))),
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
                                                             numericInput(inputId = "cost.per.trap.a", label="Fixed cost per trap ($)", value=25, width="120px")
                                                    )),
                                                
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    verbatimTextOutput("text_trap_a_cost")),
                                                
                                                div(style="display:inline-block;vertical-align:top",checkboxInput(inputId="trap_mask", label="Trapping mask", value=FALSE)),
                                                div(style="display:inline-block;vertical-align:top",conditionalPanel(
                                                  condition="input.trap_mask==1",
                                                  div(style="display:inline-block;vertical-align:top", fileInput(inputId = "trap_asc", label="Chose the trapping mask", accept=c(".asc",".tif"), multiple=TRUE, width="200px")),
                                                  div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.trap.asc", width = "450px", height="350px"))
                                                )),
                                                
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
                                                  condition="input.show_trap_b==2",
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
                                                           numericInput(inputId = "p.failure.a", label="Daily rate of failure", value=0, min=0.01, max=0.99, step=0.01, width="120px"))),
                                              
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  numericInput(inputId = "bait.g0.mean.a", label="Bait Sation Probability (g0) Mean", value=0.1, min=0.01, max=0.99, step=0.01, width="120px")),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  numericInput(inputId="bait.g0.sd.a", label='StdDev', value=.01, width="120px")),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  numericInput(inputId = "bait.g.zero.a", label="Proportion unbaitable", value=0.05, min=0.01, max=0.99, step=0.01, width="120px")),

                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="Start night of baiting.",
                                                           numericInput(inputId = "bait.start.a", label="Start night", value="35", width="120px"))),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="Number of nights baits are set for.",
                                                           numericInput(inputId = "bait.nights.a", label="Duration (nights)", value="50", width="120px"))),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(id="redtitle",title="The checking interval of the bait stations. For stations that are not cleared, set equal to Nights ",
                                                           numericInput(inputId = "bait.check.a", label="Check interval", value="25", width="120px"))),
                                              
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Number of stations checked per day",
                                                           numericInput(inputId = "bait.per.day.a", label="Bait stations checked per day", value=40, width="120px")
                                                  )),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Labour cost - day rate ($)",
                                                           numericInput(inputId = "bait.day.rate.a", label="Day rate ($)", value=400, width="135px")
                                                  )),
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Bait costs - per station per check ($)",
                                                           numericInput(inputId = "bait.cost.a", label="Bait costs ($/station/day)", value=5, width="135px")
                                                  )),
                                              
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  tags$div(title="Fixed cost ($) per bait station",
                                                           numericInput(inputId = "cost.per.bait.a", label="Fixed cost per bait station ($)", value=60, width="120px")
                                                  )),
                                              
                                              div(style="display:inline-block;vertical-align:bottom",
                                                  verbatimTextOutput("text_bait_a_cost")),
                                              
                                              
                                              
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
                                                             numericInput(inputId = "hunt.start.a", label="Start day", value="30", width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "hunt.days.a", label="Days hunted", value="10", width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "hunt.effort.a", label="Distance per day (km)", value="5", width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "hunt.rho.a", label="Kill rate", value="0.5", min = 0.01, width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "day.rate.hunt.a", label="Day rate ($)", value=2500, width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    verbatimTextOutput("text_hunt_a_cost")),
                                                

                                                
                                                div(style="display:inline-block;vertical-align:top",checkboxInput(inputId="hunt_mask", label="Hunting mask", value=FALSE)),
                                                div(style="display:inline-block;vertical-align:top",conditionalPanel(
                                                  condition="input.hunt_mask==1",
                                                  div(style="display:inline-block;vertical-align:top", fileInput(inputId = "hunt_asc", label="Chose the hunting mask", accept=c(".tif",".asc"), multiple=TRUE, width="200px")),
                                                  div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.hunt.asc", width = "450px", height="350px"))
                                                )),br(),
                                                div(style="display:inline-block;vertical-align:top",checkboxInput(inputId="hunt_effort", label="Effort vs prob of kill", value=FALSE)),
                                                div(style="display:inline-block;vertical-align:top",conditionalPanel(
                                                  condition="input.hunt_effort==1",
                                                  plotOutput(outputId = "plot_hunt",width = "450px", height="350px")))
                                                
                                              )#End of wellpanel
                                            ),#End of conditional panel 
                                            
                                            
                                            
                                            conditionalPanel(
                                              condition="input.pois_methods==1",
                                              
                                              h4("Aerial poison"),
                                              wellPanel(
                                                
                                                
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "pois.start.a", label="Start day", value="60", width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "pois.days.a", label="Operation length", value="5", width="150px"))),
                                                # div(style="display:inline-block;vertical-align:top",
                                                #     tags$div(title="help text ",
                                                #              numericInput(inputId = "pois.prop.a", label="Percentage of area poisoned", value="50", width="120px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "pois.pkill.a", label="Percent kill", value="90", min=0, max=100, step=1, width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    tags$div(title="help text ",
                                                             numericInput(inputId = "pois.per.ha.a", label="Cost per hectare ($)", value=20, width="150px"))),
                                                div(style="display:inline-block;vertical-align:top",
                                                    verbatimTextOutput("text_pois_a_cost")),
                                                div(style="display:inline-block;vertical-align:top",checkboxInput(inputId="pois_mask", label="Aerial mask", value=FALSE)),
                                                div(style="display:inline-block;vertical-align:top",conditionalPanel(
                                                  condition="input.pois_mask==1",
                                                  div(style="display:inline-block;vertical-align:top", fileInput(inputId = "pois_asc", label="Chose the aerial mask", accept=c(".asc",".tif"), multiple=TRUE, width="200px")),
                                                  div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.pois.asc", width = "450px", height="350px"))
                                                ))
                                                
                                                
                                              )#End of wellpanel
                                            )#End of conditional panel 
                                            
                                            
                                            
                                   ),

                                   tabPanel("3. Run Scenarios",
                                            # tableOutput('scenarios.table'),
                                            h4("Use this tab to check scenarios, and delete them if needed (click on the rows to delete)"), 
                                            h4("When you are happy, click 'Run Simulations' to run the simulations."),
                                            fluidRow(
                                              column(width=3,
                                                     wellPanel(
                                                       
                                                       radioButtons(inputId="sim_type", label="Simulation Type",choices=c("Individual traps"="individ"), selected="individ", width="200px"),
                                                       div(style="display:inline-block;vertical-align:bottom",
                                                           numericInput(inputId = "n.nights",label="Simulation length (days)", value=100, width="150px")),
                                                       div(style="display:inline-block;vertical-align:bottom",
                                                           numericInput(inputId = "n.its",label="Iterations per scenario", value=5, width="150px"))
                                                       
                                                     )),
                                              # column(width=1,
                                              #        # wellPanel(
                                              #          tags$style(type="text/css", "textarea {width:100%}"),
                                              #          div(style="display:inline-block;vertical-align:bottom",
                                              #              actionButton("act.btn.trapsim",strong("Run Simulations")))
                                              #          
                                              #        # )
                                              # )
                                              
                                              column(width=2,
                                                     # fluidRow(
                                                   p(),p(),
                                                       actionButton("act.btn.trapsim",strong("Run Simulations"), icon("paper-plane"),style="background-color:green")
                                                     )
                                              
                                            ),
                                            
                                            actionButton("deleteRows", strong("Delete Selected Rows")),
                                            actionButton("deleteAllRows", strong("Delete All Rows")),br(),
                                            p(),
                                            # ),
                                            
                                            DT::dataTableOutput("tableDT")
                                   ),
                                   

                                   
                                   # ),
                                   tabPanel("4. Results",
                                            fluidRow(
                                              column(width=3,
                                                     # selectInput("result_scenario","Choose Scenario to plot",choices=NULL),
                                                     plotOutput(outputId = "plot1", width = "400px", height="350px")
                                                     
                                              ),
                                              column(width=3,
                                                     plotOutput(outputId = "plot2", width = "550px", height="350px")%>% withSpinner(type=4)  #Cumulative Captures
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
                                              
                                            ),
                                            fluidRow(
                                              downloadButton("report", "Generate report") 
                                            )
                                            
                                   ),
                                   tabPanel("5. Help",
                                            fluidRow(
                                              column(width=6,
                                              h3("Instructions"),
                                              "EradSim is based on an earlier tool called TrapSim which simulated trapping only.",
                                              "For a background report on TrapSim, ",
                                              tags$a(href="https://www.pfhb.nz/assets/Document-Library/Gormley-and-Warburton-2017-TrapSim-a-decision-support-tool.pdf", "Click here.", target="_blank"),
                                              br(),
                                              "The purpose of the tool is to simulate various control scenarios in order to see which combination of methods is most likely to achieve eradication.",
                                              "The app consist of 4 tabs:",p(),
                                              h4(span("1. Area and Pest Parameters", style="color:blue")),
                                              span( "~  Area", style="color:red"),
                                              "- First set the eradication area. (The default is Mahia Peninsula, NZ). Click ",strong('Upload Shapefile')," then ",strong('Browse')," to open a window and navigate to the folder containing either the individual shapefile components (.dbf,.prj,.shp & .shx), or a .zip file which contains these components in a single file.", br(),
                                              span( "~  Pest parameters", style="color:red"),
                                              "- Set the number (and distribution) of animals, home range size (specified as sigma), and reproductive rates (annual growth rate, start day of breeding and length of breeding period) ", br(),
                                              "- Animals will be either randomly located across the landscape, or according to habitat. If you select 'Habitat Specific', then you must upload a raster of habitat that overlaps the shapefile. Rasters can be either .asc or .tif format.",br(),
                                              "- Individuals are assumed to have a circular home-range which is defined by the parameter sigma (which represents the standard deviations of the bell-shaped bivariate-normal distribution. The 95% home range diameter is approximately 5 x Sigma", br(),
                                              p(),
                                              h4(span("2. Control Methods",style="color:blue")),
                                              "- Currently there are 4 control methods to chose from - Trapping, Baitstations, Hunting and Aerial poisoning. There are various parameters for each relating to start deployment start day/length, spacing etc",
                                              radioButtons(inputId="help_control",label="Click each for more information",choices=c("Trapping"="trap","Bait-stations"="bait","Hunting"="hunt","Aerial poison"="pois"), inline=TRUE, selected=""),
                                              " Set a control scenario by selecting one or more methods - you can model methods concurrently.",
                                              " Click 'Add Scenario' to add your control to the list of scenarios you wish to model. (These will appear in Tab 3). You can build up as many Control scenarios as you like.", 
                                              " Each scenario will be given an automatically generated name, however these can be changed ", strong("before"), " clicking the 'Add Scenario' button", br(),
                                              p(),
                                              
                                              h4(span("3. Run Scenarios ",style="color:blue")),
                                              "- When you have specified the set of control scenarios you wish to model, use this tab to check  the scenarios (you can delete one or all of them).",
                                              " Enter the number of iterations for each and the simulation length in days.",br(),
                                              "Click 'Run Simulation'", p(),
                                              h4(span("4. Results ",style="color:blue")),
                                              "- Explore the results  - graphs and tables of animals remaining, Total costs etc", br(),
                                              
                                              
                                              # "For input parameters with red labels, enter values separated by a slash, e.g. 1000/500 ",
                                              p(),
                                              strong("Hover cursor over each of the input boxes for pop-up help.")
                                              
                                              ),
                                              column(width=5,
                                              

                                              
                                              conditionalPanel(condition="input.help_control=='trap'",
                                              span( "~  Trapping", style="color:red"),br(),"Traps can be either live capture traps or kill traps (single kill or self-resetting), deployed for a set period of time,",br(),
                                              strong("-Spacing")," - traps are laid out on a grid specified by trap-spacings in the E-W and N-S directions and are clipped to fit inside the shapefile.",br(),
                                              strong("-Daily bycatch")," is the daily probability of traps catching non-target animals, but can also include the daily failure rate. This probability (between 0 - 1) is specified 'per trap'. ",br(),
                                              strong("-Max catch")," is the maximum trap capacity - generally 1, unless multi-capture or self resetting traps are used.",br(),
                                              strong("-g0")," is the nightly probability of a trap catching an individual whose home range centre is at the trap location.",
                                              "This parameter is generally estimated from spatial capture-recapture studies. Published estimates are available for some species+device combinations.",br(),
                                              strong("-Proportion untrappable")," - In some populations, there will be a proportion of the population that are untrappable due to factors such as behavioural differences or trap-shyness.", br(),
                                              strong("-Start night")," is the night of the simulation that trapping starts, and ",strong("Duration")," is the length of time traps are deployed.",br(),
                                              strong("-Check interval")," is the number of nights between checking and clearing/resetting traps. For live capture traps this must be set to 1, however for kill traps this may be longer. ",br(),
                                              strong("-Traps checked perday")," is the approximate number of traps that a field tech can check each day.",br(),
                                              strong("-Day rate ($)")," is the labour cost per person, and ",strong("Fixed cost per trap ($)")," is the cost of trap purchase (including lures etc)."
                                              ),
                                              
                                              conditionalPanel(condition="input.help_control=='bait'",
                                              span( "~  Bait-stations", style="color:red"),br(),
                                              strong("-Spacing")," - bait-stations are laid out on a grid specified by spacings in the E-W and N-S directions and are clipped to fit inside the shapefile.",
                                              strong("-Daily rate of failure")," is the daily probability of a bait-station failing (for bait-stations with a mechanical opening. This probability (between 0 - 1) is specified 'per bait-station'. ",br(),
                                              strong("-g0")," is the nightly probability of a bait-station succesfully poisoning an individual whose home range centre is at the bait location.",
                                              "This parameter is generally estimated from spatial capture-recapture studies. Published estimates are available for some species+device combinations.",br(),
                                              strong("-Proportion unbaitable")," - In some populations, there will be a proportion of the population that will not interact with bait-stations due to factors such as behavioural differences or shyness.", br(),
                                              strong("-Start night")," is the night of the simulation that baiting starts, and ",strong("Duration")," is the length of time stations are deployed.",br(),
                                              strong("-Check interval")," is the number of nights between checking and refilling bait-stations. ",br(),
                                              strong("-Bait stations checked perday")," is the approximate number of bait-stations that a field tech can check each day.",br(),
                                              strong("-Day rate ($)")," is the labour cost per person, and ",strong("Fixed cost per bait-station ($)")," is the unit cost of bait-stations."
                                              
                                              ),
                                              
                                              conditionalPanel(condition="input.help_control=='hunt'",
                                                               span( "~  Hunting", style="color:red"),br(),
                                                               "This option simulates hunting across the area (or a subset of the area). It is relatively simplified in that a distance per day is specified by the user which is converted into a probability of kill based on the kill-rate parameter. This can apply to ground or aerial hunting, however each will have a different kill rate.",br(),
                                                               strong("-Start day")," - the day of the simulation that hunting starts.",br(),
                                                               strong("-Days hunted")," - the duration of the hunting period.", br(),
                                                               strong("-Distance per day (km)")," - the liear distance per day covered by hunting in kilometers.", br(),
                                                               strong("-Kill rate")," - this parameter is combined with distance per day to calculate the probability of an individual being hunted on a single day, given by:",br(),
                                                               "p.kill = 1 - exp(killrate x Effort)",br()," where Effort = distance (km) per km\U00B2 of the area ", br(),
                                                               "Click ",strong("Effort vs prob of kill"),"to see how changing the kill rate and distance changes the probability of kill.", br(),
                                                               strong("-Hunting mask")," - You can upload a raster (.asc or .tif) that defines the area that is to be hunted. The raster must consist of grid cells where 1 corresponds to areas to be hunted and 0s otherwise.",br(),
                                                               img(src="mask_eg.png", height = 250, align="center", hspace=20,vspace=10)
                                                               
                                              ),
                                              conditionalPanel(condition="input.help_control=='pois'",
                                                               span( "~  Aerial poisoning", style="color:red"),br(),
                                                               "This option simulates aerial posioning across the area (or a subset of the area). It is relatively simplified in that a distance per day is specified by the user which is converted into a probability of kill based on the kill-rate parameter.",br(),
                                                               strong("-Start day")," - the day of the simulation that poisoning starts.",br(),
                                                               strong("-Operation length")," - the duration of the poisoning period.", br(),
                                                               strong("-Percent kill")," - this is the main parameter and is the expected percentage of the population that will be killed by the aerial control operation (spread out over the operation duration)",br(),
                                                               strong("-Aerial mask")," - as with hunting, you can upload a raster (.asc or .tif) that defines the area that is to be aerial poisoned. The raster must consist of grid cells where 1 corresponds to areas to be controlled and 0s otherwise.",br(),
                                                               img(src="mask_eg.png", height = 250, align="center", hspace=20,vspace=10)
                                                               
                                              )
                                              
                                                               
                                              # Eff<-hunt.eff.a/(mydata.shp()$ha/100)
                                              # hunt.daily.pkill<-1-exp(-(hunt.rho.a*Eff))
                                            )#End of column
                                            ))
                                   
                       )
                )
              ), #End of Row
              
              
              h6("v1.9.2: May 2022"),
              h6("For help, suggested changes, or reporting issues, contact Andrew Gormley:"),
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
#1.8.1 - included aerial poisoning - and baiting and costs
#1.9 - rasters of aerial and hunting masks can be read in. and tidy up results tab

#Need to include bait costs...?
#Group costs by labour and fixed costs...?





