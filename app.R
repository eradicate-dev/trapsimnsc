#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TrapSim Tool - Multi Method - now developed for NSC project
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Code developed by Andrew Gormley, Manaaki Whenua - Landcare Research
#Originally developed for the NSC Eco-economics project and CISS
#Updated for Margaret etc

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
# library("shinyalert")

# setwd("C:\\Users\\gormleya\\OneDrive - MWLR\\Documents\\CAEM\\IslandConservation\\TrapSimFeasibility\\Shiny")

def.shp<-"WhakatipuMahia"

# cols.tst<-brewer.pal(6,"Blues")
cols.vec<-brewer.pal(8,"Paired")
cols.eff<-brewer.pal(8,"Reds")
proj4string <- "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m"
#4326 is the EPSG for WGS84...

trap_lut<-read.csv("data//trap_lut.csv", stringsAsFactors = FALSE)


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
    
  }
  return(coords)
}


#  Make the trap locations from a x/y spacing, buffer and shapefile - and if supplied, a masking raster of 0/1
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

#Error handling for Delete buttons
modal_confirm_2 <- modalDialog(
  "Are you sure you want to continue?",
  title = "This will delete all scenarios",
  footer = tagList(
    actionButton("ok_all_scen", "Delete", class = "btn btn-danger"),
    actionButton("cancel", "Cancel")
    
  )
)
modal_confirm_1 <- modalDialog(
  "Are you sure you want to continue?",
  title = "This will delete selected scenarios",
  footer = tagList(
    actionButton("ok_sel_scen", "Delete", class = "btn btn-danger"),
    actionButton("cancel", "Cancel")
    
  )
)
#~~~~~~~~~~~~~~~~~  End functions ~~~~~~~~~~~~~~~~~~~~~~


server<-function(input, output, session) {
  iv <- InputValidator$new()
  iv$add_rule("scenname", sv_required())
  iv$add_rule("numb.poss.1", sv_required())
  iv$add_rule("numb.poss.2", sv_required())
  iv$enable()
  
  
  #~~~~~~ Set up the scenarios
  scenParam <- reactiveVal()  #Can't recall why needed..
  
  #Event to delete all scenarios when the Delete All Rows is pressed
  observeEvent(input$deleteAllRows,{
    showModal(modal_confirm_2)
  })
  observeEvent(input$ok_all_scen, {  #When OK is pressed - remove em all
    t<-scenParam()
    t <- t[-(1:dim(t)[1]),]
    scenParam(t)
    removeModal()
  })
  
  #Event to delete selected scenarios when Delete Selected Rows is pressed
  observeEvent(input$deleteRows,{
    showModal(modal_confirm_1)
  })
  observeEvent(input$ok_sel_scen, {  #When OK is pressed - remove the selected
    t<-scenParam()
    if (!is.null(input$tableDT_rows_selected)) {
      t <- t[-as.numeric(input$tableDT_rows_selected),]
      rownames(t)<-1:dim(t)[1]
    }
    scenParam(t)
    removeModal()
  })
  #Close the modal - Cancel button is same for both
  observeEvent(input$cancel, {
    removeModal()
  })
  
  
  #When the dropdown is selected, update the numericInputs
  observeEvent(input$trap_type,{
    updateNumericInput(session, "g0.mean.a", value=trap_lut$g0_poss[trap_lut$Device == input$trap_type])
    updateNumericInput(session, "max.catch.a", value=trap_lut$Max[trap_lut$Device == input$trap_type])
    updateNumericInput(session, "cost.per.trap.a", value=trap_lut$Cost[trap_lut$Device == input$trap_type])
  })
  
  
  #Code to add a scenario when the Add Scenario button is pressed - with error handling
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
    trap.name = NA
    x.space.a = NA
    y.space.a = NA
    trap.start.a = NA
    trap.nights.a = NA
    check.interval.a = NA
    g.mean.a=NA
    g.zero.a=NA    #This is the proportion that is untrappable, not g0....
    trap.mask=NA
    trap.max=NA
    trap.cost=NA
    #~~~ Bait Stations ~~~
    bait.start.a = NA
    bait.nights.a = NA
    bait.check.a = NA
    bait.x.space.a = NA
    bait.y.space.a = NA
    bait.g.mean.a = NA
    bait.g.zero.a = NA
    bait.mask = NA
    
    #Add the trap stuff
    if(input$trap_methods==1){    
      trap.name=input$trap_type
      x.space.a = input$traps.x.space.a
      y.space.a = input$traps.y.space.a
      trap.start.a = input$trap.start.a
      trap.nights.a = input$trap.nights.a
      check.interval.a = input$n.check.a
      g.mean.a=input$g0.mean.a
      # g.mean.a=trap_lut$g0_poss[trap_lut$Device == input$trap_type]
      g.zero.a=input$g0.zero.a
      
      if(input$trap_mask==1){ #Work out the trapped area and save the file path to the raster that will be used.
        trap.mask = mydata.trap()$myraster
      }
      trap.max=input$max.catch.a
      checks<-ceiling(input$trap.nights.a/input$n.check.a)+1
      n.traps.a<-dim(make.trap.locs(x.space.a, y.space.a, 100, shp))[1]
      trap.cost=trap.cost.func(a=checks, b=n.traps.a, c=input$traps.per.day.a, d=input$day.rate.a, e=input$cost.per.trap.a)
    }
    
    

    if(input$bait_methods==1){
      bait.start.a = input$bait.start.a
      bait.nights.a = input$bait.nights.a
      bait.check.a = input$bait.check.a
      bait.x.space.a = input$bait.x.space.a
      bait.y.space.a = input$bait.y.space.a
      bait.g.mean.a = input$bait.g0.mean.a
      bait.g.zero.a = input$bait.g.zero.a
      if(input$bait_mask==1){ #Work out the trapped area and save the file path to the raster that will be used.
        bait.mask = mydata.bait()$myraster
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
    
    #Full scenario to add.
    to_add <- data.frame(
      scen.name = scen.name,
      trap.name = trap.name,
      x.space.a = x.space.a,
      y.space.a = y.space.a,
      trap.start.a = trap.start.a,
      trap.nights.a = trap.nights.a,
      check.interval.a = check.interval.a,
      g.mean.a=g.mean.a,
      g.zero.a=g.zero.a,
      trap.mask=trap.mask,
      trap.max=trap.max,
      
      bait.start.a = bait.start.a,
      bait.nights.a = bait.nights.a,
      bait.check.a = bait.check.a,
      bait.x.space.a = bait.x.space.a,
      bait.y.space.a = bait.y.space.a,
      bait.g.mean.a = bait.g.mean.a,
      bait.g.zero.a = bait.g.zero.a,
      bait.mask = bait.mask,
      trap.cost=trap.cost, 
      methods = scen.code
    )
    newScenParam <- rbind(scenParam(),to_add) # adding new data
    scenParam(newScenParam) # updating data
    
    #Test for duplicates and blank scenarios
    t<-scenParam()
    if(sum(duplicated(t[,2:20]))==1){
      t<-t[!duplicated(t[,2:20]), ]        #Dont add duplicates
      showModal(modalDialog(
        title = "Error: Duplicate scenario!","This control scenario has already been added",easyClose = TRUE, fade=FALSE, size="s",
        footer =  modalButton("Cancel", icon=icon("exclamation"))
      ))
    }else if (sum(rowSums(is.na(t))==dim(t)[2]-2)==1){
      
      t<-t[!rowSums(is.na(t))==dim(t)[2]-2,]  #Dont add blank scenarios
      showModal(modalDialog(
        title = "Error: Blank scenario!","No control methods were selected",easyClose = TRUE, fade=FALSE, size="s",
        footer =  modalButton("Cancel", icon=icon("exclamation"))
      ))
    }else{#Use a showNotification rather than a popup modal
      showNotification(paste0(input$scenname," added"), type="message", duration=1)
      id<-dim((t))[1]+1
      updateTextInput(session,"scenname", "Scenario name", value=sprintf("S%03d", id))
    }
    scenParam(t)
    return(list(scenParam=scenParam))
  })
  

  #Read in the trapping mask. And output the raster as well as the path to it.
  mydata.trap<-reactive({
    if(input$trap_mask==1){
      if(is.null(input$trap_asc)==FALSE){
        myraster<-input$trap_asc$datapath   #basename for filename, dirname
        ras.trap<-raster(myraster)
      }}
    return(list(ras.trap=ras.trap, myraster=myraster))
  })
  
  #Read in the bait station mask. And output the raster as well as the path to it.
  mydata.bait<-reactive({
    if(input$bait_mask==1){
      if(is.null(input$bait_asc)==FALSE){
        myraster<-input$bait_asc$datapath   #basename for filename, dirname
        ras.bait<-raster(myraster)
      }}
    return(list(ras.bait=ras.bait, myraster=myraster))
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
      }
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
    n.poss<-input$numb.poss.1#.i
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
      #Call the function
      validate(need(input$habitat_asc !="","Upload a habitat raster"))
      # need(input$traps.x.space.a != "", "Please enter a value for the X trap spacing"),
      ras.habitat<-raster(mydata.habitat()$myraster)
      shp.buff<-gBuffer(shp, width=100)
      ras.habitat<-crop(ras.habitat, shp.buff)
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
    updateTabsetPanel(session, "inTabset", selected = '4')#"4. Results")
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
      need(input$numb.poss.1 != "", "Please enter a value for the min number of animals "),
      need(input$numb.poss.2 != "", "Please enter a value for the max number of animals ")
    )
    
    #Read in all the parameter values for the scenarios
    params<-as.data.frame(scenParam())
    #Set up some places to store results 
    n.scen<-dim(params)[1]    
    pop.size.list<-vector("list",n.scen)
    trap.catch.list<-vector("list",n.scen)
    # hunt.catch.list<-vector("list",n.scen)
    bait.catch.list<-vector("list",n.scen)
    # pois.catch.list<-vector("list",n.scen)
    # pop.size.zone.mat[[ii]][,t+1]
    
    withProgress(message="Running simulation ",value=0,{
      #Go through the scenarios one by one...
      for(kk in 1:n.scen){  #For each scenario...
        incProgress(kk/n.scen, detail = paste("Doing scenario ", kk," of", n.scen))
        
        #Pass the parameters from params to a parameter name that will be used..
        # Traps 
        trap.start.a<-params$trap.start.a[kk]
        trap.nights.a<-params$trap.nights.a[kk]
        n.check.a<-params$check.interval.a[kk]
        x.space.a<-params$x.space.a[kk]
        y.space.a<-params$y.space.a[kk]
        buffer.a<-buffer
        g0.mean.a<-params$g.mean.a[kk]
        g.zero.a<-params$g.zero.a[kk]
        trap.mask<-params$trap.mask[kk]
        trap.max<-params$trap.max[kk]
        
        #Bait
        bait.start.a<-params$bait.start.a[kk]
        bait.nights.a<-params$bait.nights.a[kk]
        bait.check.a<-params$bait.check.a[kk]
        bait.x.space.a<-params$bait.x.space.a[kk]
        bait.y.space.a<-params$bait.y.space.a[kk]
        bait.buff.a<-buffer
        bait.g0.mean.a<-params$bait.g.mean.a[kk]
        bait.g.zero.a<- params$bait.g.zero.a[kk]
        bait.mask<-params$bait.mask[kk]
        
        #How long to run the simulation for.
        n.nights<-input$n.nights
        
        # #The g0 uncertainty, and sigma values for traps and bait. The sds are not in params - but should be!! 
        # g0.sd.a<-input$g0.sd.a
        # bait.g0.sd.a<-input$bait.g0.sd.a
        g0.sd.a<-0.001
        bait.g0.sd.a<-0.001
        
        hr.ha<-input$hr.ha
        sigma.mean<-sqrt(hr.ha*10000/pi)/2.45
        sigma.sd<-sigma.mean/20
        # sigma.mean<-input$sigma.mean 
        # sigma.sd<-input$sigma.sd
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
        
        #Make the bait station locations
        if(is.na(bait.start.a)==FALSE){
          if(is.na(bait.mask)==TRUE){
            baits.a<-make.trap.locs(bait.x.space.a, bait.y.space.a, bait.buff.a, shp)
          }else{
            baits.a<-make.trap.locs(bait.x.space.a, bait.y.space.a, bait.buff.a, shp, ras=raster(bait.mask))
            # traps.a<-make.trap.locs(x.space.a, y.space.a,buffer.a,shp, ras=raster(trap.mask))            
          }
          coordinates(baits.a) <- c( "X", "Y" )
          proj4string(baits.a) <- CRS(proj4string)
          baits.xy.a<-as.data.frame(baits.a)
          n.baits.a<-dim(baits.xy.a)[1]
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
        
        if(is.na(bait.start.a)==FALSE){
          bait.period.a<-seq(from=bait.start.a, to=(bait.start.a+bait.nights.a-1), by=1)
          # #This sets the checking interval. - cleared and reset on these nights only...
          bait.check.vec.a<-seq(from=bait.start.a, to=(bait.start.a+bait.nights.a), by=bait.check.a)
          p.failure.a<-0#input$p.bycatch.a
        }
        

        #Carrying capacity for the area in terms of total number of animals - should be grid based?
        K.tot<-K.poss*ha
        
        n.poss.in.1<-input$numb.poss.1 #- should this be a parameter...? Or better to leave - maybe leave cause can re-run with same params.
        n.poss.in.2<-input$numb.poss.2 #- should this be a parameter...? Or better to leave - maybe leave cause can re-run with same params.
        if(is.na(n.poss.in.1)){
          n.poss<-0
        }
        max.catch.a<-trap.max#input$max.catch.a
        #~~~~Calculate the costs~~~~
        trap.cost.sim<-params$trap.cost[kk]
        bait.cost.sim<-0

        if(is.na(bait.start.a)==FALSE){
          checks<-ceiling(input$bait.nights.a/input$bait.check.a)+1
          # bait.cost.sim<-trap.cost.func(a=bait.check.vec.a, b=n.baits.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a)
          bait.cost.sim<-bait.cost.func(a=checks, b=n.baits.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a, f=input$bait.cost.a)        
        }
        

        n_its<-input$n.its  #Number of iterations
        pop.size.mat<-matrix(NA,nrow=n_its,ncol=n.nights+1)
        trap.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of trapping - for each iteration - how many that night/day
        bait.catch.mat<-matrix(NA,nrow=n_its,ncol=n.nights)  #To keep track of baiting - for each iteration - how many that night/day

        n.poss.vec<-floor(runif(n_its, min=n.poss.in.1, max=n.poss.in.2))
        for(ii in 1:n_its){
          n.poss<-n.poss.vec[ii]
          pop.size.mat[ii,1]<-n.poss
          
          #~~~~~~~~~Make some animals~~~~~~~~~~
          if((input$ras_hab)=="Ran"){
            #1. Random locations.
            n.poss.tmp<-(runifpoint(n.poss,shp))
            animals.xy<-as.data.frame(n.poss.tmp)
          }else{
            #2. Grid specific densities
            ras.habitat<-raster(mydata.habitat()$myraster)
            shp.buff<-gBuffer(shp, width=100)
            ras.habitat<-crop(ras.habitat, shp.buff)
            animals.xy<-get.pest.locs(ras.habitat, n.poss, shp)
          }
          colnames(animals.xy)<-c("X","Y")
          #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          
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
            animals.xy$g0.a<-0
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
          
          locshp<-get.loc.shape(sigma.mean, sigma.sd)
          animals.xy$Sigma<-rlnorm(n.animals ,meanlog=locshp$location, sdlog=locshp$shape)
          
            if(is.na(trap.start.a)==FALSE){
              dist2.xy.a<-matrix(NA,n.traps.a,n.animals)
              prob.xy.a<-matrix(0,n.traps.a,n.animals)
              dist.xy.a<-dist(as.data.frame(traps.xy.a), as.data.frame(animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
              prob.xy.a<-exp(-(dist.xy.a^2)/(2*animals.xy$Sigma^2))*animals.xy$g0.a #Use the g0u for sampling...
              rm(dist.xy.a)
            }
            #Bait stations
            if(is.na(bait.start.a)==FALSE){
              dist2.xy.c<-matrix(NA,n.baits.a,n.animals)
              prob.xy.c<-matrix(0,n.baits.a,n.animals)
              dist.xy.c<-dist(as.data.frame(baits.xy.a), as.data.frame(animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
              prob.xy.c<-exp(-(dist.xy.c^2)/(2*animals.xy$Sigma^2))*animals.xy$g0.bait #Use the g0u for sampling...
              rm(dist.xy.c)
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
          if(is.na(bait.start.a)==FALSE){
            bait.catch.a<-matrix(0,n.baits.a,n.nights)  			#Trap.catch stores the captures in each trap each night
            # bait.remain.a<-rep(max.catch.b, n.traps.b)
          }

          for (t in 1:n.nights){								#For each night
            not.caught<-(1:n.animals)[animals.xy$Dead==0]			#Animals not already caught
            
            if(is.na(trap.start.a)==FALSE){
              if(t%in%trap.period.a==TRUE){
                if(t%in%check.vec.a==TRUE){#If it is a trap clearance day...then reset the traps to T *before* trappig starts!
                  # trap.remain<-rep(T,n.traps)		
                  trap.remain.a<-rep(max.catch.a,n.traps.a)
                  # if (sim_type=='grid'){
                  #   grid.traps<-grid.traps.master
                  # }
                }
                #Turn off some of the traps according to the random probability
                # trap.remain[rbinom(n.traps,1, p.bycatch)==1]<-FALSE #Nedd to only turn off those that are on...Might be okay...
                trap.remain.a<-trap.remain.a-rbinom(n.traps.a, trap.remain.a, p.bycatch.a) #This modifcation deals with multiple capture traps 
                trap.remain.a[trap.remain.a<0]<-0
                
                if(sum(not.caught)>0){
                  for (j in not.caught){ 							#For each animal not already caught
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
                  #Call the function
                  ras.habitat<-raster(mydata.habitat()$myraster)
                  shp.buff<-gBuffer(shp, width=100)
                  ras.habitat<-crop(ras.habitat, shp.buff)
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
                
                if(is.na(bait.start.a)==FALSE){
                  new.animals.xy$g0.bait<-rbeta(N.new, alpbet.bait$alpha, alpbet.bait$beta)
                  new.animals.xy$g0.bait[new.animals.xy$g0.bait<0]<-0 #Probably not needed
                }else{
                  new.animals.xy$g0.bait<-0
                }
                new.animals.xy$Sigma<-rlnorm(N.new ,meanlog=locshp$location, sdlog=locshp$shape)        
                
                  # #Calculate the probabilities for these new ones...
                  if(is.na(trap.start.a)==FALSE){
                    dist.xy.new<-matrix(NA,n.traps.a,N.new)
                    prob.xy.new<-matrix(0,n.traps.a,N.new)
                    dist.xy.new<-dist(as.data.frame(traps.xy.a), as.data.frame(new.animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
                    prob.xy.new<-exp(-(dist.xy.new^2)/(2*new.animals.xy$Sigma^2))*new.animals.xy$g0.a #Use the g0u for sampling...
                    rm(dist.xy.new)
                    prob.xy.a<-cbind(prob.xy.a, prob.xy.new)
                  }
                  
                  if(is.na(bait.start.a)==FALSE){
                    dist2.xy.new.c<-matrix(NA,n.baits.a,N.new)
                    prob.xy.new.c<-matrix(0,n.baits.a,N.new)
                    dist.xy.new.c<-dist(as.data.frame(baits.xy.a), as.data.frame(new.animals.xy)[,1:2], method="euclidean") #Distance (not squared) - faster than using outer
                    prob.xy.new.c<-exp(-(dist.xy.new.c^2)/(2*new.animals.xy$Sigma^2))*new.animals.xy$g0.bait #Use the g0u for sampling...
                    rm(dist.xy.new.c)
                    prob.xy.c<-cbind(prob.xy.c, prob.xy.new.c)
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
            if(is.na(trap.start.a)==FALSE){
              trap.catch.mat[ii,]<-colSums(trap.catch.a)
            }else{
              trap.catch.mat[ii,]<-rep(0, n.nights)
            }
            if(is.na(bait.start.a)==FALSE){
              bait.catch.mat[ii,]<-colSums(bait.catch.a)
            }else{
              bait.catch.mat[ii,]<-rep(0, n.nights)
            }
        }#End of iteration ii
        
        params$TrapCost[kk]<-trap.cost.sim
        params$BaitCost[kk]<-bait.cost.sim
        params$TotalCost[kk]<-trap.cost.sim+bait.cost.sim
        
        params$MeanPopSize[kk]<-round(mean(pop.size.mat[,n.nights+1]),2)
        
        pop.size.list[[kk]]<-pop.size.mat
        trap.catch.list[[kk]]<-trap.catch.mat
        bait.catch.list[[kk]]<-bait.catch.mat
      } #End kk - end of scenario
    })  #End of progress
    #Re-order parameters for the table
    params<-params[,c(1,21,25,24,22:23,2:19)]
    
    return(list(trap.catch.mat=trap.catch.mat, bait.catch.mat=bait.catch.mat, pop.size.mat=pop.size.mat, animals.xy=animals.xy, params=params, pop.size.list=pop.size.list, trap.catch.list=trap.catch.list, bait.catch.list=bait.catch.list))#, pop.zone.list=pop.zone.list))#, animals.done.xy=animals.xy))    
  }
  )
  
  
  
  output$mymap<-renderLeaflet({
    shp<-mydata.shp()$shp
    shp.proj<-spTransform(shp,CRS("+proj=longlat +datum=WGS84"))
    m<-leaflet() %>%
      addTiles(group="Default")%>%
      addProviderTiles("OpenTopoMap", group = "Topo")%>%
      addProviderTiles("Esri.WorldImagery", group = "Aerial")%>%
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

    nights.vec<-1:input$n.nights
    ymax.t<-max(apply(trap.catch.mat,1,cumsum))
    ymax.b<-max(apply(bait.catch.mat,1,cumsum))
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
    mtext("Cumulative kills",3, cex=1.5, line=1)
  })
  
  output$plot4<-renderPlot({
    validate(need(input$results.table_rows_selected !="","Select a results row"))
    
    idx<-as.numeric(input$results.table_rows_selected)

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
  output$plot.bait.asc<-renderPlot({
    validate(
      need(input$bait_asc != "", "upload a tif or asc file to mask the baiting area"),
    )
    ras.bait<-raster(mydata.bait()$myraster)
    plot(ras.bait)
    plot(mydata.shp()$shp, add=TRUE)
    
  }) 
  
  
  
  output$plot.habitat.asc<-renderPlot({
    # validate(need(dim(coords)[1]>0,"Raster does not overlap with the shapefile"))
    validate(
      need(input$habitat_asc != "", "Upload a tif or asc file of the relative abundance"),
    )
    ras.habitat<-raster(mydata.habitat()$myraster)
    plot(ras.habitat)
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
    hr.ha<-input$hr.ha
    sigma.mean<-sqrt(hr.ha*10000/pi)/2.45
    return(paste0("Home range sigma: ",round(sigma.mean,1), " m"))
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
    shp<-mydata.shp()$shp 
    if(input$trap_mask==1){
      validate(
        need(input$trap_asc != "", "NA"),
      )
      trap.mask<-mydata.trap()$ras.trap
      traps.a<-make.trap.locs(input$traps.x.space.a, input$traps.y.space.a, 100, shp, ras=(trap.mask))
    }else{
      traps.a<-make.trap.locs(input$traps.x.space.a, input$traps.y.space.a, 100, shp)  
    }
    n.traps.a<-dim(traps.a)[1]
    checks<-ceiling(input$trap.nights.a/input$n.check.a)+1  #The number of checks - copes with check intervals that dont fit nealy into the duration
    cost<-trap.cost.func(a=checks, b=n.traps.a, c=input$traps.per.day.a, d=input$day.rate.a, e=input$cost.per.trap.a)
    return(paste0(n.traps.a," Traps\nCost = $", cost ))
  })
  

  output$text_bait_a_cost<-renderText({
    shp<-mydata.shp()$shp   
    if(input$bait_mask==1){
      validate(
        need(input$bait_asc != "", "NA"),
      )
      bait.mask<-mydata.bait()$ras.bait
      bait.a<-make.trap.locs(input$bait.x.space.a, input$bait.y.space.a, 100, shp, ras=(bait.mask))
    }else{
      # traps.a<-make.trap.locs(input$traps.x.space.a, input$traps.y.space.a, 100, shp)  
      bait.a<-make.trap.locs(input$bait.x.space.a, input$bait.y.space.a, 100, shp)
    }
    
    n.bait.a<-dim(bait.a)[1]
    # bait.check.vec.a<-seq(from=input$bait.start.a, to=(input$bait.start.a+input$bait.nights.a), by=input$bait.check.a)
    checks<-ceiling(input$bait.nights.a/input$bait.check.a)+1
    cost<-bait.cost.func(a=checks, b=n.bait.a, c=input$bait.per.day.a, d=input$bait.day.rate.a, e=input$cost.per.bait.a, f=input$bait.cost.a)
    
    return(paste0(n.bait.a," Ground baiting \nCost = $", cost ))
  })
  
  
  
  #The number of animals per hectare shown in Tab 1 under Pest Parameters  
  output$text_density<-renderText({
    need(input$numb.poss.1 != "", "Please enter a value for the number of animals ")
    ha<-mydata.shp()$ha
    # numb<-200
    numb<-input$numb.poss.1
    return(paste0(round(numb/ha,2)," per ha"))
  })
  
  
  output$results.table<-renderDataTable(
    datab()$params, selection='single'
  )
  
  #This contains the scenarios on Tab 3: Run Scenarios
  output$tableDT <- DT::renderDataTable(
    scenParam()[c(1,21,2:11,20,12:19)]
  )
  

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
              title="TrapSim",
              tags$head(
                tags$style(HTML("
                                .shiny-output-error-validation {
                                color: blue;
                                }
    .shiny-notification {
             position:fixed;
             top: calc(50% - 150px);
             left: calc(50% - 150px);
             width: 10%;
}
                                "))),
              fluidRow(
                column(width=4,  
                       h1("TrapSim: Eradication Feasibility"),
                       strong("How much control?"),
                       "This app is intended to provide guidance as to the approximate level of control required to achieve a desired level of pest reduction.",p(),
                       
                ),
                column(8,
                       img(src="manaaki_logo.png", height = 90, align="right", hspace=20,vspace=10),
                       img(src="bhnsc.png", height = 80, align="right", hspace=20,vspace=10)
                )
              ),
              
              
              fixedRow(),      

              fixedRow(
                column(width=12,
                       tabsetPanel(id="inTabset",
                                   tabPanel(value=1, title=strong("1. Area and Pest Parameters"), 
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
                                                         numericInput(inputId = "numb.poss.1", label="Number (min)", value=200, width="140px")),
                                                     div(style="display:inline-block;vertical-align:top",title="Number of animals",
                                                         numericInput(inputId = "numb.poss.2", label="Number (max)", value=250, width="140px")
                                                     ),
                                                     div(style="display:inline-block;vertical-align:top",
                                                         verbatimTextOutput("text_density")),
                                                     
                                                     div(style="display:inline-block;vertical-align:top",title="Pests distributed randomly or according to habitat (which requires reading in a raster of relative distribution)",
                                                         radioButtons(inputId = "ras_hab", label="Pest distribution", choices=c("Random"="Ran","Habitat Specific"="Hab"), selected="Ran",width="250px")
                                                     ),
                                                     
                                                     
                                                     div(style="display:inline-block;vertical-align:top",
                                                         conditionalPanel(condition="input.ras_hab == 'Hab'",
                                                                          div(style="display:inline-block;vertical-align:top", 
                                                                              fileInput(inputId = "habitat_asc", label="Upload a raster of habitat (.asc or .,tif)", accept=c(".tif",".asc"), multiple=FALSE, width="200px")),
                                                                          div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.habitat.asc", width = "450px", height="350px"))
                                                                          # plot.habitat.asc
                                                                          # div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.hunt.asc", width = "450px", height="350px"))
                                                         )),

                                                     
                                                     br(),
                                                     div(style="display:inline-block",title="Home range where an indivudual spends 95% of its time.",
                                                         numericInput(inputId = "hr.ha", label="Home Range (hectares)", value=20, width="120px")),
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
                                   tabPanel(value=2, title=strong("2. Control Methods"),
                                            h4("Chose one or more control methods to build a control scenario. Click Add Scenario to add it to the Run Scenarios tab."),
                                            h4("Repeat this for all control scenarios you want to run"),
                                            
                                            fluidRow(
                                              column(width=2,
                                                     wellPanel(
                                                       textInput(inputId = "scenname", label=h4("Scenario name"), value="S001"))
                                              ),
                                              
                                              column(width=2,
                                                     wellPanel(
                                                       h4("Control method(s)"),
                                                       div(style="width:300px;height:20px;vertical-align:top",
                                                           title="Ground trapping including kill traps and/or capture traps",
                                                           checkboxInput(inputId="trap_methods", label="Trapping", value=FALSE)),
                                                       div(style="width:300px;height:20px;vertical-align:top",
                                                           title="Bait locations including bait stations, or point-based hand-laid baits",
                                                           checkboxInput(inputId="bait_methods", label="Bait stations", value=FALSE))
                                                     )),
                                              column(width=1,
                                                     actionButton("update", strong("Add Scenario"), icon("layer-group"),style="background-color:green")
                                              )
                                              
                                            ),
                                            
                                            conditionalPanel(
                                              condition="input.trap_methods==1",
                                              
                                              h4("Trapping"),
                                              wellPanel(
                                                
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    selectInput(inputId="trap_type", label="Trap Type",selected="",
                                                                choices=c(trap_lut$Device), width="185px")),
                                                div(id="numeric_input_container",
                                                    tags$div(title="Maximum catch per trap  ",
                                                             numericInput(inputId = "max.catch.a", label="Maximum catch", value=1, width="150px"))),
                                                div(id="numeric_input_container",
                                                    numericInput(inputId = "g0.mean.a", label="Trap Probability (g0)", value=0.1, min=0.01, max=0.99, step=0.01, width="150px")),
                                                div(id="numeric_input_container",
                                                    numericInput(inputId = "cost.per.trap.a", label="Fixed cost per trap ($)", value=25, width="150px")),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(title="Nightly probability of by-catch, false triggers etc  ",
                                                             numericInput(inputId = "p.bycatch.a", label="Daily bycatch", value=0.01, min=0.01, max=0.99, step=0.01, width="120px"))),
                                                # div(style="display:inline-block;vertical-align:bottom",
                                                #     numericInput(inputId="g0.sd.a", label='StdDev', value=.01, width="120px")),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    numericInput(inputId = "g0.zero.a", label="Proportion untrappable", value=0.05, min=0.01, max=0.99, step=0.01, width="120px")),
                                                
                                                #Tag style for boxes that are pre-populated from the look-up-table
                                                tags$head(tags$style(HTML(
                                                  "#numeric_input_container {position: relative;display:inline-block;vertical-align:bottom;}
                                      #numeric_input_container input {background-color: lightblue; width: 100%;box-sizing: border-box;}"))),

                                                p(),
                                                
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The trap spacing in the east-west direction",
                                                             numericInput(inputId = "traps.x.space.a", label="Spacing E-W (m)", value="500",width="135px"))),
                                                div(style="display:inline-block;vertical-align:bottom",
                                                    tags$div(id="redtitle",title="The trap spacing in the north-south direction",
                                                             numericInput(inputId = "traps.y.space.a", label="Spacing N-S (m)", value="500",width="135px"))),
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
                                                    verbatimTextOutput("text_trap_a_cost")),
                                                br(),
                                                div(style="display:inline-block;vertical-align:top",checkboxInput(inputId="trap_mask", label="Trapping mask", value=FALSE)),
                                                div(style="display:inline-block;vertical-align:top",conditionalPanel(
                                                  condition="input.trap_mask==1",
                                                  div(style="display:inline-block;vertical-align:top", fileInput(inputId = "trap_asc", label="Chose the trapping mask", accept=c(".asc",".tif"), multiple=TRUE, width="200px")),
                                                  div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.trap.asc", width = "450px", height="350px"))
                                                )),
                                                
                                                p(),
                                                tags$style(type="text/css", "#redtitle {color: black}"),
                                                p(),
                                              )
                                            ),
                                            
                                            
                                            conditionalPanel(
                                              condition="input.bait_methods==1",
                                              
                                              h4("Bait Stations"),
                                              wellPanel(
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
                                                # div(style="display:inline-block;vertical-align:bottom",
                                                #     numericInput(inputId="bait.g0.sd.a", label='StdDev', value=.01, width="120px")),
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
                                                
                                                br(),
                                                div(style="display:inline-block;vertical-align:top",checkboxInput(inputId="bait_mask", label="Ground bait mask", value=FALSE)),
                                                div(style="display:inline-block;vertical-align:top",conditionalPanel(
                                                  condition="input.bait_mask==1",
                                                  div(style="display:inline-block;vertical-align:top", fileInput(inputId = "bait_asc", label="Chose the ground bait mask", accept=c(".asc",".tif"), multiple=TRUE, width="200px")),
                                                  div(style="display:inline-block;vertical-align:top", plotOutput(outputId = "plot.bait.asc", width = "450px", height="350px"))
                                                )),
                                                
                                                tags$style(type="text/css", "#redtitle {color: black}")
                                                # ),
                                              )
                                            ),#End of Well Panel
                                            
                                   ),
                                   
                                   tabPanel(value=3, title=strong("3. Run Scenarios"),
                                            # tableOutput('scenarios.table'),
                                            h4("Use this tab to check scenarios, and delete them if needed (click on the rows to delete)"), 
                                            h4("When you are happy, click 'Run Simulations' to run the simulations."),
                                            fluidRow(
                                              column(width=3,
                                                     wellPanel(
                                                       
                                                       # radioButtons(inputId="sim_type", label="Simulation Type",choices=c("Individual traps"="individ"), selected="individ", width="200px"),
                                                       div(style="display:inline-block;vertical-align:bottom",
                                                           numericInput(inputId = "n.nights",label="Simulation length (days)", value=100, width="150px")),
                                                       div(style="display:inline-block;vertical-align:bottom",
                                                           numericInput(inputId = "n.its",label="Iterations per scenario", value=5, width="150px"))
                                                       
                                                     )),
                                              
                                              column(width=2,
                                                     p(),p(),
                                                     actionButton("act.btn.trapsim",strong("Run Simulations"), icon("paper-plane"),style="background-color:green")
                                              )
                                            ),
                                            actionButton("deleteRows", strong("Delete Selected Rows"), icon("eraser"),style="background-color:indianred"),
                                            actionButton("deleteAllRows", strong("Delete All Rows"), icon("trash"),style="background-color:firebrick"),br(),
                                            p(),
                                            # ),
                                            
                                            DT::dataTableOutput("tableDT")
                                   ),
                                   
                                   # ),
                                   tabPanel(value=4, title=strong("4. Results"),
                                            fluidRow(
                                              column(width=3,
                                                     plotOutput(outputId = "plot1", width = "400px", height="350px")
                                                     
                                              ),
                                              column(width=3,
                                                     plotOutput(outputId = "plot2", width = "550px", height="350px")%>% withSpinner(type=4)  #Cumulative Captures
                                              ),
                                              column(width=3,
                                                     plotOutput(outputId = "plot4", width = "550px", height="350px") #Population Size
                                              )
                                            ),
                                            fluidRow(
                                              dataTableOutput('results.table')
                                              
                                            ),
                                            fluidRow(
                                              downloadButton("report", strong("Generate report"),icon=icon("file-export")) 
                                            )
                                   ),
                                   tabPanel(value=5, title=strong("5. Help"),
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
                                                     "- Animals will be either randomly located across the landscape, or according to habitat. If you select",span("'Habitat Specific'", style="color:red"),", then you must upload a raster of habitat that overlaps the shapefile. Rasters can be either .asc or .tif format.",br(),
                                                     "- Individuals are assumed to have a circular home-range which is defined by the parameter sigma (which represents the standard deviations of the bell-shaped bivariate-normal distribution. The 95% home range diameter is approximately 5 x Sigma", br(),
                                                     p(),
                                                     h4(span("2. Control Methods",style="color:blue")),
                                                     "- Currently there are 4 control methods to chose from - Trapping, Baitstations, Hunting and Aerial poisoning. There are various parameters for each relating to start deployment start day/length, spacing etc",
                                                     radioButtons(inputId="help_control",label="Click each for more information",choices=c("Trapping"="trap","Bait-stations"="bait","Aerial poison"="pois"), inline=TRUE, selected=""),
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
              
              
              h6("v1.9.6: November 2023"),
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
#1.9.3 - mask for trapping
#1.9.5 - inlcude ability to have distribution of individuals, and tab switching works better

#Need to include bait costs...?
#Group costs by labour and fixed costs...?


#From Al Robley
#1. Save sets of parameters...
#2. Different bait densities... - yeah nah - would need a lot of data on bait interaction rates
#3. Multiple baiting runs per year...
#4. Download data behind graphs...
#5. Need to tidy up the kill rate parameter stuff for hunting...



# ras.trap<-raster(myraster)
