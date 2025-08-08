# functions for parallel individual-level estimation in batches
# on a SLURM batch system (named "greatlakes" after 
# University of Michigan's system)


setup.greatlakes<-function (hsamples, dname, froot = "s") 
{
  dir.create(dname)
  for (i in 1:length(hsamples)) {
    onam <- paste(froot, i, sep = ".")
    fnam <- paste(dname, "/", onam, ".RData", sep = "")
    assign(onam, hsamples[[i]])
    save(list = onam, file = fnam)
  }
}


run.greatlakes.dmc<-function(
        fname,model.dir,model.file,user,account,n.add, # usually 1/3 of initial size
        nservers=1, ncpus=1, # Number os servers and cpus/server
         GB=2, # GB per job 
         wall.hours=100, # Wall time in hours
         # required by grid to access R
         module_load="module load R gsl",
         froot="s",
         verbose=TRUE,
         RUN=TRUE,
         # h.RUN.dmc parameters
         force.RUN=FALSE, # dont run if file has an auto attribute 
         max.try=20,
         minN=NA,meanN=NA,
         cut.unstuck=10,
         cut.flat.location=.5,
         cut.flat.scale=.5,
         cut.converge=1.1,split=TRUE,
         # run.dmc parameters
         farjump=NA, 
         force=FALSE,
         # common parameters
         p.migrate=.05,
         report=10,
         cleanup=TRUE,
         gamma.mult=2.38
)
{
  oname <- load(paste(fname,"RData",sep="."))
  hsamples <- get(oname)
  setup.greatlakes(hsamples,fname,froot)
  # Write sh and R run files
  cat(paste("#!/bin/bash\n#SBATCH --ntasks=",ncpus,"\n",
            "#SBATCH --mem=",GB,"gb\n",
            "#SBATCH --time=",wall.hours,":00:00\n#SBATCH --account=",account,
            "\n#SBATCH --partition=standard",
            "\n#SBATCH --job-name=",fname,
            "\n",module_load,"\n","cd $SLURM_SUBMIT_DIR\nRscript run-",fname,".R\n",sep=""),
      file=paste("run-",fname,".sh",sep=""))
  if (RUN) run.call <- paste("RUN.dmc(get(tmp), cores =",ncpus,", max.try =",max.try,
                             ", minN =",minN,", meanN =",meanN,
                             ", cut.flat.location =",cut.flat.location,
                             ", cut.flat.scale =",cut.flat.scale,", n.add =",n.add,
                             ", cut.unstuck =",cut.unstuck,", p.migrate =",p.migrate,
                             ", cut.converge =",cut.converge,", split =",split,
                             ", gamma.mult =",gamma.mult,", force =",force.RUN,
                             ", verbose =",verbose,", report =",report,")") else
                               run.call <- paste("run.dmc(get(tmp),report=",report,",cores=",ncpus,
                                                 ",p.migrate=",p.migrate,",gamma.mult=",gamma.mult,
                                                 ",farjump=",farjump,",force=",force,")")
  cat(paste("fname <- \"",fname,"\"; froot <- \"",froot,"\"\n",
            "source (\"dmc/dmc.R\")\nload_model(\"",model.dir,"\",\"",model.file,"\")\n",
            "num <- Sys.getenv(\"SLURM_ARRAY_TASK_ID\")\n",
            "onam <- paste(froot,num,sep=\".\")\n",
            "fnam <- paste(fname,\"/\",onam,\".RData\",sep=\"\")\n",
            "tmp=load(fnam)\ntmp <-",run.call,"\n",
            "assign(onam,tmp)\nfnam <- paste(fname,\"/results.\",onam,\".RData\",sep=\"\")\n",
            "save(list=onam,file=fnam)",sep=""),file=paste("run-",fname,".R",sep=""))
  jnum <- system(paste("sbatch --array=1-",length(hsamples)," run-",fname,".sh",sep=""),intern=TRUE)
  jnum <- gsub("Submitted batch job ","",jnum)
  save(user,jnum,hsamples,oname,
       froot,fname,verbose,cleanup,
       file=paste0(fname,"_RUN_info.RData"))
  
}

cleanup.greatlakes.dmc<-function(rinfo){
  
  # load run information file
  
  load(rinfo)
  
  # is run finished?
    qstat <- system(paste("squeue -u",user),intern=TRUE)
    
    if ( length(qstat[grep(jnum,qstat)])==0 ) {
 
  cat("Harvesting files: ")
  for (i in 1:length(hsamples)) {
    onam <- paste(froot,i,sep=".")
    fnam <- paste(fname,"/results.",onam,".RData",sep="")
    tmp <- try(load(fnam),silent=TRUE)
    if (class(tmp)=="try-error") bad <- TRUE else {
      hsamples[[i]] <- get(tmp)
      bad <- FALSE
    }
    if (bad) attr(hsamples[[i]],"auto") <- "GRID FAIL"  
    if (verbose) {
      if (bad) cat("Read fail\n") else 
        cat(paste(onam,":",attr(hsamples[[i]],"auto"),"\n")) 
    } else {
      if (bad) cat("F") else cat(".")
    }
  }
  cat("\n")
  auto <- unlist(lapply(hsamples,function(x){attr(x,"auto")}))
  attr(hsamples,"auto") <- auto
  assign(oname,hsamples)
  save(list=oname,file=paste(fname,"RData",sep="."))
  
  if (cleanup) { # Clean up
    unlink(fname,recursive=TRUE)
    file.remove(paste("run-",fname,".R",sep=""))
    file.remove(paste("run-",fname,".sh",sep=""))
    system(paste("rm *",jnum,"*",sep=""))
    unlink(paste(fname,".o",jnum,".*",sep=""))
  }
  
}
  
}


