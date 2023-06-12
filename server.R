########GeneBuilder####
###Amhed Missael Vargas Velazquez
###avargas0lcg@gmail.com
###Server side -  Modular version


###Load libraries
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(Biostrings)
library(Cairo)
library(stringdist)
library(shinyjs)
library(DT)
library(promises)
library(future)
#Add parallel functionality
plan(multisession) 
options(future.rng.onMisuse="ignore")

###Now load data
##Codon frequencies
CAIS=read.csv("DATA/Codon_frequencies.csv",stringsAsFactors=F,header=T)
#Format data
rownames(CAIS)= toupper(as.character(CAIS$Codon))
codons=unique(CAIS$Amino)
AAtoCodF=list()
for(i in 1:length(codons)){
  AAtoCodF=append(AAtoCodF, list(CAIS[which(CAIS$Amino == codons[i]),]))
  names(AAtoCodF)[i]=as.character(codons[i])
}

##Introns
IntronSeqs=read.table("DATA/Introns.csv",sep=",",header=T, stringsAsFactors=F)
#Format data
rownames(IntronSeqs)=IntronSeqs$Name

##piRNA sequences, please note that Pies files requires renaming so it works as before
PiesFin=read.csv("DATA/piRNA_sequence.csv",header=F,stringsAsFactors=F)
#Format data
Pies=as.character(PiesFin[,2])
PiesNA=as.character(PiesFin[,1])
PiesFin=cbind(Pies,PiesNA,as.character(PiesFin[,3]))
rownames(PiesFin)=as.character(Pies)
##New piRNA data
extraPiinfo=read.table("DATA/piRNA_abundance.csv",sep=",",header=F,stringsAsFactors=F)
rownames(extraPiinfo)=as.character(extraPiinfo[,1])
##Add extra names
extraPisite=setdiff(unique(c(rownames(extraPiinfo), as.character(PiesNA))),rownames(extraPiinfo))
extraPiinfo[extraPisite,]=cbind(extraPisite,rep("n.d.",length(extraPisite)),rep("n.d.",length(extraPisite)))

##Enzymes
enzy=read.table("DATA/Enzymes.csv", sep=",", colClasses = "character",header=T)
rownames(enzy)=as.character(enzy$Enzyme)

###Use code adapted from https://github.com/dreamRs/shinyWidgets/issues/211
choice_values <- as.character(c(enzy$Enzyme))
choice_names <- lapply(
  X = choice_values,
  FUN = function(x) {
    tags$div(
      style = "width: 160px;", x
    )
  }
)

##5primeUTRS
FivepData=read.csv("DATA/5pUTRs.csv",header=T,stringsAsFactors=F)

##3primeUTRS
ThreepData=read.csv("DATA/3pUTRs.csv",header=T,stringsAsFactors=F)

##Promoters
PromoterSeqs=read.csv("DATA/Promoters.csv",header=T,stringsAsFactors=F)
####Functions
{
  ##Sample a codon
  sampcod=function(aa,list,cai){
    newcod=sample((list[[aa]])[,6],1,prob=(list[[aa]])[,cai])
    return(toupper(as.character(newcod)))
  }
  
  ##Construct new sequence based on multiple samplings
  repcds=function(x,tabibi,list,cai){
    x=toupper(x)
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if((length(vecseq) %% 3) != 0){return(c())}
    nnseq=c()
    for(i in seq(1,length(vecseq),by=3)){
      nncod=sampcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),list,cai)
      nnseq=append(nnseq,nncod)
    }
    return(paste(nnseq,sep="",collapse=""))
  }
  
  ##Sample a codon different from input
  sampnewcod=function(aa,oldcodon,list,cai){
    oldcodon=toupper(oldcodon)
    if(nrow(list[[aa]]) < 2 ){return(oldcodon)}
    oldcodon =which(as.character(rownames(list[[aa]])) == oldcodon)
    newcod=sample((list[[aa]])[-c(oldcodon),6],1,prob=(list[[aa]])[-c(oldcodon),cai])
    return(toupper(as.character(newcod)))
  }
  
  ##Construct new sequences without repeating old codons
  repnewcds=function(x,tabibi,list,cai){
    x=toupper(x)
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if((length(vecseq) %% 3) != 0){return(c())}
    nnseq=c()
    for(i in seq(1,length(vecseq),by=3)){
      nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[i:(i+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[i:(i+2)]),sep="",collapse=""),list,cai)
      nnseq=append(nnseq,nncod)
    }
    return(paste(nnseq,sep="",collapse=""))
  }
  
  ##Modify particular positions
  modbyposiz=function(x,starts,tabibi,list,cai){
    x=toupper(x)
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if((length(vecseq) %% 3) != 0){return(c())}
    if(length(starts)>0){starts=unique(starts - (starts %% 3) + 1)}
    for(pos in c(starts)){
      nncod=sampnewcod(as.character(tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]),paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),list,cai)
      vecseq[pos:(pos+2)]=unlist(strsplit(nncod,""))
    }
    return(paste(vecseq,sep="",collapse=""))
  }
  
  ##################Related to piRNA search#################################################
  ##Pis
  countpies=function(x,y){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) != 21){return(c())}
    return((stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")))
  }
  
  Strcountpies=function(x,y){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) < 21){return(c())}
    con=0
    for(i in 1:(length(vecseq)-20)){
      con=con+countpies(paste(vecseq[i:(i+20)],collapse=""),y)
    }
    return(con)
  }
  
  countmatpies=function(x,y){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) != 21){return(c())}
    return(length(which(stringdist(paste(vecseq[1:14],collapse=""),unlist(y[paste(vecseq[15:20],collapse="")]), method="h")<6)))
  }
  
  Strcountmatpies=function(x,y){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) < 21){return(c())}
    con=0
    for(i in 1:(length(vecseq)-20)){
      con=con+countmatpies(paste(vecseq[i:(i+20)],collapse=""),y)
    }
    return(con)
  }
  
  condepies=function(x,pies,mm){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) != 20){return(c())}
    return(sum(stringdist(x,pies,method="hamming") <= mm))
  }
  
  Strcondepies=function(x,y,mm){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) < 20){return(c())}
    con=0
    for(i in 1:(length(vecseq)-19)){
      con=con+condepies(paste(vecseq[i:(i+19)],collapse=""),y,mm)
    }
    return(con)
  }
  
  findpies=function(x,pies,mm){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) != 20){return(c())}
    idx=c(stringdist(x,pies,method="hamming") <= mm)
    if(sum(idx)>0){return(pies[idx])}else{return()}
  }
  
  Strfindpies=function(x,y,mm){
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if(length(vecseq) < 20){return(c())}
    con=c()
    for(i in 1:(length(vecseq)-19)){
      con=append(con,findpies(paste(vecseq[i:(i+19)],collapse=""),y,mm))
    }
    return(con)
  }
  
  ###############################Related to HTML production and sequence viewer###################
  ###HTML Gene info
  CalculateGC = function(x){
    if(!is.character(x)){return(c())}
    x=toupper(x)
    vecseq=unlist(strsplit(x,""))
    return((countPattern("C",x)+countPattern("G",x))/length(vecseq))
  }
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  CalculateCAI = function(x,tabibi,cai,list){
    x=toupper(x)
    if(!is.character(x)){return(c())}
    vecseq=unlist(strsplit(x,""))
    if((length(vecseq) %% 3) != 0){return(c())}
    if(cai == 5 ){cai=11}
    CAIvalues=c()
    for(pos in seq(1,length(vecseq),by=3)){
      am=tabibi[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),"Amino"]
      tabs=list[[as.character(am)]]
      CAIvalues=append(CAIvalues,tabs[paste(c(vecseq[pos:(pos+2)]),sep="",collapse=""),cai]/max(tabs[,cai]))
    }
    return(gm_mean(CAIvalues))
  }
  
  HigPrintSeq= function(x,tabpatcol){
    if(!is.character(x)){return(c())}
    if(is.null(nrow(tabpatcol))){return(paste(x))}
    x=toupper(x)
    vecseq=unlist(strsplit(x,""))
    for(i in 1:nrow(tabpatcol)){
      tag=paste0("<span style=\"background-color:",tabpatcol[i,3],"\">")
      vecseq[as.integer(tabpatcol[i,1])]=paste(tag,vecseq[as.integer(tabpatcol[i,1])],sep="")
      vecseq[as.integer(tabpatcol[i,2])]=paste(vecseq[as.integer(tabpatcol[i,2])],"</span>",sep="")
    }
    return(paste(vecseq,sep="",collapse=""))
  }
  
  CDSHTMLinfo = function(x,tabibi,cai,list,tabpatcol){
    if (is.null(x)) return(NULL)
    x=toupper(x)
    if(!is.character(x)){return(c())}
    
    if(is.null(nrow(tabpatcol))){colseq=x}else{
      colseq=HigPrintSeq(x,tabpatcol)
      
      ann=c()
      paco=unique(tabpatcol[,c(3,4)])
      for(j in 1:nrow(paco)){
        ann=append(ann,paste0("<span style=\"background-color:",paco[j,1],"\">",paco[j,2],"</span><br>"))
      }
      colseq=paste(c(ann,colseq),sep="",collapse="")
    }
    
    paste0("<b>Codon Adaptation Index</b> based on codon usage selected: ", CalculateCAI(x,tabibi,cai,list), 
           "<br><b>GC content</b>: ", as.integer((CalculateGC(x))*100), "%<br>",
           "<p align=\"justify\"><tt>",
           colseq,
           "</tt></p>"
    )
  }
  
  GCHTMLinfo = function(x,tabibi,cai,list){
    if (is.null(x)) return(NULL)
    x=toupper(x)
    if(!is.character(x)){return(c())}
    
    paste0("<b>Codon Adaptation Index</b>: ", round(CalculateCAI(x,tabibi,cai,list),2), 
           "<br><b>GC content</b>: ", as.integer((CalculateGC(x))*100), "%<br>"
    )
  }
  
  #####Modify function to make it work better
  NewSequenceViewer = function(title,div_id,sequence,patterns,colors,tooltips,dflegends,starseq="",endseq=""){
    if (is.null(sequence)){return(NULL)}
    if (is.null(div_id)){return(NULL)}
    if(!is.character(sequence)){return(c())}
    if(length(patterns) != length(colors)){return(c())}
    if(length(tooltips) != length(colors)){return(c())}
    
    if(length(patterns) != length(unique(patterns))){return(c())}
    
    if(starseq==""){starseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[1:10],sep="",collapse="")}
    if(endseq==""){endseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[(nchar(sequence)-10):nchar(sequence)],sep="",collapse="")}
    
    if(length(patterns)==0){
      paste0("<div id=\"",div_id,"\"/></div>",
             "<script type=\"text/javascript\">",
             "var seq",div_id," = new Sequence(\'",
             sequence,
             "\');",
             "seq",div_id,".render(\'#",div_id,"\',{",
             "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
</script>                   
                   ")
    }else{
      ##New paradigm, just add start and endseq to patterns.. not the best
      patterns=append(patterns,starseq)
      colors=append(colors,"#a9dfbf")
      tooltips=append(tooltips,"Start")
      
      patterns=append(patterns,endseq)
      colors=append(colors,"#e6b0aa")
      tooltips=append(tooltips,"End")
      
      patitos=c()
      
      for(pat in patterns){
        patitos=append(patitos,matchPattern(DNAString(as.character(pat)),DNAString(paste(sequence,sep="",collapse="")),fixed=T))
      }
      
      patitos=patitos[order(start(patitos)),]
      
      subset=patitos[1]
      subnames=as.character(patitos[1])
      
      
      if(length(patitos)>1){
        for(n in 2:length(patitos)){
          if(length(disjoin(c(subset,patitos[n]))) != (length(subset)+1)){
            
            if(length(disjoin(c(subset,patitos[n]))) == (length(subset)+2)){
              subset=disjoin(c(subset,patitos[n]))
              subnames=c(subnames,as.character(patitos[n]), paste(as.character(patitos[n-1]),"_;_",as.character(patitos[n]),sep=""))
            }
            
            if(length(disjoin(c(subset,patitos[n]))) == (length(subset))){
              subnames[n-1]=paste(as.character(subnames[n-1]),"_;_",as.character(patitos[n]),sep="")
            }
            
          }else{
            subset=c(subset,patitos[n])
            subnames=c(subnames,as.character(patitos[n]))
            
          }     
          
        }}
      
      subcol=c("")
      subtol=c("")
      
      if(length(subnames)>1){
        for(n in 1:length(subnames)){
          if(subnames[n] %in% patterns){
            subcol[n]=colors[which(subnames[n]==patterns)]
            subtol[n]=tooltips[which(subnames[n]==patterns)]
          }else{
            subcol[n]="#bfc9ca"
            #subtol[n]=subnames[n]
            vecna=unlist(strsplit(subnames[n],"_;_"))
            compna=c()
            for(ja in vecna){
              compna=append(compna,tooltips[which(ja==patterns)])
            }
            subtol[n]=paste(compna,collapse=";")
          }
        }
      }
      
      seqcoverage="{}"
      stpos=start(subset)
      edpos=end(subset)
      if(length(stpos)>0){
        for(s in 1:length(stpos)){
          seqcoverage=paste(seqcoverage,",
                                    {start: ",stpos[s]-1,", end: ",edpos[s],", color: \"white\", bgcolor: \"",subcol[s],"\", underscore: false, tooltip: \"",subtol[s],"\"}",sep="",collapse="")
        }
      }
      
      
      newdf=data.frame(Type=c("Start","End","Multiple annotations"),Color=c("#a9dfbf","#e6b0aa","#bfc9ca"))
      
      
      if(is.data.frame(dflegends)){
        newdf=rbind(newdf,dflegends)
      }
      
      newdf=unique(newdf)
      
      LegendList=paste("{name: \"",newdf[1,1],"\", color: \"",newdf[1,2],"\", underscore: false}",sep="",collapse="")
      for(n in 2:nrow(newdf)){
        LegendList=paste(LegendList,",
                                    {name: \"",newdf[n,1],"\", color: \"",newdf[n,2],"\", underscore: false}",sep="",collapse="")
      }
      
      
      paste0("<div id=\"",div_id,"\"/></div>",
             "<script type=\"text/javascript\">",
             "var seq",div_id," = new Sequence(\'",
             sequence,
             "\');",
             "seq",div_id,".render(\'#",div_id,"\',{",
             "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
                 ",
"var Sequence",div_id,"Coverage = [",seqcoverage,"];

",

"var Sequence",div_id,"Legend = [
",LegendList,"
];

",
"seq",div_id,".coverage(Sequence",div_id,"Coverage);
",

"seq",div_id,".addLegend(Sequence",div_id,"Legend);

",

"</script>")
      
    }
  }
  
  ##Nope, revised new version for proper finding of stat and stop codons
  #####Modify function to make it work better
  NewSequenceViewerDu = function(title,div_id,sequence,patterns,colors,tooltips,dflegends,starseq="",endseq=""){
    if(is.null(sequence)){return(NULL)}
    if(is.null(div_id)){return(NULL)}
    if(!is.character(sequence)){return(c())}
    if(length(patterns) != length(colors)){return(c())}
    if(length(tooltips) != length(colors)){return(c())}
    
    if(length(patterns) != length(unique(patterns))){return(c())}
    
    if(nchar(sequence) > 25){
      if(starseq==""){starseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[1:25],sep="",collapse="")}
      if(endseq==""){endseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[(nchar(sequence)-25):nchar(sequence)],sep="",collapse="")}
    }else{
      if(starseq==""){starseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[1:10],sep="",collapse="")}
      if(endseq==""){endseq=paste(unlist(strsplit(paste(sequence,sep="",collapse=""),""))[(nchar(sequence)-10):nchar(sequence)],sep="",collapse="")}
    }
    
    ##New paradigm, just add start and endseq to patterns.. not the best so better use something similar to the original function
    ##But to do it properly, define what's going on on the chunk of code below
    
    ################################
    ##Initialize vector of matching pattern with the start sequence and keep the first match, let's test if instead I can do that match and then modify the end of the sequence
    patotes=matchPattern(DNAString(endseq),DNAString(sequence),fixed=T)
    patotes=patotes[length(patotes)]
    
    start(patotes)=end(patotes)-2
    
    ####
    patitos=matchPattern(DNAString(starseq),DNAString(sequence),fixed=T)[1]
    end(patitos)=start(patitos)+2
    
    for(pat in patterns){
      patitos=c(patitos,matchPattern(DNAString(as.character(pat)),DNAString(paste(sequence,sep="",collapse="")),fixed=T))
    }
    
    subnames=c("Start")
    subset=patitos[1]
    
    partpat=patitos[-1]
    if(length(partpat)>2){
      partpat=partpat[order(start(partpat))]
      patitos=c(patitos[1],partpat)
    }
    
    #patitos=patitos[order(start(patitos)),]
    
    if(length(patitos)>1){
      for(n in 2:length(patitos)){
        if(length(disjoin(c(subset,patitos[n]))) != (length(subset)+1)){
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset)+2)){
            subset=disjoin(c(subset,patitos[n]))
            subnames=c(subnames,as.character(patitos[n]), paste(as.character(patitos[n-1]),"_;_",as.character(patitos[n]),sep=""))
          }
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset))){
            subnames[n-1]=paste(as.character(subnames[n-1]),"_;_",as.character(patitos[n]),sep="")
          }
          
        }else{
          subset=c(subset,patitos[n])
          subnames=c(subnames,as.character(patitos[n]))
          
        }     
        
      }}
    
    #subcol=c("green")
    #subtol=c("ATG")
    
    subcol=c("#a9dfbf")
    subtol=c("ATG")
    
    if(length(subnames)>1){
      for(n in 2:length(subnames)){
        if(subnames[n] %in% patterns){
          subcol[n]=colors[which(subnames[n]==patterns)]
          subtol[n]=tooltips[which(subnames[n]==patterns)]
        }else{
          subcol[n]="#ccd1d1"
          #subtol[n]=subnames[n]
          vecna=unlist(strsplit(subnames[n],"_;_"))
          compna=c()
          for(ja in vecna){
            compna=append(compna,tooltips[which(ja==patterns)])
          }
          subtol[n]=paste(compna,collapse=";")
        }
      }
    }
    ####
    toorder=order(start(subset))
    subset=subset[toorder]
    subcol=subcol[toorder]
    subtol=subtol[toorder]
    subnames=subnames[toorder]
    ##########################
    subset=append(subset,patotes)
    subcol=append(subcol,"#e6b0aa")
    subtol=append(subtol,"Stop")
    subnames=append(subnames,"Stop codon")
    
    
    ################################
    
    #seqcoverage="{start: 1, end: 2, color: \"black\", bgcolor: \"white\", underscore: false, tooltip: \"\"
    #}"
    seqcoverage="{}"
    stpos=start(subset)
    edpos=end(subset)
    if(length(stpos)>0){
      if(stpos[1] >= 3){
        seqcoverage="{start: 1, end: 2, color: \"black\", bgcolor: \"white\", underscore: false, tooltip: \"\"
        }"
      }
      for(s in 1:length(stpos)){
        seqcoverage=paste(seqcoverage,",
                                    {start: ",stpos[s]-1,", end: ",edpos[s],", color: \"white\", bgcolor: \"",subcol[s],"\", underscore: false, tooltip: \"",subtol[s],"\"}",sep="",collapse="")
      }
    }
    
    if("#ccd1d1" %in% subcol){
      Type=c("Start","End","Multiple annotations")
      Color=c("#a9dfbf","#e6b0aa","#ccd1d1")
    }else{
      Type=c("Start","Stop")
      Color=c("#a9dfbf","#e6b0aa")  
    }
    
    newdf=data.frame(Type=Type,Color=Color)
    
    
    if(is.data.frame(dflegends)){
      newdf=rbind(newdf,dflegends)
    }
    
    newdf=unique(newdf)
    
    LegendList=paste("{name: \"",newdf[1,1],"\", color: \"",newdf[1,2],"\", underscore: false}",sep="",collapse="")
    for(n in 2:nrow(newdf)){
      LegendList=paste(LegendList,",
                                    {name: \"",newdf[n,1],"\", color: \"",newdf[n,2],"\", underscore: false}",sep="",collapse="")
    }
    
    
    paste0("<div id=\"",div_id,"\"/></div>",
           "<script type=\"text/javascript\">",
           "var seq",div_id," = new Sequence(\'",
           sequence,
           "\');",
           "seq",div_id,".render(\'#",div_id,"\',{",
           "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
                 ",
"var Sequence",div_id,"Coverage = [",seqcoverage,"];

",

"var Sequence",div_id,"Legend = [
",LegendList,"
];

",
"seq",div_id,".coverage(Sequence",div_id,"Coverage);
",

"seq",div_id,".addLegend(Sequence",div_id,"Legend);

",

"</script>")
    
  }
  
  ##For Twist only, remove  start and end
  NewSequenceViewerTW = function(title,div_id,sequence,patterns,colors,tooltips,dflegends){
    if(is.null(sequence)){return(NULL)}
    if(is.null(div_id)){return(NULL)}
    if(!is.character(sequence)){return(c())}
    if(length(patterns) != length(colors)){return(c())}
    if(length(tooltips) != length(colors)){return(c())}
    
    if(length(patterns) != length(unique(patterns))){return(c())}
    
    ##New paradigm, just add start and endseq to patterns.. not the best so better use something similar to the original function
    ##But to do it properly, define what's going on on the chunk of code below
    
    ################################
    ##Initialize vector of matching pattern with the start sequence and keep the first match, let's test if instead I can do that match and then modify the end of the sequence
    patitos=c()
    if(length(patterns)>1){
      for(pat in patterns){
        patitos=c(patitos,matchPattern(DNAString(toupper(as.character(pat))),DNAString(toupper(sequence)),fixed=T))
      }}else{patitos=matchPattern(DNAString(toupper(as.character(patterns))),DNAString(toupper(sequence)),fixed=T)}
    
    subnames=patterns[1]
    subset=patitos[1]
    
    #patitos=patitos[order(start(patitos)),]
    
    if(length(patitos)>1){
      for(n in 2:length(patitos)){
        if(length(disjoin(c(subset,patitos[n]))) != (length(subset)+1)){
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset)+2)){
            subset=disjoin(c(subset,patitos[n]))
            subnames=c(subnames,as.character(patitos[n]), paste(as.character(patitos[n-1]),"_;_",as.character(patitos[n]),sep=""))
          }
          
          if(length(disjoin(c(subset,patitos[n]))) == (length(subset))){
            subnames[n-1]=paste(as.character(subnames[n-1]),"_;_",as.character(patitos[n]),sep="")
          }
          
        }else{
          subset=c(subset,patitos[n])
          subnames=c(subnames,as.character(patitos[n]))
          
        }     
        
      }}
    
    #subcol=c("green")
    #subtol=c("ATG")
    
    #subcol=colors[1]
    #subtol=tooltips[1]
    
    subcol=c()
    subtol=c()
    
    
    for(n in 1:length(subnames)){
      if(subnames[n] %in% patterns){
        subcol[n]=colors[which(subnames[n]==patterns)]
        subtol[n]=tooltips[which(subnames[n]==patterns)]
      }else{
        subcol[n]="#ccd1d1"
        #subtol[n]=subnames[n]
        vecna=unlist(strsplit(subnames[n],"_;_"))
        compna=c()
        for(ja in vecna){
          compna=append(compna,tooltips[which(ja==patterns)])
        }
        subtol[n]=paste(compna,collapse=";")
      }
    }
    
    ####
    if(length(subset)>1){
      toorder=order(start(subset))
      subset=subset[toorder]
      subcol=subcol[toorder]
      subtol=subtol[toorder]
      subnames=subnames[toorder]
    }
    ##########################
    
    
    ################################
    
    #seqcoverage="{start: 1, end: 2, color: \"black\", bgcolor: \"white\", underscore: false, tooltip: \"\"
    #}"
    seqcoverage="{}"
    stpos=start(subset)
    edpos=end(subset)
    if(length(stpos)>0){
      if(stpos[1] >= 2){
        seqcoverage="{start: 1, end: 2, color: \"black\", bgcolor: \"white\", underscore: false, tooltip: \"\"
        }"
      }
      for(s in 1:length(stpos)){
        seqcoverage=paste(seqcoverage,",
                                    {start: ",stpos[s]-1,", end: ",edpos[s],", color: \"white\", bgcolor: \"",subcol[s],"\", underscore: false, tooltip: \"",subtol[s],"\"}",sep="",collapse="")
      }
    }
    
    if("#ccd1d1" %in% subcol){
      Type=c("Multiple annotations")
      Color=c("#ccd1d1")
    }else{
      Type=c("")
      Color=c("")  
    }
    
    newdf=data.frame(Type=Type,Color=Color)
    
    
    if(is.data.frame(dflegends)){
      newdf=rbind(newdf,dflegends)
    }
    
    newdf=unique(newdf)
    
    LegendList=paste("{name: \"",newdf[1,1],"\", color: \"",newdf[1,2],"\", underscore: false}",sep="",collapse="")
    for(n in 2:nrow(newdf)){
      LegendList=paste(LegendList,",
                                    {name: \"",newdf[n,1],"\", color: \"",newdf[n,2],"\", underscore: false}",sep="",collapse="")
    }
    
    
    paste0("<div id=\"",div_id,"\"/></div>",
           "<script type=\"text/javascript\">",
           "var seq",div_id," = new Sequence(\'",
           sequence,
           "\');",
           "seq",div_id,".render(\'#",div_id,"\',{",
           "\'showLineNumbers\': true,
  \'wrapAminoAcids\': true,
  \'charsPerLine\': 100,
  \'toolbar\': false,
  \'search\': false,
  \'title\' : \"",title,"\",
  \'sequenceMaxHeight\': \"300px\",
  \'badge\': false
});
                 ",
"var Sequence",div_id,"Coverage = [",seqcoverage,"];

",

"var Sequence",div_id,"Legend = [
",LegendList,"
];

",
"seq",div_id,".coverage(Sequence",div_id,"Coverage);
",

"seq",div_id,".addLegend(Sequence",div_id,"Legend);

",

"</script>")
    
  }
  
  ############################################Make Ape file#################################
  PasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,tabibi,cai,list,PiesList,extracomments=c(),optsecnotr=""){
    if(is.null(sequence)){return(NULL)}
    if(!is.character(sequence)){return(c())}
    #if(length(patterns) < 1 ){return(c(paste(sequence)))}
    if(length(patterns) != length(FWDcolors)){return(c())}
    if(length(REVcolors) != length(FWDcolors)){return(c())}
    if(length(tooltips) != length(FWDcolors)){return(c())}
    
    if(optsecnotr == ""){optsecnotr=sequence}
    CAIS=CalculateCAI(optsecnotr,tabibi,cai,list) 
    GCp=as.integer((CalculateGC(optsecnotr))*100)
    NoPies=length(Strfindpies(optsecnotr,PiesList,4))
    
    ##Save Lines
    FileLines=c()
    FileLines=append(FileLines,paste("LOCUS",paste(locus_name,sep="",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
    FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
    FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
    FileLines=append(FileLines,paste("VERSION",".",sep="     "))
    FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
    FileLines=append(FileLines,paste("ORGANISM","C.elegans",sep="     "))
    
    FileLines=append(FileLines,paste("COMMENT",paste(paste(paste(unlist(strsplit(locus_name,"_")),sep=" ",collapse=" ")," (Coding sequence is annoted in uppercase)"),sep="     ")))
    
    FileLines=append(FileLines,paste("COMMENT",paste("Codon Adaptation Index",as.character(round(CAIS,2))),sep="     "))
    FileLines=append(FileLines,paste("COMMENT",paste("GC",as.character(GCp),"%"),sep="     "))
    FileLines=append(FileLines,paste("COMMENT",paste("piRNA binding sites in sequence:",as.character(NoPies)),sep="     "))
    
    if(length(extracomments)>0){
      for(comocomo in extracomments){
        FileLines=append(FileLines,paste("COMMENT",paste(paste(comocomo,sep="",collapse="")),sep="     "))
      }
    }
    
    posipat=c()
    if(length(patterns) > 0){
      ##Match sequences
      for(i in 1:length(patterns)){
        stpos=start(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
        edpos=end(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
        if(length(stpos)>0){
          posipat=rbind(posipat, cbind(stpos,edpos,rep(tooltips[i],length(stpos)),rep(FWDcolors[i],length(stpos)),rep(REVcolors[i],length(stpos))))
        }
      }
      
      if(!(is.null(posipat))){
        colnames(posipat)=c("start","end","label","fwdc","revc")
      }
      
      if(!(is.null(posipat))){
        for(i in 1:length(patterns)){
          FileLines=append(FileLines,paste("COMMENT",paste(as.character(tooltips[i]),as.character(patterns[i])),sep="     "))
        }
      }
    }
    FileLines=append(FileLines,paste("COMMENT","Generated using www.wormbuilder.org/transgenebuilder/",sep="     "))
    FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
    
    
    if(!(is.null(posipat))){
      FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
      for(n in 1:nrow(posipat)){
        FileLines=append(FileLines,paste("     primer_bind     ",c(posipat[n,1]),"..",c(posipat[n,2]),"",sep=""))
        FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
        FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
        FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"",c(posipat[n,4]),"\"",sep=""))
        FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"",c(posipat[n,5]),"\"",sep=""))
        FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
      }
      
    }
    FileLines=append(FileLines,paste("ORIGIN"))
    
    Compseq=unlist(strsplit(sequence,""))
    
    partseq=c()
    
    for(i in seq(1,length(Compseq),10)){
      endseq=i+9
      if(length(Compseq)-i < 9){endseq=length(Compseq)}
      partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
      
    }
    
    i=1
    for(num in seq(1,length(Compseq),60)){
      index=as.character(num)
      spaces=paste(rep(" ",6-nchar(index)),collapse="")
      endseq=i+5
      if((length(partseq)-i) < 5){endseq=length(partseq)}
      FileLines=append(FileLines , paste(spaces,index," ",paste(partseq[i:(endseq)],collapse=" "),sep=""))
      
      i=i+6
    }
    
    FileLines=append(FileLines,paste("//"))
    
    return(FileLines)
  }
  
  ############################################Make Ape file#################################
  NewPasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,FeatType,tabibi,cai,list,PiesList,extracomments=c(),optsecnotr="",maximm=4){
    #sequence: sequence to use
    #patterns: match patterns
    #col, fwd, rev
    #tootips name as feature
    #tabibi, only to calculate CAIS, as well as cai and list
    #Pie list is to pass the pies to use
    #extra comments, things added at the end; to modify
    #FeatType is new
    
    ##Verify correct inputs
    if(is.null(sequence)){return(NULL)}
    if(!is.character(sequence)){return(c())}
    #if(length(patterns) < 1 ){return(c(paste(sequence)))}
    if(length(patterns) != length(FWDcolors)){return(c())}
    if(length(REVcolors) != length(FWDcolors)){return(c())}
    if(length(tooltips) != length(FWDcolors)){return(c())}
    if(length(FeatType) != length(FWDcolors)){return(c())}
    
    ##Only to produce extra information
    if(optsecnotr == ""){optsecnotr=sequence}
    CAIS=CalculateCAI(optsecnotr,tabibi,cai,list) 
    GCp=as.integer((CalculateGC(optsecnotr))*100)
    NoPies=length(Strfindpies(optsecnotr,PiesList,maximm))
    
    ##Save Lines
    FileLines=c()
    FileLines=append(FileLines,paste("LOCUS",paste(locus_name,sep="",collapse=""),paste(nchar(sequence),"bp ds-DNA", sep=""),"linear",paste(c(unlist(strsplit(date()," ")))[c(3,2,5)],sep="",collapse="-"),sep="     "))
    FileLines=append(FileLines,paste("DEFINITION",".",sep="     "))
    FileLines=append(FileLines,paste("ACCESSION",".",sep="     "))
    FileLines=append(FileLines,paste("VERSION",".",sep="     "))
    FileLines=append(FileLines,paste("SOURCE",".",sep="     "))
    FileLines=append(FileLines,paste("ORGANISM","C.elegans",sep="     "))
    
    FileLines=append(FileLines,paste("COMMENT",paste(paste(paste(unlist(strsplit(locus_name,"_")),sep=" ",collapse=" ")," (Coding sequence is annoted in uppercase)"),sep="     ")))
    
    FileLines=append(FileLines,paste("COMMENT",paste("Codon Adaptation Index",as.character(round(CAIS,2))),sep="     "))
    FileLines=append(FileLines,paste("COMMENT",paste("GC",as.character(GCp),"%"),sep="     "))
    FileLines=append(FileLines,paste("COMMENT",paste("piRNA binding sites in sequence:",as.character(NoPies)),sep="     "))
    
    if(length(extracomments)>0){
      for(comocomo in extracomments){
        FileLines=append(FileLines,paste("COMMENT",paste(paste(comocomo,sep="",collapse="")),sep="     "))
      }
    }
    
    posipat=c()
    if(length(patterns) > 0){
      ##Match sequences
      for(i in 1:length(patterns)){
        stpos=start(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
        edpos=end(matchPattern(DNAString(as.character(patterns[i])),DNAString(sequence),fixed=T))
        if(length(stpos)>0){
          posipat=rbind(posipat, cbind(stpos,edpos,rep(tooltips[i],length(stpos)),rep(FWDcolors[i],length(stpos)),rep(REVcolors[i],length(stpos)),rep(FeatType[i],length(stpos))))
        }
      }
      
      if(!(is.null(posipat))){
        colnames(posipat)=c("start","end","label","fwdc","revc","feat")
      }
      
      ##Remove comments of patterns
      # if(!(is.null(posipat))){
      #   for(i in 1:length(patterns)){
      #     FileLines=append(FileLines,paste("COMMENT",paste(as.character(tooltips[i]),as.character(patterns[i])),sep="     "))
      #   }
      # }
    }
    FileLines=append(FileLines,paste("COMMENT","Generated using www.wormbuilder.org/transgenebuilder/",sep="     "))
    FileLines=append(FileLines,paste("COMMENT","ApEinfo:methylated:1",sep="     "))
    
    
    if(!(is.null(posipat))){
      FileLines=append(FileLines,paste("FEATURES             Location/Qualifiers",sep=""))
      for(n in 1:nrow(posipat)){
        if(posipat[n,6] == "piRNA"){
          FileLines=append(FileLines,paste("     ",posipat[n,6],"     complement(",c(posipat[n,1]),"..",c(posipat[n,2]),")",sep=""))
        }else{
          FileLines=append(FileLines,paste("     ",posipat[n,6],"     ",c(posipat[n,1]),"..",c(posipat[n,2]),"",sep=""))
        }
        FileLines=append(FileLines,paste("                     ",paste("/locus_tag=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
        FileLines=append(FileLines,paste("                     ",paste("/ApEinfo_label=","\"",c(posipat[n,3]),"\"",sep="",collapse=""),sep="     "))
        FileLines=append(FileLines,paste("                     ","/ApEinfo_fwdcolor=\"",c(posipat[n,4]),"\"",sep=""))
        FileLines=append(FileLines,paste("                     ","/ApEinfo_revcolor=\"",c(posipat[n,5]),"\"",sep=""))
        FileLines=append(FileLines,paste("                     ","/ApEinfo_graphicformat=\"arrow_data {{0 0.5 0 1 2 0 0 -1 0 -0.5} {} 0} width 5 offset 0\"",sep=""))
      }
      
    }
    FileLines=append(FileLines,paste("ORIGIN"))
    
    Compseq=unlist(strsplit(sequence,""))
    
    partseq=c()
    
    for(i in seq(1,length(Compseq),10)){
      endseq=i+9
      if(length(Compseq)-i < 9){endseq=length(Compseq)}
      partseq=append(partseq,paste(Compseq[i:endseq],collapse=""))
      
    }
    
    i=1
    for(num in seq(1,length(Compseq),60)){
      index=as.character(num)
      spaces=paste(rep(" ",6-nchar(index)),collapse="")
      endseq=i+5
      if((length(partseq)-i) < 5){endseq=length(partseq)}
      FileLines=append(FileLines , paste(spaces,index," ",paste(partseq[i:(endseq)],collapse=" "),sep=""))
      
      i=i+6
    }
    
    FileLines=append(FileLines,paste("//"))
    
    return(FileLines)
  }
  
  ##############Twist checking
  LogicRange=function(x,treshold){
    if(!is.logical(x)){return(c())}
    flag=FALSE
    st=c()
    ed=c()
    for(idx in 1:length(x)){
      if(flag){
        if(!x[idx]){
          flag=FALSE
          ed=append(ed,idx-1)
        }
      }else{
        if(x[idx]){
          st=append(st,idx)
          flag=TRUE
        }
      }
    }
    tab=c()
    if(length(st) != length(ed)){ed=append(ed,idx)}
    if(length(st) > 0){
      tab=cbind(start=st,end=ed,dist=((ed-st)+1))
    }
    res=tab[which(tab[,"dist"] >= treshold),]
    return(res)
  }
  
  LogicHomoRange=function(x,treshold){
    if(!is.logical(x)){return(c())}
    flag=FALSE
    st=c()
    ed=c()
    for(idx in 1:length(x)){
      if(flag){
        if(!x[idx]){
          flag=FALSE
          ed=append(ed,idx)
        }
      }else{
        if(x[idx]){
          st=append(st,idx)
          flag=TRUE
        }
      }
    }
    tab=c()
    if(length(st) != length(ed)){ed=append(ed,idx)}
    if(length(st) > 0){
      tab=cbind(start=st,end=ed,dist=((ed-st)+1))
    }
    res=tab[which(tab[,"dist"] >= treshold),]
    return(res)
  }
  
  MicroHomPeWin=function(seq,window){
    if(!is.character(seq)){return(c())}
    if(!is.numeric(window)){return(c())}
    
    seqs=SplitSePeWin(seq,window)
    homos=c()
    for(i in 1:(length(seqs)-window)){
      if(seqs[i] == seqs[i+window]){
        homos=append(homos,TRUE)
      }else{
        homos=append(homos,FALSE)  
      }
    }
    return(homos)
  }
  
  SplitSePeWin = function(seq,window){
    if(!is.character(seq)){return(c())}
    if(length(seq) > 1){return(c())}
    if(nchar(seq) < window){return(c(seq))}
    if(window < 1){return(c(seq))}
    
    seq=toupper(seq)
    vecseq=unlist(strsplit(seq,""))
    
    seqs=c()
    for(i in 1:(length(vecseq) - window + 1)){
      seqs=append(seqs, paste(vecseq[i:(i+window-1)],sep="",collapse=""))
    }
    return(seqs)
  }
  
  DupSeqsWin=function(seq,window){
    seqs=SplitSePeWin(seq,window)
    revseqs=paste(reverseComplement(DNAStringSet(seqs)))
    dupidx=duplicated(c(seqs,revseqs))
    dupidx=dupidx[1:length(seqs)]
    res=c()
    if(sum(dupidx) > 0){res=unique(seqs[c(dupidx)])}
    return(res)
  }
  
  CalculateTM = function(x){
    if(!is.character(x)){return(c())}
    x=toupper(x)
    vecseq=unlist(strsplit(x,""))
    return(4*(countPattern("C",x)+countPattern("G",x)) + 2*(countPattern("A",x)+countPattern("T",x)))
  }
  
  TMWindow = function(seq, window){
    if(!is.character(seq)){return(c())}
    if(length(seq) > 1){return(c())}
    if(nchar(seq) < window){return(c())}
    if(window < 1){return(c())}
    
    seq=toupper(seq)
    vecseq=unlist(strsplit(seq,""))
    
    sts=c()
    gcs=c()
    for(i in 1:(length(vecseq) - window + 1)){
      tgc=CalculateTM(paste(vecseq[i:(i+window-1)],sep="",collapse=""))
      sts=append(sts, i)
      gcs=append(gcs, tgc)
    }
    dt=cbind(start=sts,end=(sts+window-1), tm=gcs)
    return(dt)
  }
  
  GCWindow = function(seq, window){
    if(!is.character(seq)){return(c())}
    if(length(seq) > 1){return(c())}
    if(nchar(seq) < window){return(c())}
    if(window < 1){return(c())}
    
    seq=toupper(seq)
    vecseq=unlist(strsplit(seq,""))
    
    sts=c()
    gcs=c()
    for(i in 1:(length(vecseq) - window + 1)){
      tgc=CalculateGC(paste(vecseq[i:(i+window-1)],sep="",collapse=""))
      sts=append(sts, i)
      gcs=append(gcs, tgc)
    }
    dt=cbind(start=sts,end=(sts+window-1), gc=gcs)
    return(dt)
  }
  
  FlagGCWindow = function(seq, window, loB, hiB){
    if(!is.character(seq)){return(c())}
    if(length(seq) > 1){return(c())}
    if(nchar(seq) < window){return(c())}
    if(window < 1){return(c())}
    seq=toupper(seq)
    vecseq=unlist(strsplit(seq,""))
    
    sts=c()
    eds=c()
    gcs=c()
    
    for(i in 1:(length(vecseq) - window + 1)){
      tgc=CalculateGC(paste(vecseq[i:(i+window-1)],sep="",collapse=""))
      if((tgc > hiB) | (tgc < loB)){
        sts=append(sts, i)
        eds=append(eds, (i+window-1))
        gcs=append(gcs, tgc)
      }
    }
    dt=c()
    if(length(sts)>1){dt=cbind(start=sts, end=eds, gc=gcs)}
    return(dt)
  }
  
  
  ##Main function
  TwistSynthesisWithCoordinates = function(sequence){
    message="pass"
    errors=c()
    regions=data.frame(Start=integer(),End=integer(),Description=character(), Color=character(),stringsAsFactors=FALSE)
    
    #Check if at least is a character
    if(!(is.character(sequence))){
      message = "Not a sequence" 
      errors=append(errors,message)
      #return (message)
    }
    
    ##Now split into characters for DNA sequence
    seqDNA=unlist(strsplit(toupper(sequence),""))
    
    ##Check if DNA sequence
    if((nchar(gsub("A|T|C|G","",toupper(paste(seqDNA,sep="",collapse="")))) != 0)){
      message = "Strange characters in DNA sequence"
      errors=append(errors,message)
      #return (message)
    }
    
    ##Check for errors in size
    #lower than 300 bp
    if(length(seqDNA) < 300){ ##Check for errors in size
      message = "Sequence is smaller than 300 bp"
      errors=append(errors,message)
      #return (message)
    }
    
    #larger than 500 bp
    if(length(seqDNA) > 5000){ ##Check for errors in size
      message = "Sequence is larger than 5 kb"
      errors=append(errors,message)
      #return (message)
    }
    
    ###Now do base composition analysis
    ##GC content
    GCc= CalculateGC(paste(seqDNA,sep="", collapse=""))
    
    if(GCc > .65){
      message = "GC content higher than 65%"
      errors=append(errors,message)
      #return (message)
    }
    
    if(GCc < .25){
      message = "GC content lower than 25%"
      errors=append(errors,message)
      #return (message)
    }
    
    ##Now compositional analysis per window
    #Check GC content per 50bp window
    GCt= GCWindow(paste(seqDNA,sep="", collapse=""), 50)
    difgc=max(GCt[,"gc"]) - min(GCt[,"gc"])
    
    if(difgc > .52){
      message = "GC content difference between the 50bp windows is larger than 52%"
      regions[nrow(regions)+1,]=c(GCt[which.max(GCt[,"gc"])[1],1],GCt[which.max(GCt[,"gc"])[1],2],"Highest GC 50bp window","purple")
      regions[nrow(regions)+1,]=c(GCt[which.min(GCt[,"gc"])[1],1],GCt[which.min(GCt[,"gc"])[1],2],"Lowest GC 50bp window","blue")
      errors=append(errors,message)
      #return (message)
    }
    
    #Strange double evening out event but seems to be working... what?
    ##I think part of the code of twist is to smooth the gc content on some regions and then re-run analysis
    ##I think, what Amhed in the past did was to smooth by 50 windows and the check once again complex regions
    
    avrm=c()
    for(i in 1:(nrow(GCt)-49)){avrm=append(avrm,mean(GCt[i:(i+49),3]))}
    
    ##Smooth region of GC content
    if(min(avrm) < .2){
      regions[nrow(regions)+1,]=c(GCt[which.min(avrm)[1],1],GCt[which.min(avrm)[1]+49,2],"Complex region with low levels of GC","green")
      message = "Complex region with low levels of GC detected"
      errors=append(errors,message)
      #return (message)
    }
    
    if(max(avrm) > .8){
      regions[nrow(regions)+1,]=c(GCt[which.max(avrm)[1],1],GCt[which.max(avrm)[1]+49,2],"Complex region with high levels of GC","pink")
      message = "Complex region with high levels of GC detected"
      errors=append(errors,message)
      #return (message)
    }
    
    ##Now homopolymer track analysis
    ##20bp repeated sequences
    ##First we identify duplicated sequences in the full DNA sequence, by chopping into 20 mers and see if these are present once again in the sequence
    dupseqs=DupSeqsWin(paste(seqDNA,sep="", collapse=""),20)
    
    #Cuunt how many duplicated sequences
    if(length(dupseqs) > 0){
      #Aerr=append(Aerr,paste("Duplicated 20-mer:",dupseqs))
      #ErrorFlag= ErrorFlag +1
      unidupseqs=unique(dupseqs)
      for(pat in unidupseqs){
        matches=matchPattern(DNAString(pat),DNAString(paste(seqDNA,sep="", collapse="")),fixed=T)
        regions[(nrow(regions)+1):(nrow(regions)+length(matches)),]=c(start(matches),end(matches),paste("Duplicated sequence:",as.character(matches)),rep("gold",length(matches)))
      }
      message = paste("There are at least",length(dupseqs), "20-mer sequences identical in DNA sequence:\n",paste0(dupseqs,sep="\n",collapse="\n"))
      errors=append(errors,message)
      #return (message)
    }
    
    ##Now melting temperature for those regions
    ##Calculate 20-mers with TM higher than 60
    TMt=TMWindow(paste(seqDNA,sep="", collapse=""),20)
    timd=which(TMt[,"tm"] > 80)
    
    #Check if the exist
    if(length(timd)>0){
      #Aerr=append(Aerr,paste("20bp regions with high TM (> 80C):",length(timd)))
      #ErrorFlag= ErrorFlag +1
      regions[(nrow(regions)+1):(nrow(regions)+length(timd)),]=c(TMt[timd,1],TMt[timd,2],rep("Region with high melting temperature",length(timd)),rep("purple",length(timd)))
      message = paste("20bp regions with high TM (> 80C):",length(timd))
      errors=append(errors,message)
      #return (message)
    }
    
    ##Highest microhomologies on kmer analysis
    kmers=c()
    for(k in 1:20){
      kmers=append(kmers,sum(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),k)))
    }
    
    ##Check
    if(which.max(kmers) > 10){
      #Aerr=append(Aerr,paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse=""))
      #ErrorFlag= ErrorFlag +1
      message = paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse="")
      errors=append(errors,message)
      #return (message)
    }
    
    ##Homo-polymer track larger than 10bp
    t1s=c(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),1),FALSE)
    homotra=LogicHomoRange(t1s,10)
    
    ##Homopolymer tracks, finally
    if(length(homotra)>0){
      #Aerr=append(Aerr,paste("Homo-polymer tracks:",paste(homotra,sep="-", collapse="-")))
      #ErrorFlag= ErrorFlag +1
      regions[(nrow(regions)+1):(nrow(regions)+nrow(homotra)),]=c(homotra[,1],homotra[,2],rep("Micro-homology found between adjacent regions",nrow(homotra)),rep("cyan",nrow(homotra)))
      message = paste("Homo-polymer tracks:",paste(homotra,sep="-", collapse="-"))
      errors=append(errors,message)
      #return (message)
    }
    
    return(list(message,errors,regions))	
  }
  
  ##REDO
  ##Main function
  TwistSynthesisWithCoordinatesREDO = function(sequence){
    PASS=TRUE
    results=c()
    errors=c()
    regions=data.frame(Start=integer(),End=integer(),Description=character(), Color=character(),stringsAsFactors=FALSE)
    
    #Check if at least is a character
    if(!(is.character(sequence))){
      message = "Not a sequence" 
      errors=append(errors,message)
    }
    
    ##Now split into characters for DNA sequence
    seqDNA=unlist(strsplit(toupper(sequence),""))
    
    ##Check if DNA sequence
    if((nchar(gsub("A|T|C|G","",toupper(paste(seqDNA,sep="",collapse="")))) != 0)){
      message = "Strange characters in DNA sequence"
      errors=append(errors,message)
    }
    
    ##Check for errors in size
    #lower than 300 bp
    if(length(seqDNA) < 300){ ##Check for errors in size
      message = "Sequence is smaller than 300 bp"
      errors=append(errors,message)
    }else{
      results=append(results,"Sequence length is above the minimum requirement in size (300 bp)")
    }
    
    #larger than 500 bp
    if(length(seqDNA) > 5000){ ##Check for errors in size
      message = "Sequence is larger than 5 kb"
      errors=append(errors,message)
      #return (message)
    }else{
      results=append(results,"Sequence length is below the maximum requirement in size (5 kb)")
    }
    
    ###Now do base composition analysis
    ##GC content
    GCc= CalculateGC(paste(seqDNA,sep="", collapse=""))
    
    if(GCc < .25){
      message = "GC content lower than 25%"
      errors=append(errors,message)
      #return (message)
    }else{
      results=append(results,"Sequence overall GC content is above the minimum requirement (25%)")
    }
    
    if(GCc > .65){
      message = "GC content higher than 65%"
      errors=append(errors,message)
      #return (message)
    }else{
      results=append(results,"Sequence overall GC content is below the maximum requirement (65%)")
    }
    
    
    ##Now compositional analysis per window
    #Check GC content per 50bp window
    GCt= GCWindow(paste(seqDNA,sep="", collapse=""), 50)
    difgc=max(GCt[,"gc"]) - min(GCt[,"gc"])
    
    if(difgc > .52){
      message = "GC content difference between the 50bp windows is larger than 52%"
      regions[nrow(regions)+1,]=c(GCt[which.max(GCt[,"gc"])[1],1],GCt[which.max(GCt[,"gc"])[1],2],"Highest GC 50bp window","#c39bd3")
      regions[nrow(regions)+1,]=c(GCt[which.min(GCt[,"gc"])[1],1],GCt[which.min(GCt[,"gc"])[1],2],"Lowest GC 50bp window","#7fb3d5")
      errors=append(errors,message)
    }else{
      results=append(results,"There is no larger GC content difference (> 52%) between 50 bp fragments of the sequence")
    }
    
    #Strange double evening out event but seems to be working... what?
    ##I think part of the code of twist is to smooth the gc content on some regions and then re-run analysis
    ##I think, what Amhed in the past did was to smooth by 50 windows and the check once again complex regions
    
    avrm=c()
    for(i in 1:(nrow(GCt)-49)){avrm=append(avrm,mean(GCt[i:(i+49),3]))}
    
    ##Smooth region of GC content
    if(min(avrm) < .2){
      regions[nrow(regions)+1,]=c(GCt[which.min(avrm)[1],1],GCt[which.min(avrm)[1]+49,2],"Complex region with low levels of GC","#f8c471")
      message = "Complex region with low levels of GC detected"
      errors=append(errors,message)
      #return (message)
    }else{
      results=append(results,"There is no complex regions with low content of GC (20%)")
    }
    
    if(max(avrm) > .8){
      regions[nrow(regions)+1,]=c(GCt[which.max(avrm)[1],1],GCt[which.max(avrm)[1]+49,2],"Complex region with high levels of GC","#e59866")
      message = "Complex region with high levels of GC detected"
      errors=append(errors,message)
      #return (message)
    }else{
      results=append(results,"There is no complex regions with high content of GC (80%)")
    }
    
    ##Now homopolymer track analysis
    ##20bp repeated sequences
    ##First we identify duplicated sequences in the full DNA sequence, by chopping into 20 mers and see if these are present once again in the sequence
    dupseqs=DupSeqsWin(paste(seqDNA,sep="", collapse=""),20)
    
    #COunt how many duplicated sequences
    if(length(dupseqs) > 0){
      #Aerr=append(Aerr,paste("Duplicated 20-mer:",dupseqs))
      #ErrorFlag= ErrorFlag +1
      unidupseqs=unique(dupseqs)
      for(pat in unidupseqs){
        matches=matchPattern(DNAString(pat),DNAString(paste(seqDNA,sep="", collapse="")),fixed=T)
        regions[(nrow(regions)+1):(nrow(regions)+length(matches)),]=c(start(matches),end(matches),paste("Duplicated sequence:",as.character(matches)),rep("#d98880",length(matches)))
      }
      message = paste("There are",length(dupseqs), "20-mer sequences identical in DNA sequence")
      errors=append(errors,message)
      message = paste0("<p style=\"color:black\">")
      errors=append(errors,message)
      message = paste0(paste0(dupseqs,sep=""))
      errors=append(errors,message)
      message = paste0("</p>")
      errors=append(errors,message)
      
      #return (message)
    }else{
      results=append(results,"There is no identical 20 bp long regions within the sequence")
    }
    
    ##Now melting temperature for those regions
    ##Calculate 20-mers with TM higher than 60
    TMt=TMWindow(paste(seqDNA,sep="", collapse=""),20)
    timd=which(TMt[,"tm"] > 80)
    
    #Check if the exist
    if(length(timd)>0){
      #Aerr=append(Aerr,paste("20bp regions with high TM (> 80C):",length(timd)))
      #ErrorFlag= ErrorFlag +1
      regions[(nrow(regions)+1):(nrow(regions)+length(timd)),]=c(TMt[timd,1],TMt[timd,2],rep("Region with high melting temperature",length(timd)),rep("#73c6b6",length(timd)))
      message = paste("20bp regions with high TM (> 80C):",length(timd))
      errors=append(errors,message)
      #return (message)
    }else{
      results=append(results,"There is no 20 bp long regions with high melting temperature")
    }
    
    ##Highest microhomologies on kmer analysis
    kmers=c()
    for(k in 1:20){
      kmers=append(kmers,sum(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),k)))
    }
    
    ##Check
    if(which.max(kmers) > 10){
      #Aerr=append(Aerr,paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse=""))
      #ErrorFlag= ErrorFlag +1
      message = paste("Most frequent micro-homologies seen at ",which.max(kmers), "bps",sep="", collapse="")
      errors=append(errors,message)
      #return (message)
    }else{
      results=append(results,"There is no pairing sequences that can cause hairpins larger than 10 bp")
    }
    
    ##Homo-polymer track larger than 10bp
    t1s=c(MicroHomPeWin(paste(seqDNA,sep="", collapse=""),1),FALSE)
    homotra=LogicHomoRange(t1s,10)
    
    ##Homopolymer tracks, finally
    if(length(homotra)>0){
      #Aerr=append(Aerr,paste("Homo-polymer tracks:",paste(homotra,sep="-", collapse="-")))
      #ErrorFlag= ErrorFlag +1
      regions[(nrow(regions)+1):(nrow(regions)+nrow(homotra)),]=c(homotra[,1],homotra[,2],rep("Micro-homology found between adjacent regions",nrow(homotra)),rep("#bfc9ca",nrow(homotra)))
      message = paste("Homo-polymer tracks:")
      errors=append(errors,message)
      message = paste0("<p style=\"color:black\">")
      errors=append(errors,message)
      message = paste(paste(homotra,sep="-", collapse="-"))
      errors=append(errors,message)
      #return (message)
      message = paste0("</p>")
      errors=append(errors,message)
      
    }else{
      results=append(results,"There is no consecutive repeated regions above 10 bp long")
    }
    
    if(length(errors) > 0){PASS=FALSE}
    
    return(list(PASS,errors,regions,results))	
  }
  
}

###Server
shinyServer(function(input, output, session) {
  
  ##Retrieve unique ID for the session
  session_id <- session$token
  ##Create temporary folder for unique user
  system(paste("mkdir -p DATA/users/",session_id,sep=""))
  ###On exit, force the remove of directory
  ##Attention: Interactive sessions create problems because the listening of the server stops in main directory and not in sub directory
  session$onSessionEnded(function(){
    system(paste("rm -rf DATA/users/",session_id,sep=""))
  }
  )
  
  ###Create path so life becomes easier
  UserPath = paste("DATA/users/",session_id,"/",sep="")
  
  ###########################################################################################
  #########################Functions#########################################################
  ###########################################################################################
  
  ###########################################################################################
  ##############################Download handlers############################################
  ###########################################################################################
  
 
  ###Testing with input name
  output$DownSeqOut <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(SeqNameIn, "fasta", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/SeqOpop.fasta", sep=""), file)
    }
  )
  
  output$DownBlockConstruct <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(SeqNameIn, "fasta", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/BlockConstruct.fasta", sep=""), file)
    }
  )
  
  
  output$DownOriApe <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(SeqNameIn, "gb", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/Seqog.gb", sep=""), file)
    }
  )
  
  output$DownOptiApe <- downloadHandler(
    filename <- function() {
      SeqNameIn=as.character(input$nameinput)
      
      ###Workout name
      SeqNameIn=gsub(" ","", SeqNameIn)
      if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
      
      paste(paste("Optimized",SeqNameIn,sep="_"), "gb", sep=".")
    },
    
    content <- function(file) {
      file.copy(paste("DATA/users/",session_id,"/Seqpop.gb", sep=""), file)
    }
  )
  
  
  ###########################################################################################
  #########################Events############################################################
  ###########################################################################################
  
  ###Control panels##########################################################################################
  observeEvent(input$link_to_tabpanel_sequenceadaptation, {
    newvalue <- "Sequence Adaptation"
    updateTabsetPanel(session, "panels", newvalue)
  })
  
  ##########################################################################################################
  ##############################Main Event that brings the dynamic UI at its initial state##################
  ##########################################################################################################
  ###Render UI controls given a function called init
  inituiui=function(){
    output$DynamicUserInterface <- renderUI({
      #fillPage(
      fluidRow(
        HTML("<h2><i>C. elegans</i> transgene adaptation<sup style=\"color:red;font-size:0.70em;\">Beta version</sup></h2>"),
        radioButtons("intypeinput", label = HTML("<h4>Input [<a href=\"\" onclick=\"$('#explain_seq_input_any').toggle(); return false;\">info</a>]</h4>"),
                     choices = list("DNA" = 1, 
                                    "Protein" = 2), 
                     selected = 1, inline = TRUE, width='100%'),
        textAreaInput("nameinput", label = HTML(""), value = "", resize="none", placeholder= "Sequence name (optional)", rows=1),
        conditionalPanel(condition = "input.intypeinput==1",
                         fluidRow(
                           column(8,
                                  textAreaInput("seqDNA", label = HTML(""), value = "", cols= 100, rows=5, width = "600px"),
                                  HTML("")
                           )
                         )),
        conditionalPanel(condition = "input.intypeinput==2",
                         fluidRow(
                           column(8,
                                  textAreaInput("seqPROT", label = HTML(""), value = "", cols= 100, rows=5, width = "600px"),
                                  HTML("")
                           ),
                         )),
        HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_seq_input_any\">
                        Input sequences must be coding sequences with a start (ATG or M) and a stop (TAA/TAG/TGA or *). 
             The sequence must have a minimum length of 39 bp/13 amino acids and a maximum length of 10 kb / 3333 amino acids. 
                          <br></div></p>"),
        
        ###Parameters
        ##HTML("<h4>Sequence manipulation:</h4>"),
        ###Codon algorithm as option instead of  mandatory table
        
        selectizeInput("selectCAI", label = HTML("<b>Optimization method
                                                           [<a href=\"\" onclick=\"$('#explain_codon').toggle(); return false;\">info</a>]
                                                           </b>"), 
                       choices = list("No codon optimization" = 100, "--- Codon table ---" = 1, 
                                      "Highly expressed ubiquitous genes (Serizay et al., Genome Res., 2020)" = 2, 
                                      #"Germline" = 3, 
                                      #"Neuronal" = 4, 
                                      #"Somatic" = 5, 
                                      "--- Published algorithms ---" = 6, 
                                      "High expression (Redemann et al., Nature Meth., 2011)" = 7, 
                                      "Germline optimized (Fielmich et al., eLife, 2018)" = 8
                       ), 
                       selected = 100),
        HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_codon\">
             <b>Highly expressed ubiquitous genes.</b> Codon frequencies from the 500 genes with highest ubiquitous expression identified by <a href=\"https://genome.cshlp.org/content/30/12/1752\">Serizay <i>et al.</i> (2020)</a>. Note that codon sampling is probabilistic and the optimized output sequence will not be invariant.
<br><br>
<b>Max. expression.</b> This setting corresponds to optimization with Codon Adaptation Index = 1, developed by <a href=\"https://tu-dresden.de/cmcb/biotec/forschungsgruppen/bringmann\">Henrik Bringmann's group</a> and described in <a href=\"http://www.nature.com/nmeth/journal/v8/n3/full/nmeth.1565.html\">Redemann <i>et al.</i> (2011)</a>
 (see <a href=\"https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py\"><i>C. elegans</i> Codon Adapter</a>).
<br><br>
<b>Germ line optimization (GLO).</b> This algorithm was developed by <a href=\"https://www.utdickinsonlab.org\">Dan Dickinson</a> and described in <a href=\"https://elifesciences.org/articles/38198\">Fielmich <i>et al.</i> (2018)</a>. 
<br><br>
<b>No codon optimization.</b>This option does not optimize the coding sequence but can be used to add introns, append 3' UTRs, optimize the RBS, or remove restriction sites.
             </div></p>"),
checkboxInput("checkPirna", label = HTML("<b>Minimize <i>C. elegans</i> piRNA homology
                                                              [<a href=\"\" onclick=\"$('#explain_piRNA').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),

###Add conditional input
conditionalPanel(condition = "input.checkPirna==1",
                 selectizeInput("selectPiMM", label = HTML("<b>Max. number of mismatches</b>"), 
                                choices = list("0" = 0, "1" = 1, 
                                               "2" = 2,
                                               "3" = 3,
                                               "4" = 4
                                ), 
                                selected = 4)
),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_piRNA\">
      This option annotates and minimizes sequence homology to piRNAs to reduce germline silencing. The algorithm uses Hamming distance to identify piRNA homology and removes, when possible, sequences with up to four mismatches to the 20-mer binding region of all annotated class I and class 2 endogenous piRNAs (<a href=\"https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/Endogenous_piRNA_list.txt\">list</a>). We demonstrated in <a href=\"https://doi.org/10.1038/s41592-021-01369-z\">Priyardarshini <i>et al.</i> (2022)</a> that synthetic piRNAs with three mismatches significantly reduced expression and observed small effects with four mismatches.
The output annotates piRNA abundance based on Priyadarshini <i>et al.</i> (2022) and tissue-specific expression based on <a href=\"https://pubmed.ncbi.nlm.nih.gov/23516384/\">Billi <i>et al.</i> (2013)</a>.
     <br>Please note that this algorithm is computationally demanding and removing piRNAs comes at the cost of reducing the optimality of algorithms (e.g., the CAI is slightly reduced for every piRNA homology removed). <a href=\"http://www.hcleelab.org/\">Heng-Chi Lee's laboratory</a> has developed an alternative algorithm for optimizing transgenes and removing piRNAs, described in <a href=\"https://academic.oup.com/nar/article/46/W1/W43/4979435\">Wu <i>et al.</i> (2018)</a> (see <a href=\"http://cosbi4.ee.ncku.edu.tw/pirScan/\">pirScan</a>).
                    </div></p>"),
#h1(" "),
#br(),
##Next section
checkboxInput("checkboxRibo", label = HTML("<b>Optimize ribosomal binding 
                                                                [<a href=\"\" onclick=\"$('#explain_ribo').toggle(); return false;\">info</a>]</b>"), value = FALSE),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_ribo\">
                        This option minimizes the folding energy of positions -4 to +39 to increase ribosomal binding. This option also adds the consensus start sequence (aaaaATG). Note that ribosomal binding site optimization is not compatible with the germline optimization algorithm.
                    </div></p>"),
checkboxInput("checkfouras", label = HTML("<b>Start sequence with consensus start site
                                                              [<a href=\"\" onclick=\"$('#explain_consensusstart').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_consensusstart\">
                        This option adds a consensus start sequence (aaaaATG).
                    </div></p>"),

##Here use intron information instead ##NOTE

checkboxInput("checkIntron", label = HTML("<b>Add three introns
                                                               [<a href=\"\" onclick=\"$('#explain_introns').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
conditionalPanel(condition = "input.checkIntron==1",
                 radioButtons("intropt", label = HTML(""),
                              #choices = list(
                              #"Synthetic, Golden Gate compatible (BsaI, 51 bp, 33% GC)" = 1, 
                              #  "<i>rps-0<\i> (avg. length 55 bp, 15% GC)" = 2, 
                              #  "rps-5 (avg. length 65 bp, 22% GC)" = 3,
                              # "rps-20 (avg. length 62 bp, 28% GC)" = 4,
                              #  "Canonical Fire lab introns" = 5
                              #),
                              #choiceNames= list(
                              #  HTML("<i>rps-0</i> (avg. length 55 bp, 15% GC)"),
                              #  HTML("<i>rps-5</i> (avg. length 65 bp, 22% GC)"),
                              #  HTML("<i>rps-20</i> (avg. length 62 bp, 28% GC)"),
                              #  HTML("Canonical Fire lab (avg. length 51 bp, 24% GC)")
                              #),
                              choiceNames= lapply(IntronSeqs$UI_HTML,function(x){HTML(x)}),
                              choiceValues= lapply(IntronSeqs$UI_HTML,function(x){which(IntronSeqs$UI_HTML==x)}),
                              selected = 1, width='100%'),
                 radioButtons("intdistop",label = HTML("Intron placement"),
                              choices = list("Early start" = 1, 
                                             "equidistant" = 2), 
                              selected = 1, width='100%', inline = FALSE),
                 checkboxInput("checkintframe", label = HTML("Force introns in reading frame"), value = FALSE)
),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_introns\">
                        Introns can improve transgene expression (<a href=\"https://pubmed.ncbi.nlm.nih.gov/8244003/\">Okkema <i>et al.,</i> 1993</a>). 
      Inserted introns sequences are lower-case and inserted at consensus splice sites (<a href=\"https://www.ncbi.nlm.nih.gov/books/NBK20075/\"><i>C. elegans</i> II, 2. ed</a>). 
                        Introns within the first 150 basepairs are particularly efficient at improving germline expression (<a href=\"https://www.nature.com/articles/s41467-020-19898-0\">Al Johani <i>et al.,</i> 2020</a>).
      <br><br>
      <b>Fire lab kit introns.</b> These introns were initially used in the Fire lab vector kits (<a href=\"https://www.addgene.org/kits/firelab/\">Addgene documentation</a>).
      <br><br>
      <b>Ribosomal introns.</b> These introns are short endogenous introns from ribosomal genes. The higher GC content of <i>rps-20</i> can facilitate gene synthesis.
      </div></p>"),
###Addition of UTRs
checkboxInput("checkUTRs", label = HTML("<b>Append UTRs
                                                               [<a href=\"\" onclick=\"$('#explain_UTRs').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
conditionalPanel(condition = "input.checkUTRs==1",
                 radioButtons("p5UTR", label = HTML("5' UTR"),
                              #choices = list(
                              #  "None" = 1, 
                              #  "Fire lab synthetic spliced" = 2
                              #),
                              choiceNames=lapply(FivepData$UI_HTML,function(x){HTML(x)}),
                              choiceValues= lapply(FivepData$UI_HTML,function(x){which(FivepData$UI_HTML==x)}),
                              selected = 1, width='100%', inline = FALSE),
                 radioButtons("p3UTR", label = HTML("3' UTR"),
                              #choices = list(
                              #  "None" = 1, 
                              #  "rps-1 3' UTR" = 2,
                              #  "rps-4 3' UTR" = 3,
                              #  "tbb-2 3' UTR" = 4
                              #),
                              choiceNames=lapply(ThreepData$UI_HTML,function(x){HTML(x)}),
                              choiceValues= lapply(ThreepData$UI_HTML,function(x){which(ThreepData$UI_HTML==x)}),
                              selected = 1, width='100%', inline = FALSE)
),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_UTRs\">
                        <b>5' UTR.</b> The synthetic 5' UTR can stimulate expression similar to introns within coding sequence (<a href=\"https://www.addgene.org/kits/firelab/\">see documentation for Fire lab 1995 vector kit</a>).
      <br><br>
      <b>3' UTRs.</b> The <i>tbb-2</i> 3'UTR is generally permissive for gene expression (<a href=\"https://pubmed.ncbi.nlm.nih.gov/18818082/\">Merritt & Seydoux, Curr. Bio., 2008</a>) but can be difficult to synthesize. <i>rps-1</i> and <i>rps-4</i> are short 3' UTRs from highly expressed ribosomal genes and are often easy to synthesize. Please note that we have not tested the quantitative effect of <i>rps-1</i> and <i>rps-4</i> on expression but synthetic genes with these 3' UTRs are expressed and can rescue phenotypes from single copy inserts (<i>e.g., unc-119</i>).
      </div></p>"),
###Addition of promoters
checkboxInput("checkPromoters", label = HTML("<b>Add a Promoter element
                                                               [<a href=\"\" onclick=\"$('#explain_promoters').toggle(); return false;\">info</a>]
                                                               </b>"), value = FALSE, width='100%'),
conditionalPanel(condition = "input.checkPromoters==1",
                 radioButtons("oppromo", label = HTML(""),
                              #choices = list(
                              #  "None" = 1, 
                              #  "Fire lab synthetic spliced" = 2
                              #),
                              choiceNames=lapply(PromoterSeqs$UI_HTML,function(x){HTML(x)}),
                              choiceValues= lapply(PromoterSeqs$UI_HTML,function(x){which(PromoterSeqs$UI_HTML==x)}),
                              selected = 1, width='100%', inline = FALSE)
                 
),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_promoters\">
                        We provide a selection of promoter sequences that we use currently in our lab.
      </div></p>"),
checkboxInput("checkEnzySites", label = HTML("<b>Remove restriction enzyme sites</b>"), value = FALSE, width='100%'),
conditionalPanel(condition = "input.checkEnzySites==1",
                 fluidRow(
                   #HTML("<tt>"),
                   column(6,
                          
                          tags$div(
                            style = "width: 800px;",
                            prettyCheckboxGroup("Oenzymes", "",
                                                choiceNames = choice_names, choiceValues = c(enzy$Enzyme), inline=TRUE))
                          # checkboxInput("checkMetSites", label = HTML("<b>Dam/Dcm
                          #                                     [<a href=\"\" onclick=\"$('#explain_metsites').toggle(); return false;\">info</a>]
                          #                                      </b>"), value = FALSE, width='100%'),
                          #   HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_metsites\">
                          # Disable restriction sites affected by <a href=\"https://blog.addgene.org/plasmids-101-methylation-and-restriction-enzymes\">Dam/Dcm</a> methylases. </div></p>")
                   ),
                   #HTML("</tt>")
                 )),
#h1(" "),
#br(),
##Next section
###Output manipulation
#HTML("<h4>Visual output:</h4>"),
#checkboxInput("checkAnno", label = HTML("<b>Annotate sequences
#                                                            [<a href=\"\" onclick=\"$('#explain_anno').toggle(); return false;\">info</a>]
#                                                            </b>"), value = FALSE, width='100%'),
#HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_anno\">
#                        Display the location of piRNA and restriction sites on both input and output sequences.
#                    </div></p>"),
conditionalPanel(condition = "input.intypeinput==1",
                 checkboxInput("checkAnal", label = HTML("<b>Analytical mode
                                                              [<a href=\"\" onclick=\"$('#explain_anal').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),
                 HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_anal\">
                        This mode analyzes the input sequence but does not perform any optimization.
                    </div></p>")
                 # checkboxInput("checkModOnly", label = HTML("<b>Skip codon optimization routine
                 #                                             [<a href=\"\" onclick=\"$('#explain_modonly').toggle(); return false;\">info</a>]
                 #                                             </b>"), value = FALSE, width='100%'),
                 # HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_modonly\">
                 #       Skip codon optimization but run the rest of the algorithm, e.g. RBS optimization, piRNA removal and addition of extra sequences.
                 #   </div></p>"),
),
checkboxInput("checkTwisty", label = HTML("<b>Check for gene synthesis
                                                              [<a href=\"\" onclick=\"$('#explain_Twisty').toggle(); return false;\">info</a>]
                                                              </b>"), value = FALSE, width='100%'),
HTML("<p align=\"justify\"><div class=\"explain\" style=\"display: none\" id=\"explain_Twisty\">
                        This option checks the optimized sequence against gene synthesis guidelines from <a href=\"https://www.twistbioscience.com/\">Twist Biosciences</a>. Repetitive DNA structures and extreme transitions in GC content make DNA synthesis difficult. The higher GC content introns and shorter 3' UTRs can facilitate synthesis.  
Specifically, we analyze the output sequence for size (> 300 bp, < 5 kb), overall GC content (> 24%, < 66%), large differences in GC content within sequence, repeats (> 20 bp), melting temperature (Tm > 60°C), and homopolymers (> 10 bp).
                    </div></p>"),

verbatimTextOutput("ErrorMessage"),
#br(),
actionButton("actionSeq", label = "Optimize sequence")
#uiOutput("AllResults")
      )
      #)
    })
    
  }
  
  ##uiOutput("DynamicUserInterface")
  observeEvent(session$clientData,{inituiui()})
  
  ######Render UI controls once reset button is clicked
  observeEvent(input$actionRESET, {
    clearAllResults()
    inituiui()})
  
  ##Clear all results tab
  clearAllResults = function(){
    
    #uiOutput("downloadoptseq"),
    output$downloadoptseq <- renderUI({})
    
    #htmlOutput("oldsequence"),
    output$oldsequence <- renderUI({}) 
    
    #uiOutput("downloadApeOriseq"),
    output$downloadApeOriseq <- renderUI({})
    
    #tableOutput("OriPiTab"),
    output$OriPiTab <- renderTable({})
    output$TextPitab <- renderUI({})
    
    #htmlOutput("newsequence"),
    output$newsequence <- renderUI({})
    
    #uiOutput("downloadApeOptiseq"),
    output$downloadApeOptiseq <- renderUI({})
    
    #tableOutput("OptiPiTab")
    output$OptiPiTab <- renderTable({})
    output$TextOPitab <- renderUI({})
    
    #outputbutton4later
    #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
    output$button4later <- renderUI({HTML("<b></b>")})
  }
  
  #######Javascript events
  observeEvent(input$clockactivity, {
    checkfilestatus(UserPath)
  })
  
  checkfilestatus = function(path){
    StatusMessage=""
    if(file.exists(paste(path,"JobStatus.txt",sep=""))){
      StatusMessage=readLines(paste(path,"JobStatus.txt",sep=""))[1]
    }else{
      StatusMessage="Initializing..."
    }
    session$sendCustomMessage("JobStatusMessage", StatusMessage)
  }
  
  #########################################Main generator after clicking input#######################
  #####Transgene generation#####
  observeEvent(input$actionSeq, {  
    #######Acquisition of parameters and set up initial variables
    ErrorFlag=0
    
    SeqNameIn=as.character(input$nameinput)
    CodonAl=as.integer(input$selectCAI)
    
    ##Check for strange naming, particularly <, >, comas, ., and ;
    if((length(grep(";|>|<",SeqNameIn)) != 0)){ ##Check for strange characters in name
      output$ErrorMessage <- renderText({
        paste("Error: Special characters >, <, and ; are not allowed in the sequence name")
      })
      return(NULL)
    }
    
    #Ori name
    OriSeqNameIn=as.character(input$nameinput)
    if(OriSeqNameIn ==""){OriSeqNameIn="Sequence"}
    
    ###Workout name
    SeqNameIn=gsub(" ","", SeqNameIn)
    if(SeqNameIn ==""){SeqNameIn="Input_sequence"}
    
    ##Remove piRNas?
    FlaPi=input$checkPirna
    ##Max Mismatches
    if(FlaPi){
      maxMM=as.integer(input$selectPiMM)
    }else{
      maxMM=4
    }
    
    ##Add introns?
    FlaIn=input$checkIntron
    
    ##Perform RBS optimization?
    FlaRi=input$checkboxRibo
    
    ##Remove RE sites?
    FlaEnz=input$checkEnzySites
    
    ##Perform only analysis?
    FlaAna=input$checkAnal
    
    ##Trailing Aas?
    Fla5paaa=input$checkfouras
    
    ##Do not codon optimize
    FlaModOnly = FALSE
    
    
    ###Patch in wrong way to mantain analysis always active
    #FlaAno=input$checkAnno
    ##Annotate sequences?
    FlaAno=TRUE
    
    typinput=input$intypeinput
    secprousr=gsub(" |\n","", input$seqPROT)
    secdnausr=gsub(" |\n","", input$seqDNA)
    gegegenzymes=c()
    gogogonzymes=input$Oenzymes
    Inenzy=c(as.character(gegegenzymes),as.character(gogogonzymes))
    
    ###Patches to parameters
    if(CodonAl == 100){
      FlaModOnly= TRUE
      CodonAl = 2
    }
    
    if((FlaModOnly)&(as.integer(typinput) == 2)){
      output$ErrorMessage <- renderText({
        paste("Error: An optimization algorithm has to be selected when protein sequences are selected as input")
      })
      return(NULL)
    }
    
    ##Additional flag so inputs are obtain from the beginning of function
    if(FlaIn){
      typeIn=as.integer(input$intropt)
      typedistint=as.integer(input$intdistop)
      Flaframeint=input$checkintframe
    }
    
    ##Now check UTRS
    FlaTURS=input$checkUTRs
    befUTR=""
    aftUTR=""
    ##Not the best approach as I'm adding on top but good enough patch for the moment
    nam5UTR=""
    nam3UTR=""
    if(FlaTURS){
      befUTR=FivepData[as.integer(input$p5UTR),"Sequence"]
      nam5UTR=FivepData[as.integer(input$p5UTR),"Name"]
      aftUTR=ThreepData[as.integer(input$p3UTR),"Sequence"]
      nam3UTR=ThreepData[as.integer(input$p3UTR),"Name"]
    }
    
    if(befUTR != ""){Fla5paaa=TRUE}
    
    ###Now check if promoters will be added
    FlaPro=input$checkPromoters
    Proseq=""
    Pronam=""
    if(FlaPro){
      Proseq=PromoterSeqs[as.integer(input$oppromo),"Sequence"]
      Pronam=PromoterSeqs[as.integer(input$oppromo),"Name"]
      }
    ##If promoter added, add trailing A's
    if(Proseq != ""){Fla5paaa=TRUE}
    
    ##CHeck if user wants to look for issues with gene synthesis
    ##This is not anymore a separated analysis
    FlaTwisty=input$checkTwisty
    
    output$ErrorMessage <- renderText({""})
    
    output$AllResults <- renderUI({})
    
    seqDNA=""
    
    ##Status file
    statfiletowork = paste(UserPath,"JobStatus.txt",sep="")
    #Remove if present
    if(file.exists(statfiletowork)){file.remove(statfiletowork)}
    
    
    if(as.integer(typinput) == 2){###Input protein
      if((ErrorFlag == 0) & (nchar(gsub("A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|X|Y|\\*","",toupper(secprousr))) != 0)){ ##Check for strange non aminoacid characters
        output$ErrorMessage <- renderText({
          paste("Error: Unrecognized characters in sequence:",secprousr)
        })
        ErrorFlag=1
        return(NULL)
      }else{
        seqPROTO=gsub("*","X",toupper(secprousr),fixed=TRUE)
        trDNA=paste(unlist(sapply(unlist(strsplit(toupper(seqPROTO),"")),function(x){sampcod(x,AAtoCodF,1)})),sep="",collapse="")
        seqDNA=unlist(strsplit(toupper(trDNA),""))
        ###Make sure to make analysis flag off and Twist off
        FlaAna=FALSE
        #FlaTwisty=FALSE
      }
    }else{
      seqDNA=unlist(strsplit(toupper(secdnausr),""))
    }
    
    ###FLaTwisty
    ##Change error flag so we skip the other checks
    ##NOpe
    # if(FlaTwisty&(ErrorFlag==0)){
    #   ErrorFlag=5 #5, becuase... why not
    # }
    
    ##Twist supercededs analytical mode
    # if((FlaTwisty)&(FlaAna)){
    #   output$ErrorMessage <- renderText({
    #     paste("Error: Both analysis cannot be active at the same time. Analytical mode will be deactivated")
    #   })
    #   FlaAna=FALSE
    #   updateCheckboxInput(session, "checkAnal", value = FALSE)
    #   }
    
    ##Analytical mode supercedes codon table
    if(!(FlaAna) & (ErrorFlag == 0)){
      if( (CodonAl == 1) | ( CodonAl == 6)){ ##Check for error on optimization algorithm
        output$ErrorMessage <- renderText({
          paste("Select a codon table or a published algorithm to use for adaptation.")
        })
        ErrorFlag=1
        return(NULL)
      }else{
        if(CodonAl > 6){CodonAl = CodonAl - 2}else{
          CodonAl = CodonAl - 1
        }
      }
    }
    
    
    if((ErrorFlag == 0) & ((length(seqDNA) %% 3) != 0)){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is not multiple of three")
      })
      ErrorFlag=1
      return(NULL)
    }
    
    if((ErrorFlag == 0) & ((length(seqDNA) < 50) & (FlaRi))){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is too short to optimize for ribosome binding")
      })
      ErrorFlag=1
      return(NULL)
    }
    
    if((ErrorFlag == 0) & (length(seqDNA) <= 39)){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence length should be at least of 39 bp (13 aa)")
      })
      ErrorFlag=1
      return(NULL)
    }
    
    if((ErrorFlag == 0) & (length(seqDNA) > 10000)){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is larger than 10 kb.")
      })
      ErrorFlag=1
      return(NULL)
    }
    
    if((ErrorFlag == 0) & ((length(seqDNA) < 450) & (FlaIn))){ ##Check for errors in size
      output$ErrorMessage <- renderText({
        paste("Error: Sequence is too short to add ~150bp spaced introns")
      })
      #updateCheckboxInput(session, "checkIntron", value = FALSE)
      ErrorFlag=1
      return(NULL)
    }
    
    if((ErrorFlag == 0) & (nchar(gsub("A|T|C|G","",toupper(paste(seqDNA,sep="",collapse="")))) != 0)){ ##Check for strange non ATCG characters
      output$ErrorMessage <- renderText({
        paste("Error: Unrecognized characters or multiple newlines found in:",paste(seqDNA,sep="",collapse=""))
      })
      ErrorFlag=1
      return(NULL)
    }
    
    if((ErrorFlag == 0) & (paste(seqDNA[1:3],sep="",collapse="") != "ATG")){ ##Check for errors in size
      if(as.integer(typinput) == 2){
        output$ErrorMessage <- renderText({
          paste("Error: Sequence does not start with Methionine")
        })
      }else{
        output$ErrorMessage <- renderText({
          paste("Error: Sequence does not start with ATG")
        })
      }
      
      ErrorFlag=1
      return(NULL)
    } 
    
    # if((ErrorFlag == 0) & ((FlaPi)&(CodonAl == 5))){ ##Check if CAI = 1 is required
    #   output$ErrorMessage <- renderText({
    #     paste("Error: Sequence cannot remove piRNAs and mantain a Codon Adaptation Index equal to 1. piRNAs will not be removed")
    #   })
    #   FlaPi=FALSE
    #   updateCheckboxInput(session, "checkPirna", value = FALSE)
    #   return(NULL)
    # } 
    
    if((ErrorFlag == 0) & ((FlaRi)&(CodonAl == 6))){ ##Check if GLO is required
      output$ErrorMessage <- renderText({
        paste("Error: Sequence cannot be optimized for ribosome binding and Germline expression at the same time. GLO algoritmh will be run on full sequence")
      })
      FlaRi=FALSE
      #updateCheckboxInput(session, "checkboxRibo", value = FALSE)
      return(NULL)
    }
    
    if(ErrorFlag == 0){##Check stop codon
      seqiqi=paste(seqDNA,sep="",collapse="")
      stpos=c()
      stpos=start(matchPattern("*",translate(DNAString(seqiqi))))
      if(length(stpos) == 0){
        output$ErrorMessage <- renderText({
          paste("Error: Sequence does not have stop codon")
        })
        ErrorFlag=1
      }
      if(length(stpos) > 1){
        output$ErrorMessage <- renderText({
          paste("Error: Sequence have multiple stop codons")
        })
        ErrorFlag=1
      }
      if(length(stpos) == 1){
        if( (stpos) != (length(translate(DNAString(seqiqi)))) ){
          output$ErrorMessage <- renderText({
            paste("Error: Sequence has a stop codon but not at its end.")
          })
          ErrorFlag=1
        }
      }
    }
    
    ####Passed parameters
    
    eeexxttpar=c()
    if(as.integer(typinput) == 1){eeexxttpar=append(eeexxttpar,paste("Input sequence: ",secdnausr,sep="",collapse=""))}
    if(as.integer(typinput) == 2){eeexxttpar=append(eeexxttpar,paste("Input sequence: ",secprousr,sep="",collapse=""))}
    if(FlaModOnly | FlaAna){eeexxttpar=append(eeexxttpar,"Skip optimization routine: Yes")}else{
      if(CodonAl == 1){eeexxttpar=append(eeexxttpar,"Optimization routine: Ubiquitous (stochastic, frequently used codons)")}
      if(CodonAl == 5){eeexxttpar=append(eeexxttpar,"Optimization routine: Redemann et al. (2011)")}
      if(CodonAl == 6){eeexxttpar=append(eeexxttpar,"Optimization routine: Fielmich et al. (2018)")}
    }
    
    
    if(FlaRi){eeexxttpar=append(eeexxttpar,"RBS optimization: Yes")}else{eeexxttpar=append(eeexxttpar,"RBS optimization: No")}
    if(FlaPi){
      eeexxttpar=append(eeexxttpar,"Remove piRNA sites: Yes")
      eeexxttpar=append(eeexxttpar,paste("Maxmimum number of mismatches to 20-mer piRNA binding region:", maxMM,sep=" "))
    }else{eeexxttpar=append(eeexxttpar,"Remove piRNA sites: No")}
    if(FlaEnz){
      eeexxttpar=append(eeexxttpar,"Remove RE sites: Yes")
      eeexxttpar=append(eeexxttpar,"Restriction sites removed:")
      for(totoro in Inenzy){
        eeexxttpar=append(eeexxttpar,paste(totoro, " ",enzy[totoro,"Site"],"",sep=""))
      }
    }else{eeexxttpar=append(eeexxttpar,"Remove RE sites: No")}
    if(FlaIn){
      eeexxttpar=append(eeexxttpar,"Add introns: Yes")
      #for(totoro in 1:3){
      #  eeexxttpar=append(eeexxttpar,paste(paste(rownames(IntronSeqs[typeIn,])," intron",sep=""), " ",as.character(IntronSeqs[typeIn,totoro]),"",sep=""))
      #}
      eeexxttpar=append(eeexxttpar,paste(paste(rownames(IntronSeqs[typeIn,])," intron",sep=""), " ",IntronSeqs[typeIn,"Sequence_1"],",",IntronSeqs[typeIn,"Sequence_2"],",",IntronSeqs[typeIn,"Sequence_3"],sep=""))
    }else{eeexxttpar=append(eeexxttpar,"Add introns: No")}
    #UTRS
    if(FlaTURS){
      eeexxttpar=append(eeexxttpar,"Add UTRs: Yes")
      eeexxttpar=append(eeexxttpar,paste("5' UTR - ",nam5UTR," ",befUTR,sep=""))
      eeexxttpar=append(eeexxttpar,paste("3' UTR - ",nam3UTR," ",aftUTR,sep=""))
    }else{eeexxttpar=append(eeexxttpar,"Add UTRs: No")}
    #Promoter
    if(FlaPro){
      eeexxttpar=append(eeexxttpar,"Add promoter: Yes")
      eeexxttpar=append(eeexxttpar,paste("Promoter - ",Pronam," ",Proseq,sep=""))
    }else{eeexxttpar=append(eeexxttpar,"Add promoter: No")}
    
    if(Fla5paaa){
      eeexxttpar=append(eeexxttpar,"Add consensus start site: Yes")
    }else{eeexxttpar=append(eeexxttpar,"consensus start site: No")}
    
    ###Also, sequences to workon
    oriseqstart=paste(seqDNA[1:30],sep="",collapse="")
    oriseqend=paste(seqDNA[(length(seqDNA)-29):(length(seqDNA))],sep="",collapse="")
    
    ##Analytical mode supercedes main routine
    ##But in case that analytical mode and piRNAi is off the former routine will be shot down just after
    if((ErrorFlag == 0)&(FlaAna)){
      ##Load all results in dynamic ui
      output$DynamicUserInterface <- renderUI({uiOutput("AllResults")})
      
      output$AllResults <- renderUI({
        fluidRow(
          #actionButton("actionRESET", label = "RESET"),
          uiOutput("button4later"),
          htmlOutput("oldsequence"),
          br(),
          tableOutput("OriPiTab"),
          htmlOutput("TextPitab"),
          br(),
          uiOutput("downloadApeOriseq")
          #dataTableOutput("OriPiTab")
        )
      })
      
      #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
      output$button4later <- renderUI({HTML("<b></b>")})
      
      if(as.integer(typinput) == 1){ ###If original sequence
        output$oldsequence <-renderUI({
          seqiqi=toupper(seqiqi)
          
          piss=Strfindpies(seqiqi,Pies,maxMM)
          
          ##New approach
          paterns=c()
          seqpaterns=c()
          colpaterns=c()
          typpaterns=c()
          dftyp=c()
          dfcol=c()
          Neweeexxttpar=eeexxttpar
          #Ezymes
          if(length(Inenzy)>0){
            paterns=append(paterns,Inenzy)
            seqpaterns=append(seqpaterns,enzy[Inenzy,"Site"])
            colpaterns=append(colpaterns,rep("#d7bde2",length(Inenzy)))
            typpaterns=append(typpaterns,rep("RE_site",length(Inenzy)))
            dftyp=append(dftyp,"RE site")
            dfcol=append(dfcol, "#d7bde2")
          }
          
          piss=Strfindpies(seqiqi,Pies,maxMM)
          if(length(piss)>0){
            popos=c()
            papas=c()
            for(pipi in piss){
              patotes=as.character(matchPattern(DNAString(pipi),DNAString(seqiqi),max.mismatch=maxMM,fixed=T))
              popos=append(popos,patotes)
              papas=append(papas,rep(PiesFin[piss,2],length(patotes)))
            }
            
            papos=c()
            pospos=unique(popos)
            for(pipi in pospos){
              papos=append(papos,paste(papas[which(pipi == popos)], collapse=";"))  
            }
            
            paterns=append(paterns,papos)
            seqpaterns=append(seqpaterns,pospos) 
            colpaterns=append(colpaterns,rep("#f9e79f",length(pospos)))
            typpaterns=append(typpaterns,rep("piRNA",length(pospos)))
            
            for(totoro in 1:length(pospos)){
              Neweeexxttpar=append(Neweeexxttpar,paste(papos[totoro]," ",pospos[totoro],sep=""))
            }
            
            dftyp=append(dftyp,"piRNA homology")
            dfcol=append(dfcol, "#f9e79f")
          }
          #PasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,tabibi,cai,list,PiesList,extracomments=c(),optsecnotr=""){
          #NewPasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,FeatType,tabibi,cai,list,PiesList,extracomments=c(),optsecnotr=""){
          #write(paste(PasteApe(OriSeqNameIn,seqiqi,c(seqpaterns),c(colpaterns),c(colpaterns),c(paterns),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
          write(paste(NewPasteApe(OriSeqNameIn,seqiqi,c(seqpaterns),c(colpaterns),c(colpaterns),c(paterns),c(typpaterns),CAIS,5,AAtoCodF,Pies,extracomments=Neweeexxttpar,maximm=maxMM),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
          
          
          legdf=data.frame(Type=c(dftyp),Color=c(dfcol))
          
          HTML("<br><b><h3>Input ",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewerDu("","oldestseq",seqiqi,c(seqpaterns),c(colpaterns),c(paterns),legdf,oriseqstart,oriseqend),"<br>")
          
          
          
        })
        
        ##ApeButtonDownOriginal
        output$downloadApeOriseq <- renderUI({
          downloadButton('DownOriApe', 'Genbank with annotations')
        })
        
        ##PiTab
        output$OriPiTab <- renderTable({
          piesinseq=Strfindpies(toupper(seqiqi),Pies,maxMM)
          #Subsampling sequence and creating table
          pimattab=c()
          for(pipat in piesinseq){
            
            #matchpattern to find targets
            matseqpie=as.character(matchPattern(DNAString(pipat),DNAString(toupper(seqiqi)),max.mismatch=maxMM,fixed=T))
            
            for(mat in matseqpie){
              pimattab=rbind(pimattab,cbind(pipat,mat,stringdist(pipat,mat,"hamming")))
            }
          }
          if(length(piesinseq) > 0){
            pimattab=cbind(PiesFin[pimattab[,1],2],pimattab)
            #pimattab[,2]=paste("substitute(paste(italic(",PiesFin[pimattab[,2],3],")))",sep="")
            pimattab[,2]=paste("",PiesFin[pimattab[,2],3],"",sep="")
            
            #if(ncol(as.data.frame(pimattab))!=1){pimattab=pimattab[order(pimattab[,4]),]}
            
            rownames(pimattab)=1:nrow(pimattab)
            #colnames(pimattab)=c("piRNA locus","21-U reverse complement sequence","Matching sequence","Edit distance")
            #pimattab=pimattab[,-2]
            colnames(pimattab)=c("piRNA","Matching piRNA","Matching sequence","Mismatches")
            pimattab=pimattab[,-3]
          }
          
          #if(ncol(pimattab) == 1){
          # pimattab=t(pimattab)
          #  }
          #cat(pimattab,class(pimattab))
          if(ncol(as.data.frame(pimattab))==1){
            pimattab=t(as.data.frame(pimattab))
            #colnames(pimattab)=c("piRNA locus","Matching sequence","Edit distance")
            colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches")
            pimattab=cbind(pimattab,extraPiinfo[as.character(pimattab[,1]),2:3])
            colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches","piRNA abundance<sup>1</sup>","Tissue expression<sup>2</sup>")
          }else{
            if(ncol(as.data.frame(pimattab)) == 3){
              pimattab=pimattab[order(pimattab[,3]),]
              pimattab=cbind(pimattab,extraPiinfo[as.character(pimattab[,1]),2:3])
              colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches","piRNA abundance<sup>1</sup>","Tissue expression<sup>2</sup>")
            }
          }
          pimattab
          
        }, sanitize.text.function = function(x) x)
        output$TextPitab <- renderUI(
          {
            #HTML("<p style=\"font-size:12px;\"><b>1:</b> Calculated<br><b>2:</b></p>")
            HTML("<b><sup>1</sup>:</b> Average reads per million seen in short RNA sequencing libraries (<a href=\"https://doi.org/10.1038/s41592-021-01369-z\">Priyardarshini <i>et al.</i> 2022</a>).<br>
                 <b><sup>2</sup>:</b> Enrichment score as described in <a href=\"https://pubmed.ncbi.nlm.nih.gov/23516384/\">Billi <i>et al.</i> (2013)</a>.")
          })
        
      }
      output$button4later <- renderUI({actionButton("actionRESET", label = "RESET")})
      ErrorFlag = 2
    }
    
    ####Where the future lays, i.e., the place where the sequence is created
    ####################################From here change the dynamic output so a table or something appears on the meanwhile, you can reander the gif for now
    
    if(ErrorFlag == 0){ ##Main routine
      ###Here's where the dynamic interface starts
      ##Render a data table, displaying current execution time
      ##Or expected completion time and status of process
      ##Status
      #output$DynamicUserInterface <- renderUI({verbatimTextOutput("Counter")})
      
      output$DynamicUserInterface <- renderUI({
        #HTML("")
        includeHTML("www/clock.html")
      })
      
      #ElapsedTime <- c(0)
      #TimerActive <- c(TRUE)
      
      ###Send custom messages to start reading status
      session$sendCustomMessage("submitted-job", TRUE)
      
      myFuture <- future({
        #withProgress(message = 'Generating transgene', style = "notification", detail = "(~4 min per kb for piRNA optimization and GLO algorithm)", value = 0, {
        ###Internal parameters
        RetrieveTop=5000
        
        writeLines(text="Performing Ribosomal optimization",con=statfiletowork)
        
        ######1st. step: Ribosomal Optimization
        if(FlaRi){ ###Ribosomal binding optimization
          testSeq=paste(c(seqDNA[1:39]),sep="",collapse="")
          write(paste("AAAA",testSeq,sep="",collapse=""),paste("DATA/users/",session_id,"/gene.fasta", sep=""))
          for(h in 1:100){
            testSeq2=repcds(testSeq,CAIS,AAtoCodF,8)
            write(paste("AAAA",testSeq2,sep="",collapse=""),paste("DATA/users/",session_id,"/gene.fasta", sep=""),append=T)
          }
          system(paste("sh bin/rnafold.sh",paste("DATA/users/",session_id,sep=""),RetrieveTop))
          inseqs=readLines(paste("DATA/users/",session_id,"/seqswithoutaaaas.txt", sep=""))
          fwd=grep("GGTCTC",x=inseqs)
          rev=grep("GAGACC",x=inseqs)
          all=unique(c(fwd,rev))
          if((length(all) != 0)&(length(all) != length(inseqs))){inseqs=inseqs[-c(all)]}
          id=order(sapply(inseqs,function(x){Strcondepies(x,Pies,maxMM)}))[1]
          SeqStart=inseqs[id]
        }else{ ###Not optimization
          if(CodonAl != 6){
            SeqStart=repcds(paste(c(seqDNA[1:39]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl)
          }
          if(FlaModOnly){
            SeqStart=paste(c(seqDNA[1:39]),sep="",collapse="")
          }
        }
        
        
        #incProgress(3/10)
        
        ######2nd. step: Sequence Optimization
        writeLines(text="Adapting sequence by specified algorithm or codon usage",con=statfiletowork)
        
        ##CHeck box that overrides sequence optimization
        if(FlaModOnly){
          SeqtoOpt=paste(c(SeqStart,seqDNA[40:length(seqDNA)]),sep="",collapse="")
        }else{
          
          
          if((CodonAl <= 5)&(!(FlaPi))){ ###Use CAI equal to one
            SeqEnd=repcds(paste(c(seqDNA[40:length(seqDNA)]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl)
            
          }else{ 
            
            ##Check for Dans algoritmh, this should sepercede ribosomal optimization
            if(CodonAl == 6){ ###Produce start seq as originally it wasnt 
              
              write(paste(seqDNA,sep="",collapse=""),paste("DATA/users/",session_id,"/gene.fasta", sep=""))
              system(paste("perl bin/GLO_CLI_one_line.pl",paste("DATA/users/",session_id,"/gene.fasta", sep=""),">", paste("DATA/users/",session_id,"/GLO.fasta", sep="")))
              
              GLOseq=readLines(paste("DATA/users/",session_id,"/GLO.fasta", sep=""))
              GLOseq=unlist(strsplit(toupper(GLOseq),""))
              
              SeqStart=paste(c(GLOseq[1:39]),sep="",collapse="")
              SeqEnd=paste(c(GLOseq[40:length(GLOseq)]),sep="",collapse="")
              
            }else{ ###Use Codon Al as column of frequencies to mimic, and sample 100 times to reduce piRNAs   
              setrep=c()
              for(j in 1:100){
                setrep=append(setrep,repcds(paste(c(seqDNA[40:length(seqDNA)]),sep="",collapse=""),CAIS,AAtoCodF,CodonAl))
              }
              inseqs=setrep
              ###Outdated, we dont remove anymore RE sites before
              #fwd=grep("GGTCTC",x=inseqs)
              #rev=grep("GAGACC",x=inseqs)
              #all=unique(c(fwd,rev))
              #if((length(all) != 0)&(length(all) != length(inseqs))){inseqs=inseqs[-c(all)]}
              id=order(sapply(inseqs,function(x){Strcondepies(x,Pies,maxMM)}))[1]
              SeqEnd=inseqs[id]
            }
          }
          
          
          SeqtoOpt=paste(c(SeqStart,SeqEnd),sep="",collapse="")
          #############
        }
        #incProgress(3/10)
        
        ######3rd. step: Sequence maniputalion for piRNA removal
        writeLines(text="Removing piRNA sites",con=statfiletowork)
        
        ##If PiRNA removal
        if(FlaPi){
          pipipis=Strfindpies(SeqtoOpt,Pies,maxMM)
          if(length(pipipis) > 0 ){
            stpos=c()
            for(pipi in pipipis){
              stpos=append(stpos,start(matchPattern(DNAString(pipi),DNAString(SeqtoOpt),max.mismatch=maxMM,fixed=T)))
            }
            stpos=unique(c(stpos,stpos+3,stpos+6,stpos+9,stpos+12,stpos+15))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
            
            pipipis=Strfindpies(SeqtoOpt,Pies,maxMM)
            if(length(pipipis) > 0 ){
              Iter=1
              nflag=TRUE
              while((Iter < 100)&(nflag)){
                
                stpos=c()
                for(pipi in pipipis){
                  stpos=append(stpos,start(matchPattern(DNAString(pipi),DNAString(SeqtoOpt),max.mismatch=maxMM,fixed=T)))
                }
                stpos=unique(c(stpos,stpos+3,stpos+6,stpos+9,stpos+12,stpos+15))
                SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
                
                pipipis=Strfindpies(SeqtoOpt,Pies,maxMM)
                
                if(length(pipipis) > 0 ){Iter=1+Iter}else{nflag=FALSE}
              }
              ##Error Message for iterations
              if(Iter==100){output$ErrorMessage <- renderText({paste("Error: PiRNA removal did not work even after 100 iterations. Final number of piRNA sites found was: ",length(pipipis))})}
            }
          }
          
        }
        
        #incProgress(3/10)
        
        ######4th. step: Sequence maniputalion for restriction site removal
        writeLines(text="Restriction Sites removal",con=statfiletowork)
        
        #Inenzy=c(as.character(gegegenzymes),as.character(gogogonzymes))
        
        #If restriction sites
        if(FlaEnz & (length(Inenzy) > 0)){
          enpat=c(as.character(enzy[Inenzy,"Site"]))
          for(pattemp in enpat){
            enpat=append(enpat,as.character(reverseComplement(DNAString(as.character(pattemp)))))
          }
          enpat=unique(enpat)
          #First search
          all=c()
          for(patito in enpat){
            all=append(all,grep(as.character(patito),x=SeqtoOpt))
          }
          all=unique(all)
          if(length(all) > 0 ){ ##Do proper biostrings match
            stpos=c()
            for(patito in enpat){
              stpos=append(stpos,start(matchPattern(DNAString(as.character(patito)),DNAString(SeqtoOpt),fixed=F)))
            }
            stpos=unique(c(stpos))
            SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
            
            all=c()
            for(patito in enpat){
              all=append(all,grep(as.character(patito),x=SeqtoOpt))
            }
            all=unique(all)
            
            if(length(all) > 0 ){
              Iter=1
              nflag=TRUE
              while((Iter < 100)&(nflag)){
                stpos=c()
                for(patito in enpat){
                  stpos=append(stpos,start(matchPattern(DNAString(as.character(patito)),DNAString(SeqtoOpt),fixed=F)))
                }
                stpos=unique(c(stpos))
                SeqtoOpt=modbyposiz(SeqtoOpt,stpos,CAIS,AAtoCodF,8)
                
                all=c()
                for(patito in enpat){
                  all=append(all,grep(as.character(patito),x=SeqtoOpt))
                }
                all=unique(all)
                
                if(length(all) > 0 ){Iter=1+Iter}else{nflag=FALSE}
              }
              ##Error Message for iterations
              if(Iter==100){output$ErrorMessage <- renderText({paste("Error: Restriction site removal did not work even after 100 iterations. Final number of sites found was: ",length(all))})}
            }
          }
          
          
        }
        
        #incProgress(1/20)
        
        ######5th. step: Sequence maniputalion to add introns
        writeLines(text="Adding introns",con=statfiletowork)
        
        finalvec=c()
        
        if(FlaIn){
          ##Make a vector where to insert introns later
          finalvec=unlist(strsplit(toupper(SeqtoOpt),""))
          
          ##Find all AAGR
          AAGRs=start(matchPattern(DNAString("AAGR"),DNAString(SeqtoOpt),fixed=F))
          #Place in position to cut start position, i.e. AAG][R
          if(length(AAGRs)>0){AAGRs=AAGRs+2}else{AAGRs=c()}
          
          ##Find all AAGN
          AAGs=start(matchPattern(DNAString("AAG"),DNAString(SeqtoOpt),fixed=F))
          #Place in position to cut start position, i.e. AAG][N
          if(length(AAGs)>0){AAGs=AAGs+2}else{AAGs=c()}
          
          ##Find all AGR
          AGRs=start(matchPattern(DNAString("AGR"),DNAString(SeqtoOpt),fixed=F))
          #Place in position to cut start position, i.e. AG][R
          if(length(AGRs)>0){AGRs=AGRs+1}else{AGRs=c()}
          
          ##Find all AG
          AGs=start(matchPattern(DNAString("AG"),DNAString(SeqtoOpt),fixed=F))
          #Place in position to cut start position, i.e. AG][N
          if(length(AGs)>0){AGs=AGs+1}else{AGs=c()}
          
          ##Find all GR
          GRs=start(matchPattern(DNAString("GR"),DNAString(SeqtoOpt),fixed=F))
          #Place in position to cut start position, i.e. G][R
          if(length(GRs)>0){GRs=GRs}else{GRs=c()}
          
          ##Find all NG
          Gs=start(matchPattern(DNAString("G"),DNAString(SeqtoOpt),fixed=T))
          #Place in position to cut start position, i.e. G][N
          if(length(Gs)==0){inposis=c(49,149,249)}else{
            trypos=list(AAGRs,AAGs,AGRs,AGs,GRs,Gs)
            inposis=c()
            trpo=1
            while(length(inposis)!=3){
              #First statement
              ttpos=trypos[[trpo]]
              inpos=c()
              
              ##Middle check
              #Filter for intron
              if(Flaframeint){ttpos=ttpos[(ttpos %% 3 == 0)]}
              ##Now check for location
              if(length(ttpos) >= 3){
                if(typedistint == 1){##if early start
                  if(sum((ttpos > 50)&(ttpos <150))>1){
                    inpos=append(inpos,c(ttpos[(ttpos > 50)&(ttpos<150)])[1])
                    if(sum(ttpos > (inpos[1]+150))>1){
                      inpos=append(inpos,c(ttpos[(ttpos > (inpos[1]+150))])[1])
                      if(sum(ttpos > (inpos[2]+150))>1){
                        inpos=append(inpos,c(ttpos[(ttpos > (inpos[2]+150))])[1])
                        inposis=inpos[order(inpos)]
                      }
                    }
                  }
                }else{
                  inposis=as.integer(quantile(ttpos,names=F)[c(2,3,4)])
                  inposis=inposis[order(inposis)]	
                }
              }
              
              #End statement
              trpo=trpo+1
              if(trpo > length(trypos)){inposis=c(49,149,249)}
            }
          }
          SeqtoOpt=paste(c(finalvec[1:inposis[1]],as.character(IntronSeqs[typeIn,"Sequence_1"]),finalvec[(inposis[1]+1):inposis[2]],as.character(IntronSeqs[typeIn,"Sequence_2"]),finalvec[(inposis[2]+1):inposis[3]], as.character(IntronSeqs[typeIn,"Sequence_3"]),finalvec[(inposis[3]+1):length(finalvec)]),sep="",collapse="")
        }
        
        
        results=list(SeqtoOpt,finalvec)
        ####################################################################################################################
        results
        ###Consider writing SeqtoOPt in usr folder so no race is invoked
      }) 
      
      
      then(myFuture, onFulfilled = function(resultslist)
        
      {#############################################Future happening############################################################3
        #############################################################
        ######6th. step: Display of results
        SeqtoOpt=unlist(resultslist[1])
        finalvec=unlist(resultslist[2])
        #TwistResults=resultslist[3]
        ##Variables to consider
        #IF: FlaTURS
        #  befUTR,SeqtoOpt,aftUTR
        
        ##If annotation, show results
        #SeqtoOpt <- .
        ###Also, it seems future forgets about all variables created previously which makes sense
        ###This allows to add dynamic outputs
        output$DynamicUserInterface <- renderUI({uiOutput("AllResults")})
        ###Send custom messages to stop status
        session$sendCustomMessage("submitted-job", FALSE)
        
        ##Extra things coming from before
        gegegenzymes=c()
        gogogonzymes=input$Oenzymes
        
        ###Clean status file
        if(file.exists(statfiletowork)){file.remove(statfiletowork)}
        
        Inenzy=c(as.character(gegegenzymes),as.character(gogogonzymes))
        
        ##This creates these dynamic outputs
        ##A clear function will maybe introduce errors asthe dynamic output differs, similarly, if those are clear just before getting displayed, I dont think will make a big difference
        
        
        #####change it to twist to add sequence viewer
        
        if(FlaTwisty){
          output$AllResults <- renderUI({
            fluidRow(
              uiOutput("button4later"),
              br(),
              ##First optimized
              
              htmlOutput("newsequence"),
              br(),
              #uiOutput("downloadoptseq"),
              tableOutput("OptiPiTab"),
              htmlOutput("TextOPitab"),
              br(),
              uiOutput("downloadApeOptiseq"),
              ###Then Twist analysis
              htmlOutput("newsequenceTwist"),
              ##Then previous
              htmlOutput("oldsequence"),
              br(),
              tableOutput("OriPiTab"),
              htmlOutput("TextPitab"),
              br(),
              uiOutput("downloadApeOriseq")
              #DT::dataTableOutput('OptiPiTab')
            )
          })
        }else{
          output$AllResults <- renderUI({
            fluidRow(
              #actionButton("actionRESET", label = "RESET"),
              #actionButton("actionRESET", label = "RESET"),
              uiOutput("button4later"),
              br(),
              ##First optimized
              #DT::dataTableOutput('OriPiTab'),
              htmlOutput("newsequence"),
              br(),
              #uiOutput("downloadoptseq"),
              tableOutput("OptiPiTab"),
              htmlOutput("TextOPitab"),
              br(),
              uiOutput("downloadApeOptiseq"),
              ##Then previous
              htmlOutput("oldsequence"),
              br(),
              tableOutput("OriPiTab"),
              htmlOutput("TextPitab"),
              br(),
              uiOutput("downloadApeOriseq")
              #DT::dataTableOutput('OptiPiTab')
            )
          })
        }
        
        #output$button4later <- renderUI({HTML("<b>Loading results ...</b>")})
        output$button4later <- renderUI({HTML("<b></b>")})
        output$ErrorMessage <- renderText({""})
        
        
        ##Now sequence manipulation
        aaaads=""
        if((FlaRi)|(Fla5paaa)){aaaads="aaaa"}
        
        if(!(FlaIn)){finalvec=unlist(strsplit(toupper(SeqtoOpt),""))}
        if(!(FlaIn)){optsec=SeqtoOpt}else{optsec=paste(finalvec,sep="",collapse="")}
        
        splitseqopt=unlist(strsplit(optsec,""))
        newseqstart=paste(splitseqopt[1:30],sep="",collapse="")
        newseqend=paste(splitseqopt[(length(splitseqopt)-29):(length(splitseqopt))],sep="",collapse="")
        
        ###OK now make sure to put things in order, first optimized, if fla twisty, thingy for flatwisty, and finally original sequence
        ##15/08/22: Add nor promoter sequence
        SeqtoOpt=paste(Proseq,befUTR,aaaads,SeqtoOpt,aftUTR,sep="",collapse="")
        ###Graphical Output
        
        if((ErrorFlag == 0) & !is.null(SeqtoOpt)) {
          
          ####New approach
          output$newsequence <-renderUI({
            
            ##New approach
            paterns=c()
            seqpaterns=c()
            colpaterns=c()
            typpaterns=c()
            dftyp=c()
            dfcol=c()
            Neweeexxttpar=eeexxttpar
            #Ezymes
            if(length(Inenzy)>0){
              paterns=append(paterns,Inenzy)
              seqpaterns=append(seqpaterns,enzy[Inenzy,"Site"])
              colpaterns=append(colpaterns,rep("#d7bde2",length(Inenzy)))
              typpaterns=append(typpaterns,rep("RE_site",length(Inenzy)))
              dftyp=append(dftyp,"RE site")
              dfcol=append(dfcol, "#d7bde2")
            }
            #Introns
            if(FlaIn){
              paterns=append(paterns,paste(rownames(IntronSeqs[typeIn,])," intron",sep=""))
              #paterns=append(paterns,"A")
              seqpaterns=append(seqpaterns,toupper(as.character(IntronSeqs[typeIn,"Sequence_1"])))
              colpaterns=append(colpaterns,"#f5cba7")
              typpaterns=append(typpaterns,"intron")
              
              paterns=append(paterns,paste(rownames(IntronSeqs[typeIn,])," intron",sep=""))
              #paterns=append(paterns,"B")
              seqpaterns=append(seqpaterns,toupper(as.character(IntronSeqs[typeIn,"Sequence_2"])))
              colpaterns=append(colpaterns,"#f5cba7")
              typpaterns=append(typpaterns,"intron")
              
              paterns=append(paterns,paste(rownames(IntronSeqs[typeIn,])," intron",sep=""))
              #paterns=append(paterns,"C")
              seqpaterns=append(seqpaterns,toupper(as.character(IntronSeqs[typeIn,"Sequence_3"])))
              colpaterns=append(colpaterns,"#f5cba7")
              typpaterns=append(typpaterns,"intron")
              dftyp=append(dftyp,"Intron")
              dfcol=append(dfcol, "#f5cba7")
            }
            ##If UTRs
            if(FlaTURS){
              if(befUTR !=""){
                paterns=append(paterns,nam5UTR)
                seqpaterns=append(seqpaterns,toupper(befUTR))
                colpaterns=append(colpaterns,"#a3e4d7")
                typpaterns=append(typpaterns,"5'UTR")
                dftyp=append(dftyp,"5' UTR")
                dfcol=append(dfcol, "#a3e4d7")
              }
              if(aftUTR !=""){
                paterns=append(paterns,nam3UTR)
                seqpaterns=append(seqpaterns,toupper(aftUTR))
                colpaterns=append(colpaterns,"#a9cce3")
                typpaterns=append(typpaterns,"3'UTR")
                dftyp=append(dftyp,"3' UTR")
                dfcol=append(dfcol, "#a9cce3")
              }
            }
            ###If promoters
            if(FlaPro){
              if(Proseq !=""){
                paterns=append(paterns,Pronam)
                seqpaterns=append(seqpaterns,toupper(Proseq))
                colpaterns=append(colpaterns,"#fad7a0")
                typpaterns=append(typpaterns,"Promoter")
                dftyp=append(dftyp,"Promoter")
                dfcol=append(dfcol, "#fad7a0")
              }
            }
            if(FlaPi){
              piss=Strfindpies(toupper(optsec),Pies,maxMM)
              if(length(piss)>0){
                popos=c()
                papas=c()
                for(pipi in piss){
                  patotes=as.character(matchPattern(DNAString(pipi),DNAString(toupper(optsec)),max.mismatch=maxMM,fixed=T))
                  popos=append(popos,patotes)
                  papas=append(papas,rep(PiesFin[piss,2],length(patotes)))
                }
                
                papos=c()
                pospos=unique(popos)
                for(pipi in pospos){
                  papos=append(papos,paste(papas[which(pipi == popos)], collapse=";"))  
                }
                
                paterns=append(paterns,papos)
                seqpaterns=append(seqpaterns,pospos) 
                colpaterns=append(colpaterns,rep("#f9e79f",length(pospos)))
                typpaterns=append(typpaterns,rep("piRNA",length(pospos)))
                
                for(totoro in 1:length(pospos)){
                  Neweeexxttpar=append(Neweeexxttpar,paste(papos[totoro]," ",pospos[totoro],sep=""))
                }
                
                dftyp=append(dftyp,"piRNA homology")
                dfcol=append(dfcol, "#f9e79f")
              }
            }
            
            #write(paste(PasteApe(OriSeqNameIn,seqiqi,c(seqpaterns),c(colpaterns),c(colpaterns),c(paterns),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqpop.gb", sep=""))
            #NewPasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,FeatType,tabibi,cai,list,PiesList,extracomments=c(),optsecnotr=""){
            #write(paste(PasteApe(paste("Optimized_",OriSeqNameIn,sep=""),SeqtoOpt,c(seqpaterns),c(colpaterns),c(colpaterns),c(paterns),CAIS,5,AAtoCodF,Pies,extracomments=Neweeexxttpar,optsecnotr=optsec),collapse="\n"),paste("DATA/users/",session_id,"/Seqpop.gb", sep=""))
            write(paste(NewPasteApe(paste("Optimized_",OriSeqNameIn,sep=""),SeqtoOpt,c(seqpaterns),c(colpaterns),c(colpaterns),c(paterns),c(typpaterns),CAIS,5,AAtoCodF,Pies,extracomments=Neweeexxttpar,optsecnotr=optsec,maximm=maxMM),collapse="\n"),paste("DATA/users/",session_id,"/Seqpop.gb", sep=""))
            
            legdf=data.frame(Type=c(dftyp),Color=c(dfcol))
            
            #HTML("<br><b><h3>Input ",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(optsec,CAIS,5,AAtoCodF),NewSequenceViewerDu("","oldestseq",seqiqi,c(seqpaterns),c(colpaterns),c(paterns),legdf,oriseqstart,oriseqend),"<br>")
            HTML("<br><b><h3>Optimized",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(optsec,CAIS,5,AAtoCodF),NewSequenceViewerDu("","newseq",SeqtoOpt,c(seqpaterns),c(colpaterns),c(paterns),legdf,newseqstart,newseqend),"<br>")
            #HTML("<br><b><h3>Optimized",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(optsec,CAIS,5,AAtoCodF),NewSequenceViewerDu("","newseq",SeqtoOpt,c(seqpaterns),c(paterns),c(colpaterns),legdf,newseqstart,newseqend),"<br>")
            
          })
          
          
          output$downloadApeOptiseq <- renderUI({
            downloadButton('DownOptiApe', 'Genbank with annotations')
          })
          
          ##PiTab
          if(FlaPi){
            output$OptiPiTab <- renderTable({
              piesinseq=Strfindpies(toupper(optsec),Pies,maxMM)
              #Subsampling sequence and creating table
              pimattab=c()
              for(pipat in piesinseq){
                
                #matchpattern to find targets
                matseqpie=as.character(matchPattern(DNAString(pipat),DNAString(toupper(optsec)),max.mismatch=maxMM,fixed=T))
                
                for(mat in matseqpie){
                  pimattab=rbind(pimattab,cbind(pipat,mat,stringdist(pipat,mat,"hamming")))
                }
              }
              
              if(length(piesinseq) > 0){
                pimattab=cbind(PiesFin[pimattab[,1],2],pimattab)
                pimattab[,2]=paste(PiesFin[pimattab[,2],3],"",sep="")
                #if(ncol(as.data.frame(pimattab))!=1){pimattab=pimattab[order(pimattab[,4]),]}
                rownames(pimattab)=1:nrow(pimattab)
                #colnames(pimattab)=c("piRNA locus","21-U reverse complement sequence","Matching sequence","Edit distance")
                #pimattab=pimattab[,-2]
                colnames(pimattab)=c("piRNA","Matching piRNA","Matching sequence","Mismatches")
                pimattab=pimattab[,-3]
              }
              if(ncol(as.data.frame(pimattab))==1){
                pimattab=t(as.data.frame(pimattab))
                #colnames(pimattab)=c("piRNA locus","Matching sequence","Edit distance")
                colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches")
                pimattab=cbind(pimattab,extraPiinfo[as.character(pimattab[,1]),2:3])
                colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches","piRNA abundance<sup>1</sup>","Tissue expression<sup>2</sup>")
              }else{
                if(ncol(as.data.frame(pimattab)) == 3){
                  pimattab=pimattab[order(pimattab[,3]),]
                  pimattab=cbind(pimattab,extraPiinfo[as.character(pimattab[,1]),2:3])
                  colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches","piRNA abundance<sup>1</sup>","Tissue expression<sup>2</sup>")
                }
              }
              pimattab
            }, sanitize.text.function = function(x) x)
            output$TextOPitab <- renderUI(
              {
                #HTML("<p style=\"font-size:12px;\"><b>1:</b> Calculated<br><b>2:</b></p>")
                HTML("<b><sup>1</sup>:</b> Average reads per million seen in short RNA sequencing libraries (<a href=\"https://doi.org/10.1038/s41592-021-01369-z\">Priyardarshini <i>et al.</i> 2022</a>).<br>
                 <b><sup>2</sup>:</b> Enrichment score as described in <a href=\"https://pubmed.ncbi.nlm.nih.gov/23516384/\">Billi <i>et al.</i> (2013)</a>.")
              })
            
          }
          
          
          ####Now Twist viewer
          if(FlaTwisty){
            TwistResults=TwistSynthesisWithCoordinatesREDO(toupper(SeqtoOpt))
            Twistcheck=TwistResults[[1]]
            errors=TwistResults[[2]]
            regions=TwistResults[[3]]
            twistres=TwistResults[[4]]
            if(Twistcheck){
              output$newsequenceTwist <-renderUI({
                #HTML("<br><b><h3>Optimized ",OriSeqNameIn," <i>in-silico</i> synthesis results</h3></b><br><p style=\"color:green\"><b>Simple synthesis - no problems detected</b><br></p><h5><p style=\"color:green\">",paste0(paste0(twistres),sep="<br>",collapse=""),"<br></p></h5>")
                HTML("<br><b><h3>Optimized ",OriSeqNameIn," <i>in-silico</i> synthesis results</h3></b><br><p style=\"color:green\"><b>Simple synthesis - no problems detected</b><br></p>")
                #NewSequenceViewer("","twistannoseq",SeqtoOpt,c(),c(),c(),c(),"","")
              })
            }else{
              starts=as.integer(regions[,1])
              ends=as.integer(regions[,2])
              
              patseqs=as.character(extractAt(DNAString(toupper(SeqtoOpt)),IRanges(start=starts,end=ends)))
              newdf=data.frame(Seqs=patseqs,Description=(regions)[,3],Color=(regions)[,4])
              
              ##Alternative df
              altndfst=unlist(strsplit(as.character(newdf[,"Description"]),":"))[c(TRUE,FALSE)]
              if(length(altndfst)==nrow(newdf)){
                newdf=data.frame(Seqs=patseqs,Description=c(altndfst),Color=(regions)[,4])
              }
              newdf=unique(newdf)
              
              legdf=data.frame(Type=as.character(newdf[,2]),Color=as.character(newdf[,3]))
              
              output$newsequenceTwist <-renderUI({
                
                HTML("<br><b><h3>Optimized ",OriSeqNameIn," <i>in-silico</i> synthesis results</h3></b><br><b><p style=\"color:red\">Complex synthesis:</b></p><br><h5><p style=\"color:red\">",paste0(paste0(errors),sep="<br>",collapse=""),"<br></p></h5>",NewSequenceViewerTW("","twistannoseq",SeqtoOpt,c(as.character(newdf[,1])),c(as.character(newdf[,3])),c(as.character(newdf[,2])),legdf))
              })
            }
          }
          ##Now data and viewer of old sequence
          if(as.integer(typinput) == 1){ ###If original sequence
            
            output$oldsequence <-renderUI({
              seqiqi=toupper(seqiqi)
              
              piss=Strfindpies(seqiqi,Pies,maxMM)
              
              ##New approach
              paterns=c()
              seqpaterns=c()
              colpaterns=c()
              typpaterns=c()
              dftyp=c()
              dfcol=c()
              Neweeexxttpar=eeexxttpar
              #Ezymes
              if(length(Inenzy)>0){
                paterns=append(paterns,Inenzy)
                seqpaterns=append(seqpaterns,enzy[Inenzy,"Site"])
                colpaterns=append(colpaterns,rep("#d7bde2",length(Inenzy)))
                typpaterns=append(typpaterns,rep("RE_site",length(Inenzy)))
                dftyp=append(dftyp,"RE site")
                dfcol=append(dfcol, "#d7bde2")
              }
              
              if(FlaPi){
                piss=Strfindpies(seqiqi,Pies,maxMM)
                if(length(piss)>0){
                  popos=c()
                  papas=c()
                  for(pipi in piss){
                    patotes=as.character(matchPattern(DNAString(pipi),DNAString(seqiqi),max.mismatch=maxMM,fixed=T))
                    popos=append(popos,patotes)
                    papas=append(papas,rep(PiesFin[piss,2],length(patotes)))
                  }
                  
                  papos=c()
                  pospos=unique(popos)
                  for(pipi in pospos){
                    papos=append(papos,paste(papas[which(pipi == popos)], collapse=";"))  
                  }
                  
                  paterns=append(paterns,papos)
                  seqpaterns=append(seqpaterns,pospos) 
                  colpaterns=append(colpaterns,rep("#f9e79f",length(pospos)))
                  typpaterns=append(typpaterns,rep("piRNA",length(pospos)))
                  
                  for(totoro in 1:length(pospos)){
                    Neweeexxttpar=append(Neweeexxttpar,paste(papos[totoro]," ",pospos[totoro],sep=""))
                  }
                  
                  dftyp=append(dftyp,"piRNA homology")
                  dfcol=append(dfcol, "#f9e79f")
                }
              }
              #PasteApe = function(locus_name,sequence,patterns,FWDcolors,REVcolors,tooltips,tabibi,cai,list,PiesList,extracomments=c(),optsecnotr=""){
              
              #write(paste(PasteApe(OriSeqNameIn,seqiqi,c(seqpaterns),c(colpaterns),c(colpaterns),c(paterns),CAIS,5,AAtoCodF,Pies,extracomments=eeexxttpar),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
              write(paste(NewPasteApe(OriSeqNameIn,seqiqi,c(seqpaterns),c(colpaterns),c(colpaterns),c(paterns),c(typpaterns),CAIS,5,AAtoCodF,Pies,extracomments=Neweeexxttpar,maximm=maxMM),collapse="\n"),paste("DATA/users/",session_id,"/Seqog.gb", sep=""))
              
              legdf=data.frame(Type=c(dftyp),Color=c(dfcol))
              
              HTML("<br><b><h3>Input ",OriSeqNameIn,"</h3></b><br>",GCHTMLinfo(seqiqi,CAIS,5,AAtoCodF),NewSequenceViewerDu("","oldestseq",seqiqi,c(seqpaterns),c(colpaterns),c(paterns),legdf,oriseqstart,oriseqend),"<br>")
              
              
              
            })
            
            ##ApeButtonDownOriginal
            output$downloadApeOriseq <- renderUI({
              downloadButton('DownOriApe', 'Genbank with annotations')
            })
            
            ##PiTab
            if(FlaPi){
              output$OriPiTab <- renderTable({
                piesinseq=Strfindpies(toupper(seqiqi),Pies,maxMM)
                #Subsampling sequence and creating table
                pimattab=c()
                for(pipat in piesinseq){
                  
                  #matchpattern to find targets
                  matseqpie=as.character(matchPattern(DNAString(pipat),DNAString(toupper(seqiqi)),max.mismatch=maxMM,fixed=T))
                  
                  for(mat in matseqpie){
                    pimattab=rbind(pimattab,cbind(pipat,mat,stringdist(pipat,mat,"hamming")))
                  }
                }
                if(length(piesinseq) > 0){
                  pimattab=cbind(PiesFin[pimattab[,1],2],pimattab)
                  pimattab[,2]=paste(PiesFin[pimattab[,2],3],"",sep="")
                  #if(ncol(as.data.frame(pimattab))!=1){pimattab=pimattab[order(pimattab[,4]),]}
                  rownames(pimattab)=1:nrow(pimattab)
                  #if(ncol(pimattab)==4){pimattab=pimattab[order(pimattab[,4]),]}
                  #colnames(pimattab)=c("piRNA locus","21-U reverse complement sequence","Matching sequence","Edit distance")
                  #pimattab=pimattab[,-2]
                  colnames(pimattab)=c("piRNA","Matching piRNA","Matching sequence","Mismatches")
                  pimattab=pimattab[,-3]
                }
                
                if(ncol(as.data.frame(pimattab))==1){
                  pimattab=t(as.data.frame(pimattab))
                  #colnames(pimattab)=c("piRNA locus","Matching sequence","Edit distance")
                  colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches")
                  pimattab=cbind(pimattab,extraPiinfo[as.character(pimattab[,1]),2:3])
                  colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches","piRNA abundance<sup>1</sup>","Tissue expression<sup>2</sup>")
                }else{
                  if(ncol(as.data.frame(pimattab)) == 3){
                    pimattab=pimattab[order(pimattab[,3]),]
                    pimattab=cbind(pimattab,extraPiinfo[as.character(pimattab[,1]),2:3])
                    colnames(pimattab)=c("piRNA","Matching piRNA","Mismatches","piRNA abundance<sup>1</sup>","Tissue expression<sup>2</sup>")
                  }
                }
                pimattab
              }, sanitize.text.function = function(x) x)
              output$TextPitab <- renderUI(
                {
                  #HTML("<p style=\"font-size:12px;\"><b>1:</b> Calculated<br><b>2:</b></p>")
                  HTML("<b><sup>1</sup>:</b> Average reads per million seen in short RNA sequencing libraries (<a href=\"https://doi.org/10.1038/s41592-021-01369-z\">Priyardarshini <i>et al.</i> 2022</a>).<br>
                 <b><sup>2</sup>:</b> Enrichment score as described in <a href=\"https://pubmed.ncbi.nlm.nih.gov/23516384/\">Billi <i>et al.</i> (2013)</a>.")
                })
            }
          }
          
          
        }
        output$button4later <- renderUI({actionButton("actionRESET", label = "RESET")})
      },
      onRejected = function(){
        inituiui()
        ###Send custom messages to stop status
        session$sendCustomMessage("submitted-job", FALSE)
        
        ###Clean status file
        if(file.exists(statfiletowork)){file.remove(statfiletowork)}
        
        output$ErrorMessage <- renderText({
          paste("Unxexpected error, please contact us with the details of your submitted job")
        })
      }
      ) ## end of main then
      
      ###End Main sequence adaptation routine
      ##Clean clock
      
      
      #}) ### end of progresses
      return(NULL)
    }
    
  })##ENd main observer 
  
}) ####End of shiny server function

