########Gene Optimizer 2022 - Modular####
###Amhed Missael Vargas Velazquez
###avargas0lcg@gmail.com

###Short description:###
##By itself, the previous versions of the transgene builder could not allow multiple users due to the way that shiny is constructed, i.e. a single R session per app
##In theory, I could implement the package in a docker and mitigate the issues by creating a new session everytime a user connects it to. However, deploying a server only for that is a little cumbersome and not really required
##Similarly, I could use shiny proxy, but that also requires some configuration on the server's end and might obstruct the use we've implemented already.
##For now, the apparent most appropriate option is to use the package promises and start somewhat from scratch regarding the code. All right

##Also, in this final iteration, the code has been modified so every aspect of the sequences and user interface can be loaded from modifiable files.

#Load libraries
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyjs)
library(DT)

####ExternalFunction
busyIndicator <- function(text = "Loading...", wait=2000) {
  shiny::tagList(
    shiny::div(class="loadbanner",id="loadmessage",text,img(src="elegans3.gif"))
    ,shiny::tags$script(sprintf(
      " setInterval(function(){
         if ($('html').hasClass('shiny-busy')) {
          setTimeout(function() {
            if ($('html').hasClass('shiny-busy')) {
              $('div.loadbanner').show()
            }
          }, %d)          
        } else {
          $('div.loadbanner').hide()
        }
      },100)
      ",wait)
    )
  ) 
}

# Define User interface
shinyUI(
    fluidPage(
      tags$head(
        tags$link(rel="stylesheet",type = "text/css", href="bootstrap.min.css")
      ),
      ###Loading message
      tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 60px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #D3D3D3;
               z-index: 105;
             }
          ")),
    ##Custom extra styles: single sliders background and title of navbar  
    tags$style(type = 'text/css', 
               ".js-irs-none .irs-single, .js-irs-none .irs-bar-edge, .js-irs-none .irs-bar {
                          background: transparent;
                          border-top-color: transparent;
                          border-bottom-color: transparent;
                          border-left-color: transparent;
                          border-right-color: transparent}
               .navbar-default .navbar-brand:hover {color: #ffffff;}
               
    .selectize-input {
        width: 600px;
        padding-top: 5px;
      }
               "),
    tags$head(tags$script(src ="sequence-viewer.bundle.js")),
    ###Extra scripts communicating java with shiny
    tags$script("
                var activitytimer;
                
                function checkActivity() {
  Shiny.setInputValue('clockactivity', Math.random());
  console.log(\" :stat\" + Math.random());
  activitytimer = setTimeout(function(){ checkActivity() }, 1000);
};

Shiny.addCustomMessageHandler('submitted-job', function(activity) {
console.log(activity + \" :stat1\");
        if(activity){
          checkActivity();
console.log(activity + \" :stat2\");
        }else{
console.log(activity + \" :statFInal\");
          clearTimeout(activitytimer);
        }

      });

Shiny.addCustomMessageHandler('JobStatusMessage', function(JobMessage) {
  document.getElementById(\"jobstat\").innerText = JobMessage;
      });
    "),
includeHTML("www/ga.html"),
busyIndicator(),
#' ##Add font for Beta style
#' tags$style(HTML("
#'       @import url('https://fonts.googleapis.com/css2?family=Gochi+Hand&display=swap');
#'       #divBetaAnnounce {
#'         font-family: 'Gochi Hand', sans-serif;
#'         color: #c8b291;
#'       }
#'       ")),
tags$style(type = "text/css", ".navbar{padding-left:15px;
       padding-right:15px ; margin-right:auto; margin-left:auto;}"),
#Main tab pages
    navbarPage(
      #title=actionLink("link_to_tabpanel_sequenceadaptation", HTML("<b>Wormbuilder</b>")),
      title=HTML("<a href=\"https://wormbuilder.org/\"><b>Wormbuilder</b></a>"),
      windowTitle="WormBuilder transgenic tools",
        id = "panels",
        tabPanel("Sequence Adaptation",
                 mainPanel(
                     uiOutput("DynamicUserInterface")
                 , width = 12)),
      ###About
      tabPanel("About",
               mainPanel(
                 h3("Transgene builder"),
                 HTML("<p align=\"justify\">
                      This application is developed by postdoc Amhed Velazquez and is currently a beta version. Please let us know if you find any bugs so we can fix the site.
                 <br>
                 This website is generated via custom modified css/html page running in R via the shiny library. All the templates, libraries, and programs used to produce this site are under the MIT and GNU licenses.
                 </p>"),
                 br(),
                 HTML("<p align=\"justify\">
                 <i>Contact info</i>:<br>Christian-Froekjaer Jensen, Ph.D. 
                 <br>Assistant Professor of Bioscience
                 <br>Laboratory of Synthetic Genome Biology
                 <br>Email: <a href=\"mailto:cfjensen@kaust.edu.sa\">cfjensen@kaust.edu.sa</a>
                      </p>")
                 , width = 12)
      )
    ),
hr(),
    HTML("<a href=\"http://www.wormbuilder.org/\">Wormbuilder</a><br>"),
    HTML("<a href=\"mailto:amhed.velazquez@kaust.edu.sa\">Contact us!</a>")
    
)
)

