library(shiny)

library(JBrowseR)

anno_dir <- "/Users/sparthib/Desktop/JBrowseR_test/"
bam_dir <- "/Users/sparthib/Desktop/JBrowseR_test/"
bam_file <- paste0(bam_dir, "EP1-BRN3B-RO_sorted.bam") 

ui <- fluidPage(
  titlePanel("View annotation JBrowseR Example"),
  # this adds to the browser to the UI, and specifies the output ID in the server
  JBrowseROutput("browserOutput")
)

# Warning message:
#   `bindFillRole()` only works on htmltools::tag() 
# objects (e.g., div(), p(), etc.), not objects of type 'shiny.tag.list'. 

server <- function(input, output, session) {
  # create the necessary JB2 assembly configuration
  assembly <- assembly(
    paste0(anno_dir, "gencode.v44.transcripts_short_header.fa"),
    bgzip = FALSE
  )
  

  # create configuration for a JB2 GFF FeatureTrack
  alignments <- track_alignments(
    bam_file,
    assembly
  )
  
  # create the tracks array to pass to browser
  tracks <- tracks(
    alignments
  )
  
  # set up the default session for the browser
  default_session <- default_session(
    assembly,
    c(alignments)
  )
  
  theme <- theme("#5da8a3", "#333")
  
  # link the UI with the browser widget
  output$browserOutput <- renderJBrowseR(
    JBrowseR(
      "View",
      assembly = assembly,
      tracks = tracks,
      location = "all",
      defaultSession = default_session,
      theme = theme
    )
  )
}

shinyApp(ui, server)



# Questions:
# JBrowseR for automation.
# way to view only reads that have translocations or only supplementary reads etc.? 
# Saving images of coverage plot.
# Linked view of supplementary reads. 
# Dealing with big bam files. 
# Chimeric reads - how do I get information out of a SAM/BAM file 
#as to which reads are linked to each other. 

# breakpoint split view. 
# linear read vs ref. 
# SA - tag added to all split alignments, or ones have translocation.

#samtools view -d SA  - has primary and supplementary reads. 
#2048 flag is only applied to n - 1 of the reads 



