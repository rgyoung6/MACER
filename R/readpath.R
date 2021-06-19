#Roxygen2 Documentation:

#' @keywords internal
#'
#' @author Robert G. Young
#'
#' @references
#' https://github.com/rgyoung6/MACER
#' Young, R. G., Gill, R., Gillis, D., Hanner, R. H. (Submitted June 2021). Molecular Acquisition, Cleaning, and Evaluation in R (MACER) - A tool to assemble molecular marker datasets from BOLD and GenBank. Biodiversity Data Journal.
#'
#' @import utils
#'
################################ SELECT TARGET DIRECTORY AND FLAGGING THE OPERATING SYSTEM FUNCTION ###########

readpath <- function()
{

  if (interactive() && .Platform$OS.type == "windows"){
    target_dir<-choose.dir(getwd(), "Choose a suitable folder")
  }else{
    system("osascript -e 'tell app \"R\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
           intern = FALSE, ignore.stderr = TRUE)
    p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
    target_dir<-ifelse(length(p), p, NA)
  }

  return(target_dir)

}

