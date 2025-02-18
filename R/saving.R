## Function to save data.frame

require(data.table)

saveDF <- function(
    spotcalldf,
    params = get('params', envir = globalenv()),
    currentFovName = NULL
){
  if(is.null(currentFovName)){
    currentFovName = params$current_fov #Try to see if there exists current FOV in params
    if(is.null(currentFovName)){
      stop('currentFovName not specified and params$current_fov not found: unable to save!')
    }
  }
  data.table::fwrite(spotcalldf, paste0(params$out_dir, 'SPOTCALL_', currentFovName, '.csv.gz' ))
}
