## ZENODO : Downloads the requested file_name to the provided local_path, and returns the local position. 'force' forces the download despite local version already exists.
zenodo_get_data <- function(file_name = NULL, zenodo_id = NULL, local_path = NULL, force = FALSE) {
  ## Checks
  if (is.null(file_name)) stop("A file name to retrieve is required !")
  if (is.null(zenodo_id)) stop("A Zenodo ID is required !")
  if (is.null(local_path)) stop("A local path to store the file is required !")
  if (!is.logical(force)) stop("'force' should be TRUE or FALSE !")
  
  ## Build future local file path
  local_file <- paste0(local_path, "/", file_name)
  
  if (file.exists(local_file)) {
    message("The requested file to retrieve from Zenodo [", zenodo_id, "] already exists in [", local_path, "] !")
    if (!force) {
      message("'force' is set to 'FALSE' : nothing will be downloaded.")
      ## Nothing done : return local_file
      return(local_file)
    }
  }
  
  ## Retrieve data
  remote_file <- paste0("https://zenodo.org/records/", zenodo_id, "/files/", file_name)
  download.file(url = remote_file, destfile = local_file)
  
  ## Print and return local path
  message("Downloaded to [", local_file, "]")
  return(local_file)
}

## BACKUP : Copies the requested backup_file to the provided local_path, and returns the local position. 'force' forces the copy despite local version already exists.
backup_get_data <- function(backup_file = NULL, local_path = NULL, force = FALSE) {
  ## Checks
  if (is.null(backup_file)) stop("A backup file to retrieve is required !")
  if (is.null(local_path)) stop("A local path to store the file is required !")
  if (!is.logical(force)) stop("'force' should be TRUE or FALSE !")
  if (!file.exists(backup_file)) stop("Provided [", backup_file, "] does not exist !")
  
  ## Build future local file path
  local_file <- paste0(local_path, "/", basename(backup_file))
  
  if (file.exists(local_file)) {
    message("The requested file to retrieve already exists in [", local_path, "] !")
    if (!force) {
      message("'force' is set to 'FALSE' : nothing will be downloaded.")
      ## Nothing done : return local_file
      return(local_file)
    }
  }
  
  ## Retrieve data
  file.copy(from = backup_file,to = local_file)
  
  ## Print and return local path
  message("Copied to [", local_file, "]")
  return(local_file)
}

