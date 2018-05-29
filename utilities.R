dl_eviction_data <- function(state, geography, format) {
    if (!dir.exists("data")) dir.create("data")
    origin_filename <- paste("https://eviction-lab-data-downloads.s3.amazonaws.com", 
                              state, geography, sep = "/")
    origin_filename <- paste0(origin_filename, ".", format)
    dest_filename <- paste0("data/ca-", geography, ".", format)
    if (!file.exists(dest_filename))
        download.file(origin_filename,
                      destfile = dest_filename)
    
}

fix_csv_names <- function(df) {
    rlang::set_names(df, ~stringr::str_replace_all(., "-", "_"))
}