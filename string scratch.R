
plot_network2 <- function (string_ids, payload_id = NULL, required_score = NULL, 
    add_link = TRUE, add_summary = TRUE) 
{
    "\nDescription:\n  Plots an image of the STRING network with the given proteins.\n\nInput parameters:\n  \"string_ids\"        a vector of STRING identifiers\n  \"payload_id\"        an identifier of payload data on the STRING server (see method post_payload for additional informations)\n  \"required_score\"   a threshold on the score that overrides the default score_threshold, that we use only for the picture\n                      As default this option is active but we suggest to deactivate it in case one is generating many images (e.g. in a loop). \n                      Deactivating this option avoids to generate and store a lot of short-urls on our server.\n  \"add_summary\"       parameter to specify whether you want to add a summary text to the picture. This summary includes a p-value and the number of proteins/interactions.\n\nAuthor(s):\n   Andrea Franceschini\n"
    if (is.null(required_score)) 
        required_score = score_threshold
    img = get_png(string_ids, payload_id = payload_id, required_score = required_score)
    if (!is.null(img)) {
        plot(1:(dim(img)[2]), type = "n", xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "", ylim = c(1, dim(img)[1]), xlim = c(1, 
                (dim(img)[2])), asp = 1)
        if (add_summary) 
            mtext(get_summary(string_ids, required_score), cex = 0.7)
        rasterImage(img, 1, 1, (dim(img)[2]), dim(img)[1])
    }
}

string_db$get_png()

img = string_db$get_png(hits, payload_id = NULL, required_score = 200)

res2 <- readPNG(source = img[2], native = TRUE, info = TRUE)

writePNG(img, target="testpng3.png")

#------------

function (string_ids, required_score = NULL, network_flavor = "evidence", 
          file = NULL, payload_id = NULL) 
{
  "\nDescription:\n  Returns a png image of a STRING protein network with the given identifiers.\n\nInput parameters:\n  \"string_ids\"        a vector of STRING identifiers.\n  \"required_score\"    minimum STRING combined score of the interactions \n                        (if left NULL we get the combined score of the object, which is 400 by default)\n  \"network_flavor\"    specify the flavor of the network (\"evidence\", \"confidence\" or \"actions\".  default \"evidence\")\n  \"file\"              file where to save the image (must have .png extension)\n\nAuthor(s):\n   Andrea Franceschini\n"
  if (length(string_ids) > 2000) {
    cat("ERROR: We do not support lists with more than 2000 genes.\nPlease reduce the size of your input and rerun the analysis. \t")
    stop()
  }
  if (is.null(required_score)) 
  
  required_score = 200
  string_ids = unique(hits)
  string_ids = string_ids[!is.na(string_ids)]
  stable_url <- string_db$stable_url
  network_flavor <- "evidence"
  species <- string_db$species
  payload_id <- NULL
  urlStr = paste(stable_url, "/api/tsv-no-header/get_link?", sep = "")
  identifiers = ""
  for (id in string_ids) {
    identifiers = paste(identifiers, id, sep = "%0d")
  }
  params = list(required_score = required_score, required_score = required_score, 
                network_flavor = network_flavor, identifiers = identifiers, 
                species = species, caller_identity = "STRINGdb-package")
  if (!is.null(payload_id)) 
    params["internal_payload_id"] = payload_id
  img_link <- postFormSmart(urlStr, .params = params)
  if (!is.null(file)) 
    writePNG(img, file)
  return(img)
}

s

rasterImage(img, 1.2, 1.27, 1.8, 1.73, interpolate=FALSE)


res = postForm(urlStr, .params = params2, .opts = curlOptions(url = urlStr),
               curl = getCurlHandle(), style = 'HTTPPOST',
               .encoding = integer(), binary = NA, .checkParams = TRUE,
               .contentEncodeFun = curlEscape)

params2 = list(required_score = required_score, 
              network_flavor = network_flavor, identifiers = identifiers, 
              species = species, caller_identity = "STRINGdb-package")
#--------------------

postFormSmart <- function(uri, ..., .params = list(), .opts = curlOptions(url = uri),
                          curl = getCurlHandle(), style = 'HTTPPOST',
                          .encoding = integer(), binary = NA, .checkParams = TRUE,
                          .contentEncodeFun = curlEscape){
  
  res = postForm(uri, ..., .params = .params, .opts = .opts,
                 curl = curl, style = style,
                 .encoding = .encoding, binary = binary, .checkParams = .checkParams,
                 .contentEncodeFun = .contentEncodeFun)
  
  
  suppressWarnings( if(grepl("The document has moved", res)){
    
    begin <- regexpr("href",res)+6
    mys2=substr(res, begin, 10000000)
    end <- regexpr('"',mys2)-1
    uriNew = substr(mys2, 1, end)
    
    res=postForm(uriNew, ..., .params = .params, .opts = .opts,
                 curl = curl, style = style,
                 .encoding = .encoding, binary = binary, .checkParams = .checkParams,
                 .contentEncodeFun = .contentEncodeFun)
  } )
  
  return(res)
  
}
