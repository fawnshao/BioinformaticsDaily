setwd("~/Data/myworkData/WuXi/GTEx.imprint/SD-ASM/")
library(jsonlite)
# library(data.table)
# library(tidyverse)
# library(dplyr)
# library(purrr)
# mytable <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.jsonld", flatten = TRUE)
# mytable.def <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.jsonld")
# mytable.fmt <- mytable.def %>% purrr::flatten() %>% map_if(is_list, as_tibble) %>% map_if(is_tibble, list) %>% bind_cols()
# mytable.df <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.jsonld") %>% as.data.frame
# write.table(mytable, file = "AllelicEpigenome-sigOnly-AllDocs.jsonld.tsv", sep = "\t")
mytable <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.jsonld")
# mytable.df <- apply(mytable, 2, as.character)
mytable.df <- apply(jsonlite::flatten(mytable, recursive = T), 2, as.character)
# output$Title <- vapply(output$Title, paste, collapse = ", ", character(1L))
write.table(mytable.df, file = "AllelicEpigenome-sigOnly-AllDocs.jsonld.tsv", sep = "\t")
# 
# 
# mytable.def <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.jsonld")
# mytable.def.fmt <- lapply(mytable.def, function(x) {
#     x[sapply(x, is.null)] <- NA
#     unlist(x)
# })
# output <- do.call("rbind", mytable.def.fmt)
# write.table(t(output), file = "AllelicEpigenome-sigOnly-AllDocs.jsonld.tsv", sep = "\t")
# 
# mytable.def <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.json")
# mytable.def.fmt <- lapply(mytable.def, function(x) {
#     x[sapply(x, is.null)] <- NA
#     unlist(x)
# })
# output <- do.call("rbind", mytable.def.fmt)
# write.table(t(output), file = "AllelicEpigenome-sigOnly-AllDocs.json.tsv", sep = "\t")

# mytable <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.json", flatten = TRUE)
# mytable.def <- fromJSON("AllelicEpigenome-sigOnly-AllDocs.json")
# str(mytable)
# str(flatten(mytable))
# str(flatten(mytable, recursive = FALSE))
# 
# str(mytable.def)
# str(flatten(mytable.def, recursive = FALSE))
# mytable.fmt <- jsonlite::flatten(mytable.def, recursive = TRUE) %>% purrr::flatten() %>% 
#     purrr::map_if(is_list, as_tibble) %>% purrr::map_if(is_tibble, list) %>% dplyr::bind_cols()
# 
# mytable.fmt <- mytable.def %>%
#     
#     # remove classification level
#     purrr::flatten() %>%
#     
#     # turn nested lists into dataframes
#     purrr::map_if(is_list, as_tibble) %>%
#     
#     # bind_cols needs tibbles to be in lists
#     purrr::map_if(is_tibble, list) %>%
#     
#     # creates nested dataframe
#     dplyr::bind_cols()
# # write.table(mytable.fmt, file = "AllelicEpigenome-sigOnly-AllDocs.json.tsv1", sep = "\t")
# mytable.fmt <- as.data.table(mytable.fmt)
# fwrite(mytable.fmt, file = "AllelicEpigenome-sigOnly-AllDocs.json.tsv1", sep = "\t")
# # > for(i in 1:42){print(c(i,typeof(mytable[,i])))}
# # [1] "1"         "character"
# # [1] "2"         "character"
# # [1] "3"         "character"
# # [1] "4"         "character"
# # [1] "5"         "character"
# # [1] "6"         "character"
# # [1] "7"         "character"
# # [1] "8"         "character"
# # [1] "9"         "character"
# # [1] "10"        "character"
# # [1] "11"        "character"
# # [1] "12"        "character"
# # [1] "13"        "character"
# # [1] "14"        "character"
# # [1] "15"        "character"
# # [1] "16"        "character"
# # [1] "17"        "character"
# # [1] "18"        "character"
# # [1] "19"        "character"
# # [1] "20"        "character"
# # [1] "21"        "character"
# # [1] "22"        "character"
# # [1] "23"        "character"
# # [1] "24"        "character"
# # [1] "25"        "character"
# # [1] "26"        "character"
# # [1] "27"        "character"
# # [1] "28"        "character"
# # [1] "29"   "list"
# # [1] "30"      "integer"
# # [1] "31"   "list"
# # [1] "32"      "integer"
# # [1] "33"   "list"
# # [1] "34"      "integer"
# # [1] "35"        "character"
# # [1] "36"   "list"
# # [1] "37"      "integer"
# # [1] "38"        "character"
# # [1] "39"        "character"
# # [1] "40"        "character"
# # [1] "41"   "list"
# # [1] "42"      "integer"
# 
# mytable.df <- apply(mytable, 2, as.character)
# # output$Title <- vapply(output$Title, paste, collapse = ", ", character(1L))
# write.table(mytable, file = "AllelicEpigenome-sigOnly-AllDocs.json.tsv", sep = "\t")
# 
# data2 <- fromJSON("https://api.github.com/users/hadley/repos")
# data1 <- fromJSON("https://api.github.com/users/hadley/repos", simplifyMatrix = c("name","full_name"))


# 
# library(jsonlite)
# library(tidyverse)
# 
# your_df <- mytable %>%
#     
#     # remove classification level
#     purrr::flatten() %>%
#     
#     # turn nested lists into dataframes
#     map_if(is_list, as_tibble) %>%
#     
#     # bind_cols needs tibbles to be in lists
#     map_if(is_tibble, list) %>%
#     
#     # creates nested dataframe
#     bind_cols()
