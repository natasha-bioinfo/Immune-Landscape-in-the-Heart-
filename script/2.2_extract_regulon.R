type="DC"

#load packages
library(jsonlite)

#read json file
input_json=paste0(type,".json")
regulons=read_json(input_json, simplify=T)

save(regulons,file=paste0(type,"_all_regulon.rds"))
