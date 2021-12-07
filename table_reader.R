get_hs_info <- function (file_name){
  hs <- read.table(file_name,header=TRUE);
  number_of_loci<-ncol(hs)-2
  number_of_individual<-((nrow(hs))/2)
  tot_num_of_pairwise_comparison<-(number_of_loci*(number_of_loci-1)/2)
  table_info <- list("number_of_loci" = number_of_loci, "number_of_individual" = number_of_individual, "number_of_pairwise_comparison" = tot_num_of_pairwise_comparison)
  return(table_info)
}

get_etaik_info <- function (file_name){
  etaik<-read.table(file_name);
  number_of_population <- (ncol(etaik)-5)
  table_info <- list ("number_of_population" = number_of_population)
  return(table_info)
}

# hs table (example.str) #what does hs stands for?
# number_of_loci<-ncol(hs)-2
# number_of_individual<-(nrow(hs)-1/2)
# tot.Num.Of.Pairwise.Comparison<-function(n){
#   return(n*(n-1)/2)
#   
#   etaik (example.indivq)
#   first column confirm the number of individual (column 3)
#   figure out the number of populations after (:) this is your "K population"
#   figure out which population an individual belongs to 
#   
#   pkla (example.pkla)
