barNaming <- function(vec) {
retVec <- vec
for (k in 2:length(vec)) {
if (vec[k - 1] == vec[k])
retVec[k] <- ""
}
return(retVec)
}