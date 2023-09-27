combineDiploid <- function(lNormals)
{
    nlNormals <- lNormals[[1]]
    if(length(lNormals)>1)
        for(i in 2:length(lNormals))
        {
            for(j in 1:length(lNormals[[1]]))
                nlNormals[[j]]$records <- nlNormals[[j]]$records+lNormals[[i]][[j]]$records
        }
    nlNormals
}
