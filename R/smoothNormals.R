smoothNormals <- function(logr, lNormals)
{
    notalreadyinpanel <- !sapply(lNormals,function(x) cor(log2(unlist(lapply(x,function(y) y$records)) +1),logr))>.999
    normals <- sapply(lNormals,function(x) log2(unlist(lapply(x,function(y) y$records)) +1))
    cat(".")
    predicted <- lm(y ~ . - 1, data = data.frame(y = 2^logr-1,
                                                 X = 2^normals[, notalreadyinpanel])-1)$fitted.values
    predicted[predicted < 1] <- 1
    print(summary(log2(predicted+1)))
    logr <- logr-log2(predicted+1)
    logr
}
