\name{SEGtool-package}
\alias{SEGtool-package}
\alias{SEGtool}
\docType{package}
\title{
Detecting specifically expressed gene in expression profile
}
\description{
This package aims to detect specifically expressed genes in
expression profile. Such genes are significantly
up- or down-expressed in a small number of samples. This work
performed by modified FCM Algorithm, modified jaccard coeffience,
and some criteria. Expression data can be from different tissues, organs
or any other sources that you wish to compare in expression profile.
}
\details{
\tabular{ll}{
Package: \tab SEGtool\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-01-19\cr
License: \tab GPL (>= 2.0) \cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Qiong Zhang & Chun-Jie Liu

Maintainer: zhang qiong <zhangqiong@hust.edu.cn>

}
\references{
Qiu, C., Xiao, J., Yu, L., Han, L., & Iqbal, M. N. (2013). A modified interval type-2 fuzzy C-means algorithm with application in MR image segmentation. Pattern Recognition Letters, 34(12), 1329-1338.

Kadota, K., Ye, J., Nakai, Y., Terada, T., & Shimizu, K. (2006). ROKU: a novel method for identification of tissue-specific genes. BMC bioinformatics, 7(1), 294.

McLachlan, G. J., Bean, R. W., & Jones, L. B. T. (2006). A simple implementation of a normal mixture approach to differential gene expression in multiclass microarrays. Bioinformatics, 22(13), 1608-1615.
}
\keyword{ Tissue specific expression }
\seealso{
http://bioinfo.life.hust.edu.cn/SEGtool/
\code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(SEGtool)
data(EbiHumanExpression)
SEGtool(EbiHumanExpression, exp_cutoff = 5,multi_cpu = 4,detect_mod=2,result_outdir='SEGtool_result',draw_SEG_plot=FALSE,html_report=TRUE)
}
