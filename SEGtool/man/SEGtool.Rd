\name{SEGtool}
\alias{SEGtool}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detecting specifically expressed gene from expression profile
}
\description{
This package aims to detect specifically expressed genes in
expression profile. Such genes are significantly
up- or down-expressed in a small number of samples. This work
performed by modified FCM Algorithm, modified jaccard coeffience,
and some criteria. Expression data can be from different tissues, organs
or any other sources that you wish to compare in expression profile.
}
\usage{
SEGtool(X=NULL, exp_cutoff = 5, multi_cpu = 4, detect_mod = 2, result_outdir = "SEGtool_result", draw_plot = TRUE, html_report = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{expression dataframe : the numeric expression values store in data.frame format
%%     ~~Describe \code{X} here~~
}
  \item{exp_cutoff}{numeric value : expression value threshold, any vector in row which is lower than threshold will not be done with the SEG analysis
%%     ~~Describe \code{exp_cutoff} here~~
}
  \item{multi_cpu}{ integer value : how many cpus will be used for analysis
%%     ~~Describe \code{multi_cpu} here~~
}
  \item{detect_mod}{ integer value : detection mode for the analysis,the bigger value the more critical for SEG detection
%%     ~~Describe \code{detect_mod} here~~
}
  \item{result_outdir}{ result output directory 
%%     ~~Describe \code{result_outdir} here~~
}
  \item{draw_plot}{ Bool logic value : draw each SEG expression plot or not, default is TRUE. This argument is needed for creating html report and will take longer time.
%%     ~~Describe \code{draw_plot} here~~
}
  \item{html_report}{ Bool logic value : generate html report, default is TRUE.
%%     ~~Describe \code{html_report} here~~
} 
}
\details{ 
    Expression data can be obtained from different tissues, organs or any other sources 
that you wish to compare in expression profile. Any matrix-like expression value 
could be accepted for this package. log2-like value is not accepted in this package, please 
transform those value to normal style.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ 
    The object return a list of SEGs result informations.
	\item{Allsummary} all summary result of input data
	\item{SEGinSample}{0, 1, -1 matrix about SEG information in samples, 1 represents high SEG gene in the sample, -1 means low SEG and 0 is non-SEG}
}
\references{
Qiu, C., Xiao, J., Yu, L., Han, L., & Iqbal, M. N. (2013). A modified interval type-2 fuzzy C-means algorithm with application in MR image segmentation. Pattern Recognition Letters, 34(12), 1329-1338.

Kadota, K., Ye, J., Nakai, Y., Terada, T., & Shimizu, K. (2006). ROKU: a novel method for identification of tissue-specific genes. BMC bioinformatics, 7(1), 294.

McLachlan, G. J., Bean, R. W., & Jones, L. B. T. (2006). A simple implementation of a normal mixture approach to differential gene expression in multiclass microarrays. Bioinformatics, 22(13), 1608-1615.

Ye Chengyin WX: TSGA: an R package for tissue specific genes analysis 2008 [http://www.cab.zju.edu.cn/ics/faculty/zhuj/software/tsga/index.htm]

}
\author{
Qiong Zhang & Chun-Jie Liu

Maintainer: zhang qiong <zhangqiong@hust.edu.cn>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{http://bioinfo.life.hust.edu.cn/SEGtool/}}
}
\examples{
library(SEGtool)
data(EbiHumanExpression)
result<-SEGtool(EbiHumanExpression, exp_cutoff = 5,multi_cpu = 4,detect_mod=2,result_outdir='SEGtool_result',draw_plot=FALSE,html_report=TRUE)

## Should be DIRECTLY executable
## Define data, use random, or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, exp_cutoff = 5, multi_cpu = 4, detect_mod = 2, result_outdir = "SEGtool_result", 
    draw_SEG_plot = FALSE, html_report = FALSE) 
{
    .makedir("SEGtool_result")
    filter_bool_list <- apply(X, 1, .filter, exp_cutoff = exp_cutoff, 
        detect_mod = detect_mod)
    Candidate_TS_G <- X[filter_bool_list, ]
    Candidate_TS_names <- rownames(Candidate_TS_G)
    Uniform_list <- X[!filter_bool_list, ]
    tmp_list <- split(Candidate_TS_G, rownames(Candidate_TS_G))
    clusterd_id_list <- mclapply(tmp_list, .FCM_cluster_id, mc.cores = multi_cpu, 
        detect_mod = detect_mod)
    matrix_cts_names <- matrix(Candidate_TS_names)
    candidate_TS_matrix <- apply(matrix_cts_names, 1, function(x) {
        .judge_TS(Candidate_TS_G[x, ], clusterd_id_list[x][[1]][[1]], 
            clusterd_id_list[x][[1]][[2]], 2, clusterd_id_list[x][[1]][[3]])
    })
    summary_TS <- .get_TS_result(candidate_TS_matrix, Candidate_TS_G, 
        result_outdir, draw_plot = FALSE)
    all_matrix <- t(candidate_TS_matrix)
    all_df <- as.data.frame(all_matrix, row.names = rownames(Candidate_TS_G))
    colnames(all_df) <- colnames(Candidate_TS_G)
    tissue_df <- as.data.frame(candidate_TS_matrix, row.names = colnames(Candidate_TS_G))
    colnames(tissue_df) <- rownames(Candidate_TS_G)
    candicate_tissue_df <- tissue_df[apply(tissue_df, 1, function(z) {
        return(any(z != 0))
    }), ]
    TS_df <- all_df[apply(all_df, 1, function(z) {
        return(any(z != 0))
    }), ]
    TS_df_names <- names(TS_df)
    TS_gene_names <- rownames(TS_df)
    TS_exp <- Candidate_TS_G[TS_gene_names, ]
    summary_all_samples <- dim(X)
    summary_result <- c(summary_all_samples, summary_TS)
    summary_result <- data.frame(matrix(summary_result, 1))
    names(summary_result) <- c("total_genes", "total_samples", 
        "total_TS_genes", "total_TS_samples", "only_high_SEG", 
        "only_low_SEG", "gene_with_high_low_SEG")
    summary_result_file <- paste(result_outdir, "summary_result.txt", 
        sep = "/")
    write.table(summary_result, file = summary_result_file, col.names = TRUE, 
        row.names = FALSE, sep = "\t", quote = FALSE)
    if (html_report) {
        .get_visualisation_html(indir = result_outdir, page.name = "SEGtool_result", 
            page.title = "specific expressed gene analysis results", 
            outdir = result_outdir, draw_heatmap = TRUE, draw_pca = TRUE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{specially }
\keyword{expressed }% __ONLY ONE__ keyword per line
