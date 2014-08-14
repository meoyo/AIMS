\name{applyAIMS}
\alias{applyAIMS}
\title{
Function that will assign the AIMS subtype to a gene expression matrix
}
\description{
Given a gene expression matrix D where rows correspond to genes and
columns to samples and a list of Entrez gene ids, this function will
assign the breast cancer molecular subtype using Absolute PAM50 (AIMS)
(Paquet et al. under review at JNCI).
}
\usage{
applyAIMS (eset, EntrezID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
  An ExpressionSet (Biobase) or a  matrix D of gene expression. Rows are genes and columns are samples (tumor from breast cancer patients).
}
  \item{EntrezID}{
  A character vector corresponding to the Entrez Ids of the genes of matrix D.
}
}

\details{
  We defined Absolute assignment of breast cancer Molecular Intrinsic Subtypes (AIMS) to stabilize the actual PAM50. The idea of the approach is to use simple rules of the form "if gene A expression < gene B expression assigned sample to subtype X". By using those simple rules we are able to assign subtype using only the expression values of one patient. We have shown AIMS recapitulates PAM50 subtype assignments and preserved the prognostic values of the subtypes.
  This function will return the subtype assignment as well as the posterior probabilities for all the subtypes.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{cl}{Subtypes identified by the AIMS. It could be either "Basal", "Her2", "LumA", "LumB" or "Normal".}
\item{prob}{A vector corresponding to the posterior probabilities of the subtypes in cl.}
\item{all.probs}{A matrix of all the posterior probabilities for all the samples and all the subtypes}
\item{rules.matrix}{The matrix of 100 rules used by AIMS to assign the breast cancer subtypes. Rows correspond to rules and columns to samples. This is a 0,1 matrix in which 1 represents the rules is true and 0 the rule is false}
\item{data.used}{The expression values used to evaluate the simple rules.}
\item{EntrezID.used}{The list of Entrez ids used by AIMS}
}
\references{
Parker, Joel S. and Mullins, Michael and Cheang, Maggie C.U. and Leung, Samuel and Voduc, David and Vickery, Tammi and Davies, Sherri and Fauron, Christiane and He, Xiaping and Hu, Zhiyuan and Quackenbush, John F. and Stijleman, Inge J. and Palazzo, Juan and Marron, J.S. and Nobel, Andrew B. and Mardis, Elaine and Nielsen, Torsten O. and Ellis, Matthew J. and Perou, Charles M. and Bernard, Philip S. (2009) "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes", \emph{Journal of Clinical Oncology}, \bold{27}(8):1160--1167

Donald Geman, Christian d'Avignon, Daniel Q. Naiman and Raimond L. Winslow (2004) "Classifying Gene Expression Profiles from Pairwise mRNA Comparisons", \emph{Stat Appl Genet Mol Biol.}, \bold{3}: Article19.
}
\author{
Eric R. Paquet (eric.r.paquet@gmail.com)
}

\seealso{
\code{\link[AIMS]{AIMSmodel}}, \code{\link[AIMS]{mcgillExample}}
}
\examples{
## Load the McGill dataset used in the paper
data(mcgillExample)

## Convert the expression matrix to an ExpressionSet.
## could also send directly an expresion matrix
mcgillExample$D <- ExpressionSet(assayData=mcgillExample$D)

## Assigne AIMS on McGill dataset
mcgill.AIMS.subtypes.batch <- applyAIMS (mcgillExample$D,
                                         mcgillExample$EntrezID)

## Print a summary of all the subtype in the dataset
table(mcgill.AIMS.subtypes.batch$cl)

## We can do the samething using only one sample
mcgill.AIMS.subtypes.first <- applyAIMS (mcgillExample$D[,1,drop=FALSE], 
                                         mcgillExample$EntrezID)
table(mcgill.AIMS.subtypes.first$cl)

if (mcgill.AIMS.subtypes.batch$cl[1] == mcgill.AIMS.subtypes.first$cl[1]){
   message("Identical assignment batch and first sample")
}else{
   message("Different assignment batch and first sample")
}

## We can do the samething for the first 20 samples
mcgill.AIMS.subtypes.first20 <- applyAIMS (mcgillExample$D[,1:20,drop=FALSE], 
                                           mcgillExample$EntrezID)
table(mcgill.AIMS.subtypes.first20$cl)

if (all(mcgill.AIMS.subtypes.batch$cl[1:20] == mcgill.AIMS.subtypes.first20$cl)){
   message("Identical assignment batch and first 20")
}else{
   message("Different assignment batch and first 20")
}

## We can do the samething using 50 randomly selected samples, no
## set.seed needed
sel.rand.50 = sample(1:ncol(mcgillExample$D),50)
mcgill.AIMS.subtypes.rand50 <- applyAIMS (mcgillExample$D[,sel.rand.50,drop=FALSE], 
                                          mcgillExample$EntrezID)
table(mcgill.AIMS.subtypes.rand50$cl)

if (all(mcgill.AIMS.subtypes.batch$cl[sel.rand.50] == mcgill.AIMS.subtypes.rand50$cl)){
   message("Identical assignment batch and random 50")
}else{
   message("Different assignment batch and random 50")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }
