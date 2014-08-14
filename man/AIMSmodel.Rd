\name{AIMSmodel}
\alias{AIMSmodel}
\title{
The AIMS model
}
\description{
This is the model definition for AIMS. It contains the naive bayes
classifier composed of the 100 rules described in Paquet et al. "Absolute
assignment of breast cancer intrinsic molecular subtype" (under review at
JNCI).
}
\usage{
AIMSmodel
}
%- maybe also 'usage' for other objects documented here.

\details{
  This is the AIMS model define using 100 simple rules of the form gene
  A < gene B and combine within a naive bayes classifier within
  e1071. (Paquet et al. under review JNCI)

  Briefly, using a suitably large training dataset(~5000 gene breast
  cancer gene expression profiles), the approach
  identifies a small set of simple binary rules (~20) that examine the raw
  expression measurements for pairs of genes from a single breast cancer patient, and
  only that patient. The binary rules are of the form "if the expression
  of gene x is greater than gene y, then tend to assign subtype z for that
  patient". Subtypes could be : Basal, Her2, LumA, LumB, or Normal.
  The collection of binary rules is combined for a single
  estimation of a patient subtype via a single probabilistic model using
  naiveBayes in e1071. In
  this way, since only expression levels of genes with a single patient is
  considered, the method represents a promising approach to ablate the
  instability caused by relativistic approaches (Paquet et al. in review
  at JNCI).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{all.pairs}{The 100 rules in AIMS in the form EntrezID gene A < EntrezID gene B}
\item{k}{The selected number of optimal rules. For AIMS we have shown it is 20.}
\item{one.vs.all.tsp}{The Naive bayes classifier used in combination with the 100 rules}
\item{selected.pairs.list}{The list of rules sorted from the best discriminating rule to the least discriminating rules subdivided by subtype.}
}
\author{
Eric R. Paquet (eric.r.paquet@gmail.com)
}

\seealso{
\code{\link[AIMS]{applyAIMS}}, \code{\link[AIMS]{mcgillExample}}, \code{\link[e1071]{naiveBayes}}
}
\examples{
## Load a sample of the McGill dataset used in the paper
data(AIMSmodel)

## list the top-scoring rules for the individual subtypes
AIMSmodel$selected.pairs.list

## List the posterior probability tables for the 100-rules
AIMSmodel$one.vs.all.tsp[[20]]$table

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
