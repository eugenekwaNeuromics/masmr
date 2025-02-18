\name{establishCellSegModel}
\alias{establishCellSegModel}
\title{establishCellSegModel}
\description{
Prepares cell segmentation model.
}
\usage{
establishCellSegModel(
    model = 'cellpose',
    gpu = F,
    cellposePreTrainedModel = NULL,
    cellposeModelType = 'nuclei',
    stardistModelType = '2D_versatile_fluo',
    params = get('params', envir = globalenv())
)
}
\arguments{
\item{  model  }{ Character. Choice of cell segmentation model. Currently accepts \code{cellpose} or \code{stardist}. }
\item{  gpu  }{ Boolean. Determines whether GPU to be used for chosen cell segmentation model. }
\item{  cellposePreTrainedModel  }{ Character. Location of custom pre-trained model for Cellpose. If unspecified, packaged pre-trained models loaded instead. }
\item{  cellposeModelType  }{  Character. One of the default packaged pre-trained models for Cellpose: see \code{model_type} in Cellpose documentation for details. }
\item{  stardistModelType  }{  Character. One of the default packaged pre-trained models for Stardist: see \code{from_pretrained} in Stardist documentation for details. }
\item{ params }{ The parameters environment. }
}
\value{
  Loads models to the \code{cellSeg} environment.
}
\details{
}
\seealso{
\link{ Cellpose: https://cellpose.readthedocs.io/en/latest/ }
\link{ Stardist: https://stardist.net/ }
}
\examples{
}
