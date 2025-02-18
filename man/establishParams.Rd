\name{establishParams}
\alias{establishParams}
\title{establishParams}
\description{
Establishes parameters.
}
\usage{
establishParams(
    ## User to fill in
    imageDirs,
    imageFormat,
    outDir,
    codebookFileName,
    dapiDir,

    ## Will attempt automatic fill
    metaFormat,
    dapiFormat,
    pythonLocation,
    codebookColumnNames,

    ## Automatically filled in
    seed = 12345,
    resumeMode = T,

    ## Optional
    fpkmFileName = NULL,
    nProbesFileName = NULL
)
}
\arguments{
\item{ imageDirs }{ Character vector. Containing directory/ies that contain image data. If no images found in immediate directory but subdirectories with images present, subdirectories are automatically referenced instead. }
\item{ imageFormat }{ Character. Image file format - should be readable with BioFormats (e.g. '\code{.ome.tif}'), or '\code{.dax}' (binary).}
\item{ outDir }{ Character. Directory to output data in. If folder does not exist, automatically created. }
\item{ codebookFileName }{ Character. File that contains codebook. }
\item{ dapiDir }{ Character or character vector. Directory/ies containing DAPI images. If unspecified, will presently skip -- but will flag an error during nuclei segmentation.}


\item{ metaFormat }{ Character. Meta file format. Likely \code{.ome.tif} for \code{.ome.tif}, and \code{.xml} for \code{.dax}. If unspecified, attempt automatic fill-in. }
\item{ dapiFormat }{ Character. DAPI image file format. Assumed to be the same as \code{imageFormat}. If unspecified, attempt automatic fill-in. }
\item{ pythonLocation }{ Character. Location of Python environment. Used for nuclei segmentation. If unspecified, defaults to \code{reticulate::miniconda_path()}. }
\item{ codebookColumnNames }{ Character vector. Column names (specifying channel and cycle) for codebook. If unspecified, a placeholder will be provided -- that should be corrected upon \code{readImageMetaData}. }

\item{ seed }{ Integer. Seed for pseudo-random number generation. }
\item{ resumeMode }{ Boolean (default \code{TRUE}). If \code{TRUE}, resumes from last entry (should algorithm be interrupted). Otherwise, overwrites existing files. }

\item{ fpkmFileName }{ (Optional) Character. File that contains FPKM data. }
\item{ nProbesFileName }{ (Optional) Character. File that contains number of probes per gene. }
}
\value{
  An environment object containing user input parameters, to be used downstream.
}
\seealso{
\link{readImageMetaData}
}
\examples{
data_dir = 'C:/DATA/'
establishParams(
  imageDirs = paste0(data_dir, 'IMAGES/'),
  outDir = 'C:/OUT/'
  codebookFileName = paste0(data_dir, 'CODEBOOK.csv'),
  codebookColumnNames = paste0(rep(c(0:3), 4), '_', rep(c('cy3', 'alexa594', 'cy5', 'cy7'), each=4) ),
  imageFormat = '.dax'
)

## Params can be updated as needed:
# params$resumeMode <- F
}
