
//#include <Eigen/SVD>
#include <algorithm> // std::sort, std::stable_sort
#include <itkDataObject.h>
#include <itkImageRegionIterator.h>
#include <itkIndex.h>
#include <itkVariableLengthVector.h>
#include <itkVectorImage.h>
#include <m2CoreCommon.h>
#include <m2PcaImageFilter.h>
#include <m2Timer.h>
#include <mitkIOUtil.h>
#include <mitkImage.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <numeric> // std::iota
#include <vnl/vnl_matrix.h>
// #include <boost/progress.hpp>

void m2::PcaImageFilter::initMatrix()
{
  // this->GetValidIndices();
  auto input = this->GetIndexedInputs();
  auto mitkImage = dynamic_cast<mitk::Image *>(input.front().GetPointer());
  size_t pixels = 1;
  for (unsigned int i = 0; i < mitkImage->GetDimension(); ++i)
    pixels *= mitkImage->GetDimensions()[i];
  const unsigned long numberOfrow = pixels;
  const unsigned long numberOfcolumn = this->GetIndexedInputs().size();

  this->m_DataMatrix.resize(numberOfrow, numberOfcolumn);
  unsigned int c = 0;
  // boost::progress_display p(numberOfcolumn);
  
  /*Fill matrix with image values one column includes values of one image*/
  for (auto it = input.begin(); it != input.end(); ++it, ++c)
  {
    mitkImage = dynamic_cast<mitk::Image *>(it->GetPointer());
    mitk::ImagePixelReadAccessor<m2::DisplayImagePixelType, 3> access(mitkImage);
    std::copy(access.GetData(), access.GetData() + pixels, m_DataMatrix.col(c).data());
    // ++p;
  }

  // auto maxCoeffs = m_DataMatrix.colwise().maxCoeff();
  // auto minCoeffs = m_DataMatrix.colwise().minCoeff();
  // auto scalingFactors = maxCoeffs - minCoeffs;
  // m_DataMatrix = ((m_DataMatrix.rowwise() - minCoeffs).array().rowwise() / scalingFactors.array()).matrix();

  // auto means = m_DataMatrix.colwise().mean();
  // m_DataMatrix = m_DataMatrix.rowwise() - means;
  // m_DataMatrix /= (m_DataMatrix.rows() - 1);

  // for (unsigned int c = 0; c < m_DataMatrix.cols(); ++c)
  // m_DataMatrix.col(c) = m_DataMatrix.col(c) / stdDevs[c];
}

Eigen::MatrixXf m2::PcaImageFilter::GetEigenImageMatrix(){
  return m_EigenImageMatrix;
}

Eigen::VectorXf m2::PcaImageFilter::GetMeanImage(){
  return m_MeanImage;
}

void m2::PcaImageFilter::GenerateData()
{
  auto timer = m2::Timer("PCA - Generate data ...");
  this->initMatrix();

  // eigenionimages
  m_MeanImage = m_DataMatrix.rowwise().mean();
  Eigen::MatrixXf eigenionData = m_DataMatrix.colwise() - m_MeanImage;

  Eigen::JacobiSVD<Eigen::MatrixXf> svd(eigenionData, Eigen::ComputeThinU | Eigen::ComputeThinV);
  MITK_INFO << "S size: " << svd.singularValues().rows() << " " << svd.singularValues().cols();
  MITK_INFO << "U size: " << svd.matrixU().rows() << " " << svd.matrixU().cols();
  MITK_INFO << "V size: " << svd.matrixV().rows() << " " << svd.matrixV().cols();
  MITK_INFO << "m_DataMatrix: " << m_DataMatrix.rows() << " " << m_DataMatrix.cols();

  const Eigen::MatrixXf &U = svd.matrixU();
  m_EigenImageMatrix = U;

  auto eigenIonVectorImage = initializeItkVectorImage(m_NumberOfComponents);
  m2::DisplayImagePixelType *data;
  data = eigenIonVectorImage->GetBufferPointer();

  for (unsigned int c = 0; c < m_NumberOfComponents; ++c)
  {
    const auto &col = U.col(c);
    for (unsigned int p = 0; p < U.rows(); ++p)
    {
      data[p * m_NumberOfComponents + c] = col(p);
    }
  }

  // feature pca
  // Eigen::MatrixXf pcaData = m_DataMatrix.rowwise() - m_DataMatrix.colwise().mean();
  // Eigen::MatrixXf cov = (pcaData.transpose() * pcaData) / (pcaData.rows() - 1);

  // Eigen::EigenSolver<Eigen::MatrixXf> solver(cov);
  // Eigen::VectorXf values = solver.eigenvalues().real();
  // Eigen::MatrixXf vectors = solver.eigenvectors().real();
  // Eigen::MatrixXf vectorsSorted(vectors);
  // // fill indices
  // std::vector<unsigned int> indices(values.size());
  // std::iota(indices.begin(), indices.end(), 0);
  // // sort indices according to
  // std::stable_sort(indices.begin(), indices.end(), [&values](auto i1, auto i2) { return values[i1] > values[i2]; });

  // unsigned int i = 0;
  // for (auto u : indices)
  //   vectorsSorted.col(i++) = vectors.col(u);

  // Eigen::MatrixXf pc = pcaData * vectorsSorted;
  // auto pcVectorImage = initializeItkVectorImage(m_NumberOfComponents);
  // data = pcVectorImage->GetBufferPointer();

  // for (unsigned int c = 0; c < m_NumberOfComponents; ++c)
  // {
  //   const auto &col = pc.col(c);
  //   for (unsigned int p = 0; p < pc.rows(); ++p)
  //   {
  //     data[p * m_NumberOfComponents + c] = col(p);
  //   }
  // }
  

  mitk::Image::Pointer eigenIonImage = this->GetOutput(0);
  mitk::CastToMitkImage(eigenIonVectorImage, eigenIonImage);
  eigenIonImage->SetSpacing(this->GetInput()->GetGeometry()->GetSpacing());
  eigenIonImage->SetOrigin(this->GetInput()->GetGeometry()->GetOrigin());

  // mitk::Image::Pointer pcImage = this->GetOutput(1);
  // mitk::CastToMitkImage(pcVectorImage, pcImage);
  // pcImage->SetSpacing(this->GetInput()->GetGeometry()->GetSpacing());
  // pcImage->SetOrigin(this->GetInput()->GetGeometry()->GetOrigin());
 
}
