/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <m2KMeansImageFilter.h>
#include <m2SpectrumImage.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>

// OpenMP
#include <omp.h>

void m2::KMeansImageFilter::GenerateData()
{
  m_ValidIndicesMap.clear();
  m_Outputs.clear();
  MITK_INFO << "Start KMeansImageFilter ....";
  
  if (!m_Intervals.size()){
    MITK_INFO << "Intervals are not set";
    return;
  }

  for (auto [imageId, image] : m_Inputs)
  {
    auto spectrumImage = dynamic_cast<m2::ImzMLSpectrumImage *>(image.GetPointer());
    if (!spectrumImage->GetImageAccessInitialized())
      MITK_INFO << "Image access not initialized";

    if (!spectrumImage->GetMaskImage())
      MITK_INFO << "Mask image not set";


    std::vector<itk::Index<3>> validIndices;
    auto maskImage = spectrumImage->GetMaskImage();
    mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType, 3> maskAcc(maskImage);

    for (auto s : spectrumImage->GetSpectra())
      if (maskAcc.GetPixelByIndex(s.index) != 0)
        validIndices.push_back(s.index);

    m_ValidIndicesMap[imageId] = validIndices;
    MITK_INFO << "Image found with " << validIndices.size() << " valid indices";
  }

  auto N = std::accumulate(m_ValidIndicesMap.begin(),
                           m_ValidIndicesMap.end(),
                           0,
                           [](auto &acc, auto &pair) { return acc + pair.second.size(); });
  MITK_INFO << "N: " << N;

  Eigen::MatrixXd data(N, m_Intervals.size());
  

  MITK_INFO << "Start filling data matrix ....";
  size_t offset = 0;
  for (auto [imageId, image] : m_Inputs)
  {
    auto validIndices = m_ValidIndicesMap[imageId];
    auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(image.GetPointer());

    // fill the data matrix
    mitk::Image::Pointer ionImage = spectrumImage->mitk::Image::Clone();
    for (size_t col = 0; col < m_Intervals.size(); ++col)
    {
      const auto mz = m_Intervals.at(col).x.mean();
      spectrumImage->GetImage(mz, spectrumImage->ApplyTolerance(mz), spectrumImage->GetMaskImage(), ionImage);

      size_t v = offset;
      for (auto index : validIndices)
      {
        mitk::ImagePixelReadAccessor<m2::DisplayImagePixelType, 3> acc(ionImage);
        data(v++, col) = acc.GetPixelByIndex(index);
      }
    }
    offset += validIndices.size();
  }

  MITK_INFO << "Start clustering ....";
  std::vector<int> clusterAssignments(data.rows(), -1);
  DoKMeans(data, m_NumberOfClusters, clusterAssignments);

  auto classLabelIterator = clusterAssignments.begin();
  for (auto [imageId, image] : m_Inputs)
  {
    auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(image.GetPointer());

    auto clusteredImage = dynamic_cast<mitk::LabelSetImage *>(this->GetOutput(imageId).GetPointer());
    auto maskImage = spectrumImage->GetMaskImage();
    clusteredImage->Initialize(maskImage);
    
    {
      mitk::ImagePixelWriteAccessor<mitk::LabelSetImage::PixelType, 3> c_acc(clusteredImage);
      for (auto s : m_ValidIndicesMap[imageId])
        c_acc.SetPixelByIndex(s, *classLabelIterator++ + 1);
    }

    clusteredImage->InitializeByLabeledImage(clusteredImage->Clone());

    //   if(imageId > 0 && m_Outputs.find(0) != m_Outputs.end()){
    //     auto refClusterImage = dynamic_cast<mitk::LabelSetImage *>(this->GetOutput(0).GetPointer());
    //     auto refClusterImageConst = const_cast<const mitk::LabelSetImage *>(refClusterImage);
    //     mitk::LabelSetImage::ConstLabelVectorType layer = refClusterImageConst->GetLabels();
    //     clusteredImage->AddLayer(layer);

    //   }

    //   for (unsigned long i = 0; i < m_NumberOfClusters; ++i)
    //   {
    //     mitk::Label::Pointer label = mitk::Label::New();
    //     label->SetName("Cluster " + std::to_string(i));
    //     label->SetValue(i);
    //     label->SetOpacity(1);
    //     label->SetColor(m2::RandomColor());
    //     clusteredImage->AddLabel(label, 0);
    //   }

    //   clusteredImage->SetActiveLayer(0);
  }
}

// Function to calculate Euclidean distance between two points
double euclideanDistance(const Eigen::VectorXd &point1, const Eigen::VectorXd &point2)
{
  return (point1 - point2).norm();
}

// Function to calculate correlation distance between two points
double correlationDistance(const Eigen::VectorXd &point1, const Eigen::VectorXd &point2)
{
  double mean1 = point1.mean();
  double mean2 = point2.mean();
  double numerator = (point1.array() - mean1).matrix().dot((point2.array() - mean2).matrix());
  double denominator = std::sqrt((point1.array() - mean1).square().sum() * (point2.array() - mean2).square().sum());
  return 1.0 - (numerator / denominator);
}

// Function to calculate cosine similarity between two points
double cosineSimilarity(const Eigen::VectorXd &point1, const Eigen::VectorXd &point2)
{
  double dotProduct = point1.dot(point2);
  double norm1 = point1.norm();
  double norm2 = point2.norm();
  return dotProduct / (norm1 * norm2);
}

void m2::KMeansImageFilter::DoKMeans(const Eigen::MatrixXd &data, int k, std::vector<int> &clusterAssignments)
{
  // Randomly initialize centroids

  MITK_INFO << "Start KMeans: rows: " << data.rows() << " cols: " << data.cols();

  std::vector<Eigen::VectorXd> centroids(k);
  for (int i = 0; i < k; ++i)
  {
    centroids[i] = data.row(rand() % data.rows());
  }

  bool centroidsChanged = true;
  while (centroidsChanged)
  {
    centroidsChanged = false;

    // Assign each point to the nearest centroid
// #pragma omp parallel for
    for (int i = 0; i < data.rows(); ++i)
    {
      double minDistance = std::numeric_limits<double>::max();
      int closestCentroid = -1;
      for (int j = 0; j < k; ++j)
      {
        double distance = euclideanDistance(data.row(i), centroids[j]);
        if (distance < minDistance)
        {
          minDistance = distance;
          closestCentroid = j;
        }
      }
// #pragma omp critical
      {
        if (clusterAssignments[i] != closestCentroid)
        {
          clusterAssignments[i] = closestCentroid;
          centroidsChanged = true;
        }
      }
    }

    // Update centroids
// #pragma omp parallel for
    for (int i = 0; i < k; ++i)
    {
      Eigen::VectorXd newCentroid = Eigen::VectorXd::Zero(data.cols());
      int count = 0;
      for (int j = 0; j < data.rows(); ++j)
      {
        if (clusterAssignments[j] == i)
        {
          newCentroid += data.row(j);
          count++;
        }
      }
      if (count > 0)
      {
        newCentroid /= count;
      }
// #pragma omp critical
      {
        centroids[i] = newCentroid;
      }
    }
  }
}
