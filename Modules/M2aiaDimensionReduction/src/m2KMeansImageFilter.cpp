/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <m2SpectrumImage.h>
#include <m2KMeansImageFilter.h>
#include <mitkImage.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkImageCast.h>
#include <mitkLabelSetImage.h>


// OpenMP
#include <omp.h>


void m2::KMeansImageFilter::GenerateData(){

    auto image = dynamic_cast<const m2::ImzMLSpectrumImage *>(m_Input.GetPointer());

    if (!image->GetImageAccessInitialized())
        MITK_INFO << "Image access not initialized";
    

    if(!m_Intervals.size())
        MITK_INFO << "Intervals are not set";
    
    std::vector<itk::Index<3>> validIndices;
    auto maskImage = image->GetMaskImage();
    mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType, 3> maskAcc(maskImage->Clone());

    for(auto s : image->GetSpectra())
       if(maskAcc.GetPixelByIndex(s.index)!=0)
           validIndices.push_back(s.index);
       
    
        

    Eigen::MatrixXd data(validIndices.size(), m_Intervals.size());   
    mitk::Image::Pointer ionImage = image->mitk::Image::Clone();
    
    for (size_t row = 0; row < m_Intervals.size(); ++row)
    {
      const auto mz = m_Intervals.at(row).x.mean();
      image->GetImage(mz, image->ApplyTolerance(mz), image->GetMaskImage(), ionImage);
      unsigned v = 0;
      for (auto s : validIndices)
      {
        mitk::ImagePixelReadAccessor<m2::DisplayImagePixelType, 3> acc(ionImage);
        data(v, row) = acc.GetPixelByIndex(s);
        ++v;
      }
    }
    
    std::vector<int> clusterAssignments(data.rows(), -1);
    DoKMeans(data, m_NumberOfClusters, clusterAssignments);
    MITK_DEBUG << "KMeans done";
    
    
    auto clusteredImage = this->GetOutput(0);
    MITK_DEBUG << "Clustered image initialized";
    clusteredImage->Initialize(maskImage);
    MITK_DEBUG << "Clustered image initialized";
    // clusteredImage->AddLayer(mitk::LabelSet::New());
    // auto labelSet = clusteredImage->GetLabelSet(0);
    for(unsigned long i = 0; i < m_NumberOfClusters; ++i){
        mitk::Label::Pointer label = mitk::Label::New();
        label->SetName("Cluster " + std::to_string(i));
        label->SetValue(i);
        label->SetOpacity(1);
        label->SetColor(m2::RandomColor());
        clusteredImage->AddLabel(label,0);
    }

    clusteredImage->SetActiveLayer(0);
    MITK_DEBUG << "Clustered image initialized";
    {
    mitk::ImagePixelWriteAccessor<mitk::LabelSetImage::PixelType, 3> c_acc(clusteredImage);

    unsigned long v = 0;
    for(auto s: validIndices)
    {
        if( maskAcc.GetPixelByIndex(s) != 0){
            c_acc.SetPixelByIndex(s, clusterAssignments[v] + 1);
            ++v;
        }
    }
    MITK_DEBUG << "Clustered image filled";
    }
}


// Function to calculate Euclidean distance between two points
double euclideanDistance(const Eigen::VectorXd& point1, const Eigen::VectorXd& point2) {
    return (point1 - point2).norm();
}


void m2::KMeansImageFilter::DoKMeans(const Eigen::MatrixXd& data, int k, std::vector<int>& clusterAssignments){
    // Randomly initialize centroids
    std::vector<Eigen::VectorXd> centroids(k);
    for (int i = 0; i < k; ++i) {
        centroids[i] = data.row(rand() % data.rows());
    }

    bool centroidsChanged = true;
    while (centroidsChanged) {
        centroidsChanged = false;

        // Assign each point to the nearest centroid
#pragma omp parallel for
        for (int i = 0; i < data.rows(); ++i) {
            double minDistance = std::numeric_limits<double>::max();
            int closestCentroid = -1;
            for (int j = 0; j < k; ++j) {
                double distance = euclideanDistance(data.row(i), centroids[j]);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestCentroid = j;
                }
            }
#pragma omp critical
            {
                if (clusterAssignments[i] != closestCentroid) {
                    clusterAssignments[i] = closestCentroid;
                    centroidsChanged = true;
                }
            }
        }

        // Update centroids
#pragma omp parallel for
        for (int i = 0; i < k; ++i) {
            Eigen::VectorXd newCentroid = Eigen::VectorXd::Zero(data.cols());
            int count = 0;
            for (int j = 0; j < data.rows(); ++j) {
                if (clusterAssignments[j] == i) {
                    newCentroid += data.row(j);
                    count++;
                }
            }
            if (count > 0) {
                newCentroid /= count;
            }
#pragma omp critical
            {
                centroids[i] = newCentroid;
            }
        }
    }
}

