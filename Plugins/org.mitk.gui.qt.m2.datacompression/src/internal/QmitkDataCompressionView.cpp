/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDataCompressionView.h"
#include <QFileDialog>

#include <QMessageBox>
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// m2
#include <m2CoreCommon.h>
#include <m2IntervalVector.h>
#include <m2MultiSliceFilter.h>
#include <m2PcaImageFilter.h>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageHelper.h>
#include <m2TSNEImageFilter.h>
#include <m2ImzMLImageIO.h>
#include <m2SpectrumImageHelper.h>

// mitk
#include <mitkDockerHelper.h>
#include <mitkImage.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateFunction.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateProperty.h>
#include <mitkProgressBar.h>
#include <usModuleRegistry.h>

// itk
#include <itkBinaryThresholdImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkShrinkImageFilter.h>
#include <itkVectorImageToImageAdaptor.h>
#include <itkVectorImage.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkImageKmeansModelEstimator.h>
#include <itkImageRegionIterator.h>

// boost
#include <berryPlatformUI.h>
#include <boost/algorithm/string.hpp>

// Eigen
// #include "Eigen/Dense"
// OpenMP
#include <omp.h>

// Don't forget to initialize the VIEW_ID.
const std::string QmitkDataCompressionView::VIEW_ID = "org.mitk.views.m2.datacompression";
using DisplayImageType = itk::Image<m2::DisplayImagePixelType, 3>;
using VectorImageAdaptorType = itk::VectorImageToImageAdaptor<m2::DisplayImagePixelType, 3>;
using VectorImageType = itk::VectorImage<m2::DisplayImagePixelType, 3>;

void QmitkDataCompressionView::CreateQtPartControl(QWidget *parent)
{
  m_Parent = parent;
  // Setting up the UI is a true pleasure when using .ui files, isn't it?
  m_Controls.setupUi(parent);

  auto NodePredicateIsCentroidSpectrum = mitk::NodePredicateFunction::New(
    [this](const mitk::DataNode *node) -> bool
    {
      if (auto intervals = dynamic_cast<m2::IntervalVector *>(node->GetData()))
        return ((unsigned int)(intervals->GetType())) & ((unsigned int)(m2::SpectrumFormat::Centroid));
      return false;
    });

  auto NodePredicateIsActiveHelperNode = mitk::NodePredicateFunction::New(
    [this](const mitk::DataNode *node) { return node->IsOn("helper object", nullptr, false); });
  auto NodePredicateNoActiveHelper = mitk::NodePredicateNot::New(NodePredicateIsActiveHelperNode);

  m_Controls.imageSelection->SetDataStorage(GetDataStorage());
  m_Controls.imageSelection->SetNodePredicate(
    mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<m2::SpectrumImage>::New(), NodePredicateNoActiveHelper));
  m_Controls.imageSelection->SetSelectionIsOptional(true);
  m_Controls.imageSelection->SetEmptyInfo(QString("Image selection"));
  m_Controls.imageSelection->SetPopUpTitel(QString("Image"));

  m_Controls.peakListSelection->SetDataStorage(GetDataStorage());
  m_Controls.peakListSelection->SetNodePredicate(
    mitk::NodePredicateAnd::New(NodePredicateIsCentroidSpectrum, NodePredicateNoActiveHelper));
  m_Controls.peakListSelection->SetSelectionIsOptional(true);
  m_Controls.peakListSelection->SetEmptyInfo(QString("PeakList selection"));
  m_Controls.peakListSelection->SetPopUpTitel(QString("PeakList"));

  // Wire up the UI widgets with our functionality.
  // connect(m_Controls.imageSelection,
  //         &QmitkSingleNodeSelectionWidget::CurrentSelectionChanged,
  //         this,
  //         &QmitkDataCompressionView::OnImageChanged);

  connect(m_Controls.btnRunPCA, SIGNAL(clicked()), this, SLOT(OnStartPCA()));
  connect(m_Controls.btnRunKMeans, SIGNAL(clicked()), this, SLOT(OnStartKMeans()));
  connect(m_Controls.btnRunTSNE, SIGNAL(clicked()), this, SLOT(OnStartTSNE()));
 
}

void QmitkDataCompressionView::SetFocus() {}

  
// Function to calculate Euclidean distance between two points
double euclideanDistance(const Eigen::VectorXd& point1, const Eigen::VectorXd& point2) {
    return (point1 - point2).norm();
}


void QmitkDataCompressionView::DoKMeans(const Eigen::MatrixXd& data, int k, std::vector<int>& clusterAssignments){
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

// void QmitkDataCompressionView::DoKMeans(unsigned int numberOfClasses){

//   using PixelType = float;
//   constexpr unsigned int Dimension = 3;
//   using ImageType = itk::VectorImage<PixelType, Dimension>;

  
//   using EstimatorType = itk::ImageKmeansModelEstimator<ImageType>;
//   EstimatorType::Pointer estimator = EstimatorType::New();
//   estimator->SetInputImage(inputImage);
//   estimator->SetNumberOfClasses(numberOfClasses);

//     try {
//         estimator->Update();
//     } catch (itk::ExceptionObject & error) {
//         std::cerr << "Error: " << error << std::endl;
//         return EXIT_FAILURE;
//     }

//     // Assign labels to each pixel
//     using LabelImageType = itk::Image<unsigned char, Dimension>;
//     using AssignFilterType = itk::KmeansModelEstimator<ImageType>::AssignLabelsFilterType;
//     AssignFilterType::Pointer assignFilter = AssignFilterType::New();
//     assignFilter->SetInputImage(inputImage);
//     assignFilter->SetModelEstimator(estimator);
//     assignFilter->Update();

//     LabelImageType::Pointer labeledImage = assignFilter->GetOutput();
// }


void QmitkDataCompressionView::OnStartKMeans()
{
  for (auto imageNode : m_Controls.imageSelection->GetSelectedNodesStdVector())
  {
    for (auto vectorNode : m_Controls.peakListSelection->GetSelectedNodesStdVector())
    {
      auto image = dynamic_cast<const m2::SpectrumImage *>(imageNode->GetData());
      if (!image->GetImageAccessInitialized())
        return;

      auto vector = dynamic_cast<m2::IntervalVector *>(vectorNode->GetData());
      const auto &intervals = vector->GetIntervals();

      itk::Image<double, 3>::Pointer itkImage;
      mitk::CastToItkImage(image, itkImage);
      unsigned int N = image->GetDimension(0) * image->GetDimension(1) * image->GetDimension(2);

      // std::vector<mitk::Image::Pointer> temporaryImages;
      // auto inputImage = itk::Image<m2::DisplayImagePixelType, 3>::New();
      // auto region = itkImage->GetLargestPossibleRegion();
      // region.SetSize(2, intervals.size());
      // inputImage->SetRegions(region);
      // inputImage->Allocate();



      auto progressBar = mitk::ProgressBar::GetInstance();

      auto m = image->GetMaskImage();
      mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType, 3> maskAcc(m);
      auto validPixel = std::accumulate(maskAcc.GetData(), maskAcc.GetData() + N, 0);

      Eigen::MatrixXd data(validPixel, intervals.size());
      MITK_INFO << validPixel << " " << data.cols() << " " << data.rows();

      {
        auto ionImage = mitk::Image::New();
        ionImage->Initialize(image);

        progressBar->AddStepsToDo(intervals.size() + 1);
        for (size_t row = 0; row < intervals.size(); ++row)
        {
          progressBar->Progress();
          const auto mz = intervals.at(row).x.mean();
          
          image->GetImage(mz, image->ApplyTolerance(mz), image->GetMaskImage(), ionImage);
          mitk::ImagePixelReadAccessor<m2::DisplayImagePixelType, 3> acc(ionImage);
          // std::copy(acc.GetData(), acc.GetData() + N, inputImage->GetBufferPointer() + row * N);

          unsigned v = 0;
          for(unsigned long i = 0; i < N; ++i)
          {
            if( maskAcc.GetData()[i] != 0){
              data(v, row) = acc.GetData()[i];
              ++v;
            }
          }
        }
      }

      // auto result = mitk::Image::New();
      // mitk::CastToMitkImage(inputImage, result);

      // mitk::DataNode::Pointer node = mitk::DataNode::New();
      // node->SetData(result);
      // node->SetName("VectorImage");
      // this->GetDataStorage()->Add(node);
  


      std::vector<int> clusterAssignments(data.rows(), -1);
      DoKMeans(data, m_Controls.kmeans_clusters->value(), clusterAssignments);
      for (int i = 0; i < data.rows(); ++i) {
        std::cout << "Data point " << i << " assigned to cluster " << clusterAssignments[i] << std::endl;
      }
      
      auto clusteredImage = mitk::LabelSetImage::New();
      clusteredImage->Initialize(m);
      // clusteredImage->GetActiveLabelSet();
      // auto labelset = clusteredImage->GetActiveLabelSet();
      for(int i = 0; i < m_Controls.kmeans_clusters->value(); ++i){
        mitk::Label::Pointer label = mitk::Label::New();
        label->SetName("Cluster " + std::to_string(i));
        label->SetValue(i);
        label->SetOpacity(0.5);
        label->SetColor(m2::RandomColor());
        clusteredImage->AddLabel(label,0);
      }
      
      {
        mitk::ImagePixelWriteAccessor<mitk::LabelSetImage::PixelType, 3> c_acc(clusteredImage);

        unsigned long v = 0;
        for(unsigned long i = 0; i < N; ++i)
        {
          if( maskAcc.GetData()[i] != 0){
            c_acc.GetData()[i] = clusterAssignments[v];
            ++v;
          }
        }
      }

      mitk::DataNode::Pointer node2 = mitk::DataNode::New();
      node2->SetData(clusteredImage);
      node2->SetName("ClusteredImage");
      this->GetDataStorage()->Add(node2);

    }
  }
}


void QmitkDataCompressionView::OnStartPCA()
{
  for (auto imageNode : m_Controls.imageSelection->GetSelectedNodesStdVector())
  {
    for (auto vectorNode : m_Controls.peakListSelection->GetSelectedNodesStdVector())
    {
      auto image = dynamic_cast<const m2::SpectrumImage *>(imageNode->GetData());
      auto vector = dynamic_cast<m2::IntervalVector *>(vectorNode->GetData());
      const auto &intervals = vector->GetIntervals();

      if (!image->GetImageAccessInitialized())
        continue;

      auto filter = m2::PcaImageFilter::New();
      filter->SetMaskImage(image->GetMaskImage());

      std::vector<mitk::Image::Pointer> temporaryImages;
      auto progressBar = mitk::ProgressBar::GetInstance();
      progressBar->AddStepsToDo(intervals.size() + 1);
      size_t inputIdx = 0;
      for (size_t row = 0; row < intervals.size(); ++row)
      {
        progressBar->Progress();
        temporaryImages.push_back(mitk::Image::New());
        temporaryImages.back()->Initialize(image);
        const auto mz = intervals.at(row).x.mean();
        image->GetImage(mz, image->ApplyTolerance(mz), image->GetMaskImage(), temporaryImages.back());
        filter->SetInput(inputIdx, temporaryImages.back());
        ++inputIdx;
      }

      if (temporaryImages.size() <= 2)
      {
        progressBar->Progress();
        QMessageBox::warning(nullptr,
                             "Select image,s first!",
                             "Select at least three peaks!",
                             QMessageBox::StandardButton::NoButton,
                             QMessageBox::StandardButton::Ok);
        continue;
      }

      filter->SetNumberOfComponents(m_Controls.pca_dims->value());
      filter->Update();
      progressBar->Progress();

      auto outputNode = mitk::DataNode::New();
      mitk::Image::Pointer data = filter->GetOutput(0);
      outputNode->SetData(data);
      outputNode->SetName("PCA");
      this->GetDataStorage()->Add(outputNode, const_cast<mitk::DataNode *>(imageNode.GetPointer()));
    }
  }

  // const auto &peakList = m_PeakList;

  // auto outputNode2 = mitk::DataNode::New();
  // mitk::Image::Pointer data2 = filter->GetOutput(1);
  // outputNode2->SetData(data2);
  // outputNode2->SetName("pcs");
  // this->GetDataStorage()->Add(outputNode2, node.GetPointer());
}

void QmitkDataCompressionView::OnStartTSNE()
{
  for (auto node : m_Controls.imageSelection->GetSelectedNodesStdVector())
  {
    if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
    {
      
      auto child = this->GetDataStorage()->GetNamedDerivedNode("PCA", node);
      if (!child)
        return;

      auto pcaImage = dynamic_cast<mitk::Image *>(child->GetData());
      const auto pcaComponents = pcaImage->GetPixelType().GetNumberOfComponents();

      auto filter = m2::TSNEImageFilter::New();
      filter->SetPerplexity(m_Controls.tsne_perplexity->value());
      filter->SetIterations(m_Controls.tnse_iters->value());
      filter->SetTheta(m_Controls.tsne_theta->value());

      using MaskImageType = itk::Image<mitk::LabelSetImage::PixelType, 3>;
      MaskImageType::Pointer maskImageItk;
      mitk::Image::Pointer maskImage;
      mitk::CastToItkImage(image->GetMaskImage(), maskImageItk);
      auto caster = itk::ShrinkImageFilter<MaskImageType, MaskImageType>::New();
      caster->SetInput(maskImageItk);
      caster->SetShrinkFactor(0, m_Controls.tsne_shrink->value());
      caster->SetShrinkFactor(1, m_Controls.tsne_shrink->value());
      caster->SetShrinkFactor(2, 1);
      caster->Update();

      mitk::CastToMitkImage(caster->GetOutput(), maskImage);

      filter->SetMaskImage(maskImage);
      // const auto &peakList = image->GetPeaks();
      size_t index = 0;

      mitk::ImageReadAccessor racc(pcaImage);
      auto *inputData = static_cast<const typename DisplayImageType::PixelType *>(racc.GetData());

      std::vector<mitk::Image::Pointer> bufferedImages(pcaComponents);
      unsigned int n = pcaImage->GetDimensions()[0] * pcaImage->GetDimensions()[1] * pcaImage->GetDimensions()[2];
      for (auto &I : bufferedImages)
      {
        I = mitk::Image::New();
        I->Initialize(image);
        {
          mitk::ImageWriteAccessor acc(I);
          auto outCData = static_cast<typename DisplayImageType::PixelType *>(acc.GetData());
          for (unsigned int k = 0; k < n; ++k)
            *(outCData + k) = *(inputData + (k * pcaComponents) + index);
        }

        DisplayImageType::Pointer cImage;
        mitk::CastToItkImage(I, cImage);

        auto caster = itk::ShrinkImageFilter<DisplayImageType, DisplayImageType>::New();
        caster->SetInput(cImage);
        caster->SetShrinkFactor(0, m_Controls.tsne_shrink->value());
        caster->SetShrinkFactor(1, m_Controls.tsne_shrink->value());
        caster->SetShrinkFactor(2, 1);
        caster->Update();

        // Buffer the image
        mitk::CastToMitkImage(caster->GetOutput(), I);
        filter->SetInput(index, I);
        ++index;
      }
      filter->Update();

      auto outputNode = mitk::DataNode::New();
      auto data = m2::MultiSliceFilter::ConvertMitkVectorImageToRGB(ResampleVectorImage(filter->GetOutput(), image));
      outputNode->SetData(data);
      outputNode->SetName("tSNE");
      this->GetDataStorage()->Add(outputNode, const_cast<mitk::DataNode *>(node.GetPointer()));

    }
  }
}

mitk::Image::Pointer QmitkDataCompressionView::ResampleVectorImage(mitk::Image::Pointer vectorImage,
                                                                   mitk::Image::Pointer referenceImage)
{
  const unsigned int components = vectorImage->GetPixelType().GetNumberOfComponents();

  VectorImageType::Pointer vectorImageItk;
  mitk::CastToItkImage(vectorImage, vectorImageItk);

  DisplayImageType::Pointer referenceImageItk;
  mitk::CastToItkImage(referenceImage, referenceImageItk);

  auto resampledVectorImageItk = VectorImageType::New();
  resampledVectorImageItk->SetOrigin(referenceImageItk->GetOrigin());
  resampledVectorImageItk->SetDirection(referenceImageItk->GetDirection());
  resampledVectorImageItk->SetSpacing(referenceImageItk->GetSpacing());
  resampledVectorImageItk->SetRegions(referenceImageItk->GetLargestPossibleRegion());
  resampledVectorImageItk->SetNumberOfComponentsPerPixel(components);
  resampledVectorImageItk->Allocate();
  itk::VariableLengthVector<m2::DisplayImagePixelType> v(components);
  v.Fill(0);
  resampledVectorImageItk->FillBuffer(v);

  auto inAdaptor = VectorImageAdaptorType::New();
  auto outAdaptor = VectorImageAdaptorType::New();
  using LinearInterpolatorType = itk::LinearInterpolateImageFunction<VectorImageAdaptorType>;
  using TransformType = itk::IdentityTransform<m2::DisplayImagePixelType, 3>;

  for (unsigned int i = 0; i < components; ++i)
  {
    inAdaptor->SetExtractComponentIndex(i);
    inAdaptor->SetImage(vectorImageItk);
    inAdaptor->SetOrigin(vectorImageItk->GetOrigin());
    inAdaptor->SetDirection(vectorImageItk->GetDirection());
    inAdaptor->SetSpacing(vectorImageItk->GetSpacing());

    outAdaptor->SetExtractComponentIndex(i);
    outAdaptor->SetImage(resampledVectorImageItk);

    auto resampler = itk::ResampleImageFilter<VectorImageAdaptorType, DisplayImageType>::New();
    resampler->SetInput(inAdaptor);
    resampler->SetOutputParametersFromImage(referenceImageItk);
    resampler->SetInterpolator(LinearInterpolatorType::New());
    resampler->SetTransform(TransformType::New());
    resampler->Update();

    itk::ImageAlgorithm::Copy<DisplayImageType, VectorImageAdaptorType>(
      resampler->GetOutput(),
      outAdaptor,
      resampler->GetOutput()->GetLargestPossibleRegion(),
      outAdaptor->GetLargestPossibleRegion());
  }

  mitk::Image::Pointer outImage;
  mitk::CastToMitkImage(resampledVectorImageItk, outImage);

  return outImage;
}