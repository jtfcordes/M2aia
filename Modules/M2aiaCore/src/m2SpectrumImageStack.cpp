/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <array>
#include <cstdlib>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itksys/SystemTools.hxx>
#include <m2ImzMLSpectrumImage.h>
#include <m2ISpectrumImageSource.h>
#include <m2SpectrumImageStack.h>
#include <mitkExtractSliceFilter.h>
#include <mitkIOUtil.h>
#include <mitkITKImageImport.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageReadAccessor.h>
#include <mitkImageWriteAccessor.h>
#include <mitkProgressBar.h>
#include <signal/m2PeakDetection.h>

#include <mitkCoreServices.h>
#include <mitkIPreferencesService.h>
#include <mitkIPreferences.h>

namespace m2
{
  SpectrumImageStack::SpectrumImageStack(unsigned int stackSize, double spacingZ) : 
  m_StackSize(stackSize), m_SpacingZ(spacingZ)
  {
    SetPropertyValue<double>("m2aia.xs.min", std::numeric_limits<double>::max());
    SetPropertyValue<double>("m2aia.xs.max", std::numeric_limits<double>::min());

    m_SliceTransformers.resize(stackSize);
  }


  void SpectrumImageStack::Insert(unsigned int sliceId, std::shared_ptr<m2::ElxRegistrationHelper> transformer)
  {
    m_SliceTransformers[sliceId] = transformer;

    if (auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(transformer->GetMovingImage().GetPointer()))
    {
      auto newMin = spectrumImage->GetPropertyValue<double>("m2aia.xs.min");
      auto newMax = spectrumImage->GetPropertyValue<double>("m2aia.xs.max");
      auto xLabel = spectrumImage->GetSpectrumType().XAxisLabel;
      
      auto currentMin = GetPropertyValue<double>("m2aia.xs.min");
      auto currentMax = GetPropertyValue<double>("m2aia.xs.max");
      if (newMin < currentMin)
        SetPropertyValue<double>("m2aia.xs.min", newMin);
      if (newMax > currentMax)
        SetPropertyValue<double>("m2aia.xs.max", newMax);

      this->GetSpectrumType().XAxisLabel = xLabel;
    }
    else
    {
      mitkThrow() << "Spectrum image base object expected!";
    }
  }



  void SpectrumImageStack::InitializeGeometry()
  {
    unsigned int dims[3];
    mitk::Vector3D spacing;

    auto &transformer = m_SliceTransformers.front();
    auto image = transformer->GetMovingImage();
    if (!transformer->GetTransformation().empty())
    {
      image = transformer->WarpImage(image);
    }
    std::copy(image->GetDimensions(), image->GetDimensions() + 3, dims);
    spacing = image->GetGeometry()->GetSpacing();

    dims[2] = m_SliceTransformers.size();
    this->Initialize(mitk::MakeScalarPixelType<m2::DisplayImagePixelType>(), 3, dims);

    spacing[2] = m_SpacingZ;
    this->GetGeometry()->SetSpacing(spacing);

    unsigned int sliceId = 0;

    // fill with data
    for (auto &transformer : m_SliceTransformers)
    {      
      if (auto movingImage = transformer->GetMovingImage())
      {
        // create temp image and copy requested image range to the stack
        if (!transformer->GetTransformation().empty())
        {
          movingImage = transformer->WarpImage(movingImage);
          CopyWarpedImageToStackImage(movingImage, this, sliceId);
        }
        else
        {
          CopyWarpedImageToStackImage(movingImage, this, sliceId);
        }
      }
      ++sliceId;
    }
  }

  void SpectrumImageStack::InitializeProcessor()
  {
    if(m_SliceTransformers.empty()){
      MITK_ERROR("SpectrumImageStack::InitializeProcessor") << "No transformer found!";
      return;
    }

    // assign import spectrum type based on first transformer
    // this must be equally for all slices/images
    auto movingImage = m_SliceTransformers.front()->GetMovingImage().GetPointer();
    if(auto specImage =dynamic_cast<m2::SpectrumImage *>(movingImage))
      m_SpectrumType.Format = specImage->GetSpectrumType().Format;
    
    
    for (auto &transformer : m_SliceTransformers)
    {
      auto specImage = dynamic_cast<m2::SpectrumImage *>(transformer->GetMovingImage().GetPointer());
      
      if (m_SpectrumType.Format != specImage->GetSpectrumType().Format)
      {
        MITK_ERROR("SpectrumImageStack::InitializeProcessor") << "Different import modes detected";
      }
    }
  }

  void SpectrumImageStack::InitializeImageAccess()
  {
    std::list<m2::Interval> peaks;
    auto* preferencesService = mitk::CoreServices::GetPreferencesService();
    auto* preferences = preferencesService->GetSystemPreferences();
    const auto bins = preferences->GetInt("m2aia.view.spectrum.bins", 15000);
    //    const auto normalizationStrategy = GetNormalizationStrategy();

    double max = 0;
    double min = std::numeric_limits<double>::max();
    double binSize = 1;

    for (auto &transformer : m_SliceTransformers)
    {
      auto specImage = dynamic_cast<m2::SpectrumImage *>(transformer->GetMovingImage().GetPointer());
      auto xAxis = specImage->GetXAxis();
      min = std::min(min, xAxis.front());
      max = std::max(max, xAxis.back());
    }

    binSize = (max - min) / double(bins);
    
    std::vector<double> xSumVec(bins);
    std::vector<double> hits(bins);

    std::vector<double> ySumVec = GetSumSpectrum();
    ySumVec.resize(bins, 0);
    std::vector<double> yMeanVec = GetMeanSpectrum();
    yMeanVec.resize(bins, 0);
    std::vector<double> yMaxVec = GetSkylineSpectrum();
    yMaxVec.resize(bins, 0);

    for (auto &transformer : m_SliceTransformers)
    {
      auto specImage = dynamic_cast<m2::SpectrumImage *>(transformer->GetMovingImage().GetPointer());
      auto sliceXAxis = specImage->GetXAxis();
      auto sliceSumVec = specImage->GetSumSpectrum();
      auto sliceMaxVec = specImage->GetSkylineSpectrum();
      auto sliceMeanVec = specImage->GetMeanSpectrum();

      for (unsigned int k = 0; k < sliceXAxis.size(); ++k)
      {
        auto j = (long)((sliceXAxis[k] - min) / binSize);

        if (j >= bins)
          j = bins - 1;
        else if (j < 0)
          j = 0;

        xSumVec[j] += sliceXAxis[k];
        ySumVec[j] += sliceSumVec[k];
        yMeanVec[j] += sliceMeanVec[k];
        yMaxVec[j] = std::max(yMaxVec[j], double(sliceMaxVec[k]));
        hits[j]++;
      }
    }

    std::vector<double> &xVecFinal = GetXAxis();
    xVecFinal.clear();
    std::vector<double> &ySumVecFinal = GetSumSpectrum();
    ySumVecFinal.clear();
    std::vector<double> &yMeanVecFinal = GetMeanSpectrum();
    yMeanVecFinal.clear();
    std::vector<double> &yMaxVecFinal = GetSkylineSpectrum();
    yMaxVecFinal.clear();

    for (int k = 0; k < bins; ++k)
    {
      if (hits[k] > 0)
      {
        xVecFinal.push_back(xSumVec[k] / (double)hits[k]);      // mean uof sum of x values within bin range
        ySumVecFinal.push_back(ySumVec[k] / (double)hits[k]);   // mean of sums
        yMeanVecFinal.push_back(yMeanVec[k] / (double)hits[k]); // mean of means
        yMaxVecFinal.push_back(yMaxVec[k]);                     // max
      }
    }


    SetPropertyValue<double>("m2aia.xs.min", xVecFinal.front());
    SetPropertyValue<double>("m2aia.xs.max", xVecFinal.back());

    m_ImageAccessInitialized = true;
    
  }

  void SpectrumImageStack::SpectrumImageStack::CopyWarpedImageToStackImage(mitk::Image *warped,
                                                                           mitk::Image *stack,
                                                                           unsigned i) const
  {
    auto warpN = warped->GetDimensions()[0] * warped->GetDimensions()[1];
    auto stackN = stack->GetDimensions()[0] * stack->GetDimensions()[1];
    if (warpN != stackN)
      mitkThrow() << "Slice dimensions are not equal for target slice with index !" << i;
    if (i >= stack->GetDimensions()[2])
      mitkThrow() << "Stack index is invalid! Z dim is " << stack->GetDimensions()[2];

    mitk::ImageWriteAccessor stackAccess(stack);
    auto stackData = static_cast<m2::DisplayImagePixelType *>(stackAccess.GetData());

    if(m_UseSliceWiseMaximumNormalization){
      AccessByItk(warped,
                (
                  [&](auto itkImg)
                  {
                    auto warpedData = itkImg->GetBufferPointer();
                    auto max = *(std::max_element(warpedData, warpedData + stackN));
                    std::transform(warpedData, warpedData + stackN, stackData + (i * stackN), [max](auto v){return v/max;});
                  }));
    }else{


    AccessByItk(warped,
                (
                  [&](auto itkImg)
                  {
                    auto warpedData = itkImg->GetBufferPointer();
                    std::copy(warpedData, warpedData + stackN, stackData + (i * stackN));
                  }));
    }
  }

  void SpectrumImageStack::GetImage(double center, double tol, const mitk::Image * /*mask*/, mitk::Image *img) const
  {
    mitk::ProgressBar::GetInstance()->AddStepsToDo(m_SliceTransformers.size());
    int sliceId = 0;
    for (auto &transformer : m_SliceTransformers)
    {
      mitk::ProgressBar::GetInstance()->Progress();
      
      
      if (auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(transformer->GetMovingImage().GetPointer()))
      {
        // create temp image and copy requested image range to the stack
        auto imageTemp = mitk::Image::New();
        imageTemp->Initialize(spectrumImage);
        spectrumImage->GetImage(center, tol, spectrumImage->GetMaskImage(), imageTemp);
        if (!transformer->GetTransformation().empty())
        {
          imageTemp = transformer->WarpImage(imageTemp);
          CopyWarpedImageToStackImage(imageTemp, img, sliceId);
        }
        else
        {
          CopyWarpedImageToStackImage(imageTemp, img, sliceId);
        }
      }
      ++sliceId;
    }

    // current->GetImage(mz, tol, current->GetMaskImage(), current);

    // auto warped = Transform(
    //   current, transformations, [this](auto trafo) { ReplaceParameter(trafo, "ResultImagePixelType", "\"double\"");
    //   });

    // auto warpedIndexImage = TransformImageUsingNNInterpolation(i);
    // CopyWarpedImageToStack<IndexImagePixelType>(warpedIndexImage, GetIndexImage(), i);
  }

} // namespace m2
