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
#include <mitkDataNode.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkLevelWindowProperty.h>
#include <mitkLookupTableProperty.h>
#include <mitkOperation.h>
#include <signal/m2PeakDetection.h>

namespace m2
{
  itkEventMacroDefinition(InitializationFinishedEvent, itk::AnyEvent);
} // namespace m2

double m2::SpectrumImage::ApplyTolerance(double xValue) const
{
  if (this->GetUseToleranceInPPM())
    return m2::PartPerMillionToFactor(this->GetTolerance()) * xValue;
  else
    return this->GetTolerance();
}

mitk::Image::Pointer m2::SpectrumImage::GetNormalizationImage()
{
    return GetNormalizationImage(m_NormalizationStrategy);
}

mitk::Image::Pointer m2::SpectrumImage::GetNormalizationImage(m2::NormalizationStrategyType type)
{
  if (m_NormalizationImages.find(type) != m_NormalizationImages.end())
    
      return m_NormalizationImages.at(type).image;
  return nullptr;
}

mitk::Image::Pointer m2::SpectrumImage::GetNormalizationImage() const
{
    return GetNormalizationImage(m_NormalizationStrategy);
}

mitk::Image::Pointer m2::SpectrumImage::GetNormalizationImage(m2::NormalizationStrategyType type) const
{
  if (m_NormalizationImages.find(type) != m_NormalizationImages.end())
      return m_NormalizationImages.at(type).image;
  return nullptr;
};

bool m2::SpectrumImage::GetNormalizationImageStatus(m2::NormalizationStrategyType type)
{
  return m_NormalizationImages[type].isInitialized;
};

void m2::SpectrumImage::SetNormalizationImageStatus(m2::NormalizationStrategyType type, bool initialized){
  m_NormalizationImages[type].isInitialized = initialized;
}

void m2::SpectrumImage::SetNormalizationImage(mitk::Image::Pointer image, m2::NormalizationStrategyType type)
{
  m_NormalizationImages[type].image = image;
}


// void m2::SpectrumImage::Check(const std::string &key, mitk::Image *img)
// {
//   m_ImageArtifacts[key] = img;
//   const unsigned int DIMS = 3;

//   if (img->GetDimension() != DIMS)
//   {
//     mitkThrow() << "SpectrumBaseImage related image artifacts require 3 dimensions.";
//   }
//   auto aD = this->GetDimensions();
//   auto bD = img->GetDimensions();

//   if (!std::equal(aD, aD + DIMS, bD))
//   {
//     mitkThrow() << "SpectrumBaseImage related image artifacts require identical image dimensions.";
//   }

//   auto aS = this->GetGeometry()->GetSpacing();
//   auto bS = img->GetGeometry()->GetSpacing();

//   if (!std::equal(
//         aS.GetDataPointer(), aS.GetDataPointer() + DIMS, bS.GetDataPointer(), [](const auto &a, const auto &b) {
//           return itk::Math::FloatAlmostEqual(a, b);
//         }))
//   {
//     mitkThrow() << "SpectrumBaseImage related image artifacts require identical spacings.";
//   }

//   // if spacing and dimensions are equal, copy origin and vtkMatrix to the new image artifact.
//   img->SetClonedTimeGeometry(this->GetTimeGeometry());
// }

void m2::SpectrumImage::ApplyMoveOriginOperation(const mitk::Vector3D &v)
{
  auto geometry = this->GetGeometry();
  geometry->Translate(v);
  std::vector<mitk::BaseData *> imageList{m_IndexImage, m_MaskImage, m_Points};
  for(auto kv : m_NormalizationImages)
    imageList.push_back(kv.second.image);
  
  for (auto image : imageList)
  {
    if (image)
    {
      geometry = image->GetGeometry();
      geometry->Translate(v);
    }
  }
}

void m2::SpectrumImage::ApplyGeometryOperation(mitk::Operation *op)
{
  auto manipulatedGeometry = this->GetGeometry()->Clone();
  manipulatedGeometry->ExecuteOperation(op);
  this->GetGeometry()->SetIdentity();
  this->GetGeometry()->Compose(manipulatedGeometry->GetIndexToWorldTransform());

  std::vector<mitk::BaseData *> imageList{m_IndexImage, m_MaskImage, m_Points};
  
  for(auto kv : m_NormalizationImages)
    imageList.push_back(kv.second.image);
  
  for (auto image : imageList)
  {
    if (image)
    {
      auto manipulatedGeometry = image->GetGeometry()->Clone();
      manipulatedGeometry->ExecuteOperation(op);
      image->GetGeometry()->SetIdentity();
      image->GetGeometry()->Compose(manipulatedGeometry->GetIndexToWorldTransform());
    }
  }
}

std::vector<double> &m2::SpectrumImage::GetSkylineSpectrum()
{
  return m_SpectraArtifacts[(SpectrumType::Maximum)];
}

std::vector<double> &m2::SpectrumImage::GetMeanSpectrum()
{
  return m_SpectraArtifacts[(SpectrumType::Mean)];
}

std::vector<double> &m2::SpectrumImage::GetSumSpectrum()
{
  return m_SpectraArtifacts[(SpectrumType::Sum)];
}

std::vector<double> &m2::SpectrumImage::GetXAxis()
{
  return m_XAxis;
}

const std::vector<double> &m2::SpectrumImage::GetXAxis() const
{
  return m_XAxis;
}

void m2::SpectrumImage::GetImage(double, double, const mitk::Image *, mitk::Image *) const
{
  MITK_WARN("SpectrumImage") << "Get image is not implemented in derived class!";
}

m2::SpectrumImage::~SpectrumImage() {}
m2::SpectrumImage::SpectrumImage() : mitk::Image() {}
