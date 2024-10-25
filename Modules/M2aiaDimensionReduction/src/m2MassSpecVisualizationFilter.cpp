/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <m2MassSpecVisualizationFilter.h>
#include <m2ImzMLSpectrumImage.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>

void m2::MassSpecVisualizationFilter::SetMaskImage(mitk::Image::Pointer mask)
{
  m_MaskImage = mask;
}

void m2::MassSpecVisualizationFilter::GetValidIndices()
{
  if (this->GetIndexedInputs().empty())
  {
    return;
  }
  else
  {
    AccessFixedDimensionByItk(m_MaskImage,
                              [this](auto input) {
                                MITK_INFO << input->GetLargestPossibleRegion().GetSize(0) << " "
                                          << input->GetLargestPossibleRegion().GetSize(1) << " "
                                          << input->GetLargestPossibleRegion().GetSize(2);
                                for (unsigned int x = 0; x < input->GetLargestPossibleRegion().GetSize(0); x++)
                                {
                                  for (unsigned int y = 0; y < input->GetLargestPossibleRegion().GetSize(1); y++)
                                  {
                                    for (unsigned int z = 0; z < input->GetLargestPossibleRegion().GetSize(2); z++)
                                    {
                                      if (input->GetPixel({x, y, z}) > 0)
                                      {
                                        this->m_ValidIndices.push_back({x, y, z});
                                      }
                                    }
                                  }
                                }
                              },
                              3);
  }
}

void m2::MassSpecVisualizationFilter::initializeItkImage(itk::Image<RGBPixel, 3>::Pointer itkImage)
{
  itk::Image<RGBPixel, 3>::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  itk::Image<RGBPixel, 3>::SizeType size;
  auto dimensions = this->GetInput()->GetDimensions();
  size[0] = dimensions[0];
  size[1] = dimensions[1];
  size[2] = dimensions[2];
  itk::Image<RGBPixel, 3>::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  itkImage->SetRegions(region);
  itkImage->Allocate();

  itk::ImageRegionIterator<itk::Image<RGBPixel, 3>> imageIterator(itkImage, itkImage->GetRequestedRegion());
  imageIterator.GoToBegin();
  auto pixel = RGBPixel();
  pixel.SetRed(0);
  pixel.SetBlue(0);
  pixel.SetGreen(0);

  while (!imageIterator.IsAtEnd())
  {
    imageIterator.Set(pixel);
    ++imageIterator;
  }
}

itk::VectorImage<m2::DisplayImagePixelType, 3>::Pointer m2::MassSpecVisualizationFilter::initializeItkVectorImage(unsigned int components)
{
  auto vectorImage = itk::VectorImage<PixelType, 3>::New();
  vectorImage->SetVectorLength(components);

  itk::Image<PixelType, 3>::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  itk::Image<PixelType, 3>::SizeType size;
  auto dimensions = this->GetInput()->GetDimensions();
  size[0] = dimensions[0];
  size[1] = dimensions[1];
  size[2] = dimensions[2];
  itk::Image<PixelType, 3>::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);
  vectorImage->SetRegions(region);
  
  vectorImage->Allocate();

  itk::VariableLengthVector<PixelType> initialPixelValues;
  initialPixelValues.SetSize(m_NumberOfComponents);
  initialPixelValues.Fill(0.0);
  itk::ImageRegionIterator<itk::VectorImage<PixelType, 3>> imageIterator(vectorImage,
                                                                         vectorImage->GetRequestedRegion());

  while (!imageIterator.IsAtEnd())
  {
    imageIterator.Set(initialPixelValues);
    ++imageIterator;
  }

  return vectorImage;
}
