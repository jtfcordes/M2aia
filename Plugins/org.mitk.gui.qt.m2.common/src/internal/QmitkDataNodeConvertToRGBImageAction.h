/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#pragma once
#include "QmitkAbstractDataNodeAction.h"

// mitk core
#include <type_traits>

#include <itkImageRegionIterator.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <mitkDataNode.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkImageAccessByItk.h>
#include <mitkImagePixelWriteAccessor.h>
// qt
#include <QAction>


class QmitkDataNodeConvertToRGBImageAction : public QAction, public QmitkAbstractDataNodeAction
{
  Q_OBJECT

public:

  QmitkDataNodeConvertToRGBImageAction(QWidget* parent, berry::IWorkbenchPartSite::Pointer workbenchPartSite);
  QmitkDataNodeConvertToRGBImageAction(QWidget* parent = nullptr, berry::IWorkbenchPartSite* workbenchPartSite = nullptr);
  

protected:

  void InitializeAction() override;

  using RGBPixel = itk::RGBPixel<unsigned char>;

  static mitk::Image::Pointer ConvertMitkVectorImageToRGB(mitk::Image::Pointer vImage)
    {
      

      // AccessVectorFixedDimensionByItk(vImage, ([](auto itkImage){
      //   using ImageType = typename std::remove_ptr<decltype(itkImage)>::type;
      //   using PixelType = typename ImageType::PixelType;
        
      //   auto region =  itkImage->GetLargestPossibleRegion();
      //   itk::ImageRegionConstIterator<ImageType> it(itkImage, region);
        
      //   std::vector<float> minValues(3, std::numeric_limits<float>::max());
      //   std::vector<float> maxValues(3, std::numeric_limits<float>::lowest());

      //   for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
      //       PixelType pixel = it.Get();
      //       for (unsigned int i = 0; i < numComponents; ++i) {
      //           if (pixel[i] < minValues[i]) minValues[i] = pixel[i];
      //           if (pixel[i] > maxValues[i]) maxValues[i] = pixel[i];
      //       }
      //   }
      // })
      // , 3);

      
      mitk::Image::Pointer result;


      itk::VectorImage<unsigned char, 3>::Pointer vectorImage;
      mitk::CastToItkImage(vImage, vectorImage);
      mitk::CastToMitkImage(ConvertVectorImageToRGB(vectorImage), result);
      return result;
    }

    /*This function casts a itk vector Image with vector length of 3, to a RGB itk Image. The buffer of the
      vector image is copied to the RGB image.*/
    static itk::Image<itk::RGBPixel<unsigned char>, 3>::Pointer ConvertVectorImageToRGB(
      itk::VectorImage<unsigned char, 3>::Pointer vectorImage)
    {
      itk::Image<RGBPixel, 3>::Pointer rgbImage = itk::Image<RGBPixel, 3>::New();

      itk::Image<RGBPixel, 3>::IndexType start;
      start[0] = 0;
      start[1] = 0;
      start[2] = 0;

      itk::Image<RGBPixel, 3>::SizeType size;
      auto dimensions = vectorImage->GetLargestPossibleRegion().GetSize();
      size[0] = dimensions[0];
      size[1] = dimensions[1];
      size[2] = dimensions[2];
      itk::Image<RGBPixel, 3>::RegionType region;
      region.SetSize(size);
      region.SetIndex(start);

      rgbImage->SetRegions(region);
      rgbImage->Allocate();

      rgbImage->SetOrigin(vectorImage->GetOrigin());
      rgbImage->SetSpacing(vectorImage->GetSpacing());
      rgbImage->SetDirection(vectorImage->GetDirection());

      itk::ImageRegionIterator<itk::Image<RGBPixel, 3>> imageIterator(rgbImage, rgbImage->GetRequestedRegion());
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

      

      memcpy(
        rgbImage->GetBufferPointer(), vectorImage->GetBufferPointer(), sizeof(RGBPixel) * size[0] * size[1] * size[2]);


      return rgbImage;
    }

};