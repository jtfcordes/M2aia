/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#include <M2aiaCoreIOExports.h>
#include <itkTIFFImageIO.h>
#include <m2MicroscopyTiffImageIO.h>
#include <mitkAbstractFileIO.h>
#include <mitkIOMimeTypes.h>
#include <mitkImage.h>
#include <mitkImageAccessByItk.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkItkImageIO.h>

#define RGB_PIXEL_TYPE_SEQUENCE                                                                                             \
  (itk::RGBPixel<char>)(itk::RGBPixel<unsigned char>)(                                                                      \
    itk::RGBPixel<                                                                                                          \
      unsigned short>)(itk::                                                                                                \
                         RGBPixel<                                                                                          \
                           short>)(itk::                                                                                    \
                                     RGBPixel<                                                                              \
                                       float>)(itk::                                                                        \
                                                 RGBAPixel<                                                                 \
                                                   char>)(itk::                                                             \
                                                            RGBAPixel<                                                      \
                                                              unsigned char>)(itk::                                         \
                                                                                RGBAPixel<                                  \
                                                                                  unsigned short>)(itk::                    \
                                                                                                     RGBAPixel<             \
                                                                                                       short>)(itk::        \
                                                                                                                 RGBAPixel< \
                                                                                                                   float>)

namespace m2
{
  MicroscopyTiffImageIO::MicroscopyTiffImageIO()
    : AbstractFileIO(mitk::Image::GetStaticNameOfClass(), MICROSCOPYTIFF_MIMETYPE(), "ITK TIFFImageIO (m2aia)")
  {
    AbstractFileWriter::SetRanking(0);
    AbstractFileReader::SetRanking(0);

    Options defaultOptions;

    // defaultOptions[OPTION_MERGE_POINTS()] = us::Any(true);
    defaultOptions[A_IMPORT_VOLUME] = us::Any(true);
    defaultOptions[D_THICKNESS] = us::Any(10.0f);

    defaultOptions[E_USE_USER_DEFS] = us::Any(false);
    defaultOptions[F_PIXEL_SIZE_X] = us::Any(10.0f);
    defaultOptions[G_PIXEL_SIZE_Y] = us::Any(10.0f);
    this->SetDefaultReaderOptions(defaultOptions);

    this->RegisterService();
  }

  std::vector<mitk::BaseData::Pointer> MicroscopyTiffImageIO::DoRead()
  {
    // read the file
    auto io = itk::TIFFImageIO::New();
    mitk::ItkImageIO mitkItkImageIO(MICROSCOPYTIFF_MIMETYPE(), io, 0);
    mitkItkImageIO.mitk::AbstractFileReader::SetInput(this->GetInputLocation());
    auto resultVec = mitkItkImageIO.Read();

    // cast result and make shared pointer, remove from list
    mitk::Image::Pointer image = dynamic_cast<mitk::Image *>(resultVec.back().GetPointer());
    
    
    // adapt spacing
    auto spacing = image->GetGeometry()->GetSpacing();
    if (us::any_cast<bool>(GetReaderOption(E_USE_USER_DEFS)))
    {
      spacing[0] = us::any_cast<float>(GetReaderOption(F_PIXEL_SIZE_X)) / 1000.0f;
      spacing[1] = us::any_cast<float>(GetReaderOption(G_PIXEL_SIZE_Y)) / 1000.0f;
    }
    spacing[2] = us::any_cast<float>(GetReaderOption(D_THICKNESS)) / 1000.0f;
    image->SetSpacing(spacing);

    const auto *D = image->GetDimensions();
    const std::array<unsigned int, 3> dims{D[0], D[1], 1};
    const auto N = std::accumulate(D, D + image->GetDimension(), 1, std::multiplies<unsigned int>());
    const auto imagePixelType = image->GetPixelType();
    const auto dimensions = image->GetDimension();
    const auto numComponents = imagePixelType.GetNumberOfComponents();
    
    // const bool splitChannels = us::any_cast<bool>(GetReaderOption(B_SPLIT_CHANNELS));
    // const bool importVolume = us::any_cast<bool>(GetReaderOption(A_IMPORT_VOLUME));

    MITK_INFO << "Components:" << numComponents;
    MITK_INFO << "Dimensions:" << dimensions;
    MITK_INFO << "PixelType:" << imagePixelType.GetComponentTypeAsString();

    // split components
    if(numComponents >= 3){
      resultVec.pop_back();
      const auto splitRGB = [&](auto itkImage) {
        using RGBAPixelType = typename std::remove_pointer<decltype(itkImage)>::type::PixelType;
        using PixelType = typename RGBAPixelType::ComponentType;
        const auto newPixelType = mitk::MakePixelType<itk::Image<PixelType, 3>>(1);

        for (unsigned int i = 0; i < numComponents; ++i)
        {
          auto slice = mitk::Image::New();
          slice->Initialize(newPixelType, dims.size(), dims.data());
          mitk::ImagePixelWriteAccessor<PixelType, 3> sAcc(slice);
          std::transform(itkImage->GetBufferPointer(),
                          itkImage->GetBufferPointer() + N,
                          sAcc.GetData(),
                          [&i](const RGBAPixelType &d) { return d.GetElement(i); });
        
        
          slice->SetSpacing(spacing);
          resultVec.push_back(slice);
        }
      };

      AccessFixedPixelTypeByItk(image, splitRGB, RGB_PIXEL_TYPE_SEQUENCE);

      if (us::any_cast<bool>(GetReaderOption(C_REMOVE_DUPLICATES)))
      {
        decltype(resultVec) noDuplicates;
        for (unsigned int i = 0; i < numComponents; ++i)
        {
          auto tmp = resultVec[i];
          using namespace std;
          if (any_of(begin(noDuplicates), end(noDuplicates), [&tmp](auto current) {
                return mitk::Equal(*dynamic_cast<mitk::Image *>(tmp.GetPointer()),
                                   *dynamic_cast<mitk::Image *>(current.GetPointer()),
                                   mitk::eps,
                                   false);
              }))
            continue; // do no add iamge to vector if this and previous are equal eg RGB channels
          else
            noDuplicates.push_back(tmp);
        }
        resultVec = std::move(noDuplicates);
      }

    }


    else if (us::any_cast<bool>(GetReaderOption(A_IMPORT_VOLUME)) && image->GetDimension() == 2)
    {
      resultVec.pop_back();
      const auto lambda = [&](auto itkImage) {
        using PixelType = typename std::remove_pointer<decltype(itkImage)>::type::PixelType;
        const auto newPixelType = mitk::MakePixelType<itk::Image<PixelType, 3>>(1);
        const std::array<unsigned int, 3> dims{D[0], D[1], 1};

        auto slice = mitk::Image::New();
        slice->Initialize(newPixelType, 3, dims.data());
        mitk::ImagePixelWriteAccessor<PixelType, 3> sAcc(slice);
        std::copy(itkImage->GetBufferPointer(), itkImage->GetBufferPointer() + N, sAcc.GetData());

        slice->SetSpacing(spacing);
        resultVec.push_back(slice);
      };
      AccessByItk(image, (lambda));
    }
    else
    {
      // only spacing changed
      resultVec.push_back(image);
    }

    return resultVec;
  } // namespace m2

  mitk::DataStorage::SetOfObjects::Pointer MicroscopyTiffImageIO::Read(mitk::DataStorage &ds)
  {
    mitk::DataStorage::SetOfObjects::Pointer result = mitk::DataStorage::SetOfObjects::New();
    std::vector<mitk::BaseData::Pointer> data = this->mitk::AbstractFileReader::Read();
    unsigned int i = 0;
    for (auto iter = data.begin(); iter != data.end(); ++iter)
    {
      mitk::DataNode::Pointer node = mitk::DataNode::New();
      node->SetData(*iter);
      this->SetDefaultDataNodeProperties(node, this->GetInputLocation());
      node->SetName(node->GetName() + "_" + std::to_string(i));
      ++i;
      ds.Add(node);
      result->InsertElement(result->Size(), node);
    }
    return result;
  }

  void MicroscopyTiffImageIO::Write() {}

  mitk::AbstractFileIO::ConfidenceLevel MicroscopyTiffImageIO::GetReaderConfidenceLevel() const
  {
    return ConfidenceLevel::Supported;
  }

  mitk::AbstractFileIO::ConfidenceLevel MicroscopyTiffImageIO::GetWriterConfidenceLevel() const
  {
    return ConfidenceLevel::Unsupported;
  }

  MicroscopyTiffImageIO *MicroscopyTiffImageIO::IOClone() const { return new MicroscopyTiffImageIO(*this); }

} // namespace m2
