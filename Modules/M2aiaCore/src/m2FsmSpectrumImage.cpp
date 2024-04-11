/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <m2FsmSpectrumImage.h>
#include <m2Process.hpp>
#include <m2Timer.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>
#include <mitkProperties.h>
#include <signal/m2Baseline.h>
#include <signal/m2Normalization.h>
#include <signal/m2PeakDetection.h>
#include <signal/m2Pooling.h>
#include <signal/m2RunningMedian.h>
#include <signal/m2Smoothing.h>

void m2::FsmSpectrumImage::GetImage(double cmInv, double tol, const mitk::Image *mask, mitk::Image *destImage) const
{

  AccessByItk(destImage, [](auto itkImg) { itkImg->FillBuffer(0); });
  using namespace m2;
  // accessors
  mitk::ImagePixelWriteAccessor<DisplayImagePixelType, 3> imageAccess(destImage);

  std::shared_ptr<mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType, 3>> maskAccess;

  MITK_INFO("FSM") << "Image generation started!";

  if (mask)
  {
    maskAccess.reset(new mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType, 3>(mask));
    MITK_INFO << "> Use mask image";
  }
  
  GetPropertyList()->SetProperty("cm¯¹", mitk::DoubleProperty::New(cmInv));
  GetPropertyList()->SetProperty("tol", mitk::DoubleProperty::New(tol));
  

  const auto &xs = GetXAxis();
  std::vector<double> kernel;

  // Profile (continuous) spectrum

  const auto subRes = m2::Signal::Subrange(xs, cmInv - tol, cmInv + tol);
  const unsigned long n = m_Spectra.size();
  // map all spectra to several threads for processing
  const unsigned int t = GetNumberOfThreads();
  

  m2::Process::Map(n,
                   t,
                   [&](auto /*id*/, auto a, auto b)
                   {
                     for (unsigned int i = a; i < b; ++i)
                     {
                       auto &spectrum = m_Spectra[i];
                       auto &ys = spectrum.data;
                       auto s = std::next(std::begin(ys), subRes.first);
                       auto e = std::next(std::begin(ys), subRes.first + subRes.second);

                      //  if (maskAccess && maskAccess->GetPixelByIndex(spectrum.index) == 0)
                      //  {
                      //    imageAccess.SetPixelByIndex(spectrum.index, 0);
                      //    continue;
                      //  }

                      
                      imageAccess.SetPixelByIndex(spectrum.index, Signal::RangePooling<float>(s, e, GetRangePoolingStrategy()));
                     }
                   });
}

void m2::FsmSpectrumImage::InitializeProcessor()
{
  // this->m_Processor.reset((m2::ISpectrumImageSource *)new FsmProcessor(this));
}

void m2::FsmSpectrumImage::InitializeGeometry()
{
  
  std::array<itk::SizeValueType, 3> imageSize = {GetPropertyValue<unsigned>("dim_x"), // n_x
                                                 GetPropertyValue<unsigned>("dim_y"), // n_y
                                                 GetPropertyValue<unsigned>("dim_z")};

  std::array<double, 3> imageOrigin = {
    GetPropertyValue<double>("[IMS:1000053] absolute position offset x") * 0.001, // x_init
    GetPropertyValue<double>("[IMS:1000054] absolute position offset y") * 0.001, // y_init
    GetPropertyValue<double>("absolute position offset z") * 0.001};

  using ImageType = itk::Image<m2::DisplayImagePixelType, 3>;
  auto itkIonImage = ImageType::New();
  ImageType::IndexType idx;
  ImageType::SizeType size;

  idx.Fill(0);

  for (unsigned int i = 0; i < imageSize.size(); i++)
    size[i] = imageSize[i];
  ImageType::RegionType region(idx, size);
  itkIonImage->SetRegions(region);
  itkIonImage->Allocate();
  itkIonImage->FillBuffer(0);

  auto s = itkIonImage->GetSpacing();
  auto o = itkIonImage->GetOrigin();
  o[0] = imageOrigin[0];
  o[1] = imageOrigin[1];
  o[2] = imageOrigin[2];

  //
  s[0] = GetPropertyValue<double>("spacing_x"); // x_delta
  s[1] = GetPropertyValue<double>("spacing_y"); // y_delta
  s[2] = GetPropertyValue<double>("spacing_z");

  auto d = itkIonImage->GetDirection();
  d[0][0] = -1;

  itkIonImage->SetSpacing(s);
  itkIonImage->SetOrigin(o);
  itkIonImage->SetDirection(d);

  {
    using LocalImageType = itk::Image<m2::DisplayImagePixelType, 3>;
    auto caster = itk::CastImageFilter<ImageType, LocalImageType>::New();
    caster->SetInput(itkIonImage);
    caster->Update();
    InitializeByItk(caster->GetOutput());

    mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> acc(this);
    std::memset(acc.GetData(), 0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(m2::DisplayImagePixelType));
  }

  {
    using LocalImageType = itk::Image<m2::IndexImagePixelType, 3>;
    auto caster = itk::CastImageFilter<ImageType, LocalImageType>::New();
    caster->SetInput(itkIonImage);
    caster->Update();
    auto indexImage = mitk::Image::New();
    SetIndexImage(indexImage);
    indexImage->InitializeByItk(caster->GetOutput());

    mitk::ImagePixelWriteAccessor<m2::IndexImagePixelType, 3> acc(indexImage);
    std::memset(acc.GetData(), 0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(m2::IndexImagePixelType));
  }

  {
    mitk::LabelSetImage::Pointer image = mitk::LabelSetImage::New();
    SetMaskImage(image.GetPointer());

    image->Initialize((mitk::Image *)this);
    auto ls = image->GetActiveLabelSet();

    mitk::Color color;
    color.Set(0, 1, 0);
    auto label = mitk::Label::New();
    label->SetColor(color);
    label->SetName("Valid Spectrum");
    label->SetOpacity(0.0);
    label->SetLocked(true);
    label->SetValue(1);
    ls->AddLabel(label);
  }

  mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> acc(this);
  auto max_dim0 = GetDimensions()[0];
  auto max_dim1 = GetDimensions()[1];
  acc.SetPixelByIndex({0, 0, 0}, 1);
  acc.SetPixelByIndex({0, max_dim1 - 1, 0}, max_dim1 / 2);
  acc.SetPixelByIndex({max_dim0 - 1, 0, 0}, max_dim0 / 2);
  acc.SetPixelByIndex({max_dim0 - 1, max_dim1 - 1, 0}, max_dim1 + max_dim0);

  this->SetImageGeometryInitialized(true);
}

void m2::FsmSpectrumImage::InitializeImageAccess()
{
  using namespace m2;

  auto accMask = std::make_shared<mitk::ImagePixelWriteAccessor<mitk::LabelSetImage::PixelType, 3>>(GetMaskImage());
  auto accIndex = std::make_shared<mitk::ImagePixelWriteAccessor<m2::IndexImagePixelType, 3>>(GetIndexImage());
  // auto accNorm = std::make_shared<mitk::ImagePixelWriteAccessor<m2::NormImagePixelType,
  // 3>>(GetNormalizationImage());

  auto &xs = GetXAxis();

  // ----- PreProcess -----

  // if the data are available as continuous data with equivalent mz axis for all
  // spectra, we can calculate the skyline, sum and mean spectrum over the image
  std::vector<std::vector<double>> skylineT;
  std::vector<std::vector<double>> sumT;

  SetPropertyValue<unsigned>("m2aia.xs.n", xs.size());
  SetPropertyValue<double>("m2aia.xs.min", xs.front());
  SetPropertyValue<double>("m2aia.xs.max", xs.back());

  skylineT.resize(GetNumberOfThreads(), std::vector<double>(xs.size(), 0));
  sumT.resize(GetNumberOfThreads(), std::vector<double>(xs.size(), 0));

  // m2::Timer t("Initialize image");

  m2::Process::Map(
    GetSpectra().size(),
    GetNumberOfThreads(),
    [&](unsigned int t, unsigned int a, unsigned int b)
    {
      m2::Signal::SmoothingFunctor<float> Smoother;
      Smoother.Initialize(GetSmoothingStrategy(), GetSmoothingHalfWindowSize());

      m2::Signal::BaselineFunctor<float> BaselineSubtractor;
      BaselineSubtractor.Initialize(GetBaselineCorrectionStrategy(), GetBaseLineCorrectionHalfWindowSize());

      std::vector<float> baseline(xs.size());

      auto &spectra = GetSpectra();

      // const auto divides = [&val](const auto &a) { return a / val; };
      const auto maximum = [](const auto &a, const auto &b) { return a > b ? a : b; };
      const auto plus = std::plus<>();

      for (unsigned long int i = a; i < b; i++)
      {
        auto &spectrum = spectra[i];
        auto &ys = spectrum.data;

        Smoother(std::begin(ys), std::end(ys));
        BaselineSubtractor(std::begin(ys), std::end(ys), std::begin(baseline));

        std::transform(std::begin(ys), std::end(ys), sumT.at(t).begin(), sumT.at(t).begin(), plus);
        std::transform(std::begin(ys), std::end(ys), skylineT.at(t).begin(), skylineT.at(t).begin(), maximum);
      }
    });

  const auto &spectra = GetSpectra();
  m2::Process::Map(spectra.size(),
                   GetNumberOfThreads(),
                   [&](unsigned int /*t*/, unsigned int a, unsigned int b)
                   {
                     for (unsigned int i = a; i < b; i++)
                     {
                       const auto &spectrum = spectra[i];
                       accIndex->SetPixelByIndex(spectrum.index, i);
                       accMask->SetPixelByIndex(spectrum.index, 1);
                     }
                   });

  auto &skyline = GetSkylineSpectrum();
  skyline.resize(xs.size(), 0);
  for (unsigned int t = 0; t < GetNumberOfThreads(); ++t)
    std::transform(skylineT[t].begin(),
                   skylineT[t].end(),
                   skyline.begin(),
                   skyline.begin(),
                   [](auto &a, auto &b) { return a > b ? a : b; });

  auto &mean = GetMeanSpectrum();
  auto &sum = GetSumSpectrum();

  mean.resize(xs.size(), 0);
  sum.resize(xs.size(), 0);

  // accumulate valid spectra defined by mask image
  auto N = std::accumulate(accMask->GetData(),
                           accMask->GetData() + GetSpectra().size(),
                           mitk::LabelSetImage::PixelType(0),
                           [](const auto &a, const auto &b) -> mitk::LabelSetImage::PixelType { return a + (b > 0); });

  for (unsigned int t = 0; t < GetNumberOfThreads(); ++t)
    std::transform(
      sumT[t].begin(), sumT[t].end(), sum.begin(), sum.begin(), [](const auto &a, const auto &b) { return a + b; });
  std::transform(sum.begin(), sum.end(), mean.begin(), [&N](const auto &a) { return a / double(N); });


  this->SetImageAccessInitialized(true);
}


m2::FsmSpectrumImage::~FsmSpectrumImage()
{
  MITK_INFO << GetStaticNameOfClass() << " destroyed!";
}

m2::FsmSpectrumImage::FsmSpectrumImage()
{
  MITK_INFO << GetStaticNameOfClass() << " created!";

  m_SpectrumType.XAxisLabel = "cm¯¹";
  m_SpectrumType.Format = m2::SpectrumFormat::ContinuousProfile;
}
