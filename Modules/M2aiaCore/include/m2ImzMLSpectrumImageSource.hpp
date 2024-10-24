/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/
#pragma once

#include <M2aiaCoreExports.h>
#include <itkCastImageFilter.h>
#include <m2ISpectrumImageSource.h>
#include <m2ImzMLSpectrumImage.h>
#include <m2CoreCommon.h>
#include <m2Process.hpp>
#include <m2Timer.h>
#include <mitkImageAccessByItk.h>
#include <mitkCoreServices.h>
#include <mitkIPreferences.h>
#include <mitkIPreferencesService.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLabelSetImage.h>
#include <mitkProperties.h>
#include <mutex>
#include <signal/m2Baseline.h>
#include <signal/m2Morphology.h>
#include <signal/m2Normalization.h>
#include <signal/m2PeakDetection.h>
#include <signal/m2Pooling.h>
#include <signal/m2RunningMedian.h>
#include <signal/m2Smoothing.h>
#include <signal/m2Transformer.h>
namespace m2
{

  template <typename InputIterator, typename AccumulateType, typename BinaryOperation, typename Predicate>
  const AccumulateType accumulate_if(InputIterator first,
                                    const InputIterator last,
                                    AccumulateType init,
                                    BinaryOperation &&binary_op,
                                    Predicate &&predicate)
  {
    for (; first != last; ++first)
      if (predicate(*first))
        init = binary_op(init, *first);
    return init;
  }

  template <typename InputIterator, typename OutputIterator, typename UnaryOperation, typename Predicate>
  OutputIterator transform_if(
    InputIterator first, InputIterator last, OutputIterator dest_first, UnaryOperation unary_op, Predicate predicate)
  {
    for (; first != last; ++first, ++dest_first)
    {
      if (predicate(*first))
        *dest_first = unary_op(*first);
    }
    return dest_first;
  }

  template <typename InputIterator, typename Predicate>
  auto minmax_if(InputIterator first, InputIterator last, Predicate predicate)
  {
      using ValueType = typename std::iterator_traits<InputIterator>::value_type;
      bool found_any = false;
      ValueType min_value{}, max_value{};

      for (; first != last; ++first)
      {
          if (predicate(*first))
          {
              auto value = *first;
              if (!found_any)
              {
                  min_value = max_value = value;
                  found_any = true;
              }
              else
              {
                  if (value < min_value) min_value = value;
                  if (value > max_value) max_value = value;
              }
          }
      }

      if (!found_any)
          throw std::runtime_error("No elements satisfying predicate");

      return std::make_pair(min_value, max_value);
  }



  template <typename InputIterator, typename MaskIterator, typename OutputIterator>
  void StandardizeImage(InputIterator first,
                                  InputIterator last,
                                  MaskIterator first_mask,
                                  OutputIterator dest_first)
  {
    auto maskIt = first_mask;
    double N = accumulate_if(
      first, last, 0, [](auto s, auto) -> int { return s + 1; }, [&maskIt](auto) { return *maskIt++ > 0; });

    double sum = accumulate_if(first, last, 0.0, std::plus<>(), [&maskIt](auto) { return *maskIt++ > 0; });
    auto mean = sum / N;

    // reset mask iterator
    maskIt = first_mask;
    sum = accumulate_if(
      first,
      last,
      0.0,
      [mean](auto s, auto val) { return s + std::pow(val - mean, 2); },
      [&maskIt](auto) { return *maskIt++ > 0; });
    auto stddev = std::sqrt(sum / N);

    // reset mask iterator
    maskIt = first_mask;
    transform_if(
      first,
      last,
      dest_first,
      [mean, stddev](auto val) { return (val - mean) / stddev; },
      [&maskIt](auto) { return *maskIt++ > 0; });
  
  }


  template <typename InputIterator, typename MaskIterator, typename OutputIterator>
  void MinMaxNormalizeImage(InputIterator first,
                            InputIterator last,
                            MaskIterator first_mask,
                            OutputIterator dest_first)
  {
    auto maskIt = first_mask;
    const auto minMax = minmax_if(first, last, [&maskIt](auto) { return *maskIt++ > 0; });
    const auto minVal = minMax.first;
    const auto maxVal = minMax.second;

    maskIt = first_mask;
    transform_if(first, last, dest_first, [minVal, maxVal](auto val){return (val-minVal)/(maxVal-minVal);}, [&maskIt](auto) { return *maskIt++ > 0; });

  }

  struct BinaryDataAccessHelper
  {
    unsigned dataOffset;
    unsigned dataOffsetReverse;
    unsigned dataModifiedOffset;
    unsigned dataModifiedOffsetReverse;

    unsigned dataModifiedLength;
    unsigned dataPaddingLeft;
    unsigned dataPaddingRight;
  };

  /**
   * @brief Get the modifed binary data query range for spectral data access.
   * The range may be modified by requirements of signal processing methods - e.g. baseline correction, smoothing
   * filters
   */
  template <class MassAxisType>
  BinaryDataAccessHelper GetBinaryDataAccessHelper(const std::vector<MassAxisType> &xs,
                                                   double xRangeCenter,
                                                   double xRangeTol,
                                                   unsigned padding)
  {
    // MITK_INFO << "(GetBinaryDataAccessHelper) " << xs.empty() << " " << xs.size();

    // Image generation is based on range operations on the full spectrum for each pixel.
    // We are not interested in values outside of the range.
    // We can read only values of interest from the *.ibd file.
    // For this we have to manipulate the byte offset position.
    // 1) This is the full spectrum from start '|' to end '|'
    // |-----------------------------------------------------------------|

    // 2) Subrange from '(' to ')' with center 'c', offset left '>' and offset right '<'
    // |>>>>>>>>>>>>>>>(********c********)<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|
    auto subRes = m2::Signal::Subrange(xs, xRangeCenter - xRangeTol, xRangeCenter + xRangeTol);

    m2::BinaryDataAccessHelper offsetHelper;
    offsetHelper.dataOffset = subRes.first;
    // MITK_INFO << "(GetBinaryDataAccessHelper) [offsetHelper.dataOffset]" << offsetHelper.dataOffset;

    offsetHelper.dataOffsetReverse = (xs.size() - (subRes.first + subRes.second));
    // MITK_INFO << "(GetBinaryDataAccessHelper) [offsetHelper.dataOffsetReverse]" << offsetHelper.dataOffsetReverse;

    // check if padding can be applied on the left side
    // 3) Pad the range for overcoming border-problems with kernel based operations
    // |>>>>>>>>>[^^^^^(********c********)^^^^^]<<<<<<<<<<<<<<<<<<<<<<<<<|
    if (offsetHelper.dataOffset / padding >= 1)
      offsetHelper.dataPaddingLeft = padding;
    else
      offsetHelper.dataPaddingLeft = offsetHelper.dataOffset;
    // MITK_INFO << "(GetBinaryDataAccessHelper) [offsetHelper.dataPaddingLeft]" << offsetHelper.dataPaddingLeft;

    // check if padding can be applied on the right side
    if (offsetHelper.dataOffsetReverse / padding >= 1)
      offsetHelper.dataPaddingRight = padding;
    else
      offsetHelper.dataPaddingRight = offsetHelper.dataOffsetReverse;
    // MITK_INFO << "(GetBinaryDataAccessHelper) [offsetHelper.dataPaddingRight]" << offsetHelper.dataPaddingRight;

    // 4) We read data from the *ibd from '[' to ']' using a padded left offset
    // Continue at 5.
    offsetHelper.dataModifiedLength = subRes.second + offsetHelper.dataPaddingLeft + offsetHelper.dataPaddingRight;
    // MITK_INFO << "(GetBinaryDataAccessHelper) [offsetHelper.dataModifiedLength]" << offsetHelper.dataModifiedLength;

    offsetHelper.dataModifiedOffset = subRes.first - offsetHelper.dataPaddingLeft;
    // MITK_INFO << "(GetBinaryDataAccessHelper) [offsetHelper.dataModifiedOffset]" << offsetHelper.dataModifiedOffset;

    return offsetHelper;
  }

} // namespace m2

namespace m2
{
  class ImzMLSpectrumImage;

  /**
   * @class ImzMLSpectrumImageSource
   * @brief A class that represents a source of spectrum images in the ImzML format.
   *
   */
  template <class MassAxisType, class IntensityType>
  class ImzMLSpectrumImageSource : public m2::ISpectrumImageSource
  {
  private:
    m2::ImzMLSpectrumImage *p;

  public:
    explicit ImzMLSpectrumImageSource(m2::ImzMLSpectrumImage *owner) : p(owner) {}
    virtual void GetImagePrivate(double mz, double tol, const mitk::Image *mask, mitk::Image *image);

    /**
     * @brief Initialize the geometry of the image.
     * This method initializes the geometry of the image by setting the image size, origin, spacing, and direction.
     *
     */
    void InitializeGeometry() override;

    /**
     * @brief Initialize the image access.
     * This method initializes the image access by setting the normalization image, mask image, and index image.
     */
    void InitializeImageAccess() override;

    virtual void InitializeImageAccessContinuousProfile();

    /**
     * @brief Provides optimized access to centroid data.
     * No binning is applied. Normalization factors are provided.
     */
    virtual void InitializeImageAccessContinuousCentroid();

    /**
     * @brief See InitializeImageAccessProcessedData()
     */
    virtual void InitializeImageAccessProcessedProfile();

    /**
     * @brief See InitializeImageAccessProcessedData()
     */
    virtual void InitializeImageAccessProcessedCentroid();

    /**
     * @brief Provides access to processed centroid and processed profile spectra.
     * It includes binning for visualization purposes of the overview spectra
     * and invokes the calculation of normalization factors, but no further
     * signal-processing is supported here.
     */
    virtual void InitializeImageAccessProcessedData();

    /**
     * @brief Calculate and store the normalization image
     */
    void InitializeNormalizationImage(m2::NormalizationStrategyType type) override;

    /**
     * @brief Convert binary data to a vector.
     * @tparam OffsetType Type of the offset.
     * @tparam LengthType Type of the length.
     * @tparam DataType Type of the data.
     * @param f Input file stream.
     * @param offset The offset to start reading from.
     * @param length The length of data to read.
     * @param vec Pointer to the vector to store the data.
     */
    template <class OffsetType, class LengthType, class DataType>
    static void binaryDataToVector(std::ifstream &f, OffsetType offset, LengthType length, DataType *vec) noexcept
    {
      f.seekg(offset);
      f.read((char *)vec, length * sizeof(DataType));
    }

    template <class ItXFirst, class ItXLast, class ItYFirst, class ItYLast>
    static inline double GetNormalizationFactor(
      m2::NormalizationStrategyType strategy, ItXFirst xFirst, ItXLast xLast, ItYFirst yFirst, ItYLast yLast)
    {
      using namespace std;
      switch (strategy)
      {
        case m2::NormalizationStrategyType::TIC:
          return m2::Signal::TotalIonCurrent(xFirst, xLast, yFirst);
        case m2::NormalizationStrategyType::Sum:
          return accumulate(yFirst, yLast, double(0.0));
        case m2::NormalizationStrategyType::Mean:
          return accumulate(yFirst, yLast, double(0.0)) / double(std::distance(yFirst, yLast));
        case m2::NormalizationStrategyType::Max:
          return *max_element(yFirst, yLast);
        case m2::NormalizationStrategyType::RMS:
          return m2::Signal::RootMeanSquare(yFirst, yLast);
        case m2::NormalizationStrategyType::None:
        case m2::NormalizationStrategyType::Internal:
        case m2::NormalizationStrategyType::External:
        default:
          return 1;
      }
    }

    using XIteratorType = typename std::vector<MassAxisType>::iterator;
    using YIteratorType = typename std::vector<IntensityType>::iterator;
    m2::Signal::SmoothingFunctor<IntensityType> m_Smoother;
    m2::Signal::BaselineFunctor<IntensityType> m_BaselineSubtractor;
    m2::Signal::IntensityTransformationFunctor<IntensityType> m_Transformer;

    virtual void GetYValues(unsigned int id, std::vector<float> &yd) { GetYValues<float>(id, yd); }
    virtual void GetYValues(unsigned int id, std::vector<double> &yd) { GetYValues<double>(id, yd); }
    virtual void GetXValues(unsigned int id, std::vector<float> &yd) { GetXValues<float>(id, yd); }
    virtual void GetXValues(unsigned int id, std::vector<double> &yd) { GetXValues<double>(id, yd); }

  private:
    template <class OutputType>
    void GetYValues(unsigned int id, std::vector<OutputType> &yd);
    template <class OutputType>
    void GetXValues(unsigned int id, std::vector<OutputType> &xd);
  };

} // namespace m2

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeNormalizationImage(
  m2::NormalizationStrategyType type)
{
  // initialize the normalization iamge
  auto image = p->GetNormalizationImage(type);

  // create a write accessor
  using WriteAccessorType = mitk::ImagePixelWriteAccessor<NormImagePixelType, 3>;
  auto accNorm = std::make_shared<WriteAccessorType>(image);

  // get individual spectrum meta data
  auto &spectra = p->GetSpectra();
  int threads = p->GetNumberOfThreads();
  using namespace std;

  // split image in individual regions and process in parallel each spectrum
  Process::Map(spectra.size(),
               threads,
               [&](unsigned int /*thread*/, unsigned int a, unsigned int b)
               {
                 ifstream f(p->GetBinaryDataPath(), ifstream::binary);
                 vector<MassAxisType> mzs;
                 vector<IntensityType> ints;

                 for (unsigned long int i = a; i < b; i++)
                 {
                   auto &spectrum = spectra[i];

                   mzs.resize(spectrum.mzLength);
                   ints.resize(spectrum.intLength);
                   binaryDataToVector(f, spectrum.mzOffset, spectrum.mzLength, mzs.data());
                   binaryDataToVector(f, spectrum.intOffset, spectrum.intLength, ints.data());

                   double v;
                   if (type == NormalizationStrategyType::Internal)
                     v = spectrum.inFileNormalizationFactor;
                   else
                     v = GetNormalizationFactor(type, begin(mzs), end(mzs), begin(ints), end(ints));

                   accNorm->SetPixelByIndex(spectrum.index, v);
                 }
               });
}

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::GetImagePrivate(double xRangeCenter,
                                                                                double xRangeTol,
                                                                                const mitk::Image *mask,
                                                                                mitk::Image *destImage)
{
  using namespace m2;

  m_Smoother.Initialize(p->GetSmoothingStrategy(), p->GetSmoothingHalfWindowSize());
  m_BaselineSubtractor.Initialize(p->GetBaselineCorrectionStrategy(), p->GetBaseLineCorrectionHalfWindowSize());
  m_Transformer.Initialize(p->GetIntensityTransformationStrategy());

  // MITK_INFO <<" p->GetIntensityTransformationStrategy() " << (unsigned int)(p->GetIntensityTransformationStrategy());
  // MITK_INFO <<" p->GetBaselineCorrectionStrategy() " << (unsigned int)(p->GetBaselineCorrectionStrategy());
  // MITK_INFO <<" p->GetSmoothingStrategy() " << (unsigned int)(p->GetSmoothingStrategy());

  std::shared_ptr<mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType, 3>> maskAccess;
  if (mask)
    maskAccess.reset(new mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType, 3>(mask));

  using ShiftImageAccessorType = mitk::ImagePixelReadAccessor<m2::ShiftImageType, 3>;
  std::shared_ptr<ShiftImageAccessorType> accShift;
  if(p->GetShiftImage())
    accShift = std::make_shared<ShiftImageAccessorType>(p->GetShiftImage());

  // MITK_INFO << "Use Mask " << mask;

  if (!destImage)
    mitkThrow() << "Please provide an image into which the data can be written.";

  // clear the image conten
  // {
  // AccessByItk(destImage, [](auto itkImg) { itkImg->FillBuffer(0); });
  // }

  // Get the normalization type
  const auto currentType = p->GetNormalizationStrategy();

  // Check normalization strategy

  mitk::ImagePixelReadAccessor<NormImagePixelType, 3> normAccess(p->GetNormalizationImage());
  mitk::ImagePixelWriteAccessor<DisplayImagePixelType, 3> imageAccess(destImage);

  {
    auto d = destImage->GetDimensions();
    auto N = std::accumulate(d, d+destImage->GetDimension(), 1, std::multiplies<>());
    auto dataPointer = imageAccess.GetData();
    std::fill(dataPointer, dataPointer + N, 0);
  
  }
    

  // check if the normalization image for a given type
  // was already initialized. Initialize the image if necessary.
  // MITK_INFO << "Normalization image status: " << p->GetNormalizationImageStatus(currentType);
  if (!p->GetNormalizationImageStatus(currentType))
  {
    // MITK_INFO << "Normalization image type: " << m2::NormalizationStrategyTypeNames.at((unsigned int)(currentType));
    InitializeNormalizationImage(currentType);
    p->SetNormalizationImageStatus(currentType, true);
  }

  // Get the profile type
  const auto spectrumType = p->GetSpectrumType();
  const auto threads = p->GetNumberOfThreads();

  p->SetProperty("m2aia.xs.selection.center", mitk::DoubleProperty::New(xRangeCenter));
  p->SetProperty("m2aia.xs.selection.tolerance", mitk::DoubleProperty::New(xRangeTol));

  // Access each spectrum with identical binary offset and length parameters
  if (spectrumType.Format == m2::SpectrumFormat::ContinuousProfile)
  {
    unsigned padding = 0;
    if (p->GetBaselineCorrectionStrategy() != m2::BaselineCorrectionType::None)
      padding = p->GetBaseLineCorrectionHalfWindowSize();

    const auto mzs = p->GetXAxis();
    auto binaryDataAccessHelper = GetBinaryDataAccessHelper<double>(mzs, xRangeCenter, xRangeTol, padding);

    const auto &spectra = p->GetSpectra();
    m2::Process::Map(
      spectra.size(),
      threads,
      [&](auto /*id*/, auto a, auto b)
      {
        // create a input stream for the binary data file
        std::ifstream f(p->GetBinaryDataPath(), std::iostream::binary);

        // prepare data vectors for raw data and processing data
        std::vector<IntensityType> ints(binaryDataAccessHelper.dataModifiedLength);
        std::vector<IntensityType> baseline(binaryDataAccessHelper.dataModifiedLength);

        // 5) (For a specific thread), save the true range positions '(' and ')'
        // for pooling in the data vector. Continue at 6.

        // ints contains the padded data.
        // |>>>>>>>>>[^^^^^(********c********)^^^^^]<<<<<<<<<<<<<<<<<<<<<<<<<|
        // s,e are the start and end of the data without padding

        // |>>>>>>>>>[^^^^^s********c********e^^^^^]<<<<<<<<<<<<<<<<<<<<<<<<<|
        auto s = std::next(std::begin(ints), binaryDataAccessHelper.dataPaddingLeft);
        auto e = std::prev(std::end(ints), binaryDataAccessHelper.dataPaddingRight);

        for (unsigned int i = a; i < b; ++i)
        {
          const auto &spectrum = spectra[i];

          // check if outside of mask
          if (maskAccess && maskAccess->GetPixelByIndex(spectrum.index) == 0)
          {
            imageAccess.SetPixelByIndex(spectrum.index, 0);
            continue;
          }
          // 6) access the binary data in the file.
          // - use the spectrum.intOffset to find spectrum data in the binary file
          // - add the offset to find the right subrange of the spectrum data
          auto binaryFileOffset =  spectrum.intOffset + binaryDataAccessHelper.dataModifiedOffset * sizeof(IntensityType);
            
          if(accShift) binaryFileOffset += accShift->GetPixelByIndex(spectrum.index)*sizeof(IntensityType);
          // access the binary data and read a (padded) subrange of the intensities (y values)
          binaryDataToVector(f, binaryFileOffset, binaryDataAccessHelper.dataModifiedLength, ints.data());

          // ----- Normalization

          IntensityType norm = normAccess.GetPixelByIndex(spectrum.index);
          std::transform(std::begin(ints), std::end(ints), std::begin(ints), [&norm](auto &v) { return v / norm; });

          // ----- Smoothing
          m_Smoother(std::begin(ints), std::end(ints));

          // ----- Baseline Substraction
          m_BaselineSubtractor(std::begin(ints), std::end(ints), std::begin(baseline));

          // ----- Intensity Transformation
          m_Transformer(std::begin(ints), std::end(ints));

          // ----- Pool the range
          const auto val = Signal::RangePooling<IntensityType>(s, e, p->GetRangePoolingStrategy());

          // finally set the pixel value
          imageAccess.SetPixelByIndex(spectrum.index, val);
        }
      });
  }

  else if (any(spectrumType.Format & (m2::SpectrumFormat::ContinuousCentroid | m2::SpectrumFormat::ProcessedCentroid |
                                      m2::SpectrumFormat::ProcessedProfile)))
  {
    const auto &spectra = p->GetSpectra();
    m2::Process::Map(
      spectra.size(),
      threads,
      [&](auto /*id*/, auto a, auto b)
      {
        std::ifstream f(p->GetBinaryDataPath(), std::iostream::binary);
        std::vector<IntensityType> ints;
        std::vector<MassAxisType> mzs;
        for (unsigned int i = a; i < b; ++i)
        {
          auto &spectrum = spectra[i];
          if (maskAccess && maskAccess->GetPixelByIndex(spectrum.index) == 0)
          {
            imageAccess.SetPixelByIndex(spectrum.index, 0);
            continue;
          }

          mzs.resize(spectrum.mzLength);
          binaryDataToVector(
            f, spectrum.mzOffset, spectrum.mzLength, mzs.data()); // !! read mass axis for each spectrum

          auto binaryDataAccessHelper = GetBinaryDataAccessHelper<MassAxisType>(mzs, xRangeCenter, xRangeTol, 0);
          ints.resize(binaryDataAccessHelper.dataModifiedLength);

          if (binaryDataAccessHelper.dataModifiedLength == 0)
          {
            imageAccess.SetPixelByIndex(spectrum.index, 0);
            continue;
          }

          const auto binaryFileOffset =
            spectrum.intOffset + binaryDataAccessHelper.dataModifiedOffset * sizeof(IntensityType);
          binaryDataToVector(f, binaryFileOffset, binaryDataAccessHelper.dataModifiedLength, ints.data());

          // TODO: Is it useful to normalize centroid data?
          IntensityType norm = normAccess.GetPixelByIndex(spectrum.index);
          std::transform(std::begin(ints), std::end(ints), std::begin(ints), [&norm](auto &v) { return v / norm; });

          auto val =
            Signal::RangePooling<IntensityType>(std::begin(ints), std::end(ints), p->GetRangePoolingStrategy());
          imageAccess.SetPixelByIndex(spectrum.index, val);
        }
      });
  }

 // normalize ion image post generation
  const auto bufferN = std::accumulate(destImage->GetDimensions(), destImage->GetDimensions() + 3, 1, std::multiplies<>());
  if(p->GetImageNormalizationStrategy() == m2::ImageNormalizationStrategyType::zScore)
    StandardizeImage(imageAccess.GetData(), imageAccess.GetData()+bufferN, maskAccess->GetData(), imageAccess.GetData());
  else if(p->GetImageNormalizationStrategy() == m2::ImageNormalizationStrategyType::MinMax)
    MinMaxNormalizeImage(imageAccess.GetData(), imageAccess.GetData()+bufferN, maskAccess->GetData(), imageAccess.GetData());
 
}



template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeGeometry()
{
  std::array<itk::SizeValueType, 3> imageSize = {
    p->GetPropertyValue<unsigned int>("[IMS:1000042] max count of pixels x"),
    p->GetPropertyValue<unsigned int>("[IMS:1000043] max count of pixels y"),
    p->GetPropertyValue<unsigned int>("max count of pixels z")};

  std::array<double, 3> imageOrigin = {p->GetPropertyValue<double>("[IMS:1000053] absolute position offset x"),
                                       p->GetPropertyValue<double>("[IMS:1000054] absolute position offset y"),
                                       p->GetPropertyValue<double>("absolute position offset z")};

  using ImageType = itk::Image<m2::DisplayImagePixelType, 3>;
  auto itkIonImage = ImageType::New();

  itkIonImage->SetRegions({{0, 0, 0}, {imageSize[0], imageSize[1], imageSize[2]}});
  itkIonImage->Allocate();
  itkIonImage->FillBuffer(0);

  auto s = itkIonImage->GetSpacing();
  auto o = itkIonImage->GetOrigin();
  o[0] = imageOrigin[0];
  o[1] = imageOrigin[1];
  o[2] = imageOrigin[2];

  //
  s[0] = p->GetPropertyValue<double>("[IMS:1000046] pixel size x");
  s[1] = p->GetPropertyValue<double>("[IMS:1000047] pixel size y");
  s[2] = p->GetPropertyValue<double>("pixel size z");

  auto d = itkIonImage->GetDirection();

  itkIonImage->SetSpacing(s);
  itkIonImage->SetOrigin(o);
  itkIonImage->SetDirection(d);

  {
    using LocalImageType = itk::Image<m2::DisplayImagePixelType, 3>;
    auto caster = itk::CastImageFilter<ImageType, LocalImageType>::New();
    caster->SetInput(itkIonImage);
    caster->Update();
    p->InitializeByItk(caster->GetOutput());

    mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> acc(p);
    std::memset(acc.GetData(), 0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(m2::DisplayImagePixelType));
  }

  {
    using LocalImageType = itk::Image<m2::IndexImagePixelType, 3>;
    auto caster = itk::CastImageFilter<ImageType, LocalImageType>::New();
    caster->SetInput(itkIonImage);
    caster->Update();
    auto indexImage = mitk::Image::New();
    p->SetIndexImage(indexImage);
    indexImage->InitializeByItk(caster->GetOutput());

    mitk::ImagePixelWriteAccessor<m2::IndexImagePixelType, 3> acc(indexImage);
    std::memset(acc.GetData(), 0, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(m2::IndexImagePixelType));
  }

  if(p->GetMaskImage().IsNull()){
    auto image = mitk::LabelSetImage::New();
    image->SetProperty("m2aia.mask.initialization", mitk::StringProperty::New("internal"));
    p->SetMaskImage(image.GetPointer());
    image->Initialize((mitk::Image *)p);

    mitk::Color color;
    color.Set(0.0, 1, 0.0);
    auto label = mitk::Label::New();
    label->SetColor(color);
    label->SetName("Valid");
    label->SetOpacity(0.0);
    label->SetLocked(true);
    label->SetValue(1);
    image->AddLabel(label,0);
  }else{
    auto image = p->GetMaskImage();
    image->SetProperty("m2aia.mask.initialization", mitk::StringProperty::New("external"));
  }

  auto max_dim0 = p->GetDimensions()[0];
  auto max_dim1 = p->GetDimensions()[1];
  std::vector<std::thread> threads;

  for (auto type : m2::NormalizationStrategyTypeList)
  {
    m2::Timer t("Initialization of the " + m2::NormalizationStrategyTypeNames.at((unsigned int)(type)) +
                " images took");
    t.printIf = [](m2::Timer::Duration d) -> bool { return d.count() > 0.5; };

    using LocalImageType = itk::Image<m2::NormImagePixelType, 3>;
    auto caster = itk::CastImageFilter<ImageType, LocalImageType>::New();
    caster->SetInput(itkIonImage);
    caster->Update();
    auto normImage = mitk::Image::New();
    p->SetNormalizationImage(normImage, type);
    normImage->InitializeByItk(caster->GetOutput());

    { // scope for the accessor
      mitk::ImagePixelWriteAccessor<m2::NormImagePixelType, 3> acc(normImage);
      std::memset(acc.GetData(), 1, imageSize[0] * imageSize[1] * imageSize[2] * sizeof(m2::NormImagePixelType));
    }

    p->SetNormalizationImageStatus(type, false);
    // }
  }

  mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> acc(p);
  acc.SetPixelByIndex({0, 0, 0}, 1);
  acc.SetPixelByIndex({0, max_dim1 - 1, 0}, max_dim1 / 2);
  acc.SetPixelByIndex({max_dim0 - 1, 0, 0}, max_dim0 / 2);
  acc.SetPixelByIndex({max_dim0 - 1, max_dim1 - 1, 0}, max_dim1 + max_dim0);
}

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeImageAccess()
{
  p->SetImageAccessInitialized(false);

  m_Transformer.Initialize(p->GetIntensityTransformationStrategy());
  m_BaselineSubtractor.Initialize(p->GetBaselineCorrectionStrategy(), p->GetBaseLineCorrectionHalfWindowSize());
  m_Smoother.Initialize(p->GetSmoothingStrategy(), p->GetSmoothingHalfWindowSize());

  //////////---------------------------
  const auto spectrumType = p->GetSpectrumType();
  const auto currentType = p->GetNormalizationStrategy();
  if (p->GetNormalizationImageStatus(currentType) == false)
  {
    MITK_INFO << "First access to the normalization image: "
              << m2::NormalizationStrategyTypeNames.at((unsigned int)currentType);
    p->InitializeNormalizationImage(currentType);
    p->SetNormalizationImageStatus(currentType, true);
  }


  if (spectrumType.Format == m2::SpectrumFormat::ProcessedProfile)
  {
    // mitkThrow() << m2::ImzMLSpectrumImage::GetStaticNameOfClass() << R"(
    // This ImzML file seems to contain profile spectra in a processed memory order.
    // This is not supported in M2aia! If there are really individual m/z axis for
    // each spectrum, please resample the m/z axis and create one that is commonly
    // used for all spectra. Save it as continuous ImzML!)";
    InitializeImageAccessProcessedProfile();
  }
  else if (spectrumType.Format == m2::SpectrumFormat::ContinuousProfile)
    InitializeImageAccessContinuousProfile();
  else if (spectrumType.Format == m2::SpectrumFormat::ProcessedCentroid)
    InitializeImageAccessProcessedCentroid();
  else if (spectrumType.Format == m2::SpectrumFormat::ContinuousCentroid)
    InitializeImageAccessContinuousCentroid();

  // DEFAULT
  // INITIALIZE MASK, INDEX, NORMALIZATION IMAGES
  std::shared_ptr<mitk::ImagePixelWriteAccessor<mitk::LabelSetImage::PixelType, 3>> accMask;
  
  auto prop = p->GetMaskImage()->GetProperty("m2aia.mask.initialization");
  if(prop && prop->GetValueAsString() == "internal"){
    accMask =  std::make_shared<mitk::ImagePixelWriteAccessor<mitk::LabelSetImage::PixelType, 3>>(p->GetMaskImage());
  }
  auto accIndex = std::make_shared<mitk::ImagePixelWriteAccessor<m2::IndexImagePixelType, 3>>(p->GetIndexImage());
  auto accNorm = std::make_shared<mitk::ImagePixelWriteAccessor<m2::NormImagePixelType, 3>>(p->GetNormalizationImage());

  const auto &spectra = p->GetSpectra();
  m2::Process::Map(spectra.size(),
                   p->GetNumberOfThreads(),
                   [&](unsigned int /*t*/, unsigned int a, unsigned int b)
                   {
                     for (unsigned int i = a; i < b; i++)
                     {
                       const auto &spectrum = spectra[i];

                       accIndex->SetPixelByIndex(spectrum.index, i);
                       
                       if(accMask) 
                        accMask->SetPixelByIndex(spectrum.index, 1);

                       // If it is a processed file, normalization maps are set to 1 - assuming that spectra were
                       // already processed if (any(importMode & (m2::SpectrumFormatType::ProcessedCentroid |
                       // m2::SpectrumFormatType::ProcessedProfile)))
                       //   accNorm->SetPixelByIndex(spectrum.index + source.m_Offset, 1);
                     }
                   });
  p->SetNumberOfValidPixels(spectra.size());
  p->SetImageAccessInitialized(true);
}

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeImageAccessProcessedProfile()
{
  // MITK_INFO("m2::ImzMLSpectrumImage") << "Start InitializeImageAccessProcessedProfile";
  InitializeImageAccessProcessedData();
}

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeImageAccessContinuousProfile()
{
  std::vector<std::vector<double>> skylineT;
  std::vector<std::vector<double>> sumT;
  std::vector<MassAxisType> mzs;

  const unsigned int threads = p->GetNumberOfThreads();
  const auto currentType = p->GetNormalizationStrategy();
  
  using NormalizationImageAccessorType = mitk::ImagePixelReadAccessor<m2::DisplayImagePixelType, 3>;
  auto accNorm = std::make_shared<NormalizationImageAccessorType>(p->GetNormalizationImage(currentType)); 
  
  using ShiftImageAccessorType = mitk::ImagePixelReadAccessor<m2::ShiftImageType, 3>;
  std::shared_ptr<ShiftImageAccessorType> accShift;
  if(p->GetShiftImage())
    accShift = std::make_shared<ShiftImageAccessorType>(p->GetShiftImage());
  
  auto &mzAxis = p->GetXAxis();

  int maxUpShift, maxDownShift;
  // load continuous x axis
  if(p->GetShiftImage()){
    // **** Prepare Shifted Images

    auto N = std::accumulate(p->GetDimensions(),p->GetDimensions()+3, uint32_t(1), std::multiplies<>());
    auto minMaxElement = std::minmax_element(accShift->GetData(), accShift->GetData() + N);
    maxUpShift = std::min(0,*minMaxElement.first); // (-)
    maxDownShift = std::max(0,*minMaxElement.second); // (+)
  
    const auto &spectra = p->GetSpectra();
    std::ifstream f(p->GetBinaryDataPath(), std::ios::binary);
    mzs.resize(spectra[0].mzLength,0);
    binaryDataToVector(f, spectra[0].mzOffset, spectra[0].mzLength, mzs.data());

    // **** Prepare Overview mzAxis
    auto &mzAxis = p->GetXAxis();
    const int L = spectra[0].mzLength + std::abs(maxDownShift) + std::abs(maxUpShift);
    mzAxis.clear();
    mzAxis.reserve(L);

    
    auto mzAxisIt = std::back_inserter(mzAxis);
    // **** extrapolate mz values (right)
    const auto lowMzDelta = mzs[1] - mzs[0];
    for(int i = 0 ; i < std::abs(maxDownShift); ++i)
      mzAxisIt = mzs.front() - lowMzDelta*(std::abs(maxDownShift) - i);

    // **** add original
    std::copy(std::begin(mzs), std::end(mzs), mzAxisIt);

    // **** extrapolate mz values (right)
    const auto highMzDelta = mzs[mzs.size()-2] - mzs[mzs.size()-1];
    for(int i = 0 ; i < std::abs(maxUpShift); ++i)
      mzAxisIt = mzs.back() + highMzDelta*(i+1);


  }else{ 
    const auto &spectra = p->GetSpectra();
    std::ifstream f(p->GetBinaryDataPath(), std::ios::binary);
    mzs.resize(spectra[0].mzLength);
    binaryDataToVector(f, spectra[0].mzOffset, spectra[0].mzLength, mzs.data());
   
    mzAxis.clear();
    std::copy(std::begin(mzs), std::end(mzs), std::back_inserter(mzAxis));
  }

  p->SetPropertyValue<unsigned>("m2aia.xs.n", mzAxis.size());
  p->SetPropertyValue<double>("m2aia.xs.min", mzAxis.front());
  p->SetPropertyValue<double>("m2aia.xs.max", mzAxis.back());

  skylineT.resize(threads, std::vector<double>(mzAxis.size(), 0));
  sumT.resize(threads, std::vector<double>(mzAxis.size(), 0));

  const auto Maximum = [](const auto &a, const auto &b) { return a > b ? a : b; };
  const auto plus = std::plus<>();
  

  {
    auto &spectra = p->GetSpectra();

    m2::Process::Map(
      spectra.size(),
      p->GetNumberOfThreads(),
      [&](unsigned int t, unsigned int a, unsigned int b)
      {
        std::vector<IntensityType> baseline(mzs.size(), 0);
        std::vector<IntensityType> ints(mzs.size(), 0);
        std::ifstream f(p->GetBinaryDataPath(), std::ifstream::binary);

        double nFac = 1.0;
        for (unsigned long int i = a; i < b; i++)
        {
          auto &spectrum = spectra[i];
          // Read data from file ------------
          ints.resize(spectrum.intLength);
          binaryDataToVector(f, spectrum.intOffset, spectrum.intLength, ints.data());

          try
          {
            nFac = accNorm->GetPixelByIndex(spectrum.index);
            // Signal processing
            std::transform(
              std::begin(ints), std::end(ints), std::begin(ints), [&nFac](const auto &a) { return a / nFac; });
          }
          catch (const std::exception &e)
          {
            MITK_INFO << "Could not access normalization image";
          }

          m_Smoother(std::begin(ints), std::end(ints));
          m_BaselineSubtractor(std::begin(ints), std::end(ints), std::begin(baseline));
          m_Transformer(std::begin(ints), std::end(ints));


          if(p->GetShiftImage()){
            double shift = accShift->GetPixelByIndex(spectrum.index);              
            const auto insertPosition = std::abs(maxDownShift) - shift; 
            std::transform(std::begin(ints), std::end(ints), sumT.at(t).begin()+insertPosition, sumT.at(t).begin()+insertPosition, plus);
            std::transform(std::begin(ints), std::end(ints), skylineT.at(t).begin()+insertPosition, skylineT.at(t).begin()+insertPosition, Maximum);
          }else{
          std::transform(std::begin(ints), std::end(ints), sumT.at(t).begin(), sumT.at(t).begin(), plus);
          std::transform(std::begin(ints), std::end(ints), skylineT.at(t).begin(), skylineT.at(t).begin(), Maximum);
          }
        }
      });
  }


  auto &skyline = p->GetSkylineSpectrum();
  auto &mean = p->GetMeanSpectrum();
  auto &sum = p->GetSumSpectrum();

  
  skyline.clear();
  skyline.resize(mzAxis.size(), 0);
  
  for (unsigned int t = 0; t < threads; ++t)
    std::transform(skylineT[t].begin(), skylineT[t].end(), skyline.begin(), skyline.begin(), Maximum);

MITK_INFO << "mzAxis Size " << mzAxis.size();
  mean.clear();
  mean.resize(mzAxis.size(), 0);

  sum.clear();
  sum.resize(mzAxis.size(), 0);


  auto N = p->GetSpectra().size();
  for (unsigned int t = 0; t < threads; ++t)
    std::transform(sumT[t].begin(), sumT[t].end(), sum.begin(), sum.begin(), plus);
  std::transform(sum.begin(), sum.end(), mean.begin(), [&](auto &a) { return a / double(N); });

}

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeImageAccessContinuousCentroid()
{
  std::vector<MassAxisType> mzs;

  { // load continuous x axis
    const auto &spectra = p->GetSpectra();
    mzs.resize(spectra[0].mzLength);

    std::ifstream f(p->GetBinaryDataPath(), std::ios::binary);
    binaryDataToVector(f, spectra[0].mzOffset, spectra[0].mzLength, mzs.data());

    auto &massAxis = p->GetXAxis();
    massAxis.clear();
    std::copy(std::begin(mzs), std::end(mzs), std::back_inserter(massAxis));
    p->SetPropertyValue<unsigned>("m2aia.xs.n", mzs.size());
    p->SetPropertyValue<double>("m2aia.xs.min", mzs.front());
    p->SetPropertyValue<double>("m2aia.xs.max", mzs.back());
  }

  // each thread pixelwise accumulate peaks
  std::vector<std::vector<m2::Interval>> peaksT(p->GetNumberOfThreads());
  for (auto &peaks : peaksT)
    peaks.resize(mzs.size());

  // initialize normalization image accessor
  auto currentType = p->GetNormalizationStrategy();
  using NormImageReadAccess = mitk::ImagePixelReadAccessor<NormImagePixelType, 3>;
  NormImageReadAccess accNorm(p->GetNormalizationImage(currentType));

  auto &spectra = p->GetSpectra();

  m2::Process::Map(spectra.size(),
                   p->GetNumberOfThreads(),
                   [&](unsigned int t, unsigned int a, unsigned int b)
                   {
                     std::ifstream f;
                     f.open(p->GetBinaryDataPath(), std::ios::binary);

                     std::vector<IntensityType> ints;

                     for (unsigned i = a; i < b; i++)
                     {
                       auto iO = spectra[i].intOffset;
                       auto iL = spectra[i].intLength;
                       ints.resize(iL);
                       
                       binaryDataToVector(f, iO, iL, ints.data());

                       auto nFac = accNorm.GetPixelByIndex(spectra[i].index);
                       std::transform(
                         std::begin(ints), std::end(ints), std::begin(ints), [&nFac](auto &v) { return v / nFac; });

                       for (size_t i = 0; i < mzs.size(); ++i)
                       {
                         //  peaksT[t][i].index(i);
                         peaksT[t][i].x(mzs[i]);
                         peaksT[t][i].y(ints[i]);
                       }
                     }

                     f.close();
                   });

  auto &skyline = p->GetSkylineSpectrum();
  auto &sum = p->GetSumSpectrum();
  auto &mean = p->GetMeanSpectrum();
  skyline.clear();
  sum.clear();
  mean.clear();

  // merge all peaks
  std::vector<m2::Interval> finalPeaks(mzs.size());
  for (auto &peaks : peaksT)
    for (size_t i = 0; i < peaks.size(); ++i)
    {
      finalPeaks[i] += peaks[i];
    }

  for (const auto &peak : finalPeaks)
  {
    skyline.push_back(peak.y.max());
    sum.push_back(peak.y.sum());
    mean.push_back(peak.y.mean());    
  }
}

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeImageAccessProcessedCentroid()
{
  InitializeImageAccessProcessedData();
}

template <class MassAxisType, class IntensityType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::InitializeImageAccessProcessedData()
{
  // initialize normalization image accessor
  const auto currentType = p->GetNormalizationStrategy();
  if (p->GetNormalizationImageStatus(currentType) == false)
  {
    // MITK_INFO << "First access to the normalization image: "
    //           << m2::NormalizationStrategyTypeNames.at((unsigned int)currentType);
    p->InitializeNormalizationImage(currentType);
    p->SetNormalizationImageStatus(currentType, true);
  }

  using NormImageReadAccess = mitk::ImagePixelReadAccessor<NormImagePixelType, 3>;
  NormImageReadAccess accNorm(p->GetNormalizationImage(currentType));
  std::vector<std::list<m2::Interval>> peaksT(p->GetNumberOfThreads());

  int binsN;
  // int minHits;
  if (auto *preferencesService = mitk::CoreServices::GetPreferencesService())
    if (auto *preferences = preferencesService->GetSystemPreferences())
    {
      binsN = preferences->GetInt("m2aia.view.spectrum.bins", 15000);
      // minHits = preferences->GetInt("m2aia.view.spectrum.minimum.hits", 30);
      MITK_INFO << "Generating processed Centroid/Profile imzML overview spectra )";
      MITK_INFO << "Number of bins: " << binsN << " (can be changed in the preferences: Window->Preferences->M2aia)";
    }

  auto &spectra = p->GetSpectra();
  const auto &T = p->GetNumberOfThreads();
  std::vector<double> xMin(T, std::numeric_limits<double>::max());
  std::vector<double> xMax(T, std::numeric_limits<double>::min());

  std::vector<std::vector<double>> yT(T, std::vector<double>(binsN, 0));
  std::vector<std::vector<double>> yMaxT(T, std::vector<double>(binsN, 0));
  std::vector<std::vector<unsigned int>> hT(T, std::vector<unsigned int>(binsN, 0));
  std::vector<std::vector<double>> xT(T, std::vector<double>(binsN, 0));

  // Find min max x values
  m2::Process::Map(spectra.size(),
                   T,
                   [&](unsigned int t, unsigned int a, unsigned int b)
                   {
                     std::ifstream f(p->GetBinaryDataPath(), std::ios::binary);
                     std::vector<MassAxisType> mzs;
                     //  std::list<m2::Interval> peaks, tempList;
                     // find x min/max
                     for (unsigned i = a; i < b; i++)
                     {
                       const auto &mzO = spectra[i].mzOffset;
                       const auto &mzL = spectra[i].mzLength;
                       mzs.resize(mzL);
                       binaryDataToVector(f, mzO, mzL, mzs.data());
                       xMin[t] = std::min(xMin[t], (double)mzs.front());
                       xMax[t] = std::max(xMax[t], (double)mzs.back());
                     }
                   });

  // find overall min/max
  double binSize = 1;
  double max = std::numeric_limits<double>::min();
  double min = std::numeric_limits<double>::max();

  max = *std::max_element(std::begin(xMax), std::end(xMax));
  min = *std::min_element(std::begin(xMin), std::end(xMin));
  binSize = (max - min) / double(binsN);


  m2::Process::Map(spectra.size(),
                   T,
                   [&](unsigned int t, unsigned int a, unsigned int b)
                   {
                     std::ifstream f(p->GetBinaryDataPath(), std::ios::binary);
                     std::vector<MassAxisType> mzs;
                     std::vector<IntensityType> ints;

                     for (unsigned i = a; i < b; i++)
                     {
                       auto &spectrum = spectra[i];
                       const auto &mzO = spectrum.mzOffset;
                       const auto &mzL = spectrum.mzLength;
                       mzs.resize(mzL);
                       binaryDataToVector(f, mzO, mzL, mzs.data());

                       const auto &intO = spectrum.intOffset;
                       const auto &intL = spectrum.intLength;
                       ints.resize(intL);
                       binaryDataToVector(f, intO, intL, ints.data());

                       // Normalization
                       if (p->GetNormalizationStrategy() != m2::NormalizationStrategyType::None)
                       {
                         double nFac = accNorm.GetPixelByIndex(spectrum.index);
                         std::transform(
                           std::begin(ints), std::end(ints), std::begin(ints), [&nFac](auto &v) { return v / nFac; });
                       }

                       for (unsigned int k = 0; k < mzs.size(); ++k)
                       {
                         // find index of the bin for the k'th m/z value of the pixel
                         auto j = (long)((mzs[k] - min) / binSize);

                         if (j >= binsN)
                           j = binsN - 1;
                         else if (j < 0)
                           j = 0;

                         xT[t][j] += mzs[k];                                   // mass sum
                         yT[t][j] += ints[k] < 10e-256 ? 0 : ints[k];          // intensitiy sum
                         yMaxT[t][j] = std::max(yMaxT[t][j], double(ints[k])); // intensitiy max
                         hT[t][j]++;                                           // hits
                       }
                     }

                     f.close();
                   });

  // REDUCE
  for (unsigned int i = 1; i < T; ++i)
  {
    for (int k = 0; k < binsN; ++k)
    {
      xT[0][k] += xT[i][k];
      yT[0][k] += yT[i][k];
      yMaxT[0][k] = std::max(yMaxT[0][k], yMaxT[i][k]);
      hT[0][k] += hT[i][k];
    }
  }

  // for(int i = 0 ; i < binsN; ++i){
  //   if(hT[0][i] > 0){
  //     MITK_INFO << xT[0][i]/(double)hT[0][i] << " " <<  yT[0][i] << " max(" << yMaxT[0][i] << ") [" << hT[0][i] <<
  //     "]";
  //   }
  // }

  auto &mzAxis = p->GetXAxis();
  auto &sum = p->GetSumSpectrum();
  auto &mean = p->GetMeanSpectrum();
  auto &skyline = p->GetSkylineSpectrum();

  mzAxis.clear();
  sum.clear();
  mean.clear();
  skyline.clear();

  for (int k = 0; k < binsN; ++k)
  {
    if (hT[0][k] > 0)
    {
      mzAxis.push_back(xT[0][k] / (double)hT[0][k]);
      sum.push_back(yT[0][k]);
      mean.push_back(yT[0][k] / (double)hT[0][k]);
      skyline.push_back(yMaxT[0][k]);
    }
  }

  xT.clear();
  yT.clear();
  yMaxT.clear();
  hT.clear();

  p->SetPropertyValue<double>("m2aia.xs.min", mzAxis.front());
  p->SetPropertyValue<double>("m2aia.xs.max", mzAxis.back());
  p->SetPropertyValue<unsigned>("m2aia.xs.n", mzAxis.size());
}

template <class MassAxisType, class IntensityType>
template <class OutputType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::GetXValues(unsigned int id, std::vector<OutputType> &xd)
{
  std::ifstream f(p->GetBinaryDataPath(), std::ios::binary);

  const auto &spectrum = p->GetSpectra()[id];
  const auto &length = spectrum.mzLength;
  const auto &offset = spectrum.mzOffset;

  if (std::is_same<MassAxisType, OutputType>::value)
  {
    xd.resize(length);
    binaryDataToVector(f, offset, length, xd.data());
  }
  else
  {
    std::vector<MassAxisType> xs;
    xs.resize(length);
    binaryDataToVector(f, offset, length, xs.data());
    // copy and convert
    xd.resize(length);
    std::copy(std::begin(xs), std::end(xs), std::begin(xd));
  }
}

template <class MassAxisType, class IntensityType>
template <class OutputType>
void m2::ImzMLSpectrumImageSource<MassAxisType, IntensityType>::GetYValues(unsigned int id, std::vector<OutputType> &yd)
{
  std::ifstream f(p->GetBinaryDataPath(), std::ios::binary);

  const auto &spectrum = p->GetSpectra()[id];
  const auto &length = spectrum.intLength;
  const auto &offset = spectrum.intOffset;

  mitk::ImagePixelReadAccessor<m2::NormImagePixelType, 3> normAccess(p->GetNormalizationImage());

  {
    std::vector<IntensityType> ys;
    ys.resize(length);
    binaryDataToVector(f, offset, length, ys.data());

    IntensityType norm = normAccess.GetPixelByIndex(spectrum.index);
    std::transform(std::begin(ys), std::end(ys), std::begin(ys), [&norm](auto &v) { return v / norm; });

    // ----- Smoothing
    m_Smoother(std::begin(ys), std::end(ys));

    // ----- Baseline Substraction
    std::vector<IntensityType> baseline(length);
    m_BaselineSubtractor(std::begin(ys), std::end(ys), std::begin(baseline));

    // ----- Intensity Transformation
    m_Transformer(std::begin(ys), std::end(ys));

    // copy and convert
    yd.resize(length);
    std::copy(std::begin(ys), std::end(ys), std::begin(yd));
  }
}
