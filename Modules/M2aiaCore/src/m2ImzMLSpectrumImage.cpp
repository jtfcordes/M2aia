/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include <m2ImzMLSpectrumImage.h>
#include <m2ImzMLSpectrumImageSource.hpp>

#include <m2Process.hpp>
#include <m2Timer.h>
#include <mitkImageAccessByItk.h>
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

double m2::ImzMLSpectrumImage::GetXMin() const{
  return GetXAxis().front();
}

double m2::ImzMLSpectrumImage::GetXMax() const{  
  return GetXAxis().back();
}


void m2::ImzMLSpectrumImage::GetImage(double mz, double tol, const mitk::Image *mask, mitk::Image *img) const
{
  try{
    while(!m_IsIonImageGenerated);
    m_IsIonImageGenerated = false;
    m_SpectrumImageSource->GetImagePrivate(mz, tol, mask, img);
    m_IsIonImageGenerated = true;
    m_CurrentX = mz;
  }catch(std::exception & e){
    MITK_ERROR << "Ion image could not be generated! Queried range is [" << mz-tol << ", " <<mz+tol << "]\n" << e.what();
  }
}

void m2::ImzMLSpectrumImage::InitializeProcessor()
{
  m_MzGroupID = GetPropertyValue<std::string>("m2aia.imzml.mzGroupID");
  m_IntensityGroupID = GetPropertyValue<std::string>("m2aia.imzml.intensityGroupID");
  
  auto intensitiesDataTypeString = GetPropertyValue<std::string>("m2aia.imzml." + m_IntensityGroupID + ".value_type");
  auto mzValueTypeString = GetPropertyValue<std::string>("m2aia.imzml." + m_MzGroupID + ".value_type");

  if (mzValueTypeString.compare("32-bit float") == 0)
  {
    m_SpectrumType.XAxisType = m2::NumericType::Float;
    if (intensitiesDataTypeString.compare("32-bit float") == 0)
    {
      m_SpectrumType.YAxisType = m2::NumericType::Float;
      this->m_SpectrumImageSource.reset((ISpectrumImageSource *)new ImzMLSpectrumImageSource<float, float>(this));
    }
    else if (intensitiesDataTypeString.compare("64-bit float") == 0)
    {
      m_SpectrumType.YAxisType = m2::NumericType::Double;
      this->m_SpectrumImageSource.reset((ISpectrumImageSource *)new ImzMLSpectrumImageSource<float, double>(this));
    }
    else if (intensitiesDataTypeString.compare("32-bit integer") == 0)
    {
      MITK_ERROR(GetStaticNameOfClass()) << "Using 32-bit integer. Not implemented!";
      // this->m_SpectrumImageSource.reset((ProcessorBase *)new Processor<float, long int>(this));
      // SetYInputType(m2::NumericType::Double);
    }
    else if (intensitiesDataTypeString.compare("64-bit integer") == 0)
    {
      MITK_ERROR(GetStaticNameOfClass()) << "Using 64-bit integer. Not implemented!";
      // this->m_SpectrumImageSource.reset((ProcessorBase *)new Processor<float, long long int>(this));
      // SetYInputType(m2::NumericType::Double);
    }
  }
  else if (mzValueTypeString.compare("64-bit float") == 0)
  {
    m_SpectrumType.XAxisType = m2::NumericType::Double;
    if (intensitiesDataTypeString.compare("32-bit float") == 0)
    {
      m_SpectrumType.YAxisType = m2::NumericType::Float;
      this->m_SpectrumImageSource.reset((ISpectrumImageSource *)new ImzMLSpectrumImageSource<double, float>(this));
    }
    else if (intensitiesDataTypeString.compare("64-bit float") == 0)
    {
      m_SpectrumType.YAxisType = m2::NumericType::Double;
      this->m_SpectrumImageSource.reset((ISpectrumImageSource *)new ImzMLSpectrumImageSource<double, double>(this));
    }
    else if (intensitiesDataTypeString.compare("32-bit integer") == 0)
    {
      MITK_ERROR(GetStaticNameOfClass()) << "Using 32-bit integer. Not implemented!";
      // this->m_SpectrumImageSource.reset((ProcessorBase *)new Processor<double, long int>(this));
      // SetYInputType(m2::NumericType::Double);
    }
    else if (intensitiesDataTypeString.compare("64-bit integer") == 0)
    {
      // this->m_SpectrumImageSource.reset((ProcessorBase *)new Processor<double, long long int>(this));
      MITK_ERROR(GetStaticNameOfClass()) << "Using 64-bit integer. Not implemented!";
      // SetYInputType(m2::NumericType::Double);
    }
  }
}

void m2::ImzMLSpectrumImage::InitializeGeometry()
{
  this->SetImageGeometryInitialized(false);

  if (!m_SpectrumImageSource)
    this->InitializeProcessor();

  GetSpectraArtifacts().clear();
 this->m_SpectrumImageSource->InitializeGeometry();
  
  this->SetImageGeometryInitialized(true);
}

void m2::ImzMLSpectrumImage::InitializeImageAccess()
{

  // Reset Normalization strategy type if set to external but no external image was found
  if(GetNormalizationStrategy() == m2::NormalizationStrategyType::External && GetNormalizationImage(m2::NormalizationStrategyType::External).IsNull()){
    SetNormalizationStrategy(m2::NormalizationStrategyType::None);
    MITK_ERROR << "'External' Normalization strategy chosen but no External image exist.\n"
    "To use external normalization provide an image in the NRRD file format.\n"
    "This image requires the same name as the imzML file, with file ending *.norm.nrrd, pixel type double, identical pixel spacing, and image dimensions.";
  }

  this->SetImageAccessInitialized(false); 

  

  this->m_SpectrumImageSource->InitializeImageAccess();

  auto sx = this->GetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x");
  auto sy = this->GetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y");
  auto sz = this->GetPropertyValue<unsigned>("max count of pixels z");

  auto px = this->GetPropertyValue<double>("[IMS:1000046] pixel size x");
  auto py = this->GetPropertyValue<double>("[IMS:1000047] pixel size y");
  auto pz = this->GetPropertyValue<double>("pixel size z");



    MITK_INFO << "[imzML]: " + this->GetImzMLDataPath() +
              "\n\t[pixel size (mm)]: " + std::to_string(px)  + "x" + std::to_string(py) +"x" + std::to_string(pz) +
              "\n\t[image dimension]: " + std::to_string(sx) + "x" +std::to_string(sy) + "x" + std::to_string(sz)  +
              "\n\t[num spectra]: " + std::to_string(this->GetNumberOfValidPixels()) +
              "\n\t[spec. type]: " + to_string(this->GetSpectrumType().Format) +
              "\n\t[mass range]: " + std::to_string(GetXMin()) + " to " + std::to_string(GetXMax()) + " with #" + std::to_string(GetXAxis().size()) + " measurements" +
              ""+(     to_underlying(GetBaselineCorrectionStrategy()) ? ("\n\t[baseline correction]: " + m2::BaselineCorrectionTypeNames.at(to_underlying(GetBaselineCorrectionStrategy())) + "(" + std::to_string(GetBaseLineCorrectionHalfWindowSize()) + ")"):"") +
              ""+(              to_underlying(GetSmoothingStrategy()) ? ("\n\t[smoothing]: "           + m2::SmoothingTypeNames.at(to_underlying(GetSmoothingStrategy())) + "(" + std::to_string(GetSmoothingHalfWindowSize()) + ")"):"") +
              ""+(          to_underlying(GetNormalizationStrategy()) ? ("\n\t[normalization]: "       + m2::NormalizationStrategyTypeNames.at(to_underlying(GetNormalizationStrategy()))):"") +
              ""+(to_underlying(GetIntensityTransformationStrategy()) ? ("\n\t[processing]: "          + m2::IntensityTransformationTypeNames.at(to_underlying(GetIntensityTransformationStrategy()))):"");
              
  this->SetImageAccessInitialized(true); 
}


void m2::ImzMLSpectrumImage::InitializeNormalizationImage(m2::NormalizationStrategyType type){
  if(GetNormalizationImageStatus(type) == false){
   m_SpectrumImageSource->InitializeNormalizationImage(type);
   SetNormalizationImageStatus(type, true);
  }
}

// m2::ImzMLSpectrumImage::Pointer m2::ImzMLSpectrumImage::Combine(const m2::ImzMLSpectrumImage *A,
//                                                                 const m2::ImzMLSpectrumImage *B,
//                                                                 const char stackAxis)
// {
//   auto C = m2::ImzMLSpectrumImage::New();
//   C->SetSpectrumType(A->GetSpectrumType());
//   {
//     auto A_mzGroup = A->GetPropertyValue<std::string>("mzGroupName");
//     auto A_mzValueType = A->GetPropertyValue<std::string>(A_mzGroup + " value type");
//     auto B_mzGroup = B->GetPropertyValue<std::string>("mzGroupName");
//     auto B_mzValueType = B->GetPropertyValue<std::string>(B_mzGroup + " value type");
//     if (A_mzValueType.compare(B_mzValueType) != 0)
//       mitkThrow() << "m/z value types must the same!";

//     auto typeSizeInBytes = A->GetPropertyValue<unsigned>(A_mzGroup + " value type (bytes)");
//     auto refGroupID = A->GetPropertyValue<std::string>(A_mzGroup);

//     C->SetPropertyValue<std::string>(A_mzGroup, refGroupID);
//     C->SetPropertyValue<unsigned>(A_mzGroup + " value type (bytes)", typeSizeInBytes);
//     C->SetPropertyValue<std::string>(A_mzGroup + " value type", A_mzValueType);
//     C->SetPropertyValue<std::string>("mzGroupName", A_mzGroup);
//   }

//   {
//     auto A_intensityGroup = A->GetPropertyValue<std::string>("intensityGroupName");
//     auto A_intensityValueType = A->GetPropertyValue<std::string>(A_intensityGroup + " value type");
//     auto B_intensityGroup = B->GetPropertyValue<std::string>("intensityGroupName");
//     auto B_intensityValueType = B->GetPropertyValue<std::string>(B_intensityGroup + " value type");
//     if (B_intensityValueType.compare(A_intensityValueType) != 0)
//       mitkThrow() << "Intensity value types must the same!";

//     auto typeSizeInBytes = A->GetPropertyValue<unsigned>(A_intensityGroup + " value type (bytes)");
//     auto refGroupID = A->GetPropertyValue<std::string>(A_intensityGroup);

//     C->SetPropertyValue<std::string>(A_intensityGroup, refGroupID);
//     C->SetPropertyValue<unsigned>(A_intensityGroup + " value type (bytes)", typeSizeInBytes);
//     C->SetPropertyValue<std::string>(A_intensityGroup + " value type", A_intensityValueType);
//     C->SetPropertyValue<std::string>("intensityGroupName", A_intensityGroup);
//   }

//   {
//     auto A_mass = A->GetXAxis();
//     auto B_mass = B->GetXAxis();
//     auto isMassEqual = std::equal(std::begin(A_mass), std::end(A_mass), std::begin(B_mass));
//     if (!isMassEqual)
//       mitkThrow() << "Mass axis must be the same!";
//   }

//   {
//     auto A_x = A->GetPropertyValue<double>("[IMS:1000046] pixel size x");
//     auto A_y = A->GetPropertyValue<double>("[IMS:1000047] pixel size y");
//     auto A_z = A->GetPropertyValue<double>("pixel size z");

//     auto B_x = B->GetPropertyValue<double>("[IMS:1000046] pixel size x");
//     auto B_y = B->GetPropertyValue<double>("[IMS:1000047] pixel size y");
//     auto B_z = B->GetPropertyValue<double>("pixel size z");

//     if ((std::abs(A_x - B_x) > 10e-5 && std::abs(A_y - B_y) > 10e-5) || std::abs(A_z - B_z) > 10e-5)
//       mitkThrow() << "Pixel size must be the same!";

//     C->SetPropertyValue<double>("[IMS:1000046] pixel size x", A_x);
//     C->SetPropertyValue<double>("[IMS:1000047] pixel size y", A_y);
//     C->SetPropertyValue<double>("pixel size z", A_z);
//   }

//   unsigned A_x, A_y, A_z;
//   unsigned B_x, B_y, B_z;
//   A_x = A->GetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x");
//   B_x = B->GetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x");
//   A_y = A->GetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y");
//   B_y = B->GetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y");
//   A_z = A->GetPropertyValue<unsigned>("max count of pixels z");
//   B_z = B->GetPropertyValue<unsigned>("max count of pixels z");

//   switch (stackAxis)
//   {
//     case 'x':
//     {
//       // Images are stacked along the x axis
//       C->SetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x", A_x + B_x);
//       C->SetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y", std::max(A_y, B_y));
//       C->SetPropertyValue<unsigned>("max count of pixels z", std::max(A_z, B_z));
//     }
//     break;
//     case 'y':
//     {
//       C->SetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x", std::max(A_x, B_x));
//       C->SetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y", A_y + B_y);
//       C->SetPropertyValue<unsigned>("max count of pixels z", std::max(A_z, B_z));
//     }
//     break;
//     case 'z':
//     {
//       C->SetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x", std::max(A_x, B_x));
//       C->SetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y", std::max(A_y, B_y));
//       C->SetPropertyValue<unsigned>("max count of pixels z", A_z + B_z);
//     }
//     break;
//   }

//   {
//     C->SetPropertyValue<double>("[IMS:1000053] absolute position offset x", 0.0);
//     C->SetPropertyValue<double>("[IMS:1000054] absolute position offset y", 0.0);
//     C->SetPropertyValue<double>("absolute position offset z", 0.0);
//   }

//   C->SetPropertyValue<bool>("continuous", true);

//   auto &C_sources = C->GetImzMLSpectrumImageSourceList();

//   const auto &A_sources = A->GetImzMLSpectrumImageSourceList();
//   const auto &B_sources = B->GetImzMLSpectrumImageSourceList();

//   C_sources = A_sources;
//   // MITK_INFO(m2::ImzMLSpectrumImage::GetStaticNameOfClass()) << "Merge sources";

//   std::transform(std::begin(B_sources), std::end(B_sources), std::back_inserter(C_sources), [&](ImzMLImageSource s) {
//     if (stackAxis == 'x')
//       s.m_Offset[0] += A_x; // shift along x-axis
//     if (stackAxis == 'y')
//       s.m_Offset[1] += A_y; // shift along y-axis
//     if (stackAxis == 'z')
//       s.m_Offset[2] += A_z; // shift along z-axis
//     return s;
//   });

//   C->InitializeGeometry();
//   std::set<m2::Interval> set;
//   set.insert(std::begin(A->GetIntervals()), std::end(A->GetIntervals()));
//   set.insert(std::begin(B->GetIntervals()), std::end(B->GetIntervals()));

//   C->GetIntervals().insert(std::end(C->GetIntervals()), std::begin(set), std::end(set));
//   return C;
// }


void m2::ImzMLSpectrumImage::GetSpectrumFloat(unsigned int id,
                                         std::vector<float> &xs,
                                         std::vector<float> &ys) const
{
  m_SpectrumImageSource->GetXValues(id, xs);
  m_SpectrumImageSource->GetYValues(id, ys);
}

void m2::ImzMLSpectrumImage::GetSpectrum(unsigned int id,
                                         std::vector<double> &xs,
                                         std::vector<double> &ys) const
{
  m_SpectrumImageSource->GetXValues(id, xs);
  m_SpectrumImageSource->GetYValues(id, ys);
}

void m2::ImzMLSpectrumImage::GetIntensitiesFloat(unsigned int id, std::vector<float> &ys) const
{
  m_SpectrumImageSource->GetYValues(id, ys);
}

void m2::ImzMLSpectrumImage::GetIntensities(unsigned int id, std::vector<double> &ys) const
{
  m_SpectrumImageSource->GetYValues(id, ys);
}

// void m2::ImzMLSpectrumImage::GetIntensities(unsigned int id,
//                                             std::vector<m2::Interval> &I,
//                                             std::vector<float> &pys,
//                                             unsigned int source) const
// {
//   using namespace std;
//   vector<float> xs, ys;
//   this->GetSpectrum(id,xs,ys, source);
//   auto it = begin(xs);
//   auto xsBegin = begin(xs);
//   auto xsEnd = end(xs);
//   pys.clear();
//   auto inserter =std::back_inserter(pys);
//   for (unsigned int p = 0; p < I.size(); ++p)
//   {
//     it = lower_bound(it,xsEnd,I[p].x.min());
//     auto lower = distance(xsBegin, it);
//     it = upper_bound(it,xsEnd,I[p].x.max());
//     auto upper = distance(xsBegin, it);
    
//     inserter = Signal::RangePooling<float>(begin(ys)+lower, begin(ys)+upper, GetRangePoolingStrategy());

//   }
// }


// template <class MassAxisType, class IntensityType>
// void m2::ImzMLSpectrumImage::ImzMLImageProcessor<MassAxisType, IntensityType>::GetSpectrumPrivate(
//   unsigned int id, std::vector<float> &xd, std::vector<float> &ydId)
// {
//   mitk::ImagePixelReadAccessor<m2::NormImagePixelType, 3> normAcc(p->GetNormalizationImage());

//   const auto &source = p->m_SourcesList[sourceId];
//   std::ifstream f;
//   f.open(source.m_BinaryDataPath, std::ios::binary);

//   const auto &spectrumType = p->GetSpectrumType();
//   const auto &spectrum = source.m_Spectra[id];
//   const auto &length = spectrum.intLength;
//   const auto &ysOffset = spectrum.intOffset;

//   std::vector<IntensityType> ys;
//   ys.resize(length);
//   binaryDataToVector(f, ysOffset, length, ys.data());

//   // No preprocessing is applied for single spectrum export
//   yd.resize(length);
//   xd.resize(length);

//   std::copy(std::begin(ys), std::end(ys), std::begin(yd));

//   // continuous data -> mass axis already exist and can be used directly
//   if (any(spectrumType.Format & m2::SpectrumFormatType::ContinuousCentroid |
//   m2::SpectrumFormatType::ContinuousProfile))
//   {
//     auto &xAxis = p->GetXAxis();
//     std::copy(std::begin(xAxis), std::end(xAxis), std::begin(xd));
//   }
//   else
//   {
//     // processed data -> mass axis for individual pixels
//     const auto &xsOffset = spectrum.mzOffset;
//     if (spectrumType.XAxisType == m2::NumericType::Float)
//     {
//       // float must be converted to output type "double"
//       std::vector<MassAxisType> xs;
//       xs.resize(length);
//       binaryDataToVector(f, xsOffset, length, xs.data());
//       std::copy(std::begin(xs), std::end(xs), std::begin(xd));
//     }
//     else if (spectrumType.XAxisType == m2::NumericType::Double)
//     {
//       // copy directly
//       xd.resize(length);
//       binaryDataToVector(f, xsOffset, length, xd.data());
//     }
//   }
// }

// template <class MassAxisType, class IntensityType>
// void m2::ImzMLSpectrumImage::ImzMLImageProcessor<MassAxisType, IntensityType>::GrabMassPrivate(
//   unsigned long int index, std::vector<double> &mzsIndex) const
// {
//   const auto &source = p->GetImzMLSpectrumImageSourceList()[sourceIndex];
//   const auto &spectrum = source.m_Spectra[index];
//   const auto &mzo = spectrum.mzOffset;
//   const auto &mzl = spectrum.mzLength;

//   std::ifstream f;
//   f.open(source.m_BinaryDataPath, std::ios::binary);

//   std::vector<MassAxisType> mzs_get;
//   binaryDataToVector(f, mzo, mzl, mzs_get);
//   mzs.clear();

//   std::copy(std::begin(mzs_get), std::end(mzs_get), std::back_inserter(mzs));
// }

m2::ImzMLSpectrumImage::~ImzMLSpectrumImage()
{
  // MITK_INFO(m2::ImzMLSpectrumImage::GetStaticNameOfClass()) << " destroyed!";
}

m2::ImzMLSpectrumImage::ImzMLSpectrumImage() : m2::SpectrumImage()
{
  // MITK_INFO(m2::ImzMLSpectrumImage::GetStaticNameOfClass()) << " created!";
  m_SpectrumType.XAxisLabel = "m/z";
  m_ExportSpectrumType.XAxisLabel = "m/z";
}
