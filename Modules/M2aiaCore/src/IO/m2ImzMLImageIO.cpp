/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/
#define BOOST_TIMER_ENABLE_DEPRECATED
#include <Poco/SHA1Engine.h>
#include <boost/progress.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <itkMath.h>
#include <itksys/SystemTools.hxx>
#include <m2CoreCommon.h>
#include <m2ImzMLEngine.h>
#include <m2ImzMLImageIO.h>
#include <m2ImzMLParser.h>
#include <m2Timer.h>
#include <mitkIOUtil.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkLocaleSwitch.h>
#include <signal/m2PeakDetection.h>
#include <signal/m2Pooling.h>


template <class ConversionType, class ItFirst, class ItLast, class OStreamType>
void writeData(ItFirst itFirst, ItLast itLast, OStreamType &os)
{
  std::for_each(itFirst,
                itLast,
                [&os](const auto &v)
                {
                  const ConversionType c(v);
                  os.write((char *)(&c), sizeof(ConversionType));
                });
}

namespace m2
{
  ImzMLImageIO::ImzMLImageIO() : AbstractFileIO(mitk::Image::GetStaticNameOfClass(), IMZML_MIMETYPE(), "imzML Image")
  {
    AbstractFileWriter::SetRanking(10);
    AbstractFileReader::SetRanking(10);
    this->RegisterService();
  }

  mitk::IFileIO::ConfidenceLevel ImzMLImageIO::GetWriterConfidenceLevel() const
  {
    if (AbstractFileIO::GetWriterConfidenceLevel() == Unsupported)
      return Unsupported;

    auto input1 = dynamic_cast<const m2::ImzMLSpectrumImage *>(this->GetInput());
    auto input2 = dynamic_cast<const m2::SpectrumImageStack *>(this->GetInput());
    if (input1 || input2)
    {
      return Supported;
    }
    else
      return Unsupported;
  }

  std::string ImzMLImageIO::GetIBDOutputPath() const
  {
    auto outFile = GetImzMLOutputPath();
    std::string ibdOutPath = itksys::SystemTools::GetParentDirectory(outFile) + "/" +
                             itksys::SystemTools::GetFilenameWithoutExtension(outFile) + ".ibd";
    return itksys::SystemTools::ConvertToOutputPath(ibdOutPath);
  }

  std::string ImzMLImageIO::GetImzMLOutputPath() const
  {
    auto outFile = this->GetOutputLocation();
    return itksys::SystemTools::ConvertToOutputPath(outFile);
  }

  /*if (any(input->GetExportMode() & ImzMLFormatType::ContinuousCentroid))
  {
    std::vector<double> is;
    auto mzIt = std::cbegin(mzs); // mzs contains already masked mzs
    for (const auto &i : maskIndices)
    {
      double mz = *mzIt++;
      if (input->GetMassPickingTolerance() != 0.0)
      {
        auto iLower = std::lower_bound(mzs.begin(), mzs.end(), mz - input->GetMassPickingTolerance() * 10e-6);
        auto iUpper = std::upper_bound(mzs.begin(), mzs.end(), mz + input->GetMassPickingTolerance() * 10e-6);
        auto s = std::next(std::begin(ints), std::distance(mzs.begin(), iLower));
        auto e = std::next(s, std::distance(iLower, iUpper));
        double val = 0;
        switch (input->GetIonImageGrabStrategy())
        {
        case m2::IonImageGrabStrategyType::Sum:
          val = std::accumulate(s, e, double(0));
          break;
        case m2::IonImageGrabStrategyType::Mean:
          val = std::accumulate(s, e, double(0)) / double(std::distance(s, e));
          break;
        case m2::IonImageGrabStrategyType::Maximum:
          val = *std::max_element(s, e);
          break;
        case m2::IonImageGrabStrategyType::Median:
        {
          const unsigned int _N = std::distance(s, e);
          double median = 0;
          if ((_N % 2) == 0)
          {
            std::nth_element(s, s + _N * 0.5, e);
            std::nth_element(s, s + _N * 0.5 + 1, e);
            median = 0.5 * (*std::next(s, _N * 0.5) + *std::next(s, _N * 0.5 + 1));
          }
          else
          {
            std::nth_element(s, s + ((_N + 1) * 0.5), e);
            median = 0.5 * (*std::next(s, _N * 0.5) + *std::next(s, _N * 0.5 + 1));
          }
          val = median;
          break;
        }
        }
        is.push_back(val);
      }
      else
      {
        is.push_back(ints[i + input->GetConstantOffset()]);
      }
    }
    ints = std::move(is);
  }
  */

  void ImzMLImageIO::WriteContinuousProfile(m2::ImzMLSpectrumImage::SpectrumVectorType &spectra) const
  {
    const auto *input = static_cast<const m2::ImzMLSpectrumImage *>(this->GetInput());

    std::vector<float> mzs;
    std::vector<float> ints;

    std::ifstream ibd(input->GetBinaryDataPath(), std::ofstream::binary);
    std::ofstream b(GetIBDOutputPath(), std::ofstream::binary | std::ofstream::app);

    unsigned sourceId = 0;
    unsigned long long offset = 16;
    unsigned long long offsetDelta = 0;

    boost::progress_display show_progress(spectra.size() + 1);

    // write mzs
    {
      input->GetSpectrumFloat(0, mzs, ints); // get x axis

      spectra[0].mzOffset = 16;
      spectra[0].mzLength = mzs.size();

      auto start = std::begin(mzs);
      auto end = std::end(mzs);
      switch (m_DataTypeXAxis)
      {
        case m2::NumericType::Float:
          // input->GetProcessor()
          writeData<float>(start, end, b);
          offsetDelta = mzs.size() * sizeof(float);
          break;
        case m2::NumericType::Double:
          writeData<double>(start, end, b);
          offsetDelta = mzs.size() * sizeof(double);
          break;
        case m2::NumericType::None:
          mitkThrow() << "m2::NumericType of yAxisOutput not set";
      }
      offset += offsetDelta;
      ++show_progress;
    }

    MITK_INFO("ImzMLImageIO::WriteContinuousProfile") << "Write x axis done!";

    for (size_t id = 0; id < spectra.size(); ++id)
    {
      auto &s = spectra[id];
      // update mz axis info
      s.mzLength = spectra[0].mzLength;
      s.mzOffset = spectra[0].mzOffset;

      // write ints
      {
        input->GetSpectrumFloat(id, mzs, ints);

        s.intOffset = offset;
        s.intLength = ints.size();

        auto start = std::begin(ints);
        auto end = std::end(ints);

        switch (m_DataTypeYAxis)
        {
          case m2::NumericType::Float:
            writeData<float>(start, end, b);
            offsetDelta = ints.size() * sizeof(float);
            break;
          case m2::NumericType::Double:
            writeData<double>(start, end, b);
            offsetDelta = ints.size() * sizeof(double);
            break;
          case m2::NumericType::None:
            mitkThrow() << "m2::NumericType of yAxisOutput not set";
        }

        offset += offsetDelta;
      }
      b.flush();
      ++show_progress;
    }
    ++sourceId;
  }

  void ImzMLImageIO::SetIntervalVector(m2::IntervalVector::Pointer intervals)
  {
    m_Intervals = intervals;
  }

  void ImzMLImageIO::WriteContinuousCentroid(m2::ImzMLSpectrumImage::SpectrumVectorType &spectra) const
  {
    if (m_Intervals.IsNull() || m_Intervals->GetIntervals().empty())
      mitkThrow() << "No intervals provided!";

    std::vector<float> mzs, mzsMasked;
    std::vector<float> ints, intsMasked;

    // Output file stream for ibd
    const auto *input = static_cast<const m2::ImzMLSpectrumImage *>(this->GetInput());
    // auto nonConst_input = const_cast<m2::ImzMLSpectrumImage *>(input);
    mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType> macc(input->GetMaskImage());

    auto numPixels = std::accumulate(input->GetDimensions(), input->GetDimensions() + 3, 1, std::multiplies<unsigned int>());
    auto numMaskedPixel = 0 ;
    for( auto it = macc.GetData(); it <= macc.GetData() + numPixels;++it){
      if(*it > 0)
        numMaskedPixel++;
    }


    std::ofstream b(GetIBDOutputPath(), std::ofstream::binary | std::ofstream::app);

    unsigned long long offset = 16;
    unsigned long long offsetDelta = 0;
    b.seekp(offset);

    boost::progress_display show_progress(numMaskedPixel + 1);
    // write mzs
    {
      // convert to float
      const auto &xs = m_Intervals->GetXMean();
      std::copy(std::begin(xs), std::end(xs), std::back_inserter(mzsMasked));
      
      // update source spectra meta data to its actual values
      spectra[0].mzOffset = 16;
      spectra[0].mzLength = mzsMasked.size();

      switch (m_DataTypeXAxis)
      {
        case NumericType::Float:
          writeData<float>(std::begin(mzsMasked), std::end(mzsMasked), b);
          offsetDelta = spectra[0].mzLength * sizeof(float);
          break;
        case NumericType::Double:
          writeData<double>(std::begin(mzsMasked), std::end(mzsMasked), b);
          offsetDelta = spectra[0].mzLength * sizeof(double);
          break;
        case NumericType::None:
          mitkThrow() << "m2::NumericType of xAxisOutput not set";
      }
      offset += offsetDelta;
      ++show_progress;
    }

    for (size_t id = 0; id < spectra.size(); ++id)
    {
      auto &s = spectra[id];

      if (macc.GetPixelByIndex(s.index) > 0){
        // update mz axis info
        s.mzLength = spectra[0].mzLength;
        s.mzOffset = spectra[0].mzOffset;

        // write ints
        {
          input->GetSpectrumFloat(id, mzs, ints);
          intsMasked.clear();

          for (const Interval &I : m_Intervals->GetIntervals())
          {
            double val = 0;
            const auto tol = input->ApplyTolerance(I.x.mean());
            // MITK_INFO << "Tolerance: " << tol;
            const auto subRes = m2::Signal::Subrange(mzs, I.x.mean() - tol, I.x.mean() + tol);
            const auto s = std::next(std::begin(ints), subRes.first);
            const auto e = std::next(s, subRes.second);
            val = Signal::RangePooling<double>(s, e, input->GetRangePoolingStrategy());

            intsMasked.push_back(val);
          }

          s.intOffset = offset;
          s.intLength = intsMasked.size();

          switch (m_DataTypeYAxis)
          {
            case m2::NumericType::Float:
              writeData<float>(std::begin(intsMasked), std::end(intsMasked), b);
              offsetDelta = s.intLength * sizeof(float);
              break;
            case m2::NumericType::Double:
              writeData<double>(std::begin(intsMasked), std::end(intsMasked), b);
              offsetDelta = s.intLength * sizeof(double);
              break;
            case m2::NumericType::None:
              mitkThrow() << "m2::NumericType of yAxisOutput not set";
          }

          offset += offsetDelta;
        }
        b.flush();
        ++show_progress;
      }
    }
  }

  void ImzMLImageIO::WriteProcessedProfile(m2::ImzMLSpectrumImage::SpectrumVectorType & /*spectra*/) const
  {
    mitkThrow() << "Not implemented";
  }

  void ImzMLImageIO::WriteProcessedCentroid(m2::ImzMLSpectrumImage::SpectrumVectorType & /*spectra*/) const
  {
    mitkThrow() << "Not implemented";
    
  }

  void ImzMLImageIO::Write()
  {
    mitk::LocaleSwitch localeSwitch("C");
    ValidateOutputLocation();

    const auto *input = static_cast<const m2::ImzMLSpectrumImage *>(this->GetInput());
    input->SaveModeOn();

    std::string uuidString;
    {
      // Write UUID string to ibd
      std::ofstream b(GetIBDOutputPath(), std::ofstream::binary);
      boost::uuids::basic_random_generator<boost::mt19937> gen;
      boost::uuids::uuid u = gen();
      uuidString = boost::uuids::to_string(u);
      b.seekp(0);
      b.write((char *)(u.data), u.static_size());
      b.flush();
      MITK_INFO << u.data;
    }

    try
    {
      using m2::SpectrumFormat;

      // copy the whole spectra meta data container; spectra information is
      // updated based on the save mode.
      // mz and ints meta data is manipulated to write a correct imzML xml structure
      // copy of sources is discared after writing
      m2::ImzMLSpectrumImage::SpectrumVectorType spectraCopy(input->GetSpectra());


      std::string sha1;
      switch (m_SpectrumFormat)
      {
        case SpectrumFormat::ContinuousProfile:
          this->WriteContinuousProfile(spectraCopy);
          break;
        case SpectrumFormat::ProcessedCentroid:
          mitkThrow() << "ProcessedCentroid export type is not supported!";
          // this->WriteProcessedCentroid(spectraCopy);
          break;
        case SpectrumFormat::ContinuousCentroid:
          this->WriteContinuousCentroid(spectraCopy);
          break;
        case SpectrumFormat::ProcessedProfile:
          mitkThrow() << "ProcessedProfile export type is not supported!";
          break;
        default:
          break;
      }

      Poco::SHA1Engine generator;
      std::ifstream bf(GetIBDOutputPath(), std::fstream::binary);
      char *c = new char[2048];
      unsigned long long l = 0;
      while (bf)
      {
        bf.read(c, 2048);
        generator.update(c, bf.gcount());
        l += bf.gcount();
      }

      MITK_INFO << "bytes " << l;

      std::map<std::string, std::string> context;

      switch (m_SpectrumFormat)
      {
        case SpectrumFormat::None:
          mitkThrow() << "SpectrumFormatType::None type is not supported!";
          break;
        case SpectrumFormat::ContinuousCentroid:
        case SpectrumFormat::Centroid:
          context["spectrumtype"] = "centroid spectrum";
          context["mode"] = "continuous";
          break;
        case SpectrumFormat::ContinuousProfile:
          context["spectrumtype"] = "profile spectrum";
          context["mode"] = "continuous";
          break;
        case SpectrumFormat::ProcessedCentroid:
          context["spectrumtype"] = "centroid spectrum";
          context["mode"] = "processed";
          break;
        case SpectrumFormat::ProcessedProfile:
          context["spectrumtype"] = "profile spectrum";
          context["mode"] = "processed";
          break;
        default:
          break;
      }
      std::string sha1string = Poco::SHA1Engine::digestToHex(generator.digest());
      MITK_INFO << "[ibd SHA1] " << sha1string << "\n";
      MITK_INFO << "[uuid] " << uuidString << "\n";
      // context["mode"] = "[IMS:1000030] continuous";
      context["uuid"] = uuidString;
      context["sha1sum"] = sha1string;
      unsigned mzBytes = 0;
      unsigned intBytes = 0;

      switch (m_DataTypeXAxis)
      {
        case m2::NumericType::Double:
          context["mz_data_type"] = "64-bit float";
          mzBytes = 8;
          break;
        case m2::NumericType::Float:
          context["mz_data_type"] = "32-bit float";
          mzBytes = 4;
          break;
        case m2::NumericType::None:
          mitkThrow() << "m2::NumericType of xAxisOutput not set";
      }

      switch (m_DataTypeYAxis)
      {
        case m2::NumericType::Double:
          context["int_data_type"] = "64-bit float";
          intBytes = 8;
          break;
        case m2::NumericType::Float:
          context["int_data_type"] = "32-bit float";
          intBytes = 4;
          break;
        case m2::NumericType::None:
          mitkThrow() << "m2::NumericType of yAxisOutput not set";
      }

      auto nonConst_input = const_cast<m2::ImzMLSpectrumImage *>(input);
      mitk::ImagePixelReadAccessor<mitk::LabelSetImage::PixelType> macc(nonConst_input->GetMaskImage());

      context["mz_data_type_code"] = TextToCodeMap[context["mz_data_type"]];
      context["int_data_type_code"] = TextToCodeMap[context["int_data_type"]];
      context["mz_compression"] = "no compression";
      context["int_compression"] = "no compression";

      context["mode_code"] = TextToCodeMap[context["mode"]];
      context["spectrumtype_code"] = TextToCodeMap[context["spectrumtype"]];

      context["size_x"] = std::to_string(input->GetDimensions()[0]);
      context["size_y"] = std::to_string(input->GetDimensions()[1]);
      context["size_z"] = std::to_string(input->GetDimensions()[2]);

      auto xs = m2::MilliMeterToMicroMeter(input->GetGeometry()->GetSpacing()[0]);
      auto ys = m2::MilliMeterToMicroMeter(input->GetGeometry()->GetSpacing()[1]);
      auto zs = m2::MilliMeterToMicroMeter(input->GetGeometry()->GetSpacing()[2]);

      context["max dimension x"] = std::to_string(unsigned(input->GetDimensions()[0] * xs));
      context["max dimension y"] = std::to_string(unsigned(input->GetDimensions()[1] * ys));
      context["max dimension z"] = std::to_string(unsigned(input->GetDimensions()[2] * zs));

      context["pixel size x"] = std::to_string(xs);
      context["pixel size y"] = std::to_string(ys);
      context["pixel size z"] = std::to_string(zs);

      auto xo = m2::MilliMeterToMicroMeter(input->GetGeometry()->GetOrigin()[0]);
      auto yo = m2::MilliMeterToMicroMeter(input->GetGeometry()->GetOrigin()[1]);
      auto zo = m2::MilliMeterToMicroMeter(input->GetGeometry()->GetOrigin()[2]);

      context["origin x"] = std::to_string(xo);
      context["origin y"] = std::to_string(yo);
      context["origin z"] = std::to_string(zo);

      context["run_id"] = std::to_string(0);
      
      auto numPixels = std::accumulate(input->GetDimensions(), input->GetDimensions() + 3, 1, std::multiplies<unsigned int>());
      auto numMaskedPixel = 0 ;
      for( auto it = macc.GetData(); it <= macc.GetData() + numPixels;++it){
        if(*it > 0)
          numMaskedPixel++;
      }
      MITK_INFO << "numMaskedPixel: " << numMaskedPixel;
      // auto N = spectraCopy.size();
      context["num_spectra"] = std::to_string(numMaskedPixel);

      context["int_compression_code"] = TextToCodeMap[context["int_compression"]];
      context["mz_compression_code"] = TextToCodeMap[context["mz_compression"]];

      // Output file stream for imzML
      std::ofstream f(GetImzMLOutputPath(), std::ofstream::binary);

      std::string view = IMZML_TEMPLATE_START;
      f << m2::TemplateEngine::render(view, context);

      view = IMZML_SPECTRUM_TEMPLATE;
      unsigned long id = 1;

      MITK_INFO << "Write imzML data ...";
      boost::progress_display show_progress(numMaskedPixel);
      mitk::ImagePixelReadAccessor<m2::NormImagePixelType> nacc(nonConst_input->GetNormalizationImage());
      

      for (auto &s : spectraCopy)
      {
        if (macc.GetPixelByIndex(s.index) > 0){
          auto x = s.index[0] + 1; // start by 1
          auto y = s.index[1] + 1; // start by 1
          auto z = s.index[2] + 1; // start by 1

          context = {{"index", std::to_string(id++)},
                    {"x", std::to_string(x)},
                    {"y", std::to_string(y)},
                    {"z", std::to_string(z)},
                    {"mz_len", std::to_string(s.mzLength)},
                    {"mz_enc_len", std::to_string(s.mzLength * mzBytes)},
                    {"mz_offset", std::to_string(s.mzOffset)},
                    {"int_len", std::to_string(s.intLength)},
                    {"int_enc_len", std::to_string(s.intLength * intBytes)},
                    {"int_offset", std::to_string(s.intOffset)}};

          // if (nacc.GetPixelByIndex(s.index) != 1)
          // {
          //   context["tic"] = std::to_string(nacc.GetPixelByIndex(s.index));
          // }
          f << m2::TemplateEngine::render(view, context);
          f.flush();
          ++show_progress;
        }
      }

      MITK_INFO << "numMaskedPixel: " << numMaskedPixel << " id: " << id;


      f << IMZML_TEMPLATE_END;
      f.close();
    }
    catch (std::exception &e)
    {
      input->SaveModeOff();
      MITK_ERROR << "Writing ImzML Faild! " + std::string(e.what());
    }

    input->SaveModeOff();
  } // namespace mitk

  mitk::IFileIO::ConfidenceLevel ImzMLImageIO::GetReaderConfidenceLevel() const
  {
    if (AbstractFileIO::GetReaderConfidenceLevel() == Unsupported)
      return Unsupported;
    return Supported;
  }

  std::string ImzMLImageIO::RemoveExtensionFromPath(std::string path)
  {
      itksys::SystemTools::ReplaceString(path, ".imzML", "");
      itksys::SystemTools::ReplaceString(path, ".imzml", "");
      return path;  
             
      mitkThrow() << "The given file location requires an valid file ending.";
    }


  std::vector<mitk::BaseData::Pointer> ImzMLImageIO::DoRead()
  {
    std::string mzGroupId, intGroupId;
    m2::ImzMLSpectrumImage::Pointer object = m2::ImzMLSpectrumImage::New();

    auto pathWithoutExtension = RemoveExtensionFromPath(GetInputLocation());   
    if (!itksys::SystemTools::FileExists(pathWithoutExtension + ".ibd"))
      mitkThrow() << "No such file " << pathWithoutExtension;

    // m2::ImzMLSpectrumImage::ImzMLImageSource source;
    object->SetImzMLDataPath(this->GetInputLocation());
    object->SetBinaryDataPath(pathWithoutExtension + ".ibd");
    // object->GetImzMLSpectrumImageSourceList().emplace_back(source);
    object->GetSpectrumType().Format = SpectrumFormat::None;
    object->GetSpectrumType().XAxisLabel = "m/z";
    object->GetSpectrumType().YAxisLabel = "Intensity";

    {
      m2::ImzMLParser::ReadImageMetaData(object);
      m2::ImzMLParser::ReadImageSpectrumMetaData(object);
    }
    {
      object->InitializeGeometry();
    }
    {
      EvaluateSpectrumFormatType(object);
      LoadAssociatedData(object);
    }
    
    object->SetProperty("binary", mitk::BoolProperty::New(false));

    return {object.GetPointer()};
  }


 

  void ImzMLImageIO::EvaluateSpectrumFormatType(m2::SpectrumImage *object)
  {

    auto imzMLData = dynamic_cast<m2::ImzMLSpectrumImage *>(object);
    auto spectra = imzMLData->GetSpectra();
    
    bool hasTheSameMZAxisOffset = true;
    const auto mzAxisOffset = spectra[0].mzOffset;

    for(auto s : spectra){
      if(s.mzOffset != mzAxisOffset){
        hasTheSameMZAxisOffset = false;
      }

    }
    
    // check if the format type is set in the imzML fileContent element
    std::string formatType = object->GetPropertyValue<std::string>("m2aia.imzml.format_type");

    if(formatType.empty()){
      MITK_ERROR << "No format type found in the ImzML fileContent element. [IMS:1000030, IMS:1000031]";
      if(hasTheSameMZAxisOffset){
        MITK_WARN << "All spectra have the same mz offset: " << hasTheSameMZAxisOffset;
      }
      formatType = "continuous";
      MITK_WARN << "Format type set to: '" << formatType << "'";

      object->SetPropertyValue<std::string>("m2aia.imzml.format_type", formatType);
      
    }


    // check if the spectrum type is set in the imzML fileContent element
    std::string spectrumType = object->GetPropertyValue<std::string>("m2aia.imzml.spectrum_type");
    
    if(spectrumType.empty()){
      MITK_ERROR << "No spectrum type found in the ImzML fileContent element. [MS:1000127, MS:1000128]";
      MITK_WARN << "Could not load the dataset due to invalid imzML file! Make sure to properly set the format type and spectrum type in the imzML fileContent element.";
      MITK_WARN << "Fallback to centroid spectrum (MS:1000128)";
    }
      
    // MITK_WARN << "formatType: " << formatTypeProp << formatType;
    // MITK_WARN << "spectrumType: " << spectrumTypeProp << spectrumType;
    bool isProcessedSpectrum = formatType == "processed";
    bool isContinuousSpectrum = formatType == "continuous";
    bool isProfileSpectrum = spectrumType == "profile spectrum";
    bool isCentroidSpectrum = spectrumType == "centroid spectrum";

    if (!isProcessedSpectrum && !isContinuousSpectrum)
    {
      MITK_ERROR << "Set the continuos IMS:1000030) or processed (IMS:1000031) property in the ImzML "
                    "fileContent element.";
      MITK_ERROR << "Fallback to continuos (IMS:1000030)";

      object->GetSpectrumType().Format = m2::SpectrumFormat::ContinuousProfile;
    }

    if (!isCentroidSpectrum && !isProfileSpectrum)
    {
      MITK_ERROR << "Set the profile spectrum (MS:1000127) or centroid spectrum (MS:1000128) property in the ImzML "
                    "fileContent element.";
      if (isProcessedSpectrum)
      {
        MITK_ERROR << "Fallback to processed (IMS:1000031) centroid spectrum (MS:1000127)";
        object->GetSpectrumType().Format = m2::SpectrumFormat::ProcessedCentroid;
        object->SetPropertyValue("[MS:1000127] centroid spectrum", std::string());
        object->SetPropertyValue("m2aia.imzml.spectrum_type", "centroid spectrum"s);
      }

      if (isContinuousSpectrum)
      {
        MITK_ERROR << "Fallback to continuos (IMS:1000030) profile spectrum (MS:1000128)";
        object->GetSpectrumType().Format = m2::SpectrumFormat::ContinuousProfile;
        object->SetPropertyValue("[MS:1000128] profile spectrum", std::string());
        object->SetPropertyValue("m2aia.imzml.spectrum_type", "profile spectrum"s);
      }
    }

    if (object->GetProperty("[MS:1000128] profile spectrum") && object->GetProperty("[IMS:1000031] processed"))
    {
      MITK_WARN << "Processed profile spectrum is not fully supported! Check the ImzML file.";
    }


    if (isContinuousSpectrum && isProfileSpectrum)
      object->GetSpectrumType().Format = m2::SpectrumFormat::ContinuousProfile;
    else if (isProcessedSpectrum && isProfileSpectrum)
      object->GetSpectrumType().Format = m2::SpectrumFormat::ProcessedProfile;
    else if (isContinuousSpectrum && isCentroidSpectrum)
      object->GetSpectrumType().Format = m2::SpectrumFormat::ContinuousCentroid;
    else if (isProcessedSpectrum && isCentroidSpectrum)
      object->GetSpectrumType().Format = m2::SpectrumFormat::ProcessedCentroid;

  }

  void ImzMLImageIO::LoadAssociatedData(m2::ImzMLSpectrumImage *object)
  {
    auto pathWithoutExtension = RemoveExtensionFromPath(GetInputLocation());   

    auto maskPath = pathWithoutExtension + ".mask.nrrd";
    MITK_INFO << maskPath;
    if (itksys::SystemTools::FileExists(maskPath))
    {
      auto data = mitk::IOUtil::Load(maskPath).at(0);
      object->SetMaskImage(dynamic_cast<mitk::Image *>(data.GetPointer()));
    }

    
    auto shiftImagePath = pathWithoutExtension + ".index_shift.nrrd";
    MITK_INFO << shiftImagePath;
    if (itksys::SystemTools::FileExists(shiftImagePath))
    {
      auto data = mitk::IOUtil::Load(shiftImagePath).at(0);
      object->SetShiftImage(dynamic_cast<mitk::Image *>(data.GetPointer()));
    }

    auto normPath = pathWithoutExtension + ".norm.nrrd";
    if (itksys::SystemTools::FileExists(normPath))
    {
      auto data = mitk::IOUtil::Load(normPath).at(0);

      if(mitk::Equal(*object->GetGeometry(), *data->GetGeometry())){
        auto type = m2::NormalizationStrategyType::External;
        object->SetNormalizationImage(dynamic_cast<mitk::Image *>(data.GetPointer()), type);
        object->SetNormalizationImageStatus(type, true);
      }else{
        MITK_ERROR << "External normalization image geometry is not equal to the loaded data";
      }

    }


    auto pointsPath = pathWithoutExtension + ".mps";
    if (itksys::SystemTools::FileExists(pointsPath))
    {
      auto data = mitk::IOUtil::Load(pointsPath).at(0);
      object->SetPoints(dynamic_cast<mitk::PointSet *>(data.GetPointer()));
    }
  }

  ImzMLImageIO *ImzMLImageIO::IOClone() const
  {
    return new ImzMLImageIO(*this);
  }
} // namespace m2


#undef BOOST_TIMER_ENABLE_DEPRECATED