/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/
#include <future>
#include <iterator>
#include <m2ImzMLParser.h>
#include <m2Timer.h>
#include <math.h>
#include <mitkCoreServices.h>
#include <mitkIPreferences.h>
#include <mitkIPreferencesService.h>
#include <numeric>
#include <unordered_map>

auto m2::ImzMLParser::findLine(std::ifstream &f, std::string name, std::string start_tag, bool eol)
  -> unsigned long long
{
  auto line = f.tellg();
  std::string str;
  unsigned int i = 0;

  std::string bb;
  bb += start_tag[0];
  while (std::getline(f, str))
  {
    if (eol)
      line = f.tellg();
    auto pos2 = str.find(start_tag[0]);
    if (pos2 != std::string::npos)
    {
      bool match = true;
      for (i = 1; i < start_tag.size() && (pos2 + i) < str.size(); ++i)
      {
        bb += str[pos2 + i];
        if (str[pos2 + i] != start_tag[i])
        {
          match = false;
          break;
        }
      }

      if (match)
      {
        str = str.substr(pos2 + start_tag.size(), name.size());
        if (str.compare(name) == 0)
        {
          return (unsigned long long)(line);
        }
      }
    }

    line = f.tellg(); // point to beginning of the new line
  }
  return -1;
}

void m2::ImzMLParser::ReadImageMetaData(m2::ImzMLSpectrumImage::Pointer data)
{
  std::ifstream f;
  std::vector<std::string> stack, context_stack;
  std::string line, context, tag, name, value, accession;
  tag.reserve(60);
  name.reserve(60);
  value.reserve(60);
  context.reserve(60);
  accession.reserve(60);

  struct m2ImzMLParserAccessionData
  {
    FunctionType f;
    bool isRequired = false;
    bool wasFound = false;
  };

  std::unordered_map<std::string, FunctionType> accession_map;
  std::unordered_map<std::string, FunctionType> context_map;

  f.open(data->GetImzMLDataPath(), std::ios_base::binary);

  std::map<std::string, unsigned> precisionDict = {{"32-bit float", sizeof(float)},
                                                   {"64-bit float", sizeof(double)},
                                                   {"32-bit integer", sizeof(int32_t)},
                                                   {"64-bit integer", sizeof(int64_t)}};

  const auto ContextValueToStringProperty = [&](const std::string &line)
  {
    attributeValue(line, "name", name);
    attributeValue(line, "value", value);
    attributeValue(line, "accession", accession);

    if(value.empty())
      value = "true";

    if (!context.empty()){
      data->SetPropertyValue("["+accession+"] " + context + "." + name, value);
    }
    else{
      data->SetPropertyValue("["+accession+"] " + name, value);
    }
  };

  // const auto ValueToUnsignedIntProperty = [&](const std::string &line, auto converter, const std::string &
  // default_name = {})
  // {
  //   attributeValue(line, "value", value);
  //   if(default_name.empty()){
  //     attributeValue(line, "name", name);
  //     data->SetPropertyValue<unsigned>(name, std::stoul(value));
  //   }else{
  //     data->SetPropertyValue<unsigned>(default_name, std::stoul(value));
  //   }

  // };

  const auto ValueToProperty = [&](const std::string &line, auto converter, const std::string &default_name = {})
  {
    try
    {
      attributeValue(line, "accession", accession);
      attributeValue(line, "value", value);
      auto converted_value = converter(value);
      using PropertyType = decltype(converted_value);
      if (default_name.empty())
      {
        attributeValue(line, "name", name);
        data->SetPropertyValue<PropertyType>("["+accession+"] " + name, converted_value);
      }
      else
      {
        data->SetPropertyValue<PropertyType>("["+accession+"] " + default_name, converted_value);
      }
    }
    catch (std::exception &e)
    {
      MITK_INFO << "Check this line of your imzML file:\n" << line << "\nPossible incorrect imzML detected!";
    }
  };

  {
    // ------ Process a referenceableParameter Group
    context_map["referenceableParamGroup"] = [&](std::string line)
    {
      std::string id, dataType, spectrumType, spectrumTypeAccession;
      id.reserve(60);
      attributeValue(line, "id", id);

      const auto npos = std::string::npos;
      while (!f.eof())
      {
        std::getline(f, line); // read the next line

        if (line.find("</") != std::string::npos)
          break;

        attributeValue(line, "accession", accession);
        attributeValue(line, "name", name);
        attributeValue(line, "value", value);

        if(value.empty()) value = "true";

        data->SetPropertyValue<std::string>("["+accession+"] "+id+"."+name, value);
        
        if(line.find("MS:1000520") != npos ||
           line.find("MS:1000521") != npos ||
           line.find("MS:1000522") != npos ||
           line.find("MS:1000523") != npos){
          dataType = name;
           }
        
        if(line.find("MS:1000127") != npos || 
           line.find("MS:1000128") != npos){
          spectrumType = name;
          spectrumTypeAccession = accession;
        }

        if (line.find("MS:1000514") != npos){ 
          // m/z array https://github.com/m2aia/psi-ms-CV/blob/master/psi-ms.obo#L3695
          data->SetPropertyValue<std::string>("m2aia.imzml.mzGroupID", id);
        }

        if (line.find("MS:1000515") != npos){
            // intensity array https://github.com/m2aia/psi-ms-CV/blob/master/psi-ms.obo#L3704
          data->SetPropertyValue<std::string>("m2aia.imzml.intensityGroupID", id);
        }
      }
      
      if(!dataType.empty()){
        data->SetPropertyValue<unsigned>("m2aia.imzml." + id + ".value_type_in_bytes", precisionDict[dataType]);
        data->SetPropertyValue<std::string>("m2aia.imzml." + id + ".value_type", dataType);
      }
      if(!spectrumType.empty()){
        data->SetPropertyValue<std::string>("m2aia.imzml.spectrum_type", spectrumType);
        data->SetPropertyValue<std::string>("["+spectrumTypeAccession+"] "+name, std::string("true"));
      }
      
      return "";
    };

    context_map["software"] = [&](const std::string &line)
    {
      attributeValue(line, "id", name);
      attributeValue(line, "version", value);
      context = name + " " + value;
    };
    context_map["scanSettings"] = [&](const std::string &line) { attributeValue(line, "id", context); };

    using namespace std::string_literals;
    using namespace std::placeholders;

    auto stoui = [](const std::string &s) -> unsigned int { return std::stoul(s); };
    auto stod = [](const std::string &s) -> double { return std::stod(s); };
    auto stol = [](const std::string &s) -> long { return std::stol(s); };

    // Attention: we need to provide default values for the name to avoid frequently observed typos(pixel instead of
    // pixels). https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L136
    accession_map["IMS:1000042"] = [&](auto line) { ValueToProperty(line, stoui, "max count of pixels x"s); };

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L144
    accession_map["IMS:1000043"] = [&](auto line) { ValueToProperty(line, stoui, "max count of pixels y"s); };

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L153
    accession_map["IMS:1000044"] = [&](auto line) { ValueToProperty(line, stod); }; // "max dimension x"s

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L160
    accession_map["IMS:1000045"] = [&](auto line) { ValueToProperty(line, stod); }; // "max dimension y"s

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L113
    accession_map["IMS:1000030"] = [&](auto line) { // "continuous"
      ContextValueToStringProperty(line); 
    data->SetPropertyValue<std::string>("m2aia.imzml.format_type", name);
    }; 

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L119
    accession_map["IMS:1000031"] = [&](auto line) { // "processed"
      ContextValueToStringProperty(line); 
    data->SetPropertyValue<std::string>("m2aia.imzml.format_type", name);
    }; 


        // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L113
    accession_map["MS:1000127"] = [&](auto line) { // "continuous"
      ContextValueToStringProperty(line);
    data->SetPropertyValue<std::string>("m2aia.imzml.spectrum_type", name);
    }; 

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L119
    accession_map["MS:1000128"] = [&](auto line) { // "processed"
      ContextValueToStringProperty(line); 
    data->SetPropertyValue<std::string>("m2aia.imzml.spectrum_type", name);
    }; 


    // Attention: we should not specify a default value, for the name here.
    // Describes the length of a pixel in the x dimension. If no pixel size y (IMS:1000047) is explicitly specified,
    // then this also describes the length of a pixel in the y dimension."
    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L169
    accession_map["IMS:1000046"] = [&](auto line) { ValueToProperty(line, stod); }; // pixel size (x)

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L177
    accession_map["IMS:1000047"] = [&](auto line) { ValueToProperty(line, stod); }; // pixel size y

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L220
    accession_map["IMS:1000053"] = [&](auto line)
    {
      ValueToProperty(line, std::bind(m2::MicroMeterToMilliMeter, std::bind(stol, _1)));
    }; //, "[IMS:1000053] absolute position offset x"s

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L228
    accession_map["IMS:1000054"] = [&](auto line)
    {
      ValueToProperty(line, std::bind(m2::MicroMeterToMilliMeter, std::bind(stol, _1)));
    }; //, "[IMS:1000054] absolute position offset y"s

    // // Not supported by IMS obo
    // accession_map["IMS:X3"] = accession_map["M2:0000003"] = [&](const std::string &line) { // absolute position
    // offset y
    //   attributeValue(line, "value", value);
    //   data->SetPropertyValue<double>("absolute position offset z", m2::MicroMeterToMilliMeter(std::stol(value)));
    // };

    // origin size z
    // scanSettings

    context_map["source"] = [&](const std::string &) { context = "source"; };
    context_map["analyzer"] = [&](const std::string &) { context = "analyzer"; };
    context_map["detector"] = [&](const std::string &) { context = "detector"; };

    context_map["instrumentConfiguration"] = [&](const std::string &line) { attributeValue(line, "id", context); };
    context_map["dataProcessing"] = [&](const std::string &line) { attributeValue(line, "id", context); };

    std::string dataProcessingId;
    std::string processingMethodOrder;
    std::string processingMethodSoftwareRef;
    context_map["dataProcessing"] = [&](const std::string &line){
      attributeValue(line, "id", dataProcessingId);
      // context = context + ".dataProcessing (" + dataProcessingId + ")";
      // MITK_INFO << dataProcessingId;
    };

    context_map["processingMethod"] = [&](const std::string &line)
    {
      attributeValue(line, "order", processingMethodOrder);
      attributeValue(line, "softwareRef", processingMethodSoftwareRef);
      // MITK_INFO << processingMethodOrder << " " << processingMethodSoftwareRef;
    };

    context_map["userParam"] = [&](const std::string &line)
    {
      std::string description;
      std::string method_value;
      attributeValue(line, "name", description);
      attributeValue(line, "value", method_value);
      data->SetPropertyValue<std::string>("m2aia.imzml."+dataProcessingId+"." +processingMethodSoftwareRef+"."+processingMethodOrder+".description", description);
      data->SetPropertyValue<std::string>("m2aia.imzml."+dataProcessingId+"." +processingMethodSoftwareRef+"."+processingMethodOrder+".value", method_value);

      // MITK_INFO << description << " " << method_value;
    };

    // default values
    data->SetPropertyValue<unsigned>("max count of pixels z", 1);
    
    data->SetPropertyValue<double>("[IMS:1000046] pixel size x", -1);
    data->SetPropertyValue<double>("[IMS:1000047] pixel size y", -1);
    data->SetPropertyValue<double>("pixel size z", -1);

    data->SetPropertyValue<double>("[IMS:1000053] absolute position offset x", 0.0);
    data->SetPropertyValue<double>("[IMS:1000054] absolute position offset y", 0.0);
    data->SetPropertyValue<double>("absolute position offset z", 0.0);

    // -------- PROCESS FILE --------
    {
      while (!f.eof())
      {
        std::getline(f, line); // read the next line

        // assuming spectra meta data after ImzML meta data.
        if (line.find("<run") != std::string::npos)
          break; // indicates begin of spectra meta data

        // Check for end-tag. If true, pop the last context
        // element from the stack and clear the context string.
        if (line.find("</") != std::string::npos) // close context
        {
          context_stack.pop_back();
          context.clear();
        }

        // Check for empty-element tag. Can be cvParam or userParam (e.g. used by ScilsLab).
        else if (line.rfind("/>") != std::string::npos) // element
        {
          attributeValue(line, "accession", accession);
          if (!accession.empty())
          {
            // if (!context.empty())
            // {
            //   // Call context specific accession.
            //   if (evaluateAccession(line, accession + "[" + context + "]", accession_map)){
            //     MITK_INFO << "Context";
            //     continue; // if evaluation succeeded continue loop and go to next line
            //   }
            // }
            // Fallback if context specific accession is not found.
            auto status = evaluateAccession(line, accession, accession_map);
            accession.clear();
            if (!status) // Default: name + value is added to IMS data property list
              ContextValueToStringProperty(line);
          }else{
            GetElementName(line, tag);
            context_stack.push_back(tag);
            EvaluateContext(line, tag, context_map);
          }
  
        }

        else // open context
        {    // it is probably a start-tag. Could be a comment ...
          // .. but in a clean file no comments should be available.

          GetElementName(line, tag);
          context_stack.push_back(tag);
          EvaluateContext(line, tag, context_map);
        }
      }
    }

    // -------- POST PROCESS META DATA --------
    // convert micrometer in millimeter for compatibility with MITK

    // if(data->GetProperty("max dimension x") &&
    //    data->GetProperty("max dimension y"))
    // {
    //   auto dx = data->GetPropertyValue<double>("max dimension x");
    //   auto dy = data->GetPropertyValue<double>("max dimension y");
    //   dx = m2::MicroMeterToMilliMeter(dx);
    //   dy = m2::MicroMeterToMilliMeter(dy);
    //   data->SetPropertyValue<double>("max dimension y", dx);
    //   data->SetPropertyValue<double>("max dimension y", dy);

    // }

    if (data->GetProperty("[IMS:1000046] pixel size"))
    {
      auto v = data->GetPropertyValue<double>("[IMS:1000046] pixel size");
      data->SetPropertyValue<double>("[IMS:1000046] pixel size x", v);
    }

    if (data->GetPropertyValue<double>("[IMS:1000047] pixel size y") == -1 && data->GetPropertyValue<double>("[IMS:1000046] pixel size x") > 0)
    {
      // Only IMS:1000046 was set (should be standard)
      auto v = data->GetPropertyValue<double>("[IMS:1000046] pixel size x");
      double sqs = std::sqrt(v);
      data->SetPropertyValue<double>("[IMS:1000046] pixel size x", m2::MicroMeterToMilliMeter(sqs));
      data->SetPropertyValue<double>("[IMS:1000047] pixel size y", m2::MicroMeterToMilliMeter(sqs));
      data->SetPropertyValue<double>("squared pixel size", v);
    }
    else if (data->GetPropertyValue<double>("[IMS:1000047] pixel size y") > 0 && data->GetPropertyValue<double>("[IMS:1000046] pixel size x") > 0)
    {
      // IMS:1000046 and IMS:1000047 was set
      // Note IMS:1000047 is used in a false way
      // -> in this case assumed to take the y value
      // if pixel size y is used > do not sqrt
      auto pixelSizeXInMilliMetre = m2::MicroMeterToMilliMeter(data->GetPropertyValue<double>("[IMS:1000046] pixel size x"));
      auto pixelSizeYInMilliMetre = m2::MicroMeterToMilliMeter(data->GetPropertyValue<double>("[IMS:1000047] pixel size y"));
      data->SetPropertyValue("[IMS:1000046] pixel size x", pixelSizeXInMilliMetre);
      data->SetPropertyValue("[IMS:1000047] pixel size y", pixelSizeYInMilliMetre);
    }

    if (data->GetPropertyValue<double>("[IMS:1000047] pixel size y") <= 0 && data->GetPropertyValue<double>("[IMS:1000046] pixel size x") <= 0)
    {
      data->SetPropertyValue("[IMS:1000046] pixel size x", m2::MicroMeterToMilliMeter(50));
      data->SetPropertyValue("[IMS:1000047] pixel size y", m2::MicroMeterToMilliMeter(50));
      data->SetPropertyValue(
        "pixel size info",
        std::string("Pixel size x and y are default values, due to missing imzTags IMS:1000046 and IMS:1000047!"));
      MITK_WARN << "No pixel size found, set x and y spacing to 50 microns!";
    }

    if (data->GetPropertyValue<double>("pixel size z") < 0)
    {
      data->SetPropertyValue("pixel size z", m2::MicroMeterToMilliMeter(10));
      data->SetPropertyValue<unsigned>("max count of pixels z", 1);
    }
    else
    {
      auto pixelSizeZInMilliMetre = m2::MicroMeterToMilliMeter(data->GetPropertyValue<double>("pixel size z"));
      data->SetPropertyValue("pixel size z", pixelSizeZInMilliMetre);
    }
    // -------- END META DATA --------
  }
}

void m2::ImzMLParser::ReadImageSpectrumMetaData(m2::ImzMLSpectrumImage::Pointer data)
{
  std::ifstream f;
  std::vector<std::string> stack, context_stack;
  std::string line, context, tag, name, value, accession;
  tag.reserve(60);
  name.reserve(60);
  value.reserve(60);
  context.reserve(60);
  accession.reserve(60);

  std::unordered_map<std::string, std::function<void(const std::string &)>> accession_map;
  std::unordered_map<std::string, std::function<void(const std::string &)>> context_map;

  f.open((data->GetImzMLDataPath()), std::ios_base::binary);

  std::map<std::string, unsigned> precisionDict = {{"32-bit float", sizeof(float)},
                                                   {"64-bit float", sizeof(double)},
                                                   {"32-bit integer", sizeof(int32_t)},
                                                   {"64-bit integer", sizeof(int64_t)}};

  {
    accession_map.clear();
    context_map.clear();

    auto &spectra = data->GetSpectra();

    context_map["spectrumList"] = [&](auto line)
    {
      unsigned count = std::stoul(attributeValue(line, "count", value));
      data->SetPropertyValue<unsigned>("number of measurements", count);
      spectra.resize(count);
    };

    size_t spectrumIndexReference = 0;

    context_map["spectrum"] = [&](auto line)
    {
      attributeValue(line, "index", value);
      spectra[spectrumIndexReference].index.SetElement(2, 0);
      // TODO: Backtracability to spectrum in file may be lost.
      //  Maybe we can add the id string to the BinarySpectrumMetaData definition.
      //  e.g. spectrum.tag = attributValue(line, "id", value);
    };

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L196
    accession_map["IMS:1000050"] = [&](auto line)
    { spectra[spectrumIndexReference].index.SetElement(0, std::stol(attributeValue(line, "value", value)) - 1); };

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L204
    accession_map["IMS:1000051"] = [&](auto line)
    { spectra[spectrumIndexReference].index.SetElement(1, std::stol(attributeValue(line, "value", value)) - 1); };

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L213
    accession_map["IMS:1000052"] = [&](auto line)
    { spectra[spectrumIndexReference].index.SetElement(2, std::stol(attributeValue(line, "value", value)) - 1); };

    // Foreign user tags
    accession_map["3DPositionX"] = [&](auto line)
    { spectra[spectrumIndexReference].world.x = std::stod(attributeValue(line, "value", value)); };
    accession_map["3DPositionY"] = [&](auto line)
    { spectra[spectrumIndexReference].world.y = std::stod(attributeValue(line, "value", value)); };
    accession_map["3DPositionZ"] = [&](auto line)
    { spectra[spectrumIndexReference].world.z = std::stod(attributeValue(line, "value", value)); };

    accession_map["MS:1000285"] = [&](auto line)
    { spectra[spectrumIndexReference].inFileNormalizationFactor = std::stod(attributeValue(line, "value", value)); };

    context_map["referenceableParamGroupRef"] = [&](auto line) { attributeValue(line, "ref", context); };

    auto mzArrayRefName = data->GetPropertyValue<std::string>("m2aia.imzml.mzGroupID");
    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L303
    accession_map[mzArrayRefName + ".IMS:1000102"] = [&](auto line)
    { spectra[spectrumIndexReference].mzOffset = std::stoull(attributeValue(line, "value", value)); };

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L311
    accession_map[mzArrayRefName + ".IMS:1000103"] = [&](auto line)
    { spectra[spectrumIndexReference].mzLength = std::stoull(attributeValue(line, "value", value)); };

    auto intensityArrayRefName = data->GetPropertyValue<std::string>("m2aia.imzml.intensityGroupID");
    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L303
    accession_map[intensityArrayRefName + ".IMS:1000102"] = [&](auto line)
    { spectra[spectrumIndexReference].intOffset = std::stoull(attributeValue(line, "value", value)); };

    // https://github.com/m2aia/imzML/blob/master/imagingMS.obo#L311
    accession_map[intensityArrayRefName + ".IMS:1000103"] = [&](auto line)
    { spectra[spectrumIndexReference].intLength = std::stoull(attributeValue(line, "value", value)); };

    std::vector<char> buff;
    std::list<std::thread> threads;
    bool _ScilsTag3DCoordinateUsed = false;
    {
      while (!f.eof())
      {
        std::getline(f, line); // read the next line from the file
        if (line.find("</") != std::string::npos)
        {
          if (line.find("spectrum") != std::string::npos)
            ++spectrumIndexReference;
          continue;
        }
        if (line.rfind("/>") != std::string::npos) // check for end-tag
        {                                          // check for empty element-tag
          attributeValue(line, "accession", accession);
          if (!accession.empty())
          {
            // if a context is given
            if (!context.empty())
            {
              // call context specific accession
              if (evaluateAccession(line, context + "." + accession, accession_map))
                continue;
            }
            // fall back to context-less accession
            if (evaluateAccession(line, accession, accession_map))
              continue;
          }
          else
          {
            // e.g. support old 3D imzML Data (SciLs specific tags)
            attributeValue(line, "name", name);
            if (name.compare("3DPositionZ") == 0)
            {
              evaluateAccession(line, name, accession_map);
              _ScilsTag3DCoordinateUsed = true;
            }
          }
        }

        GetElementName(line, tag);
        EvaluateContext(line, tag, context_map);
      }
    }

    data->SetPropertyValue<double>("pixel size z", m2::MicroMeterToMilliMeter(10));

    std::set<int> uniques;
    if (_ScilsTag3DCoordinateUsed)
    { // check z world uniques
      // MITK_INFO << "SciLs 3D tag found";
      for (auto &s : spectra)
        uniques.insert(s.world.z);

      // MITK_INFO << "\t" << uniques.size() << " unique z positions found:";
      // std::copy(std::begin(uniques), std::end(uniques), std::ostream_iterator<unsigned int>{std::cout, ", "});

      std::map<unsigned, unsigned> worldToIndexMap;
      unsigned i = 0;
      for (const auto &u : uniques)
        worldToIndexMap[u] = i++;

      if (uniques.size() > 1)
      {
        unsigned zSpacing;
        unsigned zCount;
        std::list<unsigned> diffs, diffs_uniques;
        std::adjacent_difference(std::begin(uniques), std::end(uniques), std::back_inserter(diffs));
        diffs_uniques = diffs;
        diffs_uniques.sort();
        diffs_uniques.erase(std::unique(std::begin(diffs_uniques), std::end(diffs_uniques)), std::end(diffs_uniques));
        // MITK_INFO << "\t" << diffs_uniques.size() << " unique z distances found:";
        // std::copy(
        //   std::begin(diffs_uniques), std::end(diffs_uniques), std::ostream_iterator<unsigned>{std::cout, ", "});
        unsigned maxCount = 0;
        for (auto &uDiff : diffs_uniques)
        {
          unsigned count = std::count(std::begin(diffs), std::end(diffs), uDiff);
          // MITK_INFO << "Different values " << uDiff << " were found " << count << " times.";
          if (maxCount < count)
          {
            maxCount = count;
            zSpacing = uDiff;
          }
        }

        auto a = uniques.begin();
        auto b = std::next(uniques.begin(), 1);
        zSpacing = (*b) - (*a);
        zCount = uniques.size();
        //          bool startsByZero = uniques.find(0) != uniques.end();

        // Transform z physical coordinates to index coordinates
        for (auto &s : spectra)
          s.index.SetElement(2, worldToIndexMap[static_cast<unsigned>(s.world.z)]);

        data->SetPropertyValue<unsigned>("max count of pixels z", zCount);
        data->SetPropertyValue<double>("pixel size z", m2::MicroMeterToMilliMeter(zSpacing));
      }
    }

    bool minimalArea = true;  
    if(auto *preferencesService = mitk::CoreServices::GetPreferencesService())
      if(auto *preferences = preferencesService->GetSystemPreferences())
        minimalArea = preferences->GetBool("m2aia.view.image.minimal_area", true);

    if (minimalArea)
    {
      // Load only minimal area where valid spectra exist.
      // https://github.com/m2aia/m2aia/issues/45
      itk::IndexValueType xsMin, ysMin, zsMin;
      itk::IndexValueType xsMax, ysMax, zsMax;
      xsMin = ysMin = zsMin = std::numeric_limits<itk::IndexValueType>::max();
      xsMax = ysMax = zsMax = 0;

      for (auto &s : spectra)
      {
        xsMin = std::min(s.index[0], xsMin);
        xsMax = std::max(s.index[0], xsMax);

        ysMin = std::min(s.index[1], ysMin);
        ysMax = std::max(s.index[1], ysMax);

        zsMin = std::min(s.index[2], zsMin);
        zsMax = std::max(s.index[2], zsMax);
      }

      for (auto &s : spectra)
      {
        s.index[0] -= xsMin;
        s.index[1] -= ysMin;
        s.index[2] -= zsMin;
      }

      bool requireCorrectionX = xsMin;
      bool requireCorrectionY = ysMin;
      bool requireCorrectionZ = zsMin;

      auto imzMLSizeX = data->GetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x");
      auto imzMLSizeY = data->GetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y");
      auto imzMLSizeZ = data->GetPropertyValue<unsigned>("max count of pixels z");

      auto newSizeX = unsigned(xsMax - xsMin + 1);
      auto newSizeY = unsigned(ysMax - ysMin + 1);
      auto newSizeZ = unsigned(zsMax - zsMin + 1);

      data->SetPropertyValue<unsigned>("(original imzML value) [IMS:1000042] max count of pixels x", imzMLSizeX);
      data->SetPropertyValue<unsigned>("(original imzML value) [IMS:1000043] max count of pixels y", imzMLSizeY);
      data->SetPropertyValue<unsigned>("(original imzML value) max count of pixels z", imzMLSizeZ);

      data->SetPropertyValue<unsigned>("[IMS:1000042] max count of pixels x", newSizeX);
      data->SetPropertyValue<unsigned>("[IMS:1000043] max count of pixels y", newSizeY);
      data->SetPropertyValue<unsigned>("max count of pixels z", newSizeZ);

      if (requireCorrectionX)
        MITK_WARN << "The x coordinate index was shifted by " << xsMin << " towards 0.";
      if (requireCorrectionY)
        MITK_WARN << "The y coordinate index was shifted by " << ysMin << " towards 0.";
      if (requireCorrectionZ)
        MITK_WARN << "The z coordinate index was shifted by " << zsMin << " towards 0.";
      if (imzMLSizeX != newSizeX)
        MITK_WARN << "The max count of pixels x was adjusted from " << imzMLSizeX << " to " << newSizeX << "";
      if (imzMLSizeX != newSizeX)
        MITK_WARN << "The max count of pixels y was adjusted from " << imzMLSizeY << " to " << newSizeY << "";
      if (imzMLSizeX != newSizeX)
        MITK_WARN << "The max count of pixels z was adjusted from " << imzMLSizeZ << " to " << newSizeZ << "";
      
      // std::ofstream ofs("/home/jtfc/m2aia_imzml_correction.txt");
      // for(auto s : data->GetSpectra()){
      //   ofs << s.index[0] << " " << s.index[1] << " " << s.index[2] << std::endl;
      //   ofs << s.intLength << " " << s.intOffset << " " << s.mzLength << " " << s.mzOffset << std::endl;
      //   ofs << s.world.x << " " << s.world.y << " " << s.world.z << std::endl;
      //   ofs << "----------------" << std::endl;
      //   ofs.flush();
      // }
    }
  }
}