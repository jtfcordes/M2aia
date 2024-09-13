/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#include <m2IntervalVectorIO.h>
#include <mitkIOUtil.h>
#include <mitkLocaleSwitch.h>
#include <itksys/SystemTools.hxx>
#include <boost/algorithm/string/join.hpp>

namespace m2
{
  IntervalVectorIO::IntervalVectorIO() : AbstractFileIO(m2::IntervalVector::GetStaticNameOfClass(), INTERVALVECTOR_MIMETYPE(), "Centroid Data (M²aia)")
  {
    AbstractFileWriter::SetRanking(10);
    AbstractFileReader::SetRanking(10);
    this->RegisterService();
  }

  mitk::IFileIO::ConfidenceLevel IntervalVectorIO::GetWriterConfidenceLevel() const
  {
    MITK_INFO << this->GetInput() <<  dynamic_cast<const m2::IntervalVector *>(this->GetInput());
    const auto *input = dynamic_cast<const m2::IntervalVector *>(this->GetInput());
    if (input)
      return Supported;
    else
      return Unsupported;
  }

  void IntervalVectorIO::Write()
  {
    mitk::LocaleSwitch localeSwitch("C");
    ValidateOutputLocation();

    const auto *input = static_cast<const m2::IntervalVector *>(this->GetInput());
    std::ofstream file(GetOutputLocation());
    
    if (!file.is_open()) {
        std::cout << "Failed to open file: " << GetOutputLocation() << std::endl;
        return;
    }

    std::array<std::string, 3> position;
    if(auto prop = input->GetProperty("m2aia.spectrum.position.x"))
      position[0] = prop->GetValueAsString();
    if(auto prop = input->GetProperty("m2aia.spectrum.position.y"))
      position[1] = prop->GetValueAsString();
    if(auto prop = input->GetProperty("m2aia.spectrum.position.z"))
      position[2] = prop->GetValueAsString();

    // file << "[spectrum position] " << boost::algorithm::join(position, " ") << "\n"; 
    file << "center,max,min,mean\n";
    const auto intervals = input->GetIntervals();
    for(size_t i = 0 ; i < intervals.size(); ++i){
      if(i==0)
        file << intervals[i].x.max()+0.005 << "," << intervals[i].y.max() << "," << intervals[i].y.min() << "," << intervals[i].y.mean() << "\n";
      else if(i == intervals.size() -1)
        file << intervals[i].x.min()-0.005 << "," << intervals[i].y.max() << "," << intervals[i].y.min() << "," << intervals[i].y.mean() << "\n";
      else
        file << intervals[i].x.mean() << "," << intervals[i].y.max() << "," << intervals[i].y.min() << "," << intervals[i].y.mean() << "\n";
    }

    file.close();
  } // namespace mitk

  mitk::IFileIO::ConfidenceLevel IntervalVectorIO::GetReaderConfidenceLevel() const
  {
    if (AbstractFileIO::GetReaderConfidenceLevel() == Unsupported)
      return Unsupported;
    return Supported;
  }

  std::vector<mitk::BaseData::Pointer> IntervalVectorIO::DoRead()
  {  
    auto path = this->GetInputLocation();

    std::ifstream file(path);

    if (!file.is_open()) {
        std::cout << "Failed to open file: " << path << std::endl;
        return {};
    }

    std::string line;
    std::vector<std::string> header;

    if(std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            header.push_back(cell);
        }
    }

    
    std::vector<std::vector<std::string>> data;

    MITK_INFO << "Read centroid list [" << path << "] with " << header.size() << " columns";
    std::cout << boost::algorithm::join(header," ,") << std::endl;
    int row = 0;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::cout << "[" << std::to_string(row) << "] ";
        data.push_back({});

        while (std::getline(ss, cell, ',')) {
            data[row].push_back(cell);
        }
        std::cout << boost::algorithm::join(data[row]," ,") << std::endl;
        ++row;
    }

    file.close();
    
    auto Strip = [](const std::string & str){
      size_t start = 0;
      size_t end = str.length();

      // Find the index of the first non-whitespace character
      while (start < end && std::isspace(str[start])) {
          ++start;
      }

      // Find the index of the last non-whitespace character
      while (end > start && std::isspace(str[end - 1])) {
          --end;
      }

      // Return the stripped substring
      return str.substr(start, end - start);
      
      };

    

    
    
    auto targetData = m2::IntervalVector::New();
    targetData->SetType(m2::SpectrumFormat::Centroid);
    targetData->SetInfo("Load from file: " + path);
    std::transform(std::begin(data), std::end(data), std::back_inserter(targetData->GetIntervals()), [Strip](const std::vector<std::string> & row){
      auto x = std::stod(Strip(row[0]));
      auto y = std::stod(Strip(row[1]));
      m2::Interval I(x,y);
      return I;
    });

    // check for labels. If labels exist split into multiple lists
    auto labelIt = std::find(std::begin(header), std::end(header), "label");
    if(labelIt != std::end(header)){
      auto colIdx = std::distance(std::begin(header), labelIt);

      std::vector<unsigned int> labels;
      std::transform(std::begin(data), std::end(data), std::back_inserter(labels), [colIdx](std::vector<std::string> & v){return std::stoul(v[colIdx]);}); 
      
      auto eraseIt = std::unique(std::begin(labels), std::end(labels));
      labels.erase(eraseIt, std::end(labels));
      
      std::map<std::string, m2::IntervalVector::Pointer> labeledResults;
      auto intervalIt = std::begin(targetData->GetIntervals());
        
      for(const std::vector<std::string> &row : data){
        if(labeledResults[row[colIdx]].IsNull()){
          labeledResults[row[colIdx]] = m2::IntervalVector::New();
        }
        labeledResults[row[colIdx]]->GetIntervals().push_back(*intervalIt);
        ++intervalIt;
      }

      std::vector<mitk::BaseData::Pointer> res;
      for(auto kv : labeledResults)
        res.push_back(kv.second);

      return res;
    }

    

    return {targetData};
  }


  IntervalVectorIO *IntervalVectorIO::IOClone() const { return new IntervalVectorIO(*this); }
} // namespace m2
