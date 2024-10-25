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
#include <array>
#include <map>
#include <string>


#include <type_traits>

template <typename E>
constexpr auto to_underlying(E e) noexcept
{
    return static_cast<std::underlying_type_t<E>>(e);
}

namespace m2
{
  enum class SmoothingType : unsigned int
  {
    None = 0,
    SavitzkyGolay = 1,
    Gaussian = 2
  };

  const std::array<std::string, 3> SmoothingTypeNames = {"None", "SavitzkyGolay", "Gaussian"};

  enum class RangePoolingStrategyType : unsigned int
  {
    None = 0,
    Mean = 1,
    Median = 2,
    Maximum = 3,
    Sum = 4
  };

  const std::array<std::string, 5> RangePoolingStrategyTypeNames = {"None", "Mean", "Median", "Maximum", "Sum"};

  // NormalizationStrategyType
  enum class NormalizationStrategyType : unsigned int
  {
    None = 0, TIC ,Sum, Mean, Max, RMS, Internal, External
  };

  const std::array<NormalizationStrategyType, 8> NormalizationStrategyTypeList = {
    NormalizationStrategyType::None,
    NormalizationStrategyType::TIC,
    NormalizationStrategyType::Sum,
    NormalizationStrategyType::Mean,
    NormalizationStrategyType::Max,
    NormalizationStrategyType::RMS,
    NormalizationStrategyType::Internal,
    NormalizationStrategyType::External
  };

  const std::array<std::string, 8> NormalizationStrategyTypeNames = {
    "None", "TIC", "Sum", "Mean", "Max", "RMS", "Internal", "External"};


  inline std::string to_string(m2::NormalizationStrategyType type){
    return NormalizationStrategyTypeNames.at(to_underlying(type));
  }


  // ImageNormalizationStrategyType
  enum class ImageNormalizationStrategyType : unsigned int
  {
    None = 0, MinMax, zScore, zScore2Sigma
  };


  const std::array<ImageNormalizationStrategyType, 4> ImageNormalizationStrategyTypeList = {
    ImageNormalizationStrategyType::None,
    ImageNormalizationStrategyType::MinMax,
    ImageNormalizationStrategyType::zScore,
    ImageNormalizationStrategyType::zScore2Sigma
  };


  const std::array<std::string, 4> ImageNormalizationStrategyTypeNames = {
    "None", "MinMax", "zScore", "zScore2Sigma"};

  inline std::string to_string(m2::ImageNormalizationStrategyType type){
    return ImageNormalizationStrategyTypeNames.at(to_underlying(type));
  }

  // ImageSmoothingStrategyType
  enum class ImageSmoothingStrategyType : unsigned int
  {
    None = 0, Median, Gaussian
  };

  const std::array<ImageSmoothingStrategyType, 3> ImageSmoothingStrategyTypeList = {
    ImageSmoothingStrategyType::None,
    ImageSmoothingStrategyType::Median,
    ImageSmoothingStrategyType::Gaussian
  };

  const std::array<std::string, 3> ImageSmoothingStrategyTypeNames = {
    "None", "Median", "Gaussian"};

  inline std::string to_string(m2::ImageSmoothingStrategyType type){
    return ImageSmoothingStrategyTypeNames.at(to_underlying(type));
  }

  // BaselineCorrectionType
  enum class BaselineCorrectionType : unsigned int
  {
    None = 0,
    TopHat = 1,
    Median = 2
  };

  const std::array<std::string, 3> BaselineCorrectionTypeNames = {"None", "TopHat", "Median"};

  enum class IntensityTransformationType : unsigned int
  {
    None = 0,
    Log2 = 1,
    Log10 = 2,
    SquareRoot = 3
  };

  const std::array<std::string, 4> IntensityTransformationTypeNames = {"None", "Log2", "Log10", "SquareRoot"};

  //////////////////////////////////////////////////////////////////////
  /////////////////// ATTENTION ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  // If you have added a new enum value to the list.
  // Ensure that all enum values are well defined in the m2::SIGNAL_MAPPINGS.
  // m2::SIGNAL_MAPPINGS is defined in <m2SignalCommon.h>.
  // m2::SIGNAL_MAPPINGS is used to transform strings into unsigned int values.
  // A compiler error will be thrown if you enabled Testing and not added the
  // new enum value to the test. <m2SignalMappingsTest.cpp>
  //
  //////////////////////////////////////////////////////////////////////
  /////////////////// ATTENTION ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  const std::map<const std::string, unsigned int> SMOOTHING_MAPPINGS{
    {"None", 0}, {"SavitzkyGolay", 1}, {"Gaussian", 2}};

  const std::map<const std::string, unsigned int> BASECOR_MAPPINGS{{"None", 0}, {"TopHat", 1}, {"Median", 2}};

  const std::map<const std::string, unsigned int> NORMALIZATION_MAPPINGS{
     {"None", 0}, {"TIC", 1}, {"Sum", 2}, {"Mean", 3}, {"Max",4}, {"RMS",5}, {"Internal",6}, {"External", 7}};

  const std::map<const std::string, unsigned int> IMAGENORMALIZATION_MAPPINGS{
     {"None", 0}, {"MinMax", 1}, {"zScore", 2}};

  const std::map<const std::string, unsigned int> POOLING_MAPPINGS{
    {"None", 0}, {"Mean", 1}, {"Median", 2}, {"Maximum", 3}, {"Sum", 4}};

  const std::map<const std::string, unsigned int> INTENSITYTRANSFORMATION_MAPPINGS{
    {"None", 0}, {"Log2", 1}, {"Log10", 2}, {"SquareRoot", 3}};

  template<typename ContainerType>
  std::vector<int> argsort(const ContainerType &v)
  {
    // Initialize vector of indices
    std::vector<int> indices(v.size());
    for (unsigned int i = 0; i < v.size(); ++i)
    {
      indices[i] = i;
    }
    // Sort indices based on values in v
    std::stable_sort(indices.begin(), indices.end(), [&](int i, int j) { return v[i] < v[j]; });
    return indices;
  }


  template<typename ContainerType, typename IndexContainer>
  ContainerType argsortApply(const ContainerType &v, const IndexContainer &indices)
  {
    ContainerType c;
    c.reserve(indices.size());
    auto it = std::back_inserter(c);
    for(auto i : indices)
      it = v[i];
    return c;
  }

} // namespace m2

