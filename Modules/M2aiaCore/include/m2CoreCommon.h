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

#include <mitkDataNode.h>
#include <mitkLabelSetImage.h>
#include <type_traits>

namespace m2
{

  enum class SpectrumFormat : unsigned int
  {
    None = 0,
    Profile = 1,
    Centroid = 2,
    Continuous = 4,
    Processed = 8,
    ContinuousProfile = 16 + 4 + 1,
    ProcessedProfile = 32 + 8 + 1,
    ContinuousCentroid = 64 + 4 + 2,
    ProcessedCentroid = 128 + 8 + 2
  };

  enum class ImzMLFormat : unsigned int
  {
    None = 0,
    Continuous = 1,
    Processed = 2

  };

  inline std::string to_string(const SpectrumFormat &type) noexcept
  {
    switch (type)
    {
      case SpectrumFormat::Centroid:
        return "Centroid";
      case SpectrumFormat::Profile:
        return "Profile";
      case SpectrumFormat::Continuous:
        return "Continuous";
      case SpectrumFormat::Processed:
        return "Processed";
      case SpectrumFormat::ContinuousProfile:
        return "ContinuousProfile";
      case SpectrumFormat::ProcessedProfile:
        return "ProcessedProfile";
      case SpectrumFormat::ContinuousCentroid:
        return "ContinuousCentroid";
      case SpectrumFormat::ProcessedCentroid:
        return "ProcessedCentroid";
      case SpectrumFormat::None:
        return "None";
    }
    return "";
  }

  enum class SpectrumType : unsigned int
  {
    None = 0,
    Mean = 1,
    Median = 2,
    Maximum = 3,
    Sum = 4,
    Variance = 5,
  };

  inline std::string to_string(const SpectrumType &type) noexcept
  {
    switch (type)
    {
      case SpectrumType::None:
        return "None";
      case SpectrumType::Mean:
        return "Mean";
      case SpectrumType::Median:
        return "Median";
      case SpectrumType::Maximum:
        return "Maximum";
      case SpectrumType::Sum:
        return "Sum";
      case SpectrumType::Variance:
        return "Variance";
    }
    return "";
  }

  enum class NumericType : unsigned int
  {
    None = 0,
    Float = 1,
    Double = 2
  };

  inline std::string to_string(const NumericType &type) noexcept
  {
    switch (type)
    {
      case NumericType::None:
        return "None";
      case NumericType::Float:
        return "Float";
      case NumericType::Double:
        return "Double";
    }
    return "";
  }

  inline unsigned int to_bytes(const NumericType &type) noexcept
  {
    switch (type)
    {
      case NumericType::None:
        return 0;
      case NumericType::Float:
        return sizeof(float);
      case NumericType::Double:
        return sizeof(double);
    }
    return 0;
  }

  inline double MicroMeterToMilliMeter(double x)
  {
    return x * 10e-4;
  }
  inline double MilliMeterToMicroMeter(double x)
  {
    return x * 10e2;
  }
  inline double PartPerMillionToFactor(double x)
  {
    return x * 10e-6;
  }

  //////////////////////////////////////////////////////////////////////
  /////////////////// ATTENTION ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  // If you have added a new enum value to the list.
  // Ensure that all enum values are well defined in the m2::CORE_MAPPINGS.
  // m2::CORE_MAPPINGS is defined in <m2CoreCommon.h>.
  // m2::CORE_MAPPINGS is used to transform strings into unsigned int values.
  // A compiler error will be thrown if you enabled Testing and not added the
  // new enum value to the test. <m2CoreMappingsTest.cpp>
  //
  //////////////////////////////////////////////////////////////////////
  /////////////////// ATTENTION ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  const std::map<const std::string, unsigned int> CORE_MAPPINGS{{"None", 0},
                                                                {"ContinuousProfile", 1},
                                                                {"ProcessedProfile", 2},
                                                                {"ContinuousCentroid", 4},
                                                                {"ProcessedCentroid", 8},
                                                                {"Mean", 1},
                                                                {"Median", 2},
                                                                {"Maximum", 3},
                                                                {"Sum", 4},
                                                                {"Variance", 5},
                                                                {"PeakIndicators", 6},
                                                                {"Float", 0},
                                                                {"Double", 1}};

  using DisplayImagePixelType = double;
  using NormImagePixelType = double;
  using IndexType = unsigned int;
  using WorldCoordinateType = float;
  using ShiftImageType = int;
  using IndexImagePixelType = IndexType;

  enum class TransformationMethod
  {
    None,
    Transformix
  };

  /// m2Utils

  const auto Find = [](const auto &str, const auto &searchString, auto defaultValue, auto &map)
  {
    auto p = str.find(searchString);
    decltype(defaultValue) converted;

    if (p != std::string::npos)
    {
      char begin = ' ';
      char end = ')';
      auto s = str.find(begin, p);

      auto e = str.find(end, p);
      auto val = str.substr(s + 1, e - s - 1);
      map[searchString] = val;
      std::stringstream buffer(val);

      if (std::is_same<decltype(defaultValue), bool>::value)
        buffer >> std::boolalpha >> converted;
      else
        buffer >> converted;

      return converted;
    }

    std::stringstream buffer;
    if (std::is_same<decltype(defaultValue), bool>::value)
      buffer << std::boolalpha << defaultValue;
    else
      buffer << defaultValue;

    map[searchString] = buffer.str();
    return defaultValue;
  };

} // namespace m2

inline m2::SpectrumFormat operator|(m2::SpectrumFormat lhs, m2::SpectrumFormat rhs)
{
  return static_cast<m2::SpectrumFormat>(static_cast<std::underlying_type<m2::SpectrumFormat>::type>(lhs) |
                                         static_cast<std::underlying_type<m2::SpectrumFormat>::type>(rhs));
}

inline m2::SpectrumFormat operator|=(m2::SpectrumFormat lhs, m2::SpectrumFormat rhs)
{
  return lhs | rhs;
}

inline m2::SpectrumFormat operator&(m2::SpectrumFormat lhs, m2::SpectrumFormat rhs)
{
  return static_cast<m2::SpectrumFormat>(static_cast<std::underlying_type<m2::SpectrumFormat>::type>(lhs) &
                                         static_cast<std::underlying_type<m2::SpectrumFormat>::type>(rhs));
}

inline bool any(m2::SpectrumFormat lhs)
{
  return static_cast<std::underlying_type<m2::SpectrumFormat>::type>(lhs) != 0;
}
