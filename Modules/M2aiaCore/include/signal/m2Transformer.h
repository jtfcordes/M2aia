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
#include <functional>
#include <math.h>
#include <signal/m2SignalCommon.h>

namespace m2
{
  namespace Signal
  {
    template <class ItValueType>
    class IntensityTransformationFunctor
    {
    private:
      IntensityTransformationType m_strategy;

    public:
      void Initialize(IntensityTransformationType strategy)
      {
        m_strategy = strategy;
      }

      void operator()(typename std::vector<ItValueType>::iterator start,
                      typename std::vector<ItValueType>::iterator end)
      {
        switch (m_strategy)
        {
          case m2::IntensityTransformationType::Log10:
            std::transform(start, end, start, [](const auto & a){return std::log10(a+1);});
            break;
          case m2::IntensityTransformationType::Log2:
            std::transform(start, end, start, [](const auto & a){return std::log2(a+1);});
            break;
          case m2::IntensityTransformationType::SquareRoot:
            std::transform(start, end, start, [](const auto & a){return std::sqrt(a);});
            break;
          case m2::IntensityTransformationType::None:
            break;
        }
      }
    };
  } // namespace Signal
} // namespace m2
