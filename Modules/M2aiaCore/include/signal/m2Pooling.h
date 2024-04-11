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
#include <signal/m2Normalization.h>
#include <signal/m2SignalCommon.h>

namespace m2
{


  namespace Signal
  {
    template <typename PoolingValue, class ItFirst, class ItLast>
    PoolingValue RangePooling(ItFirst first, ItLast last, RangePoolingStrategyType strategy)
    {
      PoolingValue val = 0;
      if(std::distance(first,last) == 0)
        return val;
      switch (strategy)
      {
        case RangePoolingStrategyType::None:
          break;
        case RangePoolingStrategyType::Sum:
          val = std::accumulate(first, last, PoolingValue(0));
          break;
        case RangePoolingStrategyType::Mean:
          val = std::accumulate(first, last, PoolingValue(0)) / PoolingValue(std::distance(first, last));
          break;
        case RangePoolingStrategyType::Maximum:
          val = *std::max_element(first, last);
          break;
        case RangePoolingStrategyType::Median:
          std::vector<typename ItFirst::value_type> v(first, last);
          val = m2::Signal::Median(v.begin(),v.end());
          break;
      }

      return val;
    }
  }
} // namespace m2
