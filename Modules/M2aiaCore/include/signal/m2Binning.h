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
#include "m2SignalCommon.h"
#include <algorithm>
#include <cmath>
#include <m2IntervalVector.h>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <m2IntervalVector.h>

namespace m2
{
  namespace Signal
  {
    template <class PeakItFirst, class PeakItLast, class OutIterType>
    inline void toleranceBasedBinning(
      PeakItFirst s, PeakItLast e, OutIterType output, double tolerance, bool absoluteDistance = false)
    {
      if (s == e)
        return;

      unsigned int size = 0;
      const auto accPeak = std::accumulate(s, e, m2::Interval());
      const auto meanX = accPeak.GetX();

      bool dirty = false;
      decltype(s) i = s;
      while (i != e)
      {
        if (std::abs(i->GetX() - meanX) >= ((meanX * tolerance) * (!absoluteDistance) + (absoluteDistance * tolerance)))
        {
          dirty = true;
          break;
        }
        ++i;
      }

      if (dirty && size >= 2)
      {
        decltype(s) pivot;
        {
          decltype(s) s0 = s;
          decltype(s) s1 = std::next(s);
          pivot = s1;
          double max = 0; // max distance of mz neighbors is split position
          while (s1 != e)
          {
            if (s1->GetX() - s0->GetX() > max)
            {
              max = s1->GetX() - s0->GetX();
              pivot = s1;
            }
            ++s1;
            ++s0;
          }
        }

        toleranceBasedBinning(s, pivot, output, tolerance, absoluteDistance);
        toleranceBasedBinning(pivot, e, output, tolerance, absoluteDistance);
      }
      else
      {
        // output a new peak
        m2::Interval newPeak = *s;
        ++s;

        std::for_each(s, e, [&newPeak](const m2::Interval &p) { newPeak += p; });

        (*output) = newPeak;
        ++output;
      }
    }

    template <class XIt, class YIt, class OutIterType>
    inline void equidistantBinning(XIt xFirst, XIt xEnd, YIt yFirst, OutIterType output, unsigned int bins)
    {
      if (xFirst == xEnd)
        return;

      // bin size
      double range = *(xEnd - 1) - *xFirst;
      double binSize = range / double(bins);

      unsigned int currentBin = 0;
      unsigned int currentIndex = 0;

      auto yIt = yFirst;
      auto xIt = xFirst;
      auto oIt = output;
      while (xIt != xEnd)
      {
        m2::Interval p;

        // accumulate
        while (xIt != xEnd && *xIt < (currentBin + 1) * binSize)
        {
          // p.index(currentIndex);
          p.x.add(*xIt);
          p.y.add(*yIt);
          ++xIt;
          ++yIt;
          ++currentIndex;
        }

        if (p.x.count())
        {           // avoid empty peaks
          *oIt = p; // assign peak
          ++output;
        }
        ++currentBin;
      }
    }

    template <typename Iter>
    bool isUniqueRange(Iter first, Iter last) {
      std::set<typename std::iterator_traits<Iter>::value_type> s;
      return std::all_of(first, last, [&](const auto& x) {
          return s.insert(x).second;
      });
    }


    template <typename TypeBeginX, typename TypeEndX, typename TypeBeginY, typename TypeBeginS>
    inline double grouperStrict(
      TypeBeginX sx, TypeEndX ex, TypeBeginY /*sy*/, TypeBeginS ss, double tolerance)
    {
      using namespace std;
      auto n = distance(sx, ex);
      if (!isUniqueRange(ss, ss+n))
        return 0;
      double mean_xs = std::accumulate(sx, ex, 0.0) / double(n);
      if (std::any_of(sx, ex, [=](double x) { return (std::abs(x - mean_xs) / mean_xs) > tolerance; }))
        return 0;
      return mean_xs;
    }
    

    template <typename TypeBeginX, typename TypeEndX, typename TypeBeginY, typename TypeBeginS>
    inline double grouperRelaxed(
      TypeBeginX sx, TypeEndX ex, TypeBeginY /*sy*/, TypeBeginS /*ss*/, double tolerance)
    {
      using namespace std;
      // if (!isUniqueRange(ss, ss+n))
      //   return 0;
      auto n = distance(sx, ex);
      double mean_xs = std::accumulate(sx, ex, 0.0) / double(n);
      if (std::any_of(sx, ex, [=](double x) { return (std::abs(x - mean_xs) / mean_xs) > tolerance; }))
        return 0;
      return mean_xs;
    }

    template <typename XType, typename YType, typename SourceType, typename Functor>
    inline std::tuple<std::vector<m2::Interval>, std::vector<std::vector<SourceType>>> groupBinning(
      const std::vector<XType> &xs,
      const std::vector<YType> &ys,
      const std::vector<SourceType> &sources,
      // const std::vector<SourceType> &indices,
      Functor f,
      double tolerance)
    {
      using namespace std;
      using xvec = vector<XType>;
      // using yvec = vector<YType>;
      using svec = vector<SourceType>;

      xvec d;
      d.resize(xs.size());
      adjacent_difference(xs.begin(), xs.end(), begin(d));
      d[0] = -1;


      std::vector<m2::Interval> I;
      
      // svec bin_assignments;
      // bin_assignments.resize(xs.size(), -1);

      auto n = static_cast<int>(xs.size());
      vector<pair<int, int>> boundary{{0, n}};
      boundary.reserve(xs.size());
      
      vector<svec> bin_counts;
      bin_counts.reserve(xs.size());
            
      std::vector<m2::Interval> intervals;
      intervals.reserve(xs.size());
      
      int current_id = 0;
      while (!boundary.empty())
      {
        auto left = boundary.back().first;
        auto right = boundary.back().second;
        boundary.pop_back();


        int gapIdx = distance(begin(d), max_element(begin(d)+left+1,begin(d)+right));
        // MITK_INFO << " ---- [" << left << "/" << gapIdx << "/" << right << "]";
        

        auto l = f(begin(xs) + left, begin(xs) + gapIdx, begin(ys) + left, begin(sources) + left, tolerance);
        if (l == 0){
          boundary.push_back({left, gapIdx});
          // MITK_INFO << "left : " << left << " " << gapIdx;
        }
        else
        {
          m2::Interval i;
          auto yIt = begin(ys) + left;
          auto xIt = begin(xs) + left;
          for(; xIt != begin(xs) + gapIdx; ++xIt, ++yIt){
            i.x.add(*xIt);
            i.y.add(*yIt);
          }

          
          // fill(begin(bin_assignments) + left, begin(bin_assignments) + gapIdx, current_id);
          intervals.emplace_back(i);
          bin_counts.emplace_back(begin(sources) + left, begin(sources) + gapIdx);
          current_id += 1;
          
        }

        auto r = f(begin(xs) + gapIdx, begin(xs) + right, begin(ys) + gapIdx, begin(sources) + gapIdx, tolerance);
        if (r == 0){
          boundary.push_back({gapIdx, right});
          // MITK_INFO << "right : " << gapIdx << " " << right;
        }
        else
        {
          m2::Interval i;
          auto yIt = begin(ys) + gapIdx;
          auto xIt = begin(xs) + gapIdx;
          for(; xIt != begin(xs) + right; ++xIt, ++yIt){
            i.x.add(*xIt);
            i.y.add(*yIt);
          }

          // fill(begin(bin_assignments) + gapIdx, begin(bin_assignments) + right, current_id);
          // bin_xs.push_back(r);
          intervals.emplace_back(i);
          bin_counts.emplace_back(begin(sources) + gapIdx, begin(sources) + right);
          current_id += 1;
          // MITK_INFO << "right [pure]: " << gapIdx << " " << right << " <" << r << ">";
        }
      }

      // int binId = 0;
      // int curr = bin_assignments[0];
      // transform(begin(bin_assignments), end(bin_assignments), begin(bin_assignments), [&binId, &curr](auto & val){
      //   if(val != curr){
      //     ++binId;
      //     curr = val;
      //   }
      //   return binId;
      // });

      intervals.shrink_to_fit();
     

      // vector<svec> sorted_bin_counts;
      // sorted_bin_xs.resize(bin_xs.size());
      // sorted_bin_counts.resize(bin_counts.size());

      auto indices = m2::argsort(intervals);
      auto sorted_bin_xs = m2::argsortApply(intervals, indices);
      auto sorted_bin_counts = m2::argsortApply(bin_counts, indices);


    return std::make_tuple(sorted_bin_xs, sorted_bin_counts);

    }

    //     auto l = grouper(xs.cbegin() + left, xs.cbegin() + gapIdx, ys.cbegin() + left, ss.cbegin() + left,
    //     tolerance); if (l.first.empty()) {
    //         boundary.emplace_back(left, gapIdx);
    //     } else {
    //         std::fill_n(bin_assignments.begin() + left, gapIdx - left, current_id);
    //         bin_xs.push_back(l.first[0]);
    //         bin_counts.push_back(l.second);
    //         current_id += 1;
    //     }

    //     auto r = grouper(xs.cbegin() + gapIdx, xs.cbegin() + right, ys.cbegin() + gapIdx, ss.cbegin() + gapIdx,
    //     tolerance); if (r.first.empty()) {
    //         boundary.emplace_back(gapIdx, right);
    //     } else {
    //         std::fill_n(bin_assignments.begin() + gapIdx, right - gapIdx, current_id);
    //         bin_xs.push_back(r.first[0]);
    //         bin_counts.push_back(r.second);
    //         current_id += 1;
    //     }
    // }

    // }

  } // namespace Signal
} // namespace m2