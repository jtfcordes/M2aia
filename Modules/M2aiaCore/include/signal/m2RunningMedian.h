/*===================================================================
MIT License

Copyright (c) 2018 CSEM SA

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
===================================================================*/
#pragma once

#include <M2aiaCoreExports.h>
#include <algorithm>
#include <cstring>
#include <stdlib.h>
#include <vector>

namespace m2
{
  class M2AIACORE_EXPORT RunMedian
  {
  public:
    template <class IteratorType>
    static void apply(IteratorType start, IteratorType end, unsigned int s, IteratorType baseline_start) noexcept
    {
      MedfiltData data;
      s = s * 2 + 1;
      MedfiltNode *nodes = new MedfiltNode[s];

      medfilt_init(&data, nodes, s, *start);

      auto oit = baseline_start;
      
      for (auto it = start; it != end; ++it)
      {
        double min, mid, max;
        medfilt(&data, *it, &mid, &min, &max);
        if (mid < 0)
          *oit = max;
        else
          *oit = mid;
        ++oit;
      }

      delete[] nodes;
    }
  protected:
    typedef struct MedfiltNode
    {
      double value;
      size_t index; // Node index in the sorted table
      struct MedfiltNode *parent;
      struct MedfiltNode *sorted;
    } MedfiltNode;

    typedef struct MedfiltData
    {
      MedfiltNode *kernel; // Working filter memory
      MedfiltNode *oldest; // Reference to the oldest value in the kernel
      size_t length;       // Number of nodes
    } MedfiltData;


    static void swap(MedfiltNode **a, MedfiltNode **b) noexcept
    {
      // Swap two node references in the sorted table.
      MedfiltNode *temp = *a;
      *a = *b;
      *b = temp;

      // Preserve index. Used to retrive the node position in the sorted table.
      size_t index = (*a)->index;
      (*a)->index = (*b)->index;
      (*b)->index = index;
    }

    static void medfilt(MedfiltData *data, double input, double *median, double *min, double *max) noexcept
    {
      // New value replaces the oldest
      MedfiltNode *n = data->kernel;
      MedfiltNode *node = data->oldest;
      node->value = input;
      data->oldest = node->parent;

#define VAL(x) (n[x].sorted->value)

      // Sort the kernel

      // Check neigbhour to the right
      for (size_t i = node->index; i < data->length - 1 && VAL(i) > VAL(i + 1); i++)
      {
        swap(&n[i].sorted, &n[i + 1].sorted);
      }

      // Check neigbhour to the left
      for (size_t i = node->index; i > 0 && VAL(i) < VAL(i - 1); i--)
      {
        swap(&n[i].sorted, &n[i - 1].sorted);
      }

      // Get kernel information from sorted entries
      *min = VAL(0);
      *max = VAL(data->length - 1);
      *median = VAL(data->length / 2);

#undef VAL
    }

    static void medfilt_init(MedfiltData *data, MedfiltNode *nodes, size_t length, double init) noexcept
    {
      data->kernel = nodes;
      data->length = length;

      // Linked list initialization
      data->oldest = &data->kernel[length - 1];
      for (size_t i = 0; i < length; i++)
      {
        data->kernel[i].value = init;
        data->kernel[i].parent = data->oldest;
        data->kernel[i].index = i;
        data->kernel[i].sorted = &data->kernel[i];

        data->oldest = &data->kernel[i];
      }
    }
  };
} // namespace m2