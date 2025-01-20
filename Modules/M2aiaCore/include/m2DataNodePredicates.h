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

#include "M2aiaCoreExports.h"
#include <m2IntervalVector.h>
#include <m2SpectrumImage.h>
#include <m2SpectrumImageStack.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateFunction.h>
#include <mitkNodePredicateNot.h>

namespace m2
{
  namespace DataNodePredicates
  {

    mitk::NodePredicateFunction::Pointer IsProfileOrProcessedCentroidSpectrumImage = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node) -> bool
      {
        if (auto spectrumImage = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
          return ((unsigned int)(spectrumImage->GetSpectrumType().Format)) &
                   ((unsigned int)(m2::SpectrumFormat::Profile)) ||
                 ((unsigned int)(spectrumImage->GetSpectrumType().Format)) &
                   ((unsigned int)(m2::SpectrumFormat::ProcessedCentroid));
        return false;
      });

    mitk::NodePredicateFunction::Pointer IsSpectrumImage = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node) -> bool { return dynamic_cast<m2::SpectrumImage *>(node->GetData()); });

    mitk::NodePredicateFunction::Pointer IsSpectrumImageStack = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node) -> bool { return dynamic_cast<m2::SpectrumImageStack *>(node->GetData()); });

    mitk::NodePredicateFunction::Pointer IsVisible = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node) { return node->IsVisible(nullptr, "visible", false); });

    mitk::NodePredicateFunction::Pointer IsActiveHelperNode = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node) { return node->IsOn("helper object", nullptr, false); });

    mitk::NodePredicateFunction::Pointer IsCentroidSpectrum = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node) -> bool
      {
        if (auto intervals = dynamic_cast<m2::IntervalVector *>(node->GetData()))
          return ((unsigned int)(intervals->GetType())) & ((unsigned int)(m2::SpectrumFormat::Centroid));
        return false;
      });

    mitk::NodePredicateFunction::Pointer IsProfileSpectrum = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node) -> bool
      {
        if (auto intervals = dynamic_cast<m2::IntervalVector *>(node->GetData()))
          return ((unsigned int)(intervals->GetType())) & ((unsigned int)(m2::SpectrumFormat::Profile));
        return false;
      });

    mitk::NodePredicateFunction::Pointer IsOverviewSpectrum = mitk::NodePredicateFunction::New(
      [](const mitk::DataNode *node)
      {
        if (auto intervals = dynamic_cast<m2::IntervalVector *>(node->GetData()))
          return intervals->GetInfo().find("overview") != std::string::npos;
        return false;
      });

    mitk::NodePredicateNot::Pointer NoActiveHelper = mitk::NodePredicateNot::New(IsActiveHelperNode);

  } // namespace DataNodePredicates
} // namespace m2
