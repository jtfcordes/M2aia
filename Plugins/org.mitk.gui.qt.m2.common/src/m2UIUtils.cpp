/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#include "m2UIUtils.h"
#include <mitkDataStorage.h>
#include <mitkNodePredicateProperty.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>

m2::UIUtils::NodesVectorType::Pointer m2::UIUtils::AllNodes(const mitk::DataStorage * dataStorage)
{
  auto nodesToProcess = itk::VectorContainer<unsigned int, mitk::DataNode::Pointer>::New();
  // is checked: all nodes that fit the predicate are processed
  auto predicate =
    mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<m2::SpectrumImage>::New(),
                                mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object")));

  auto processableNodes = dataStorage->GetSubset(predicate);


  for (auto node : processableNodes->CastToSTLConstContainer())
  {
    bool visibility;
    if (node->GetBoolProperty("visible", visibility))
      if (visibility)
        nodesToProcess->push_back(node);
  }

  return nodesToProcess;
}