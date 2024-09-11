/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/
#include "QmitkDataNodeExportComponentAction.h"

#include <itkVectorIndexSelectionCastImageFilter.h>

#include "mitkRenderingManager.h"
#include <mitkImage.h>
#include <mitkImageAccessByItk.h>
#include <mitkImageReadAccessor.h>
#include <mitkImageWriteAccessor.h>
#include <mitkImageCast.h>
// needed for qApp
#include <qcoreapplication.h>
#include <QWidget>
#include <QAction>

QmitkDataNodeExportComponentAction::QmitkDataNodeExportComponentAction(QWidget *parent, berry::IWorkbenchPartSite::Pointer workbenchPartSite)
  : QAction(parent),
    QmitkAbstractDataNodeAction(workbenchPartSite)
{
  InitializeAction();
}

QmitkDataNodeExportComponentAction::QmitkDataNodeExportComponentAction(QWidget *parent, berry::IWorkbenchPartSite *workbenchPartSite)
  : QAction(parent), 
    QmitkAbstractDataNodeAction(workbenchPartSite)
{
  InitializeAction();
}

void QmitkDataNodeExportComponentAction::InitializeAction()
{
  setText(tr("Export component(s)"));
  setToolTip(tr("Generate individual images from a multi component image."
                "\n e.g. a RGB imag is split into R, G, and B."
                "\n e.g. extract the selected component of an image with number of components >3."));
  connect(this, &QAction::triggered, this, &QmitkDataNodeExportComponentAction::OnActionTriggered);
}

mitk::Image::Pointer QmitkDataNodeExportComponentAction::ExportComponentImage(mitk::Image *img,
                                                                              unsigned int component)
{
  const auto inputPixelType = img->GetPixelType();
  // const auto nComponents = inputPixelType.GetNumberOfComponents();
  // const auto &sVector = inputPixelType.GetSize();
  // const auto sElement = sVector / nComponents;
  // const auto nPixels =
  //   std::accumulate(img->GetDimensions(), img->GetDimensions() + img->GetDimension(), 1, std::multiplies<>());
  // const auto inGeometry = img->GetGeometry();

  mitk::Image::Pointer output;
  
  AccessVectorPixelTypeByItk(img, ([&](auto itkImage){
    using SourceImageType = typename std::remove_pointer<decltype(itkImage)>::type;
    using ScalarImageType = itk::Image<typename SourceImageType::ValueType::ComponentType, SourceImageType::ImageDimension>;
    using IndexSelectionType = itk::VectorIndexSelectionCastImageFilter<SourceImageType, ScalarImageType>;
    auto indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(component);
    indexSelectionFilter->SetInput(itkImage);
    indexSelectionFilter->Update();
    mitk::CastToMitkImage(indexSelectionFilter->GetOutput(), output); 
  }));

  
  // auto newPixelType = mitk::MakePixelType(inputPixelType, 1);
  // output->GetGeometry()->SetSpacing(inGeometry->GetSpacing());
  // output->GetGeometry()->SetOrigin(inGeometry->GetOrigin());
  // output->GetGeometry()->SetIndexToWorldTransformWithoutChangingSpacing(inGeometry->GetIndexToWorldTransform());

  

  // mitk::ImageReadAccessor acc(img);
  // mitk::ImageWriteAccessor occ(output);
  // char *oData = (char *)occ.GetData();
  // for (char *iData = (char *)acc.GetData(); iData != (char *)acc.GetData() + (sVector * nPixels); iData += sVector)
  // {
  //   std::copy(iData + (component * sElement), iData + ((component + 1) * sElement), oData);
  //   oData += sElement;
  // }

  return output;
}

void QmitkDataNodeExportComponentAction::OnActionTriggered(bool)
{
  auto selectedNodes = this->GetSelectedNodes();
  for (const auto &referenceNode : selectedNodes)
  {
    if (referenceNode.IsNotNull())
    {
      mitk::Image *img = dynamic_cast<mitk::Image *>(referenceNode->GetData());
      if (img)
      {
        const auto inputPixelType = img->GetPixelType();
        const auto nComponents = inputPixelType.GetNumberOfComponents();
        const auto ioPixelType = inputPixelType.GetPixelType();
        if (nComponents == 1)
          continue;

        if (!m_DataStorage)
          continue;

        if (ioPixelType == mitk::PixelType::ItkIOPixelType::RGB)
        {
          const std::string channelNames = "RGB";
          for (unsigned int i = 0; i < channelNames.size(); ++i)
          {
            auto output = ExportComponentImage(img, i);
            auto newNode = mitk::DataNode::New();
            newNode->SetData(output);
            newNode->SetName(referenceNode->GetName() + "_" + channelNames.at(i));
            m_DataStorage.Lock()->Add(newNode, referenceNode);
          }
        }

        else if (ioPixelType == mitk::PixelType::ItkIOPixelType::RGBA)
        {
          const std::string channelNames = "RGBA";
          for (unsigned int i = 0; i < channelNames.size(); ++i)
          {
            auto output = ExportComponentImage(img, i);
            auto newNode = mitk::DataNode::New();
            newNode->SetData(output);
            newNode->SetName(referenceNode->GetName() + "_" + channelNames.at(i));
            m_DataStorage.Lock()->Add(newNode, referenceNode);
          }
        }

        else if (nComponents <= 6)
        {
          for (unsigned int i = 0; i < nComponents; ++i)
          {
            auto output = ExportComponentImage(img, i);
            auto newNode = mitk::DataNode::New();
            newNode->SetData(output);
            newNode->SetName(referenceNode->GetName() + "_" + std::to_string(i));
            m_DataStorage.Lock()->Add(newNode, referenceNode);
          }
        }
        else
        {
          mitk::IntProperty *componentProperty =
            dynamic_cast<mitk::IntProperty *>(referenceNode->GetProperty("Image.Displayed Component"));
          const auto currentComponent = componentProperty->GetValue();
          auto output = ExportComponentImage(img, currentComponent);
          auto newNode = mitk::DataNode::New();
          newNode->SetData(output);
          newNode->SetName(referenceNode->GetName() + "_" + std::to_string(currentComponent));
          m_DataStorage.Lock()->Add(newNode, referenceNode);
        }
      }
    }
  }
}