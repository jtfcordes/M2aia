/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDataNodeConvertToRGBImageAction.h"

#include <mitkImageAccessByItk.h>
#include <mitkImageCast.h>
#include <mitkLookupTable.h>
#include <mitkLookupTableProperty.h>

#include <itkRescaleIntensityImageFilter.h>
#include <itkImageIOBase.h>

#include <QString>
#include <QAction>
#include <QMenu>

QmitkDataNodeConvertToRGBImageAction::QmitkDataNodeConvertToRGBImageAction(QWidget* parent, berry::IWorkbenchPartSite::Pointer workbenchPartSite): 
  QAction(parent),
  QmitkAbstractDataNodeAction(workbenchPartSite)
{
  InitializeAction();
}

QmitkDataNodeConvertToRGBImageAction::QmitkDataNodeConvertToRGBImageAction(QWidget* parent, berry::IWorkbenchPartSite* workbenchPartSite): 
  QAction(parent),
  QmitkAbstractDataNodeAction(berry::IWorkbenchPartSite::Pointer(workbenchPartSite))
{
  InitializeAction();
}

void QmitkDataNodeConvertToRGBImageAction::InitializeAction()
{
  setText(tr("Convert to RGB Image"));
  setToolTip(tr("It the image is a ."));

  
    connect(this, &QAction::triggered, this, [this]() {
      auto selectedNodes = this->GetSelectedNodes();
      for (auto referenceNode : selectedNodes)
      {
        if (referenceNode.IsNotNull())
        {
          if (mitk::Image::Pointer image = dynamic_cast<mitk::Image *>(referenceNode->GetData()))
          {       
            if(image->GetPixelType().GetNumberOfComponents()==3){
              auto newImage = ConvertMitkVectorImageToRGB(image);
              auto dataNodeNew = mitk::DataNode::New();
              dataNodeNew->SetData(newImage);
              dataNodeNew->SetVisibility(true);
              auto dummyName = referenceNode->GetName() + "_RGB";
              dataNodeNew->SetName(dummyName);
              this->m_DataStorage.Lock()->Add(dataNodeNew);
            }
          }
        }
      }
    });
}

