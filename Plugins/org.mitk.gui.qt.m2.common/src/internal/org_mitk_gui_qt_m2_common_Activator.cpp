/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#include "org_mitk_gui_qt_m2_common_Activator.h"

#include "m2BrowserPreferencesPage.h"
#include <QmitkNodeDescriptorManager.h>

#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateProperty.h>

// #include "QmitkDataNodeColorAction.h"
#include "QmitkDataNodePlotColorAction.h"
#include "QmitkDataNodeColorMapAction.h"
#include "QmitkDataNodeTextureInterpolationAction.h"
#include "QmitkDataNodeExportComponentAction.h"
#include "QmitkDataNodeConvertPixelTypeAction.h"
#include "QmitkDataNodeSliceWiseNormalizationAction.h"
#include "QmitkDataNodeReimportImageAction.h"

#include <m2UIUtils.h>
#include <m2IntervalVector.h>
#include <m2ImzMLSpectrumImage.h>
#include <m2SpectrumImageStack.h>

#include <usModuleInitialization.h>

US_INITIALIZE_MODULE



void org_mitk_gui_qt_m2_common_Activator::start(ctkPluginContext *context)
{
  BERRY_REGISTER_EXTENSION_CLASS(m2BrowserPreferencesPage, context)
  // BERRY_REGISTER_EXTENSION_CLASS(QmitkDataNodeExportComponentActionProvider, context)
  // BERRY_REGISTER_EXTENSION_CLASS(QmitkDataNodeColorMapActionProvider, context)
  // BERRY_REGISTER_EXTENSION_CLASS(QmitkDataNodeColorActionProvider, context)
  // BERRY_REGISTER_EXTENSION_CLASS(QmitkDataNodeTextureInterpolationActionProvider, context)
  // BERRY_REGISTER_EXTENSION_CLASS(QmitkDataNodeConvertPixelTypeActionProvider, context)

  auto descriptorManager = QmitkNodeDescriptorManager::GetInstance();

  auto imageDataType = mitk::TNodePredicateDataType<mitk::Image>::New();
  auto a = mitk::TNodePredicateDataType<m2::SpectrumImageStack>::New();
  auto b = mitk::TNodePredicateDataType<m2::ImzMLSpectrumImage>::New();
  auto spectrumImageDataType = mitk::NodePredicateOr::New(a, b);
  auto spectrumImageDescriptorPredicate = mitk::NodePredicateAnd::New(spectrumImageDataType, imageDataType);
  auto spectrumImageStackDescriptorPredicate = mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<m2::SpectrumImageStack>::New(), imageDataType);
  
  descriptorManager->AddDescriptor(new QmitkNodeDescriptor(
    tr("SpectrumImage"), QString(":/QmitkM2aiaCore/SpectrumImage_48.png"), spectrumImageDescriptorPredicate, this));

  descriptorManager->AddDescriptor(new QmitkNodeDescriptor(
    tr("SpectrumImageStack"), QString(":/QmitkM2aiaCore/SpectrumImage_48.png"), spectrumImageStackDescriptorPredicate, this));

  auto desc = descriptorManager->GetDescriptor("Image");
  desc->AddAction(new QmitkDataNodeConvertPixelTypeAction(), false);

  desc = descriptorManager->GetDescriptor("SpectrumImage");
  desc->AddAction(new QmitkDataNodeConvertPixelTypeAction(), false);
  desc->AddAction(new QmitkDataNodePlotColorAction(), false);
  desc->AddAction(new QmitkDataNodeReimportImageAction(), false);

  desc = descriptorManager->GetDescriptor("SpectrumImageStack");
  desc->AddAction(new QmitkDataNodeSliceWiseNormalizationAction(), false);
  

  desc = descriptorManager->GetDescriptor("MultiComponentImage");
  desc->AddAction(new QmitkDataNodeExportComponentAction(), false);

  
  descriptorManager->AddDescriptor(new QmitkNodeDescriptor(
    tr("IntervalVector"), QString(":/QmitkM2aiaCore/Spectrum_48.png"), mitk::TNodePredicateDataType<m2::IntervalVector>::New(), this));
  desc = descriptorManager->GetDescriptor("IntervalVector");
  desc->AddAction(new QmitkDataNodePlotColorAction(), false);
  


}

void org_mitk_gui_qt_m2_common_Activator::stop(ctkPluginContext *)
{
  m2::UIUtils::Instance()->disconnect();
}
