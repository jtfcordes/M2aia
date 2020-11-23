/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#include "org_mitk_m2_core_helper_Activator.h"

#include "m2BrowserPreferencesPage.h"

#include "QmitkNodeDescriptorManager.h"
#include "QmitkStyleManager.h"
#include "mitkNodePredicateDataType.h"




ctkPluginContext* org_mitk_m2_core_helper_Activator::m_Context = nullptr;
org_mitk_m2_core_helper_Activator*
org_mitk_m2_core_helper_Activator::m_Instance = nullptr;

org_mitk_m2_core_helper_Activator::org_mitk_m2_core_helper_Activator()
{
    m_Instance = this;
}

org_mitk_m2_core_helper_Activator::~org_mitk_m2_core_helper_Activator()
{
    m_Instance = nullptr;
}

void org_mitk_m2_core_helper_Activator::start(ctkPluginContext* context)
{
    BERRY_REGISTER_EXTENSION_CLASS(m2BrowserPreferencesPage, context)

        this->m_Context = context;

    QmitkNodeDescriptorManager* manager = QmitkNodeDescriptorManager::GetInstance();

    mitk::NodePredicateDataType::Pointer isMITKRegistrationWrapper =
        mitk::NodePredicateDataType::New("MAPRegistrationWrapper");
    auto desc = new QmitkNodeDescriptor(QObject::tr("MAPRegistrationWrapper"),
      QmitkStyleManager::ThemeIcon(QStringLiteral(":/QmitkMatchPointCore/MAPRegData.svg")), isMITKRegistrationWrapper, manager);

    manager->AddDescriptor(desc);

	// declare
	
}

void org_mitk_m2_core_helper_Activator::stop(ctkPluginContext* context)
{
    Q_UNUSED(context)

        this->m_Context = nullptr;
}

ctkPluginContext* org_mitk_m2_core_helper_Activator::getContext()
{
    return m_Context;
}

org_mitk_m2_core_helper_Activator*
org_mitk_m2_core_helper_Activator::getDefault()
{
    return m_Instance;
}
