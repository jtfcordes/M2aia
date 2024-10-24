/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt for details.

===================================================================*/

#include "m2Spectrum.h"

#include <iostream>
#include <memory>

// Qt includes
// #include <QApplication>
#include <QColor>
#include <QComboBox>
#include <QFileDialog>
#include <QGraphicsSimpleTextItem>
#include <QLabel>
#include <QLineSeries>
#include <QMenu>
#include <QPushButton>
#include <QShortcut>
#include <QValueAxis>
#include <QWidgetAction>
#include <QXYSeries>

// MITK includes
#include <QmitkIOUtil.h>
#include <mitkColorProperty.h>
#include <mitkCoreServices.h>
#include <mitkIPreferences.h>
#include <mitkIPreferencesService.h>
#include <mitkLookupTableProperty.h>
#include <mitkStatusBar.h>

// M2aia includes
#include <m2ImzMLSpectrumImage.h>
#include <m2SpectrumImage.h>
#include <m2UIUtils.h>
#include <signal/m2PeakDetection.h>

// internal includes
#include "m2SeriesDataProvider.h"

const std::string m2Spectrum::VIEW_ID = "org.mitk.views.m2.spectrum";

void m2Spectrum::DrawSelectedArea()
{
  if (!m_xAxis)
    return;
  if (!m_yAxis)
    return;

  auto chart = m_Controls.chartView->chart();
  if (m_SelectedArea[1] == nullptr)
  {
    m_SelectedArea[1] = new QLineSeries();
    chart->addSeries(m_SelectedArea[1]);
    m_SelectedArea[1]->attachAxis(m_xAxis);
    m_SelectedArea[1]->attachAxis(m_yAxis);

    const auto markers = chart->legend()->markers();
    for (auto *marker : markers)
    {
      QObject::connect(marker, SIGNAL(clicked()), this, SLOT(OnLegnedHandleMarkerClicked()), Qt::UniqueConnection);
    }

    QPen pen(QColor(0, 155, 0, 120));
    pen.setWidthF(2);
    pen.setMiterLimit(0);
    m_SelectedArea[1]->setPen(pen);
    m_SelectedArea[1]->setName("Selection");

    m_SelectedArea[2] = new QLineSeries();
    chart->addSeries(m_SelectedArea[2]);
    m_SelectedArea[2]->attachAxis(m_xAxis);
    m_SelectedArea[2]->attachAxis(m_yAxis);

    QPen borderpen(QColor(0, 155, 0, 255));
    borderpen.setWidthF(0.75);
    m_SelectedArea[2]->setPen(borderpen);

    chart->legend()->markers(m_SelectedArea[2])[0]->setVisible(false);

    m_SelectedArea[0] = new QLineSeries();
    chart->addSeries(m_SelectedArea[0]);

    m_SelectedArea[0]->attachAxis(m_xAxis);
    m_SelectedArea[0]->attachAxis(m_yAxis);

    m_SelectedArea[0]->setPen(borderpen);
    chart->legend()->markers(m_SelectedArea[0])[0]->setVisible(false);
  }

  UpdateSelectedArea();
}

void m2Spectrum::UpdateSelectedArea()
{
  if (m_Chart->series().empty())
    return;

  if (m_SelectedArea[0] && m_SelectedArea[1] && m_SelectedArea[2])
  {
    m_SelectedArea[1]->clear();
    *m_SelectedArea[1] << QPointF(m_SelectedAreaX[0], m_yAxis->min()) << QPointF(m_SelectedAreaX[1], m_yAxis->min());

    m_SelectedArea[0]->clear();
    *m_SelectedArea[0] << QPointF(m_SelectedAreaX[0], m_yAxis->min()) << QPointF(m_SelectedAreaX[0], m_yAxis->max());

    m_SelectedArea[2]->clear();
    *m_SelectedArea[2] << QPointF(m_SelectedAreaX[1], m_yAxis->min()) << QPointF(m_SelectedAreaX[1], m_yAxis->max());

    m_Controls.chartView->repaint();

    bool showSamplingPoints = m_M2aiaPreferences->GetBool("m2aia.view.spectrum.showSamplingPoints", false);
    
    if (showSamplingPoints)
    {
      using namespace std;
      double size = 5;
      for (auto kv : m_DataProvider)
      {
        try
        {
          for (auto item : m_NodeRelatedGraphicItems[kv.first]->childItems())
          {
            m_NodeRelatedGraphicItems[kv.first]->removeFromGroup(item);
            delete item;
          }

          if (kv.first->IsOn("helper object", nullptr)) // if it is On it is not shown in the data manager
            continue;

          if (!kv.first->IsVisible(nullptr)) // if it is shown in the DataManger make sure it is set to visible
            continue;

          auto p = kv.second;
          auto lowerX = lower_bound(begin(p->xs()), end(p->xs()), std::max(m_SelectedAreaX[0], m_xAxis->min()));
          auto upperX = upper_bound(lowerX, end(p->xs()), std::min(m_SelectedAreaX[1], m_xAxis->max()));
          auto lowerY = begin(p->ys()) + distance(begin(p->xs()), lowerX);

          for (; lowerX < upperX; ++lowerX, ++lowerY)
          {
            auto rect = new QGraphicsRectItem();
            if (auto prop = kv.first->GetProperty("spectrum.marker.size"))
            {
              if (auto sizeProp = dynamic_cast<mitk::IntProperty *>(prop))
              {
                size = sizeProp->GetValue();
              }
            }
            auto itemPos = m_Chart->mapToPosition(QPointF(*lowerX, *lowerY));
            rect->setRect(itemPos.x() - size / 2.0, itemPos.y() - size / 2.0, size, size);

            if (auto prop = kv.first->GetProperty("spectrum.marker.color"))
            {
              if (auto colorProp = dynamic_cast<mitk::ColorProperty *>(prop))
              {
                auto mc = colorProp->GetColor();
                QColor c;
                c.setRgbF(mc.GetRed(), mc.GetGreen(), mc.GetBlue());
                rect->setBrush(QBrush(c));
                rect->setPen(QPen(c));
              }
            }

            rect->setZValue(20);
            // m_NodeRealtedGraphicItems[kv.first].push_back(rect);
            // QGraphicsItemGroup
            m_NodeRelatedGraphicItems[kv.first]->addToGroup(rect);
          }
        }
        catch (std::exception &e)
        {
          MITK_ERROR << e.what();
        }
      }
    }
  }
}

void m2Spectrum::OnMousePress(QPointF pos, qreal mz, qreal intValue, Qt::MouseButton button, Qt::KeyboardModifiers mod)
{
  Q_UNUSED(pos)
  if (mod & Qt::AltModifier && !m_RangeSelectionStarted)
  {
    m_Controls.chartView->setRubberBand(QChartView::RubberBand::HorizontalRubberBand);
    m_SelectedAreaX[0] = mz;
    m_RangeSelectionStarted = true;
  }

  if (button == Qt::MouseButton::MiddleButton)
  {
    auto c = m_Controls.chartView->chart();
    const auto *xAx = dynamic_cast<QValueAxis *>(c->axes(Qt::Horizontal).front());
    const auto *yAx = dynamic_cast<QValueAxis *>(c->axes(Qt::Vertical).front());
    if (mz > xAx->min() && mz < xAx->max() && intValue > yAx->min() && intValue < yAx->max())
    {
      m_MouseDragCenterPos = mz;
      m_MouseDragLowerDelta = mz - xAx->min();
      m_MouseDragUpperDelta = xAx->max() - mz;
      m_DraggingActive = true;
    }
  }

  UpdateSelectedArea();
  // MITK_INFO("m2Spectrum::OnMousePress") << "Mouse Pressed";
}

void m2Spectrum::OnMassRangeChanged(qreal x, qreal tol)
{
  SetSelectedAreaStartX(x - tol);
  SetSelectedAreaEndX(x + tol);
  DrawSelectedArea();
}

// void m2Spectrum::OnSpectrumArtifactChanged(const mitk::DataNode * node, m2::SpectrumType type){
//   if(dynamic_cast<m2::SpectrumImage *>(node->GetData())){
//     if(type == m2::SpectrumType::None){
//       if(m_DataProvider[node]->Exists("Peaks")){
//         m_DataProvider[node]->UpdateGroup(pea)
//       }
//     }
//   }
// }

void m2Spectrum::OnDataModified(const itk::Object *caller, const itk::EventObject &)
{
  if (auto node = dynamic_cast<const mitk::DataNode *>(caller))
  {
    if (dynamic_cast<m2::IntervalVector *>(node->GetData()))
    {
      auto provider = m_DataProvider[node];
      provider->Update();
      provider->UpdateBoundaries(m_LocalMinimumX, m_LocalMaximumX);
      UpdateCurrentMinMaxY();
      UpdateGlobalMinMaxValues();
      DrawSelectedArea();
      AutoZoomUseLocalExtremaY();
    }
  }
}

void m2Spectrum::OnPropertyListChanged(const itk::Object *caller, const itk::EventObject &)
{
  if (auto node = dynamic_cast<const mitk::DataNode *>(caller))
  {
    if (dynamic_cast<m2::IntervalVector *>(node->GetData()))
    {
      if (auto provider = m_DataProvider[node])
      {
        provider->GetSeries()->setVisible(node->IsVisible(nullptr));

        if (node->IsOn("helper object", nullptr, false))
          provider->GetSeries()->setVisible(false);

        if (auto genericProperty = node->GetProperty("spectrum.plot.color"))
        {
          if (auto colorProperty = dynamic_cast<mitk::ColorProperty *>(genericProperty))
          {
            auto c = colorProperty->GetColor();
            float a = 0;
            // MITK_INFO << "opacity " << a;
            node->GetOpacity(a, nullptr);
            provider->SetColor(c.GetRed(), c.GetGreen(), c.GetBlue(), 1 - a);
          }
        }

        // UpdateCurrentMinMaxY();
        // UpdateGlobalMinMaxValues();

        // UpdateSelectedArea();

        if (auto genericProperty = node->GetProperty("spectrum.marker.color"))
        {
          if (auto colorProperty = dynamic_cast<mitk::ColorProperty *>(genericProperty))
          {
            auto mc = colorProperty->GetColor();

            for (auto item : m_NodeRelatedGraphicItems[node]->childItems())
            {
              QColor c;
              c.setRgbF(mc.GetRed(), mc.GetGreen(), mc.GetBlue());
              auto rect = dynamic_cast<QGraphicsRectItem *>(item);
              rect->setBrush(QBrush(c));
              rect->setPen(QPen(c));
            }
          }
        }

        if (auto prop = node->GetProperty("spectrum.marker.size"))
        {
          if (auto sizeProp = dynamic_cast<mitk::IntProperty *>(prop))
          {
            auto size = sizeProp->GetValue();
            for (auto item : m_NodeRelatedGraphicItems[node]->childItems())
            {
              auto rect = dynamic_cast<QGraphicsRectItem *>(item);
              auto itemPos = rect->rect().center();
              rect->setRect(itemPos.x() - size / 2.0, itemPos.y() - size / 2.0, size, size);
            }
          }
        }
      }
    }
  }
}

void m2Spectrum::NodeAdded(const mitk::DataNode *node)
{
  if (!node)
    return;

  // MITK_INFO << "NodeAdded: " << node->GetName();

  m_Chart = m_Controls.chartView->chart();

  if (auto intervals = dynamic_cast<m2::IntervalVector *>(node->GetData()))
  {
    if (m_DataProvider.empty())
    {
      m_GlobalMinimumX = m_LocalMinimumX = intervals->GetXMean().front();
      m_GlobalMaximumX = m_LocalMaximumX = intervals->GetXMean().back();
    }

    bool isVisible;
    isVisible = node->IsVisible(nullptr);

    if (node->IsOn("helper object", nullptr, false))
      isVisible = false;

    auto provider = std::make_shared<m2::SeriesDataProvider>();
    provider->Initialize(intervals);
    provider->GetSeries()->setName(node->GetName().c_str());
    provider->GetSeries()->setVisible(isVisible);

    m_DataProvider[node] = provider;
    m_NodeRelatedGraphicItems[node] = new QGraphicsItemGroup();
    m_Chart->scene()->addItem(m_NodeRelatedGraphicItems[node]);
    m_Chart->addSeries(provider->GetSeries());

    m_Chart->createDefaultAxes();
    m_xAxis = static_cast<QValueAxis *>(m_Chart->axes(Qt::Horizontal).front());
    m_yAxis = static_cast<QValueAxis *>(m_Chart->axes(Qt::Vertical).front());
    QObject::connect(m_xAxis, SIGNAL(rangeChanged(qreal, qreal)), this, SLOT(OnRangeChangedAxisX(qreal, qreal)));
    QObject::connect(m_yAxis, SIGNAL(rangeChanged(qreal, qreal)), this, SLOT(OnRangeChangedAxisY(qreal, qreal)));

    auto onPropertyListModifiedCommand = itk::NodeMemberCommand<m2Spectrum>::New();
    onPropertyListModifiedCommand->SetNode(node);
    onPropertyListModifiedCommand->SetCallbackFunction(this, &m2Spectrum::OnPropertyListChanged);
    node->GetPropertyList()->AddObserver(itk::ModifiedEvent(), onPropertyListModifiedCommand);

    auto onDataModifiedCommand = itk::MemberCommand<m2Spectrum>::New();
    onDataModifiedCommand->SetCallbackFunction(this, &m2Spectrum::OnDataModified);
    node->AddObserver(itk::ModifiedEvent(), onDataModifiedCommand);

    OnAxisXTicksChanged(m_AxisTicks[0]);
    OnAxisYTicksChanged(m_AxisTicks[1]);

    UpdateAllSeries();

    UpdateGlobalMinMaxValues();
    UpdateCurrentMinMaxY();

    // OnResetView();
    // AutoZoomUseLocalExtremaY();
    if (m_xAxis)
      m_xAxis->setRange(m_LocalMinimumX, m_LocalMaximumX);
    if (m_yAxis)
      m_yAxis->setRange(m_LocalMinimumY, m_LocalMaximumY * 1.1);
    m_Controls.chartView->repaint();

    // MITK_INFO << "SpectrumView: " << node->GetName();
    // m2::DefaultNodeProperties(node, false);
    node->GetPropertyList()->Modified();
    node->Modified();

    // UpdateTitles(node);
  }
}

void m2Spectrum::OnMouseMove(
  QPointF pos, qreal mz, qreal intValue, Qt::MouseButton /*button*/, Qt::KeyboardModifiers /*mod*/)
{
  const auto chart = m_Controls.chartView->chart();
  const auto chartPos = chart->mapToPosition(QPoint(mz, intValue));

  if (m_DraggingActive)
  {
    const auto xPos = chartPos.x();
    const auto xLastPos = chart->mapToPosition(QPoint(m_MouseDragCenterPos, 0)).x();
    const auto delta = xLastPos - xPos;
    chart->scroll(delta, 0);
  }

  // inside draw area
  if (!chart->series().empty())
  {
    if (((intValue > m_yAxis->min()) && (intValue < m_yAxis->max()) && (mz > m_xAxis->min()) && (mz < m_xAxis->max())))
    {
      m_Controls.chartView->repaint();
      m_Crosshair->setText(QString("%2 @ %1").arg(mz).arg(intValue, 0, 'e', 3));
      m_Crosshair->setPos(pos + QPoint(10, -30));

      m_Crosshair->setZValue(11);
      // m_Crosshair->updateGeometry();
      m_Crosshair->show();
    }
    else
    {
      m_Crosshair->hide();
      m_Controls.chartView->repaint();
    }
  }

  // UpdateSelectedArea();

  // MITK_INFO("m2Spectrum::OnMouseMove") << "Mouse Pressed";
}

void m2Spectrum::OnMouseRelease(
  QPointF pos, qreal mz, qreal intValue, Qt::MouseButton button, Qt::KeyboardModifiers mod)
{
  Q_UNUSED(pos)
  Q_UNUSED(intValue)
  if (mod & Qt::AltModifier && m_RangeSelectionStarted)
  {
    m_SelectedAreaX[1] = mz;
  }
  if (m_RangeSelectionStarted)
  {
    m_Controls.chartView->setRubberBand(QChartView::RubberBand::NoRubberBand);
    DrawSelectedArea();
    const auto mz = (m_SelectedAreaX[1] + m_SelectedAreaX[0]) * 0.5;
    const auto tol = std::abs(m_SelectedAreaX[1] - m_SelectedAreaX[0]) * 0.5;

    emit m2::UIUtils::Instance()->UpdateImage(mz, tol);
  }
  m_RangeSelectionStarted = false;
  m_Controls.chartView->setRubberBand(QChartView::RubberBand::NoRubberBand);

  if (button == Qt::MouseButton::MiddleButton)
  {
    m_DraggingActive = false;
  }
}

void m2Spectrum::OnMouseDoubleClick(
  QPointF pos, qreal xVal, qreal yVal, Qt::MouseButton /*button*/, Qt::KeyboardModifiers /*mod*/)
{
  Q_UNUSED(pos)

  if (m_DataProvider.empty())
    return;

  // reset axis by double click below/left of the chart plot area
  if (m_xAxis && m_yAxis)
  {
    if (yVal < m_yAxis->min()) // rest x
      m_xAxis->setRange(m_GlobalMinimumX, m_GlobalMaximumX);
    else if (xVal < m_xAxis->min())
      m_yAxis->setRange(0, m_GlobalMaximumY * 1.1); // reset y
    else if (xVal > m_xAxis->max())
    {
    } // do nothing
    else
      emit m2::UIUtils::Instance()->UpdateImage(xVal, -1);
  }
}

void m2Spectrum::OnMouseWheel(QPointF pos, qreal x, qreal y, int angle, Qt::KeyboardModifiers mod)
{
  Q_UNUSED(pos)

  if (m_DataProvider.empty())
  {
    return;
  }

  const auto modifiedMin = m_yAxis->min();
  const auto modifiedMax = m_yAxis->max() * 1.1;
  bool bothAxes = ((y > modifiedMin) && (y < modifiedMax) && (x > m_xAxis->min()) && (x < m_xAxis->max()));

  if (angle != 0)
  {
    double zoomFactor = 0.1;
    if (mod & Qt::KeyboardModifier::ControlModifier)
    {
      zoomFactor = 0.35;
    }

    if ((y < modifiedMin) || bothAxes)
    {
      const int zoomDirection = angle / std::abs(angle);
      const auto d0 = std::abs(x - m_xAxis->min()) * zoomFactor * zoomDirection;
      const auto d1 = std::abs(x - m_xAxis->max()) * zoomFactor * zoomDirection;
      if (m_xAxis)
        m_xAxis->setRange(m_xAxis->min() + d0, m_xAxis->max() - d1);
    }

    if ((x < m_xAxis->min()) || (bothAxes && (mod == Qt::ShiftModifier)))
    {
      const int zoomDirection = angle / std::abs(angle);
      const auto d1 = m_yAxis->max() * (1 - zoomDirection * zoomFactor);
      if (m_yAxis)
        m_yAxis->setMax(d1);
    }
  }

  UpdateSelectedArea();
}

void m2Spectrum::UpdateCurrentMinMaxY()
{
  if (m_DataProvider.empty())
    return;

  m_LocalMinimumY = std::numeric_limits<double>::max();
  m_LocalMaximumY = std::numeric_limits<double>::min();

  for (auto &kv : m_DataProvider)
  {
    if (!kv.first->IsVisible(nullptr))
      continue;

    const auto points = kv.second->GetSeries()->points();
    if (points.empty())
      continue;

    const auto minmax =
      std::minmax_element(points.begin(), points.end(), [](const auto &a, const auto &b) { return a.y() < b.y(); });

    m_LocalMinimumY = std::min(m_LocalMinimumY, minmax.first->y());
    m_LocalMaximumY = std::max(m_LocalMaximumY, minmax.second->y());
  }

  // MITK_INFO("UpdateLocalMinMaxValues") << minX << " " << maxX;
  // MITK_INFO("UpdateLocalMinMaxValues") << "New local minmax Y: " << m_LocalMinimumY << " " << m_LocalMaximumY;
}

void m2Spectrum::UpdateGlobalMinMaxValues()
{
  if (m_DataProvider.empty())
    return;

  m_GlobalMinimumY = std::numeric_limits<double>::max();
  m_GlobalMaximumY = std::numeric_limits<double>::min();

  for (auto kv : m_DataProvider)
  {
    if (!kv.first->IsVisible(nullptr))
      continue;

    const auto points = kv.second->GetSeries()->points();
    if (points.empty())
      continue;

    if (kv.second->xs().empty())
      continue;

    m_GlobalMinimumX = std::min(m_GlobalMinimumX, kv.second->xs().front());
    m_GlobalMaximumX = std::max(m_GlobalMaximumX, kv.second->xs().back());

    const auto cmp = [](const auto &a, const auto &b) { return a.y() < b.y(); };
    const auto minmaxY = std::minmax_element(std::begin(points), std::end(points), cmp);
    m_GlobalMinimumY = std::min(std::max(0.0, m_GlobalMinimumY), std::max(0.0, minmaxY.first->y()));
    m_GlobalMaximumY = std::max(m_GlobalMaximumY, minmaxY.second->y());
  }

  auto useMinIntensity = m_M2aiaPreferences->GetBool("m2aia.view.spectrum.useMinIntensity", true);
  if (!useMinIntensity)
    m_GlobalMinimumY = 0;

  // m_LocalMinimumY = m_GlobalMinimumY;
  // m_LocalMaximumY = m_GlobalMaximumY;

  // MITK_INFO("UpdateGlobalMinMaxValues") << "New global minmax X: " << m_GlobalMinimumX << " " << m_GlobalMaximumX;
  // MITK_INFO("UpdateGlobalMinMaxValues") << "New global minmax Y: " << m_GlobalMinimumY << " " << m_GlobalMaximumY;
  // MITK_INFO("UpdateGlobalMinMaxValues") << "New local minmax Y: " << m_LocalMinimumY << " " << m_LocalMaximumY;
}

void m2Spectrum::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);

  QShortcut *shortcutLeft = new QShortcut(QKeySequence(Qt::Key_Left), parent);
  QShortcut *shortcutRight = new QShortcut(QKeySequence(Qt::Key_Right), parent);
  QShortcut *shortcutShiftLeft = new QShortcut(QKeySequence(Qt::SHIFT | Qt::Key_Left), parent);
  QShortcut *shortcutShiftRight = new QShortcut(QKeySequence(Qt::SHIFT | Qt::Key_Right), parent);
  QShortcut *shortcutUp = new QShortcut(QKeySequence(Qt::Key_Up), parent);
  QShortcut *shortcutDown = new QShortcut(QKeySequence(Qt::Key_Down), parent);

  auto UIUtilsObject = m2::UIUtils::Instance();
  connect(UIUtilsObject, SIGNAL(RangeChanged(qreal, qreal)), this, SLOT(OnMassRangeChanged(qreal, qreal)));

  connect(shortcutRight, SIGNAL(activated()), UIUtilsObject, SIGNAL(NextImage()));
  connect(shortcutLeft, SIGNAL(activated()), UIUtilsObject, SIGNAL(PreviousImage()));
  connect(shortcutShiftRight, SIGNAL(activated()), UIUtilsObject, SIGNAL(NextPeakImage()));
  connect(shortcutShiftLeft, SIGNAL(activated()), UIUtilsObject, SIGNAL(PreviousPeakImage()));
  connect(shortcutUp, SIGNAL(activated()), UIUtilsObject, SIGNAL(IncreaseTolerance()));
  connect(shortcutDown, SIGNAL(activated()), UIUtilsObject, SIGNAL(DecreaseTolerance()));

  CreateQChartView();
  CreateQChartViewMenu();

  auto *preferencesService = mitk::CoreServices::GetPreferencesService();
	m_M2aiaPreferences = preferencesService->GetSystemPreferences();

  m_Controls.chartView->chart()->legend()->setVisible(false);
  
  m_Controls.comboBox->addItem("BlueIcy", QChart::ChartThemeBlueIcy);
  m_Controls.comboBox->addItem("BlueCerulean", QChart::ChartThemeBlueCerulean);
  m_Controls.comboBox->addItem("BlueNcs", QChart::ChartThemeBlueNcs);
  m_Controls.comboBox->addItem("BrownSand", QChart::ChartThemeBrownSand);
  m_Controls.comboBox->addItem("Dark", QChart::ChartThemeDark);
  m_Controls.comboBox->addItem("HighContrast", QChart::ChartThemeHighContrast);
  m_Controls.comboBox->addItem("Light", QChart::ChartThemeLight);
  m_Controls.comboBox->addItem("Qt", QChart::ChartThemeQt);

  connect(m_Controls.comboBox,
          qOverload<int>(&QComboBox::highlighted),
          this,
          [&](int i)
          {
            auto theme = static_cast<QChart::ChartTheme>(m_Controls.comboBox->itemData(i).toInt());          

            if (!m_Chart->series().empty())
            {
              QPen pen;
              if (auto xySeriesPre = dynamic_cast<QXYSeries *>(m_Chart->series().front()))
              {
                pen = xySeriesPre->pen();
               // MITK_INFO << xySeriesPre->pen().widthF();
              }
              
              m_Chart->setTheme(theme);

              foreach (QAbstractSeries *series, m_Chart->series())
              {
                if (auto xySeriesPost = dynamic_cast<QXYSeries *>(series))
                {
                  auto pen = xySeriesPost->pen();
                  pen.setWidthF(0.75);
                  xySeriesPost->setPen(pen);
                }
              }
            }
           });

  connect(m_Controls.comboBox,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          [&](int i)
          {
            auto theme = static_cast<QChart::ChartTheme>(m_Controls.comboBox->itemData(i).toInt());
            m_Chart->setTheme(theme);
           });

  m_Controls.comboBox->setCurrentIndex(4);

  connect(m_Controls.saveButton,
          &QPushButton::clicked,
          this,
          [&]()
          {
            auto name = QFileDialog::getSaveFileName(
              m_Controls.chartView, tr("Save File"), "/home/jana/untitled.png", tr("Images (*.png *.xpm *.jpg)"));
            auto r = m_Controls.chartView->rect();
            r.setX(r.x() + 10);
            r.setY(r.y() + 10);

            r.setWidth(r.width() - 10);
            r.setHeight(r.height() - 10);

            m_Controls.chartView->grab(r).save(name);
           });
}

void m2Spectrum::CreateQChartViewMenu()
{
  m_Menu = new QMenu(m_Controls.chartView);

  // m_FocusMenu = new QMenu(m_Menu);
  // m_FocusMenu->setTitle("Focus on ...");
  // m_Menu->addMenu(m_FocusMenu);

  // m_SpectrumSkyline = new QAction("Skyline spectrum", m_Controls.chartView);
  // m_SpectrumSkyline->setCheckable(true);

  // m_SpectrumMean = new QAction("Mean spectrum", m_Controls.chartView);
  // m_SpectrumMean->setCheckable(true);
  // m_SpectrumMean->setChecked(true);

  // m_ShowLegend = new QAction("Show legend", m_Controls.chartView);
  // m_ShowLegend->setCheckable(false);

  // m_SpectrumSelectionGroup = new QActionGroup(m_Controls.chartView);
  // m_SpectrumSelectionGroup->setExclusive(true);
  // m_SpectrumSelectionGroup->addAction(m_SpectrumSkyline);
  // m_SpectrumSelectionGroup->addAction(m_SpectrumMean);

  // connect(m_SpectrumSkyline,
  //         &QAction::toggled,
  //         this,
  //         [&](bool v)
  //         {
  //           if (!v)
  //             return;

  //           for (auto &kv : this->m_DataProvider)
  //             kv.second->SetActiveProviderTo("Max");

  //           UpdateGlobalMinMaxValues();
  //           UpdateCurrentMinMaxY();
  //           emit m_xAxis->rangeChanged(m_xAxis->min(), m_xAxis->max());
  //         });

  // connect(m_SpectrumMean,
  //         &QAction::toggled,
  //         this,
  //         [&](bool v)
  //         {
  //           if (!v)
  //             return;

  //           for (auto &kv : this->m_DataProvider)
  //             kv.second->SetActiveProviderTo("Mean");

  //           UpdateGlobalMinMaxValues();
  //           UpdateCurrentMinMaxY();
  //           emit m_xAxis->rangeChanged(m_xAxis->min(), m_xAxis->max());

  // m_Menu->addAction(m_ShowLegend);
  m_ShowAxesTitles = new QAction("Show axes title", m_Controls.chartView);
  m_ShowAxesTitles->setCheckable(true);
  m_Menu->addAction(m_ShowAxesTitles);
  m_ShowAxesTitles->setChecked(false);

  connect(m_ShowAxesTitles,
          &QAction::toggled,
          this,
          [this](bool s)
          {
            m_yAxis->setTitleText("Intensity");
            m_xAxis->setTitleText("m/z");
            m_yAxis->setTitleVisible(s);
            m_xAxis->setTitleVisible(s);
          });

  // connect(
  //   m_ShowLegend, &QAction::toggled, this, [this](bool s) { m_Controls.chartView->chart()->legend()->setVisible(s);
  //   });

  // m_ShowLegend->setChecked(true);

  m_Menu->addSection("Overview spectrum type");
  // m_Menu->addAction(m_SpectrumSkyline);
  // m_Menu->addAction(m_SpectrumMean);

  // auto sec = m_Menu->addSection("Tick control");
  m_TickCountX = new QSlider(m_Controls.chartView);
  m_TickCountX->setMinimum(3);
  m_TickCountX->setMaximum(50);
  m_TickCountX->setValue(9);
  m_TickCountX->setOrientation(Qt::Horizontal);

  m_TickCountY = new QSlider(m_Controls.chartView);
  m_TickCountY->setMinimum(3);
  m_TickCountY->setMaximum(10);
  m_TickCountY->setValue(4);
  m_TickCountY->setOrientation(Qt::Horizontal);

  auto wActionX = new QWidgetAction(m_Controls.chartView);
  wActionX->setDefaultWidget(m_TickCountX);
  connect(m_TickCountX, &QSlider::valueChanged, this, &m2Spectrum::OnAxisXTicksChanged);

  auto wActionXLabel = new QWidgetAction(m_Controls.chartView);
  wActionXLabel->setDefaultWidget(new QLabel("X axis ticks", m_Controls.chartView));

  auto wActionY = new QWidgetAction(m_Controls.chartView);
  wActionY->setDefaultWidget(m_TickCountY);
  connect(m_TickCountY, &QSlider::valueChanged, this, &m2Spectrum::OnAxisYTicksChanged);
  auto wActionYLabel = new QWidgetAction(m_Controls.chartView);
  wActionYLabel->setDefaultWidget(new QLabel("Y axis ticks", m_Controls.chartView));

  wActionY->setDefaultWidget(m_TickCountY);

  m_Menu->addAction(wActionXLabel);
  m_Menu->addAction(wActionX);

  m_Menu->addAction(wActionYLabel);
  m_Menu->addAction(wActionY);

  m_Controls.chartView->setContextMenuPolicy(Qt::CustomContextMenu);
  connect(m_Controls.chartView,
          &QChartView::customContextMenuRequested,
          this,
          [&](const QPoint &pos) { m_Menu->exec(m_Controls.chartView->viewport()->mapToGlobal(pos)); });
}

void m2Spectrum::CreateQChartView()
{
  auto chart = new QChart();
  m_Controls.chartView->setChart(chart);
  m_Chart = chart;

  chart->legend()->setAlignment(Qt::AlignRight);
  chart->legend()->setShowToolTips(true);
  chart->setAnimationOptions(QChart::NoAnimation);
  chart->setTheme(QChart::ChartThemeBlueIcy);

  m_Crosshair = new QGraphicsSimpleTextItem("", chart);
  auto b = m_Crosshair->brush();
  b.setColor(QColor::fromRgbF(1, 1, 1));
  m_Crosshair->setBrush(b);
  m_Crosshair->hide();

  connect(m_Controls.chartView, &m2::ChartView::mouseDoubleClick, this, &m2Spectrum::OnMouseDoubleClick);
  connect(m_Controls.chartView, &m2::ChartView::mousePress, this, &m2Spectrum::OnMousePress);
  connect(m_Controls.chartView, &m2::ChartView::mouseMove, this, &m2Spectrum::OnMouseMove);
  connect(m_Controls.chartView, &m2::ChartView::mouseRelease, this, &m2Spectrum::OnMouseRelease);
  connect(m_Controls.chartView, &m2::ChartView::mouseWheel, this, &m2Spectrum::OnMouseWheel);
}

void m2Spectrum::AutoZoomUseLocalExtremaY()
{
  if (m_DataProvider.empty())
    return;

  bool useMaxIntensity = m_M2aiaPreferences->GetBool("m2aia.view.spectrum.useMaxIntensity", true);
  bool useMinIntensity = m_M2aiaPreferences->GetBool("m2aia.view.spectrum.useMinIntensity", true);

  if (useMaxIntensity || useMinIntensity)
  {
    if (useMaxIntensity)
    {
      m_yAxis->setMax(m_LocalMaximumY * 1.1);
    } // else use just current zoom, no adaptions.

    if (useMinIntensity)
    {
      m_yAxis->setMin(m_LocalMinimumY);
    }
    else
    {
      m_GlobalMinimumY = 0;
      m_LocalMinimumY = 0;
      m_yAxis->setMin(0);
    }
  }
}

void m2Spectrum::OnSeriesFocused(const mitk::DataNode *node)
{
  double xMin = m_GlobalMaximumX;
  double xMax = m_GlobalMinimumX;

  if (auto provider = m_DataProvider[node])
  {
    xMin = std::min(xMin, provider->xs().front());
    xMax = std::max(xMax, provider->xs().back());
  }
  if (m_xAxis)
    m_xAxis->setRange(xMin, xMax);
  // AutoZoomUseLocalExtremaY();
}

void m2Spectrum::NodeRemoved(const mitk::DataNode *node)
{
  if (dynamic_cast<m2::SpectrumImage *>(node->GetData()) || dynamic_cast<m2::IntervalVector *>(node->GetData()))
  {
    if (m_DataProvider.find(node) == m_DataProvider.end())
      return;

    // remove observers
    const_cast<mitk::DataNode *>(node)->RemoveAllObservers();
    node->GetData()->RemoveAllObservers();

    // remove scene
    m_Chart->scene()->removeItem(m_NodeRelatedGraphicItems[node]);
    delete m_NodeRelatedGraphicItems[node];
    m_NodeRelatedGraphicItems.erase(node);

    // remove provider
    if (auto provider = m_DataProvider[node])
      m_Controls.chartView->chart()->removeSeries(provider->GetSeries());
    m_DataProvider.erase(node);

    // UpdateAxisLabels(node, true);
  if (m_DataProvider.empty())
  {
    if (m_xAxis)
      m_xAxis->setRange(0, 1);
    if (m_yAxis)
    {
      m_yAxis->setRange(0, 1);
      m_yAxis->setLabelFormat("%.6f");
      }
    }
  }

  // OnResetView();
}

void m2Spectrum::SetSeriesVisible(QAbstractSeries *series, bool visibility)
{
  auto markers = m_Chart->legend()->markers();
  std::vector<QLegendMarker *> series_marker;
  for (auto marker : markers)
  {
    if (marker->series() == series)
    {
      // switch (marker->type())
      marker->series()->setVisible(visibility);
      marker->setVisible(true);

      // Dim the marker, if series is not visible
      qreal alpha = 1.0;

      if (!marker->series()->isVisible())
        alpha = 0.5;

      QColor color;
      QBrush brush = marker->labelBrush();
      color = brush.color();
      color.setAlphaF(alpha);
      brush.setColor(color);
      marker->setLabelBrush(brush);

      brush = marker->brush();
      color = brush.color();
      color.setAlphaF(alpha);
      brush.setColor(color);
      marker->setBrush(brush);

      QPen pen = marker->pen();
      color = pen.color();
      color.setAlphaF(alpha);
      pen.setColor(color);
      marker->setPen(pen);
      return;
    }
  }
}

void m2Spectrum::OnLegnedHandleMarkerClicked()
{
  auto *marker = qobject_cast<QLegendMarker *>(sender());
  Q_ASSERT(marker);
  SetSeriesVisible(marker->series(), !marker->series()->isVisible());
}

void m2Spectrum::OnResetView()
{
  if (m_DataProvider.empty())
    return;
}

void m2Spectrum::OnAxisXTicksChanged(int v)
{
  m_AxisTicks[0] = v;
  if (m_xAxis)
  {
    m_xAxis->setTickCount(m_AxisTicks[0]);
    m_xAxis->setMinorTickCount(1);
    m_xAxis->setLabelsAngle(0);
    auto col = m_xAxis->gridLineColor();
    col.setAlphaF(0.3);
    m_xAxis->setGridLineColor(col);
    col.setAlphaF(0.15);
    m_xAxis->setMinorGridLineColor(col);
  }
}

void m2Spectrum::OnAxisYTicksChanged(int v)
{
  m_AxisTicks[1] = v;
  if (m_yAxis)
  {
    m_yAxis->setTickCount(m_AxisTicks[1]);
    m_yAxis->setMinorTickCount(1);
    auto col = m_yAxis->gridLineColor();
    col.setAlphaF(0.3);
    m_yAxis->setGridLineColor(col);
    col.setAlphaF(0.15);
    m_yAxis->setMinorGridLineColor(col);
    // yAxis->setTitleText("Intensity");
  }
}

void m2Spectrum::OnRangeChangedAxisX(qreal xMin, qreal xMax)
{
  if (m_DataProvider.empty())
    return;

  m_LocalMaximumX = xMax;
  m_LocalMinimumX = xMin;

  auto axis = m_xAxis;
  if (xMax > m_GlobalMaximumX)
  {
    // MITK_INFO("OnRangeChangedAxisX") << xMax << " > " << m_GlobalMaximumX;
    this->m_MouseDragUpperDelta = std::abs(m_GlobalMaximumX - this->m_MouseDragCenterPos);
    axis->blockSignals(true);
    axis->setMax(m_GlobalMaximumX);
    axis->blockSignals(false);
  }

  if (xMin < m_GlobalMinimumX)
  {
    // MITK_INFO("OnRangeChangedAxisX") << xMin << " > " << m_GlobalMinimumX;
    this->m_MouseDragLowerDelta = std::abs(m_GlobalMinimumX - this->m_MouseDragCenterPos);
    axis->blockSignals(true);
    axis->setMin(m_GlobalMinimumX);
    axis->blockSignals(false);
  }

  UpdateCurrentMinMaxY();
  UpdateAllSeries();
  AutoZoomUseLocalExtremaY();
}

void m2Spectrum::UpdateAllSeries()
{
  for (auto &kv : m_DataProvider)
    kv.second->UpdateBoundaries(m_LocalMinimumX, m_LocalMaximumX);
}

void m2Spectrum::OnRangeChangedAxisY(qreal min, qreal max)
{
  if (max > m_GlobalMaximumY * 1.1)
  {
    m_yAxis->blockSignals(true);
    m_yAxis->setMax(m_GlobalMaximumY * 1.1);
    m_yAxis->blockSignals(false);
  }

  if (min < m_GlobalMinimumY)
  {
    m_yAxis->blockSignals(true);
    m_yAxis->setMin(m_GlobalMinimumY);
    m_yAxis->blockSignals(false);
  }
  UpdateSelectedArea();
}

void m2Spectrum::UpdateTitles(const mitk::DataNode *node)
{
  if (!m_xAxis || !m_yAxis)
    return;

  if (auto image = dynamic_cast<m2::SpectrumImage *>(node->GetData()))
  {
    auto xlabel = image->GetSpectrumType().XAxisLabel;
    auto ylabel = image->GetSpectrumType().YAxisLabel;
    m_xAxis->setTitleText(xlabel.c_str());
    m_yAxis->setTitleText(ylabel.c_str());
  }

  // switch (m_CurrentOverviewSpectrumType)
  // {
  //   case m2::SpectrumType::Maximum:
  //     m_yAxis->setTitleText("Intensity (maximum)");
  //     break;
  //   case m2::SpectrumType::Sum:
  //     m_yAxis->setTitleText("Intensity (sum)");
  //     break;
  //   case m2::SpectrumType::Mean:
  //     m_yAxis->setTitleText("Intensity (mean)");
  //     break;
  //   default:
  //     break;
  // }
}
