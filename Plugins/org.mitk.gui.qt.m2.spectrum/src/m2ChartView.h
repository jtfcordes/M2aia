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

#include <QAreaSeries>
#include <QFuture>
#include <QLineSeries>
#include <QObject>
#include <QScatterSeries>
#include <QtCharts/QChartView>
#include <QtCharts/QChart>

#include <QtWidgets/QRubberBand>
#include <iostream>
#include <QmitkChartWidget.h>

namespace m2
{
  class ChartView : public QtCharts::QChartView
  {
    Q_OBJECT

  public:
    ChartView(QWidget *parent = 0);
    ChartView(QtCharts::QChart *chart, QWidget *parent = 0);
    //void OnResize();
//    int GetTickCountX() { return m_AxisXTicks; }

  //  int GetTickCountY() { return m_AxisYTicks; }

  
  signals:
    void mousePress(QPoint pos, qreal mz, qreal intValue, Qt::MouseButton button, Qt::KeyboardModifiers);
    void mouseMove(QPoint pos, qreal mz, qreal intValue, Qt::MouseButton button, Qt::KeyboardModifiers);
    void mouseRelease(QPoint pos, qreal mz, qreal intValue, Qt::MouseButton button, Qt::KeyboardModifiers);
    void mouseDoubleClick(QPoint pos, qreal mz, qreal intValue, Qt::MouseButton button, Qt::KeyboardModifiers);
    void mouseWheel(QPoint pos, qreal mz, qreal intValue, int angle, Qt::KeyboardModifiers);

  protected:

    void resizeEvent(QResizeEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    //void keyPressEvent(QKeyEvent *event) override;

  private:
    bool m_isTouching;
  };
} // namespace m2
