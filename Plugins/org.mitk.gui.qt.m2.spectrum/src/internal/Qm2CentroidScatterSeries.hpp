#pragma once

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QScatterSeries>
#include <QtGui/QPainter>
#include <QtWidgets/QGraphicsItem>
#include <QtCharts/QValueAxis>

QT_CHARTS_USE_NAMESPACE

class Qm2CentroidMarker : public QGraphicsItem
{
public:
    Qm2CentroidMarker(qreal x, qreal y, qreal height)
        : m_x(x), m_y(y), m_height(height) {}

    QRectF boundingRect() const override {
        return QRectF(m_x - 1, 0, 2, m_height);
    }

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) override {
        Q_UNUSED(option);
        Q_UNUSED(widget);
        painter->setPen(Qt::black);
        painter->drawLine(QPointF(m_x, m_y), QPointF(m_x, 0));
    }

private:
    qreal m_x;
    qreal m_y;
    qreal m_height;
};

class Qm2CentroidScatterSeries : public QScatterSeries
{
    Q_OBJECT

public:
    Qm2CentroidScatterSeries(QObject *parent = nullptr) : QScatterSeries(parent) {
        connect(this, &QXYSeries::pointsReplaced, this, &Qm2CentroidScatterSeries::onPointsChanged);
        connect(this, &QXYSeries::pointAdded, this, &Qm2CentroidScatterSeries::onPointsChanged);
        connect(this, &QXYSeries::pointRemoved, this, &Qm2CentroidScatterSeries::onPointsChanged);
    }

    void attachMarkersToChart(QChart *chart) {
        m_chart = chart;
        onPointsChanged();
    }

private slots:
    void onPointsChanged() {
        if (!m_chart) return;

        for (auto item : m_markers) {
            m_chart->scene()->removeItem(item);
            delete item;
        }
        m_markers.clear();

        for (const QPointF &point : points()) {
            Qm2CentroidMarker *marker = new Qm2CentroidMarker(point.x(), point.y(), point.y());
            m_chart->scene()->addItem(marker);
            m_markers.append(marker);
        }
    }

private:
    QChart *m_chart = nullptr;
    QList<QGraphicsItem *> m_markers;
};