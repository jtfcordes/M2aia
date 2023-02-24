/*============================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center (DKFZ)
All rights reserved.

Use of this source code is governed by a 3-clause BSD license that can be
found in the LICENSE file.

============================================================================*/

#include "QmitkDataCompressionView.h"

#include <QMessageBox>
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// m2
#include <m2CoreCommon.h>
#include <m2IntervalVector.h>
#include <m2PcaImageFilter.h>
#include <m2SpectrumImageBase.h>
#include <m2TSNEImageFilter.h>
#include <m2MultiSliceFilter.h>

// mitk
#include <mitkDockerHelper.h>
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkImageAccessByItk.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateFunction.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateProperty.h>
#include <mitkProgressBar.h>
#include <usModuleRegistry.h>

// itk
#include <itkShrinkImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkVectorImageToImageAdaptor.h>

// boost
#include <boost/algorithm/string.hpp>

#include <berryPlatformUI.h>

// Don't forget to initialize the VIEW_ID.
const std::string QmitkDataCompressionView::VIEW_ID = "org.mitk.views.m2.docker.datacompression";
using DisplayImageType = itk::Image<m2::DisplayImagePixelType, 3>;
using VectorImageAdaptorType = itk::VectorImageToImageAdaptor<m2::DisplayImagePixelType, 3>;
using VectorImageType = itk::VectorImage<m2::DisplayImagePixelType, 3>;

void QmitkDataCompressionView::CreateQtPartControl(QWidget *parent)
{
  // Setting up the UI is a true pleasure when using .ui files, isn't it?
  m_Controls.setupUi(parent);

  m_Controls.imageSelection->SetDataStorage(GetDataStorage());
  m_Controls.imageSelection->SetAutoSelectNewNodes(true);
  m_Controls.imageSelection->SetNodePredicate(
    mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<m2::SpectrumImageBase>::New(),
                                mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object"))));
  m_Controls.imageSelection->SetSelectionIsOptional(true);
  m_Controls.imageSelection->SetEmptyInfo(QString("Image selection"));
  m_Controls.imageSelection->SetPopUpTitel(QString("Image"));

  auto isChildOfImageNode = [this](const mitk::DataNode *node) {
    if (auto imageNode = m_Controls.imageSelection->GetSelectedNode())
    {
      auto childNodes = GetDataStorage()->GetDerivations(imageNode);
      using namespace std;
      return find(begin(*childNodes), end(*childNodes), node) != end(*childNodes);
    }
    return false;
  };

  m_Controls.peakListSelection->SetDataStorage(GetDataStorage());
  m_Controls.peakListSelection->SetNodePredicate(
    mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<m2::IntervalVector>::New(),
                                mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object")),
                                mitk::NodePredicateFunction::New(isChildOfImageNode)));

  m_Controls.peakListSelection->SetAutoSelectNewNodes(true);
  m_Controls.peakListSelection->SetSelectionIsOptional(true);
  m_Controls.peakListSelection->SetEmptyInfo(QString("PeakList selection"));
  m_Controls.peakListSelection->SetPopUpTitel(QString("PeakList"));

  // Wire up the UI widgets with our functionality.
  connect(m_Controls.imageSelection,
          &QmitkSingleNodeSelectionWidget::CurrentSelectionChanged,
          this,
          &QmitkDataCompressionView::OnImageChanged);

  connect(m_Controls.btnRunPCA, SIGNAL(clicked()), this, SLOT(OnStartPCA()));
  connect(m_Controls.btnRunTSNE, SIGNAL(clicked()), this, SLOT(OnStartTSNE()));
  connect(m_Controls.btnRunSparsePCA, SIGNAL(clicked()), this, SLOT(OnStartSparsePCA()));
  connect(m_Controls.btnRunUmap, SIGNAL(clicked()), this, SLOT(OnStartUMAP()));

  connect(m_Controls.btnExport, &QPushButton::clicked, this, []() {
    try
    {
      if (auto platform = berry::PlatformUI::GetWorkbench())
        if (auto workbench = platform->GetActiveWorkbenchWindow())
          if (auto page = workbench->GetActivePage())
            if (page.IsNotNull())
              page->ShowView("org.mitk.views.m2.imzmlexport", "", 1);
    }
    catch (berry::PartInitException &e)
    {
      BERRY_ERROR << "Error: " << e.what() << std::endl;
    }
  });


  connect(m_Controls.btnPicking, &QPushButton::clicked, this, []() {
    try
    {
      if (auto platform = berry::PlatformUI::GetWorkbench())
        if (auto workbench = platform->GetActiveWorkbenchWindow())
          if (auto page = workbench->GetActivePage())
            if (page.IsNotNull())
              page->ShowView("org.mitk.views.m2.PeakPicking", "", 1);
    }
    catch (berry::PartInitException &e)
    {
      BERRY_ERROR << "Error: " << e.what() << std::endl;
    }
  });
  

  // Make sure to have a consistent UI state at the very beginning.
  this->OnImageChanged(m_Controls.imageSelection->GetSelectedNodes());

  m_Controls.lassoLabelImageSelection->SetDataStorage(GetDataStorage());
  m_Controls.lassoLabelImageSelection->SetAutoSelectNewNodes(true);
  m_Controls.lassoLabelImageSelection->SetNodePredicate(
    mitk::NodePredicateAnd::New(mitk::TNodePredicateDataType<mitk::LabelSetImage>::New(),
                                mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object"))));
  m_Controls.lassoLabelImageSelection->SetSelectionIsOptional(true);
  m_Controls.lassoLabelImageSelection->SetEmptyInfo(QString("LabelImage selection"));
  m_Controls.lassoLabelImageSelection->SetPopUpTitel(QString("LabelImage"));
}

void QmitkDataCompressionView::SetFocus() {}

void QmitkDataCompressionView::OnImageChanged(const QmitkSingleNodeSelectionWidget::NodeList &)
{
  this->EnableWidgets(m_Controls.imageSelection->GetSelectedNode().IsNotNull());
}

void QmitkDataCompressionView::OnPeakListChanged(const QmitkSingleNodeSelectionWidget::NodeList &)
{
  this->EnableWidgets(m_Controls.imageSelection->GetSelectedNode().IsNotNull());
}

void QmitkDataCompressionView::EnableWidgets(bool /*enable*/)
{
  // m_Controls
}

void QmitkDataCompressionView::OnStartUMAP()
{
  if (auto node = m_Controls.imageSelection->GetSelectedNode())
  {
    if (auto image = dynamic_cast<m2::SpectrumImageBase *>(node->GetData()))
    {
      if (!image->GetIsDataAccessInitialized())
        return;
      if (image->GetSpectrumType().Format == m2::SpectrumFormat::ContinuousCentroid)
      {
        if (mitk::DockerHelper::CheckDocker())
        {
          mitk::DockerHelper helper("m2aia/extensions:umap");
          helper.AddData(image, "--imzml", "*.imzML");
          helper.AddAutoLoadOutput("--image", "umap_image.nrrd");
          helper.AddApplicationArgument("--num_comp", std::to_string(m_Controls.spnBxComponents->value()));
          helper.AddApplicationArgument("--num_neighbors", std::to_string(m_Controls.spnBxNeighbors->value()));
          helper.AddApplicationArgument("--min_dist", std::to_string(m_Controls.spnBxMinDistance->value()));
          helper.AddApplicationArgument("--metric", m_Controls.cmbBxMetric->currentText().toStdString());
          const auto results = helper.GetResults();

          // convert 3D image to 3D vector image
          const auto image = dynamic_cast<mitk::Image *>(results[0].GetPointer());
          const auto dims = image->GetDimensions();
          const auto components = dims[2];
          auto vimage = mitk::HelperUtils::GetVectorImage3D({dims[0], dims[1], 1}, components);
          mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> outAcc(vimage);
          mitk::ImagePixelReadAccessor<float, 3> inAcc(image);

          for (unsigned int pixelInSlice = 0; pixelInSlice < dims[0] * dims[1]; ++pixelInSlice)
            for (unsigned int i = 0; i < components; ++i)
              outAcc.GetData()[pixelInSlice * components + i] = inAcc.GetData()[i * dims[0] * dims[1] + pixelInSlice];

          auto newNode = mitk::DataNode::New();
          newNode->SetData(vimage);
          newNode->SetName(node->GetName() + "_umap");
          GetDataStorage()->Add(newNode, node);
        }
      }
    }
  }
}

void QmitkDataCompressionView::OnStartSparsePCA()
{
  if (auto node = m_Controls.imageSelection->GetSelectedNode())
  {
    if (auto image = dynamic_cast<m2::SpectrumImageBase *>(node->GetData()))
    {
      if (!image->GetIsDataAccessInitialized())
        return;
      if (image->GetSpectrumType().Format == m2::SpectrumFormat::ContinuousCentroid)
      {
        if (mitk::DockerHelper::CheckDocker())
        {
          mitk::DockerHelper helper("sparse_pca");
          helper.AddData(image, "--imzml", "*.imzML");

          auto iter = helper.AddLoadLaterOutput("--csv", "pca_data.csv");
          helper.AddAutoLoadOutput("--image", "pca_data.nrrd");
          auto results = helper.GetResults();

          // convert 3D image to 3D vector image
          const auto image = dynamic_cast<mitk::Image *>(results[0].GetPointer());
          const auto dims = image->GetDimensions();
          const auto components = dims[2];
          auto vimage = mitk::HelperUtils::GetVectorImage3D({dims[0], dims[1], 1}, components);
          mitk::ImagePixelWriteAccessor<m2::DisplayImagePixelType, 3> outAcc(vimage);
          mitk::ImagePixelReadAccessor<float, 3> inAcc(image);
          auto targetGeom = vimage->GetGeometry();
          auto sourceGeom = image->GetGeometry();
          targetGeom->SetSpacing(sourceGeom->GetSpacing());
          targetGeom->SetOrigin(sourceGeom->GetOrigin());

          for (unsigned int pixelInSlice = 0; pixelInSlice < dims[0] * dims[1]; ++pixelInSlice)
            for (unsigned int i = 0; i < components; ++i)
              outAcc.GetData()[pixelInSlice * components + i] = inAcc.GetData()[i * dims[0] * dims[1] + pixelInSlice];

          // load csv
          const auto filePath = helper.GetFilePath(iter->path);
          std::ifstream t(filePath);
          if (t.good())
          {
            std::string text((std::istreambuf_iterator<char>(t)), std::istreambuf_iterator<char>());
            std::vector<std::string> lines;
            boost::split(lines, text, [](char c) { return c == '\n'; });
            lines.erase(lines.begin());
            lines.erase(lines.rbegin().base());

            auto *table = m_Controls.sparsePcaTable;
            table->clearContents();
            table->setRowCount(lines.size());
            table->setColumnCount(3);
            unsigned int row = 0;
            table->blockSignals(true);
            table->setHorizontalHeaderLabels({"m/z", "pc", "value"});

            for (std::string line : lines)
            {
              boost::erase_all(line, "\"");
              std::vector<std::string> words;
              boost::split(words, line, [](char c) { return c == ','; });

              auto item = new QTableWidgetItem(words[1].c_str());
              table->setItem(row, 0, item);

              item = new QTableWidgetItem(words[2].c_str());
              table->setItem(row, 1, item);

              item = new QTableWidgetItem(words[3].c_str());
              table->setItem(row, 2, item);
              ++row;
            }
            table->sortByColumn(1);
            table->blockSignals(false);
          }

          auto newNode = mitk::DataNode::New();
          newNode->SetData(vimage);
          newNode->SetName(node->GetName() + "_sparse_pca");
          GetDataStorage()->Add(newNode, node);
        }
      }
    }
  }
}

void QmitkDataCompressionView::OnStartPCA()
{
  auto imageNode = m_Controls.imageSelection->GetSelectedNode();
  auto vectorNode = m_Controls.peakListSelection->GetSelectedNode();

  if (!(imageNode || vectorNode))
    return;

  auto image = dynamic_cast<m2::SpectrumImageBase *>(imageNode->GetData());
  auto vector = dynamic_cast<m2::IntervalVector *>(vectorNode->GetData());

  if (!(image || vector))
    return;

  if (!image->GetIsDataAccessInitialized())
    return;

  auto filter = m2::PcaImageFilter::New();
  filter->SetMaskImage(image->GetMaskImage());

  // const auto &peakList = m_PeakList;

  std::vector<mitk::Image::Pointer> temporaryImages;
  auto progressBar = mitk::ProgressBar::GetInstance();
  progressBar->AddStepsToDo(vector->GetIntervals().size() + 1);
  size_t inputIdx = 0;
  for (size_t row = 0; row < vector->GetIntervals().size(); ++row)
  {
    progressBar->Progress();
    temporaryImages.push_back(mitk::Image::New());
    temporaryImages.back()->Initialize(image);
    const auto mz = vector->GetIntervals().at(row).x.mean();
    image->GetImage(mz, image->ApplyTolerance(mz), image->GetMaskImage(), temporaryImages.back());
    filter->SetInput(inputIdx, temporaryImages.back());
    ++inputIdx;
  }

  if (temporaryImages.size() <= 2)
  {
    progressBar->Progress();
    QMessageBox::warning(nullptr,
                         "Select images first!",
                         "Select at least three peaks!",
                         QMessageBox::StandardButton::NoButton,
                         QMessageBox::StandardButton::Ok);
    return;
  }

  filter->SetNumberOfComponents(m_Controls.pca_dims->value());
  filter->Update();
  progressBar->Progress();

  auto outputNode = mitk::DataNode::New();
  mitk::Image::Pointer data = filter->GetOutput(0);
  outputNode->SetData(data);
  outputNode->SetName("PCA");
  this->GetDataStorage()->Add(outputNode, imageNode.GetPointer());
  image->InsertImageArtifact("PCA", data);

  // auto outputNode2 = mitk::DataNode::New();
  // mitk::Image::Pointer data2 = filter->GetOutput(1);
  // outputNode2->SetData(data2);
  // outputNode2->SetName("pcs");
  // this->GetDataStorage()->Add(outputNode2, node.GetPointer());
}

void QmitkDataCompressionView::OnStartTSNE()
{
  if (auto node = m_Controls.imageSelection->GetSelectedNode())
  {
    if (auto image = dynamic_cast<m2::SpectrumImageBase *>(node->GetData()))
    {
      auto p = node->GetProperty("name")->Clone();
      static_cast<mitk::StringProperty *>(p.GetPointer())->SetValue("PCA");
      auto derivations = this->GetDataStorage()->GetDerivations(node, mitk::NodePredicateProperty::New("name", p));
      if (derivations->size() == 0)
      {
        QMessageBox::warning(nullptr,
                             "PCA required!",
                             "The t-SNE uses the top five features of the PCA transformed images! Start a PCA"
                             "first.", QMessageBox::StandardButton::NoButton, QMessageBox::StandardButton::Ok);
        return;
      }

      auto pcaImage = dynamic_cast<mitk::Image *>(derivations->front()->GetData());
      const auto pcaComponents = pcaImage->GetPixelType().GetNumberOfComponents();

      auto filter = m2::TSNEImageFilter::New();
      filter->SetPerplexity(m_Controls.tsne_perplexity->value());
      filter->SetIterations(m_Controls.tnse_iters->value());
      filter->SetTheta(m_Controls.tsne_theta->value());

      using MaskImageType = itk::Image<mitk::LabelSetImage::PixelType, 3>;
      MaskImageType::Pointer maskImageItk;
      mitk::Image::Pointer maskImage;
      mitk::CastToItkImage(image->GetMaskImage(), maskImageItk);
      auto caster = itk::ShrinkImageFilter<MaskImageType, MaskImageType>::New();
      caster->SetInput(maskImageItk);
      caster->SetShrinkFactor(0, m_Controls.tsne_shrink->value());
      caster->SetShrinkFactor(1, m_Controls.tsne_shrink->value());
      caster->SetShrinkFactor(2, 1);
      caster->Update();

      mitk::CastToMitkImage(caster->GetOutput(), maskImage);

      filter->SetMaskImage(maskImage);
      // const auto &peakList = image->GetPeaks();
      size_t index = 0;

      mitk::ImageReadAccessor racc(pcaImage);
      auto *inputData = static_cast<const typename DisplayImageType::PixelType *>(racc.GetData());

      std::vector<mitk::Image::Pointer> bufferedImages(pcaComponents);
      unsigned int n = pcaImage->GetDimensions()[0] * pcaImage->GetDimensions()[1] * pcaImage->GetDimensions()[2];
      for (auto &I : bufferedImages)
      {
        I = mitk::Image::New();
        I->Initialize(image);
        {
          mitk::ImageWriteAccessor acc(I);
          auto outCData = static_cast<typename DisplayImageType::PixelType *>(acc.GetData());
          for (unsigned int k = 0; k < n; ++k)
            *(outCData + k) = *(inputData + (k * pcaComponents) + index);
        }
        
        DisplayImageType::Pointer cImage;
        mitk::CastToItkImage(I, cImage);

        auto caster = itk::ShrinkImageFilter<DisplayImageType, DisplayImageType>::New();
        caster->SetInput(cImage);
        caster->SetShrinkFactor(0, m_Controls.tsne_shrink->value());
        caster->SetShrinkFactor(1, m_Controls.tsne_shrink->value());
        caster->SetShrinkFactor(2, 1);
        caster->Update();

        // Buffer the image
        mitk::CastToMitkImage(caster->GetOutput(), I);
        filter->SetInput(index, I);
        ++index;
      }
      filter->Update();

      auto outputNode = mitk::DataNode::New();
      auto data = m2::MultiSliceFilter::ConvertMitkVectorImageToRGB(ResampleVectorImage(filter->GetOutput(), image));
      outputNode->SetData(data);
      outputNode->SetName("tSNE");
      this->GetDataStorage()->Add(outputNode, node.GetPointer());
      image->InsertImageArtifact("tSNE", data);
    }
  }
}

mitk::Image::Pointer QmitkDataCompressionView::ResampleVectorImage(mitk::Image::Pointer vectorImage,
                                                            mitk::Image::Pointer referenceImage)
{
  const unsigned int components = vectorImage->GetPixelType().GetNumberOfComponents();

  VectorImageType::Pointer vectorImageItk;
  mitk::CastToItkImage(vectorImage, vectorImageItk);

  DisplayImageType::Pointer referenceImageItk;
  mitk::CastToItkImage(referenceImage, referenceImageItk);

  auto resampledVectorImageItk = VectorImageType::New();
  resampledVectorImageItk->SetOrigin(referenceImageItk->GetOrigin());
  resampledVectorImageItk->SetDirection(referenceImageItk->GetDirection());
  resampledVectorImageItk->SetSpacing(referenceImageItk->GetSpacing());
  resampledVectorImageItk->SetRegions(referenceImageItk->GetLargestPossibleRegion());
  resampledVectorImageItk->SetNumberOfComponentsPerPixel(components);
  resampledVectorImageItk->Allocate();
  itk::VariableLengthVector<m2::DisplayImagePixelType> v(components);
  v.Fill(0);
  resampledVectorImageItk->FillBuffer(v);

  auto inAdaptor = VectorImageAdaptorType::New();
  auto outAdaptor = VectorImageAdaptorType::New();
  using LinearInterpolatorType = itk::LinearInterpolateImageFunction<VectorImageAdaptorType>;
  using TransformType = itk::IdentityTransform<m2::DisplayImagePixelType, 3>;

  for (unsigned int i = 0; i < components; ++i)
  {
    inAdaptor->SetExtractComponentIndex(i);
    inAdaptor->SetImage(vectorImageItk);
    inAdaptor->SetOrigin(vectorImageItk->GetOrigin());
    inAdaptor->SetDirection(vectorImageItk->GetDirection());
    inAdaptor->SetSpacing(vectorImageItk->GetSpacing());

    outAdaptor->SetExtractComponentIndex(i);
    outAdaptor->SetImage(resampledVectorImageItk);

    auto resampler = itk::ResampleImageFilter<VectorImageAdaptorType, DisplayImageType>::New();
    resampler->SetInput(inAdaptor);
    resampler->SetOutputParametersFromImage(referenceImageItk);
    resampler->SetInterpolator(LinearInterpolatorType::New());
    resampler->SetTransform(TransformType::New());
    resampler->Update();

    itk::ImageAlgorithm::Copy<DisplayImageType, VectorImageAdaptorType>(
      resampler->GetOutput(),
      outAdaptor,
      resampler->GetOutput()->GetLargestPossibleRegion(),
      outAdaptor->GetLargestPossibleRegion());
  }

  mitk::Image::Pointer outImage;
  mitk::CastToMitkImage(resampledVectorImageItk, outImage);

  return outImage;
}