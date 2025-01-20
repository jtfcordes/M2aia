/*===================================================================

MSI applications for interactive analysis in MITK (M2aia)

Copyright (c) Jonas Cordes.

All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or https://www.github.com/jtfcordes/m2aia for details.

===================================================================*/

#pragma once

#include <M2aiaCoreExports.h>

#include <mitkAbstractFileIO.h>
#include <mitkIOMimeTypes.h>
#include <mitkImage.h>
#include <mitkItkImageIO.h>

#include <m2IntervalVector.h>
#include <m2ImzMLSpectrumImage.h>
#include <m2SpectrumImageStack.h>


namespace m2
{
  /**
   * Writes/Reads ImzML images
   * @ingroup Process
   */

  using namespace std::string_literals;
  class M2AIACORE_EXPORT ImzMLImageIO : public mitk::AbstractFileIO
  {
  public:
    // typedef mitk::MassSpectrometryImage InputType;
    ImzMLImageIO();

    std::string IMZML_MIMETYPE_NAME()
    {
      static std::string name = mitk::IOMimeTypes::DEFAULT_BASE_NAME() + ".image.imzML";
      return name;
    }

    mitk::CustomMimeType IMZML_MIMETYPE()
    {
      mitk::CustomMimeType mimeType(IMZML_MIMETYPE_NAME());
      mimeType.AddExtension("imzML");
      mimeType.AddExtension("imzml");
      mimeType.SetCategory("Images");
      mimeType.SetComment("imzML");
      return mimeType;
    }

    // -------------- AbstractFileReader -------------

    /**
     * @brief Reads a number of m2::ImzMLSpectrumImage from the file system
     * @return a vector of m2::ImzMLSpectrumImage
     */
    std::vector<mitk::BaseData::Pointer> DoRead() override;
    ConfidenceLevel GetReaderConfidenceLevel() const override;

    // -------------- AbstractFileWriter -------------
    /**
     * @brief Writes a m2::ImzMLSpectrumImage  as its *.imzML and *.ibd file to disk
     *
     */
    void Write() override;
    void SetIntervalVector(m2::IntervalVector::Pointer data);
    void SetDataTypeXAxis(m2::NumericType type){m_DataTypeXAxis = type;}
    void SetDataTypeYAxis(m2::NumericType type){m_DataTypeYAxis = type;}
    void SetSpectrumFormat(m2::SpectrumFormat type){m_SpectrumFormat = type;}

    ConfidenceLevel GetWriterConfidenceLevel() const override;
    std::string GetIBDOutputPath() const;
    std::string GetImzMLOutputPath() const;
    void WriteContinuousProfile(m2::ImzMLSpectrumImage::SpectrumVectorType & spectra) const;
    void WriteContinuousCentroid(m2::ImzMLSpectrumImage::SpectrumVectorType & spectra) const;
    void WriteContinuousCentroid3DStack(const m2::SpectrumImageStack *) const;
    void WriteProcessedProfile(m2::ImzMLSpectrumImage::SpectrumVectorType & spectra) const;
    void WriteProcessedCentroid(m2::ImzMLSpectrumImage::SpectrumVectorType & spectra) const;

  private:
    void EvaluateSpectrumFormatType(m2::SpectrumImage *);
    void LoadAssociatedData(m2::ImzMLSpectrumImage *);
    ImzMLImageIO *IOClone() const override;

    std::string RemoveExtensionFromPath(std::string path);

    m2::IntervalVector::Pointer m_Intervals;
    m2::NumericType m_DataTypeXAxis = m2::NumericType::Float;
    m2::NumericType m_DataTypeYAxis = m2::NumericType::Float;
    m2::SpectrumFormat m_SpectrumFormat = m2::SpectrumFormat::None;

  
    std::map<std::string, std::string> TextToCodeMap = {{"16-bit float"s, "1000520"s},
                                                        {"32-bit integer"s, "1000519"s},
                                                        {"32-bit float"s, "1000521"s},
                                                        {"64-bit integer"s, "1000522"s},
                                                        {"64-bit float"s, "1000523"s},
                                                        {"continuous"s, "1000030"s},
                                                        {"processed"s, "1000031"s},
                                                        {"zlib compression"s, "1000574"s},
                                                        {"no compression"s, "1000576"s},
                                                        {"positive scan"s, "1000130"s},
                                                        {"negative scan"s, "1000129"s},
                                                        {"centroid spectrum"s, "1000127"s},
                                                        {"profile spectrum"s, "1000128"s}};

    // mode_code, mode --> continuous | processed
    // sha1sum
    // uuid
    // polarity_code, polarity --> negative scan | positive scan
    // mz_data_type_code, mz_data_type --> 16-bit float | 32-bit float | 64-bit float | 32-bit integer | 64-bit integer
    // mz_compression_code, mz_compression --> zlib compression | no compression
    // int_data_type_code, int_data_type --> 16-bit float | 32-bit float | 64-bit float | 32-bit integer | 64-bit
    // integer int_compression_code, int_compression --> zlib compression | no compression max count of pixels y, max
    // count of pixels run_id num_spectra recursive_append_spectrum_tag --> "

    const std::string IMZML_TEMPLATE_START =
      //@require(uuid, sha1sum, mz_data_type, int_data_type, run_id, spectra, mode, obo_codes, mz_compression,
      // int_compression, polarity)
      "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
      "xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd\" "
      "version=\"1.1\">\n"
      "<cvList count=\"3\">\n"
      "<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"1.3.1\" "
      "URI=\"http://psidev.info/ms/mzML/psi-ms.obo\"/>\n"
      "<cv id=\"UO\" fullName=\"Unit Ontology\" version=\"1.15\" "
      "URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"/>\n"
      "<cv id=\"IMS\" fullName=\"Imaging MS Ontology\" version=\"0.9.1\" "
      "URI=\"http://www.maldi-msi.org/download/imzml/imagingMS.obo\"/>\n"
      "</cvList>\n"
      "<fileDescription>\n"
      "<fileContent>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000294\" name=\"mass spectrum\" value=\"\"/>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:{spectrumtype_code}\" name=\"{spectrumtype}\" value=\"\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:{mode_code}\" name=\"{mode}\" value=\"\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000080\" name=\"universally unique identifier\" value=\"{uuid}\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000091\" name=\"ibd SHA-1\" value=\"{sha1sum}\"/>\n"
      "</fileContent>\n"
      "</fileDescription>\n"
      "<referenceableParamGroupList count=\"3\">\n"
      "<referenceableParamGroup id=\"spectrum\">\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000294\" name=\"mass spectrum\" value=\"\"/>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:{spectrumtype_code}\" name=\"{spectrumtype}\" value=\"\"/>\n"
      "{#polarity}<cvParam cvRef=\"MS\" accession=\"MS:{polarity_code}\" name=\"{polarity}\" value=\"\"/>\n{/polarity}"
      "</referenceableParamGroup>\n"
      "<referenceableParamGroup id=\"mzArray\">\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\"/>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:{mz_data_type_code}\" name=\"{mz_data_type}\" value=\"\"/>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:{mz_compression_code}\" name=\"{mz_compression}\" value=\"\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>\n"
      "</referenceableParamGroup>\n"
      "<referenceableParamGroup id=\"intensityArray\">\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\"/>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:{int_data_type_code}\" name=\"{int_data_type}\" value=\"\"/>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:{int_compression_code}\" name=\"{int_compression}\" value=\"\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000101\" name=\"external data\" value=\"true\"/>\n"
      "</referenceableParamGroup>\n"
      "</referenceableParamGroupList>\n"
      "<softwareList count=\"1\">\n"
      "<software id=\"M2aia\" version=\"1.0alpha\">\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000532\" name=\"M2aia\" value=\"\"/>\n"
      "</software>\n"
      "</softwareList>\n"
      "<scanSettingsList count=\"1\">\n"
      "<scanSettings id=\"scanSettings1\">\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000042\" name=\"max count of pixels x\" value=\"{size_x}\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000043\" name=\"max count of pixels y\" value=\"{size_y}\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000044\" name=\"max dimension x\" value=\"{max dimension x}\" "
      "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000045\" name=\"max dimension y\" value=\"{max dimension y}\" "
      "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000046\" name=\"pixel size x\" value=\"{pixel size x}\" "
      "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000047\" name=\"pixel size y\" value=\"{pixel size y}\" "
      "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000053\" name=\"absolute position offset x\" value=\"{origin x}\" "
      "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n"
      "<cvParam cvRef=\"IMS\" accession=\"IMS:1000054\" name=\"absolute position offset y\" value=\"{origin y}\" "
      "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"/>\n"
      "</scanSettings>\n"
      "</scanSettingsList>\n"
      "<instrumentConfigurationList count=\"1\">\n"
      "<instrumentConfiguration id=\"IC1\">\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"Unknown\" value=\"\"/>\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000463\" name=\"Unknown\" value=\"\"/>\n"
      "</instrumentConfiguration>\n"
      "</instrumentConfigurationList>\n"
      "<dataProcessingList count=\"1\">\n"
      "<dataProcessing id=\"M2aiaProcessing\">\n"
      "<processingMethod order=\"1\" softwareRef=\"M2aia\">\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000543\" name=\"None\" value=\"\"/>\n"
      "</processingMethod>\n"
      "</dataProcessing>\n"
      "</dataProcessingList>\n"
      "<run defaultInstrumentConfigurationRef=\"IC1\" id=\"Experiment{run_id}\">\n"
      "<spectrumList count=\"{num_spectra}\" defaultDataProcessingRef=\"M2aiaProcessing\">\n";

    const std::string IMZML_TEMPLATE_END = "</spectrumList>\n"
                                           "</run>\n"
                                           "</mzML>\n";

    // index
    // x,y,z
    // use_z
    // mz_len, mz_offset
    // int_len, int_offset
    // use_encoded_length

    const std::string IMZML_SPECTRUM_TEMPLATE =

      "<spectrum defaultArrayLength=\"0\" id=\"spectrum={index}\" index=\"{index}\">\n"
      "<referenceableParamGroupRef ref=\"spectrum\"/>\n"
      "{#tic}<cvParam accession=\"MS:1000285\" cvRef=\"MS\" name=\"total ion current\" value=\"{tic}\"/>\n{/tic}"
      "<scanList count=\"1\">\n"
      "<cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>\n"
      "<scan instrumentConfigurationRef=\"IC1\">\n"
      "<cvParam accession=\"IMS:1000050\" cvRef=\"IMS\" name=\"position x\" value=\"{x}\"/>\n"
      "<cvParam accession=\"IMS:1000051\" cvRef=\"IMS\" name=\"position y\" value=\"{y}\"/>\n"
      "{#z}<cvParam accession=\"IMS:1000052\" cvRef=\"IMS\" name=\"position z\" value=\"{z}\"/>\n{/z}"
      "</scan>\n"
      "</scanList>\n"
      "<binaryDataArrayList count=\"2\">\n"
      "<binaryDataArray encodedLength=\"0\">\n"
      "<referenceableParamGroupRef ref=\"mzArray\"/>\n"
      "<cvParam accession=\"IMS:1000103\" cvRef=\"IMS\" name=\"external array length\" value=\"{mz_len}\"/>\n"
      "<cvParam accession=\"IMS:1000104\" cvRef=\"IMS\" name=\"external encoded length\" value=\"{mz_enc_len}\"/>\n"
      "<cvParam accession=\"IMS:1000102\" cvRef=\"IMS\" name=\"external offset\" value=\"{mz_offset}\"/>\n"
      "<binary/>\n"
      "</binaryDataArray>\n"
      "<binaryDataArray encodedLength=\"0\">\n"
      "<referenceableParamGroupRef ref=\"intensityArray\"/>\n"
      "<cvParam accession=\"IMS:1000103\" cvRef=\"IMS\" name=\"external array length\" value=\"{int_len}\"/>\n"
      "<cvParam accession=\"IMS:1000104\" cvRef=\"IMS\" name=\"external encoded length\" value=\"{int_enc_len}\"/>\n"
      "<cvParam accession=\"IMS:1000102\" cvRef=\"IMS\" name=\"external offset\" value=\"{int_offset}\"/>\n"
      "<binary/>\n"
      "</binaryDataArray>\n"
      "</binaryDataArrayList>\n"
      "</spectrum>\n";

  }; 

} // namespace m2


