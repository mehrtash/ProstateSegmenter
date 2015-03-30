#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkPluginUtilities.h"
#include "itkLabelGeometryImageFilter.h"

#include "ProstateSegmenterCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  typedef unsigned char                           LabelVolumePixelType;
  typedef itk::Image<LabelVolumePixelType, 3>     LabelVolumeType;
  typedef itk::ImageFileReader<LabelVolumeType>   LabelVolumeReaderType;
  typedef itk::ImageFileWriter< LabelVolumeType>  LabelVolumeWriterType;

  typedef itk::ImageRegionConstIterator<LabelVolumeType> LabelVolumeIteratorType;

  // Read label
  LabelVolumeType::Pointer inputLabelVolume;
  if(inputLabelName!= "")
  {
    LabelVolumeReaderType::Pointer labelReader = LabelVolumeReaderType::New();
    labelReader->SetFileName(inputLabelName.c_str());
    labelReader->Update();
    inputLabelVolume = labelReader->GetOutput();

  typedef itk::LabelGeometryImageFilter< LabelVolumeType> LabelGeometryImageFilterType;
  LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
  labelGeometryImageFilter->SetInput( inputLabelVolume);
  //labelGeometryImageFilter->SetIntensityInput( intensityImage );
  // These generate optional outputs.
  labelGeometryImageFilter->CalculatePixelIndicesOn();
  labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
  labelGeometryImageFilter->CalculateOrientedLabelRegionsOn();

  labelGeometryImageFilter->Update();

  LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
  LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
  std::cout << "Number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
  std::cout << std::endl;
                                               
  for( allLabelsIt = allLabels.begin(); allLabelsIt != allLabels.end(); allLabelsIt++ )
  {
  LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
  std::cout << "\tLabel: " << (int)labelValue << std::endl;
  std::cout << "\tVolume: " << labelGeometryImageFilter->GetVolume(labelValue) << std::endl;
  //std::cout << "\tIntegrated Intensity: " << labelGeometryImageFilter->GetIntegratedIntensity(labelValue) << std::endl;
  std::cout << "\tCentroid: " << labelGeometryImageFilter->GetCentroid(labelValue) << std::endl;
  std::cout << "\tWeighted Centroid: " << labelGeometryImageFilter->GetWeightedCentroid(labelValue) << std::endl;
  std::cout << "\tAxes Length: " << labelGeometryImageFilter->GetAxesLength(labelValue) << std::endl;
  std::cout << "\tMajorAxisLength: " << labelGeometryImageFilter->GetMajorAxisLength(labelValue) << std::endl;
  std::cout << "\tMinorAxisLength: " << labelGeometryImageFilter->GetMinorAxisLength(labelValue) << std::endl;
  std::cout << "\tEccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue) << std::endl;
  std::cout << "\tElongation: " << labelGeometryImageFilter->GetElongation(labelValue) << std::endl;
  std::cout << "\tOrientation: " << labelGeometryImageFilter->GetOrientation(labelValue) << std::endl;
  std::cout << "\tBounding box: " << labelGeometryImageFilter->GetBoundingBox(labelValue) << std::endl;
  std::cout << std::endl << std::endl;
  }
  LabelVolumeType::Pointer outputLabelVolume= LabelVolumeType::New();
  outputLabelVolume->SetRegions(inputLabelVolume->GetLargestPossibleRegion());
  outputLabelVolume->Allocate();
  outputLabelVolume->FillBuffer(0);
  outputLabelVolume->CopyInformation(inputLabelVolume);
  outputLabelVolume->FillBuffer(0);

  if(outputLabelName.size()){
      LabelVolumeWriterType::Pointer writer = LabelVolumeWriterType::New();
      writer->SetInput(outputLabelVolume);
      writer->SetFileName(outputLabelName.c_str());
      writer->SetUseCompression(1);
      writer->Update();
    }
}
  try
    {

    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
