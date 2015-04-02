#include "algorithm"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkBoundingBox.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
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

  typedef itk::LabelGeometryImageFilter< LabelVolumeType> LabelGeometryImageFilterType;
  typedef itk::Matrix<double, 3, 3> MatrixType;
  // Read label
  LabelVolumeType::Pointer inputLabelVolume;
  if(inputLabelName!= "")
  {
    LabelVolumeReaderType::Pointer labelReader = LabelVolumeReaderType::New();
    labelReader->SetFileName(inputLabelName.c_str());
    labelReader->Update();
    inputLabelVolume = labelReader->GetOutput();
  //-------------------------------------------------------------------------------------------------------
  // Resampling
  //-------------------------------------------------------------------------------------------------------

    typedef itk::ResampleImageFilter<LabelVolumeType,LabelVolumeType >
      ResampleFilterType;
    typedef itk::IdentityTransform<double, 3>
      TransformType;
    typename TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();

    const typename LabelVolumeType::SpacingType& inputSpacing =
      inputLabelVolume->GetSpacing();


    const typename LabelVolumeType::RegionType& inputRegion =
      inputLabelVolume->GetLargestPossibleRegion();
    const typename LabelVolumeType::SizeType& inputSize =
      inputRegion.GetSize();

    // Compute the size of the output. The user specifies a spacing on
    // the command line. If the spacing is 0, the input spacing will be
    // used. The size (#of pixels) in the output is recomputed using
    // the ratio of the input and output sizes.
    typename LabelVolumeType::SpacingType outputSpacing;
    outputSpacing[0] = 1;
    outputSpacing[1] = 1;
    outputSpacing[2] = 1;

    typename LabelVolumeType::SizeType   outputSize;
    typedef typename LabelVolumeType::SizeType::SizeValueType SizeValueType;
    outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
    outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
    outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);

    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    itk::PluginFilterWatcher watcher(resampler, "Resample Volume",
                                     CLPProcessInformation);

    resampler->SetInput( inputLabelVolume );
    resampler->SetTransform( transform );
    resampler->SetOutputOrigin( inputLabelVolume->GetOrigin() );
    resampler->SetOutputSpacing( outputSpacing );
    resampler->SetOutputDirection( inputLabelVolume->GetDirection() );
    resampler->SetSize( outputSize );
    resampler->Update();
    LabelVolumeType::Pointer inputLabelVolumeResampled= LabelVolumeType::New();
    inputLabelVolumeResampled = resampler->GetOutput();
  //-------------------------------------------------------------------------------------------------------
  LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
  labelGeometryImageFilter->SetInput(inputLabelVolumeResampled);
  //labelGeometryImageFilter->SetIntensityInput( inputLabelVolume);
  // These generate optional outputs.
  labelGeometryImageFilter->CalculatePixelIndicesOn();
  labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
  labelGeometryImageFilter->CalculateOrientedLabelRegionsOn();

  itk::PluginFilterWatcher watcher2(labelGeometryImageFilter, "Geometry filter",
                                     CLPProcessInformation);
  labelGeometryImageFilter->Update();

  LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
  LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
  std::cout << "Number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
  std::cout << std::endl;
                                               
  for( allLabelsIt = allLabels.begin()+1; allLabelsIt != allLabels.end(); allLabelsIt++ )
  {
  LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
  //std::cout << "\tLabel: " << (int)labelValue << std::endl;
  //std::cout << "\tVolume: " << labelGeometryImageFilter->GetVolume(labelValue) << std::endl;
  //std::cout << "\tIntegrated Intensity: " << labelGeometryImageFilter->GetIntegratedIntensity(labelValue) << std::endl;
  //std::cout << "\tCentroid: " << labelGeometryImageFilter->GetCentroid(labelValue) << std::endl;
  //std::cout<<"\tEigenvectors: "<< labelGeometryImageFilter->GetEigenvectors(labelValue)<<std::endl;
  //std::cout << "\tWeighted Centroid: " << labelGeometryImageFilter->GetWeightedCentroid(labelValue) << std::endl;
  //std::cout << "\tAxes Length: " << labelGeometryImageFilter->GetAxesLength(labelValue) << std::endl;
  //std::cout << "\tMajorAxisLength: " << labelGeometryImageFilter->GetMajorAxisLength(labelValue) << std::endl;
  //std::cout << "\tMinorAxisLength: " << labelGeometryImageFilter->GetMinorAxisLength(labelValue) << std::endl;
  //std::cout << "\tEccentricity: " << labelGeometryImageFilter->GetEccentricity(labelValue) << std::endl;
  //std::cout << "\tElongation: " << labelGeometryImageFilter->GetElongation(labelValue) << std::endl;
  //std::cout << "\tOrientation: " << labelGeometryImageFilter->GetOrientation(labelValue) << std::endl;
  //std::cout << "\tBounding box: " << labelGeometryImageFilter->GetBoundingBox(labelValue) << std::endl;
  //std::cout << std::endl << std::endl;
  std::cout << "\tBounding box size: " << labelGeometryImageFilter->GetBoundingBoxSize(labelValue) << std::endl;
  if (labelValue == 1){
  // get centroid
  LabelVolumeType::PointType centroid =  labelGeometryImageFilter->GetCentroid(labelValue);
  // get major axis length
  double majorAxisLength = labelGeometryImageFilter->GetMajorAxisLength(labelValue);
  // get eigenvalues


  std::vector<double> eigenValues = labelGeometryImageFilter->GetEigenvalues(labelValue);
  // get eigenvectors
  MatrixType eigenVectors = labelGeometryImageFilter->GetEigenvectors(labelValue);
  itk::FixedArray<double, 6> boundingBox;
  boundingBox = labelGeometryImageFilter->GetBoundingBox(labelValue);
  //
  std::cout << "\tCentroid: " << centroid << std::endl;
  std::cout << "\tBounding box: " << boundingBox << std::endl;
  std::cout <<"\t major axis length: "<<majorAxisLength<<std::endl;
  std::cout << "\tEigenvalues: " << eigenValues[0]<<"\t"<<eigenValues[1]<<"\t"<<eigenValues[2]<< std::endl;
  std::cout << "\tEigenvectors: " << eigenVectors<< std::endl;

  int maxIndex = std::distance(eigenValues.begin(), std::max_element(eigenValues.begin(),eigenValues.end()));
  std::cout<<"\tmaxindex: "<<maxIndex<<std::endl;

  // calculate point 1
  // TODO: major axis length is wrong here. as spacing is not taken into account
  LabelVolumeType::PointType point1;
  point1[0] = centroid[0] + eigenVectors(0,0)*majorAxisLength/6;
  point1[1] = centroid[1] + eigenVectors(1,0)*majorAxisLength/6;
  point1[2] = centroid[2] + eigenVectors(2,0)*majorAxisLength/6;
  std::cout << "\tPoint1: " << point1<< std::endl;

  // calculate point 2
  LabelVolumeType::PointType point2;
  point2[0] = centroid[0] - eigenVectors(0,0)*majorAxisLength/6;
  point2[1] = centroid[1] - eigenVectors(1,0)*majorAxisLength/6;
  point2[2] = centroid[2] - eigenVectors(2,0)*majorAxisLength/6;
  std::cout << "\tPoint2: " << point2<< std::endl;

  LabelVolumeType::Pointer outputLabelVolume= LabelVolumeType::New();
  outputLabelVolume->SetRegions(inputLabelVolume->GetLargestPossibleRegion());
  outputLabelVolume->Allocate();
  outputLabelVolume->FillBuffer(0);
  outputLabelVolume->CopyInformation(inputLabelVolume);

  LabelVolumeType::IndexType index;


  std::cout<<"bounding box points: "<<boundingBox<<std::endl;
  LabelVolumeType::SizeType regionSize;
  regionSize[0] = boundingBox[1]-boundingBox[0];
  regionSize[1] = boundingBox[3]-boundingBox[2];
  regionSize[2] = boundingBox[5]-boundingBox[4];

  LabelVolumeType::IndexType regionIndex;
  regionIndex[0] = boundingBox[0];
  regionIndex[1] = boundingBox[2];
  regionIndex[2] = boundingBox[4];

  LabelVolumeType::RegionType region;
  region.SetSize(regionSize);
  region.SetIndex(regionIndex);

  itk::ImageRegionIteratorWithIndex<LabelVolumeType> imageIterator(outputLabelVolume,outputLabelVolume->GetLargestPossibleRegion());
  //itk::ImageRegionIteratorWithIndex<LabelVolumeType> imageIterator(outputLabelVolume,region);
  itk::ImageRegionIteratorWithIndex<LabelVolumeType> inputImageIterator(inputLabelVolume,inputLabelVolume->GetLargestPossibleRegion());
  //itk::ImageRegionIteratorWithIndex<LabelVolumeType> imageIterator(outputLabelVolume,region);

  LabelVolumeType::PointType clps;
  itk::ContinuousIndex<double,3> centroidIndex;
  centroidIndex[0] = centroid[0];
  centroidIndex[1] = centroid[1];
  centroidIndex[2] = centroid[2];

  std::cout<<"centroid: "<<centroidIndex<<std::endl;
  inputLabelVolumeResampled->TransformContinuousIndexToPhysicalPoint(centroidIndex,clps);
  std::cout<<"centroid lps: "<<clps<<std::endl;
  inputLabelVolume->TransformPhysicalPointToContinuousIndex(clps,centroidIndex);
  std::cout<<"centroid: "<<centroidIndex<<std::endl;

  int i = 0;
  while(!imageIterator.IsAtEnd())
    {
      i++;
      //if (true)
      if (inputImageIterator.Value() == labelValue)
      {
        index = imageIterator.GetIndex();
        //LabelVolumeType::PointType lps;
        //outputLabelVolume->TransformIndexToPhysicalPoint(index,lps);
        // change with max index
        double d1 = eigenVectors(0,0) * centroidIndex[0] + eigenVectors(1,0) * centroidIndex[1] + eigenVectors(2,0)* centroidIndex[2];
        double val1 = eigenVectors(0,0) * index[0] + eigenVectors(1,0) * index[1] + eigenVectors(2,0) * index[2] - d1;
        // change!!!!
        double d2 = eigenVectors(0,1) * centroidIndex[0] + eigenVectors(1,1) * centroidIndex[1] + eigenVectors(2,1)* centroidIndex[2];
        double val2 = eigenVectors(0,1) * index[0] + eigenVectors(1,1) * index[1] + eigenVectors(2,1) * index[2] - d2;
        //
        double d3 = eigenVectors(0,2) * centroidIndex[0] + eigenVectors(1,2) * centroidIndex[1] + eigenVectors(2,2)* centroidIndex[2];
        double val3 = eigenVectors(0,2) * index[0] + eigenVectors(1,2) * index[1] + eigenVectors(2,2) * index[2] - d3;
        ///*
        // change with max index
       // */
        /*
        double d3 = eigenVectors(0,0) * point1[0] + eigenVectors(1,0) * point1[1] + eigenVectors(2,0)* point1[2];
        double val3 = eigenVectors(0,0) * index[0] + eigenVectors(1,0) * index[1] + eigenVectors(2,0) * index[2] - d3;
        //
        double d4 = eigenVectors(0,0) * point2[0] + eigenVectors(1,0) * point2[1] + eigenVectors(2,0)* point2[2];
        double val4 = eigenVectors(0,0) * index[0] + eigenVectors(1,0) * index[1] + eigenVectors(2,0) * index[2] - d4;
        */
        if (i%1000 == 0){
        std::cout<<"\n( "<<index[0]<<", "<<index[1]<<", "<<index[2]<<" )";
        std::cout<<" val1: "<<val1<<" val2: "<<val2<<" val3: "<<val3;
        }
        if (val1 >0 && val2>0 && val3>0 )
        {
            imageIterator.Set(1);
        }
        if (val1 >0 && val2>0 && val3<0 )
        {
            imageIterator.Set(2);
        }
        if (val1 >0 && val2<0 && val3>0 )
        {
            imageIterator.Set(3);
        }
        if (val1 >0 && val2<0 && val3<0 )
        {
            imageIterator.Set(4);
        }
        if (val1 <0 && val2>0 && val3>0 )
        {
            imageIterator.Set(5);
        }
        if (val1 <0 && val2>0 && val3<0 )
        {
            imageIterator.Set(6);
        }
        if (val1 <0 && val2<0 && val3>0 )
        {
            imageIterator.Set(7);
        }
        if (val1 <0 && val2<0 && val3<0 )
        {
            imageIterator.Set(8);
        }
      }

    ++imageIterator;
    ++inputImageIterator;
    }
   //-------------------------------------------------------------------------------------------------------
  // Resampling again
  //-------------------------------------------------------------------------------------------------------

  /*
    const typename LabelVolumeType::SpacingType& inputSpacing =
              outputLabelVolumeResampled->GetSpacing();
    const typename LabelVolumeType::SpacingType& outputSpacing=
      inputLabelVolume->GetSpacing();
    const typename LabelVolumeType::RegionType& inputRegion =
      outputLabelVolumeResampled->GetLargestPossibleRegion();
    const typename LabelVolumeType::SizeType& inputSize =
      inputRegion.GetSize();

    outputSize[0] = static_cast<SizeValueType>(inputSize[0] * inputSpacing[0] / outputSpacing[0] + .5);
    outputSize[1] = static_cast<SizeValueType>(inputSize[1] * inputSpacing[1] / outputSpacing[1] + .5);
    outputSize[2] = static_cast<SizeValueType>(inputSize[2] * inputSpacing[2] / outputSpacing[2] + .5);

    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput( outputLabelVolumeResampled);
    resampler->SetTransform( transform );
    resampler->SetOutputOrigin( outputLabelVolumeResampled->GetOrigin() );
    resampler->SetOutputSpacing( outputSpacing );
    resampler->SetOutputDirection( outputLabelVolumeResampled->GetDirection() );
    resampler->SetSize( outputSize );
    resampler->UpdateLargestPossibleRegion();
    resampler->Update();
    LabelVolumeType::Pointer outputLabelVolume = LabelVolumeType::New();
    outputLabelVolume = resampler->GetOutput();
    itk::PluginFilterWatcher watcher(resampler, "Resample Volume", CLPProcessInformation);
    */
  //-------------------------------------------------------------------------------------------------------

  if(outputLabelName.size()){
      LabelVolumeWriterType::Pointer writer = LabelVolumeWriterType::New();
      writer->SetInput(outputLabelVolume);
      writer->SetFileName(outputLabelName.c_str());
      writer->SetUseCompression(1);
      writer->Update();
  }
 }
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
