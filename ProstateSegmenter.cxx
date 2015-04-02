/*==============================================================================
  Program: Prostate segmentor
  Copyright (c) Brigham and Women's Hospital
  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
==============================================================================*/

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
  try
    {
    // Read label
    LabelVolumeType::Pointer inputLabelVolume;
    if(inputLabelName!= "")
      {
      LabelVolumeReaderType::Pointer labelReader = LabelVolumeReaderType::New();
      labelReader->SetFileName(inputLabelName.c_str());
      labelReader->Update();
      inputLabelVolume = labelReader->GetOutput();
      typedef itk::NearestNeighborInterpolateImageFunction<LabelVolumeType,
                  double> NearestNeighborInterpolatorType;
      typename NearestNeighborInterpolatorType::Pointer 
      nearestNeighborInterpolator = NearestNeighborInterpolatorType::New();
    //--------------------------------------------------------------------------
    // Resampling
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

      // set isotropic spacing
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
      itk::PluginFilterWatcher watcher(resampler, "Resample Volume",CLPProcessInformation);

      resampler->SetInterpolator( nearestNeighborInterpolator);
      resampler->SetInput( inputLabelVolume );
      resampler->SetTransform( transform );
      resampler->SetOutputOrigin( inputLabelVolume->GetOrigin() );
      resampler->SetOutputSpacing( outputSpacing );
      resampler->SetOutputDirection( inputLabelVolume->GetDirection() );
      resampler->SetSize( outputSize );
      resampler->Update();
      LabelVolumeType::Pointer inputLabelVolumeResampled= LabelVolumeType::New();
      inputLabelVolumeResampled = resampler->GetOutput();

    //----------------------------------------------------------------------------
    // Label Geometry filter
    LabelGeometryImageFilterType::Pointer labelGeometryImageFilter = LabelGeometryImageFilterType::New();
    labelGeometryImageFilter->SetInput(inputLabelVolumeResampled);
    labelGeometryImageFilter->CalculatePixelIndicesOn();
    labelGeometryImageFilter->CalculateOrientedLabelRegionsOn();
    labelGeometryImageFilter->CalculateOrientedBoundingBoxOn();
    labelGeometryImageFilter->Update();

    LabelGeometryImageFilterType::LabelsType allLabels = labelGeometryImageFilter->GetLabels();
    LabelGeometryImageFilterType::LabelsType::iterator allLabelsIt;
    std::cout << "Number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
    std::cout << std::endl;

    for( allLabelsIt = allLabels.begin()+1; allLabelsIt != allLabels.end(); allLabelsIt++ )
      {
      LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
      // we are supposing the whole gland is label 1
      if (labelValue == 1)
        {
        // get centroid
        LabelVolumeType::PointType centroid =  labelGeometryImageFilter->GetCentroid(labelValue);
        // get major axis length
        double majorAxisLength = labelGeometryImageFilter->GetMajorAxisLength(labelValue);
        // get eigenvectors
        MatrixType eigenVectors = labelGeometryImageFilter->GetEigenvectors(labelValue);
        itk::FixedArray<double, 6> boundingBox;
        boundingBox = labelGeometryImageFilter->GetBoundingBox(labelValue);
        //
        std::cout << "\tCentroid: " << centroid << std::endl;
        std::cout <<"\t major axis length: "<<majorAxisLength<<std::endl;
        std::cout << "\tEigenvectors: " << eigenVectors<< std::endl;
        std::cout << "\tBoundingBox: " << boundingBox<< std::endl;

        // calculate point 1
        LabelVolumeType::PointType point1;
        point1[0] = centroid[0] + eigenVectors(0,2)*majorAxisLength/6;
        point1[1] = centroid[1] + eigenVectors(1,2)*majorAxisLength/6;
        point1[2] = centroid[2] + eigenVectors(2,2)*majorAxisLength/6;
        std::cout << "\tPoint1: " << point1<< std::endl;

        // calculate point 2
        LabelVolumeType::PointType point2;
        point2[0] = centroid[0] - eigenVectors(0,2)*majorAxisLength/6;
        point2[1] = centroid[1] - eigenVectors(1,2)*majorAxisLength/6;
        point2[2] = centroid[2] - eigenVectors(2,2)*majorAxisLength/6;
        std::cout << "\tPoint2: " << point2<< std::endl;

        LabelVolumeType::Pointer outputLabelVolumeResampled= LabelVolumeType::New();
        outputLabelVolumeResampled->SetRegions(inputLabelVolumeResampled->GetLargestPossibleRegion());
        outputLabelVolumeResampled->Allocate();
        outputLabelVolumeResampled->FillBuffer(0);
        outputLabelVolumeResampled->CopyInformation(inputLabelVolumeResampled);

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

        itk::ImageRegionIteratorWithIndex<LabelVolumeType> imageIterator(outputLabelVolumeResampled,region);
        itk::ImageRegionIteratorWithIndex<LabelVolumeType> inputImageIterator(inputLabelVolumeResampled,region);
        LabelVolumeType::IndexType index;

        while(!imageIterator.IsAtEnd())
          {
            if (inputImageIterator.Value() == labelValue)
            {
              index = imageIterator.GetIndex();
              double d1 = eigenVectors(0,0) * centroid[0] + eigenVectors(1,0) * centroid[1] + eigenVectors(2,0)* centroid[2];
              double val1 = eigenVectors(0,0) * index[0] + eigenVectors(1,0) * index[1] + eigenVectors(2,0) * index[2] - d1;

              double d2 = eigenVectors(0,1) * centroid[0] + eigenVectors(1,1) * centroid[1] + eigenVectors(2,1)* centroid[2];
              double val2 = eigenVectors(0,1) * index[0] + eigenVectors(1,1) * index[1] + eigenVectors(2,1) * index[2] - d2;

              double d3 = eigenVectors(0,2) * point1[0] + eigenVectors(1,2) * point1[1] + eigenVectors(2,2)* point1[2];
              double val3 = eigenVectors(0,2) * index[0] + eigenVectors(1,2) * index[1] + eigenVectors(2,2) * index[2] - d3;

              double d4 = eigenVectors(0,2) * point2[0] + eigenVectors(1,2) * point2[1] + eigenVectors(2,2)* point2[2];
              double val4 = eigenVectors(0,2) * index[0] + eigenVectors(1,2) * index[1] + eigenVectors(2,2) * index[2] - d4;
             //
              if (val1 >0 && val2>0 && val4<0)
              {
                  imageIterator.Set(1);
              }

              else if (val1 >0 && val2<0 && val4<0)
              {
                  imageIterator.Set(2);
              }
              else if (val1 <0 && val2>0 && val4<0)
              {
                  imageIterator.Set(3);
              }
              else if (val1 <0 && val2<0 && val4<0)
              {
                  imageIterator.Set(4);
              }
               //
              else if (val1 >0 && val2>0 && val3<0 && val4>0)
              {
                  imageIterator.Set(5);
              }

              else if (val1 >0 && val2<0 && val3<0 && val4>0)
              {
                  imageIterator.Set(6);
              }
              else if (val1 <0 && val2>0 && val3<0 && val4>0)
              {
                  imageIterator.Set(7);
              }
              else if (val1 <0 && val2<0 && val3<0 && val4>0)
              {
                  imageIterator.Set(8);
              }
             //
             else  if (val1 >0 && val2>0 && val3>0)
              {
                  imageIterator.Set(9);
              }
              else if (val1 >0 && val2<0 && val3>0)
              {
                  imageIterator.Set(10);
              }
              else if (val1 <0 && val2>0 && val3>0)
              {
                  imageIterator.Set(11);
              }
              else if (val1 <0 && val2<0 && val3>0)
              {
                  imageIterator.Set(12);
              }
            }
          ++imageIterator;
          ++inputImageIterator;
          }
        //------------------------------------------------------------------------
        // Resampling back to the original spacing
        //------------------------------------------------------------------------
        const typename LabelVolumeType::SpacingType& outputSpacing=
            inputLabelVolume->GetSpacing();

        const typename LabelVolumeType::RegionType& outputRegion =
            inputLabelVolume->GetLargestPossibleRegion();

        const typename LabelVolumeType::SizeType& outputSize =
            outputRegion.GetSize();

        typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
        resampler->SetInterpolator( nearestNeighborInterpolator);
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
        itk::PluginFilterWatcher watcher2(resampler, "Resample Volume",CLPProcessInformation);
        //------------------------------------------------------------------------
        if(outputLabelName.size())
          {
          LabelVolumeWriterType::Pointer writer = LabelVolumeWriterType::New();
          writer->SetInput(outputLabelVolumeResampled);
          writer->SetFileName(outputLabelName.c_str());
          writer->SetUseCompression(1);
          writer->Update();
          }
        }
      }
    }
  }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
