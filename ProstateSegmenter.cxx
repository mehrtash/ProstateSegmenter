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

  typedef itk::Matrix<double, 3, 3> MatrixType;
  typedef itk::Vector<int, 3> VectorIntType;
  typedef itk::Vector<bool, 3> VectorBoolType;

void identifyEigenVectorsDirections(MatrixType & eigenVectors,MatrixType & labelDirection,VectorIntType& lps , VectorBoolType& flipped)
{
  // This functions finds eignenvector indices which are along the right-left (lIndex), anterior-posterior
  // (pIndex) and superior-inferior (sIndex) directions of prostate gland.

  typedef itk::Vector<double, 3> VectorType;
  // r is the Right vector (The first row of the direction cosines matrix of the label).
  VectorType r;
  r[0] = labelDirection(0,0);
  r[1] = labelDirection(0,1);
  r[2] = labelDirection(0,2);
  // a is the Anterior vector (The second row of the direction cosines matrix of the label).
  VectorType a;
  a[0] = labelDirection(1,0);
  a[1] = labelDirection(1,1);
  a[2] = labelDirection(1,2);
  // s is the Superior vector (The third row of the direction cosines matrix of the label).
  VectorType s;
  s[0] = labelDirection(2,0);
  s[1] = labelDirection(2,1);
  s[2] = labelDirection(2,2);
  // eTranspose is the eigenvectors matrix transpose
  MatrixType eTranspose(eigenVectors.GetTranspose());
  // we calculate the angle cosines between the eigenvectors and the r, a and s vectors
  VectorType rCosines(eTranspose*r);
  VectorType aCosines(eTranspose*a);
  VectorType sCosines(eTranspose*s);
  // find the maximum indices of the cosines and assign to lIndex, pIndex and sIndex
  double R[3] = {fabs(rCosines[0]), fabs(rCosines[1]), fabs(rCosines[2])};
  int lIndex, pIndex, sIndex;
  lIndex = std::distance(R, std::max_element(R, R+ 3));
  lps[0] = lIndex;
  if (rCosines[lIndex] < 0)
  {
    flipped[lIndex] = true;
  }
  double A[3] = {fabs(aCosines[0]), fabs(aCosines[1]), fabs(aCosines[2])};
  pIndex = std::distance(A, std::max_element(A, A+ 3));
  lps[1] = pIndex;
  if (aCosines[pIndex] < 0)
  {
    flipped[pIndex] = true;
  }
  double S[3] = {fabs(sCosines[0]), fabs(sCosines[1]), fabs(sCosines[2])};
  sIndex = std::distance(S, std::max_element(S, S+ 3));
  lps[2] = sIndex;
  if (sCosines[sIndex] < 0)
  {
    flipped[sIndex] = true;
  }
  std::cout<<"Left eigenVector index is: "<<lIndex<<std::endl;
  std::cout<<"Posterior eigenVector index is: "<<pIndex<<std::endl;
  std::cout<<"Superior eigenVector index is: "<<sIndex<<std::endl;
}
int smallest(int x, int y, int z){
    return std::min(std::min(x, y), z);
}
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
      LabelVolumeType::DirectionType labelDirection = inputLabelVolume->GetDirection();
      //std::cout<<"direction matrix is: "<<std::endl;
      //std::cout<<labelDirection<<std::endl;
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
      //itk::PluginFilterWatcher watcher(resampler, "Resample Volume",CLPProcessInformation);

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
    //std::cout << "Number of labels: " << labelGeometryImageFilter->GetNumberOfLabels() << std::endl;
    //std::cout << std::endl;

    for( allLabelsIt = allLabels.begin()+1; allLabelsIt != allLabels.end(); allLabelsIt++ )
      {
      LabelGeometryImageFilterType::LabelPixelType labelValue = *allLabelsIt;
      // we suppose that prostate gland label value is 1
      if (labelValue == 1)
        {
        // get centroid
        LabelVolumeType::PointType centroid =  labelGeometryImageFilter->GetCentroid(labelValue);
        //std::cout<<"centroid is: "<<centroid<<std::endl;
        // get major axis length
        itk::FixedArray<double, 3> axesLength;
        axesLength = labelGeometryImageFilter->GetAxesLength(labelValue);
        //std::cout<<"axis length is: "<<axesLength<<std::endl;
        // get eigenvectors
        MatrixType eigenVectors = labelGeometryImageFilter->GetEigenvectors(labelValue);
        std::cout<<"eigenvectors are : "<<std::endl;
        std::cout<<eigenVectors<<std::endl;
        itk::FixedArray<double, 6> boundingBox;
        boundingBox = labelGeometryImageFilter->GetBoundingBox(labelValue);
        //
        VectorIntType lps;
        VectorBoolType flipped;
        flipped[0] = false;
        flipped[1] = false;
        flipped[2] = false;

        identifyEigenVectorsDirections(eigenVectors,labelDirection,lps,flipped);
        //std::cout<<"flipped: "<<flipped<<std::endl;
        int lIndex = lps[0];
        int pIndex = lps[1];
        int sIndex = lps[2];


        double axisLength = axesLength[sIndex];
        //std::cout << "\tBoundingBox: " << boundingBox<< std::endl;

        // calculate point 1
        LabelVolumeType::PointType point1;
        point1[0] = centroid[0] + eigenVectors(0,sIndex)*axisLength/6;
        point1[1] = centroid[1] + eigenVectors(1,sIndex)*axisLength/6;
        point1[2] = centroid[2] + eigenVectors(2,sIndex)*axisLength/6;
        //std::cout << "\tPoint1: " << point1<< std::endl;

        // calculate point 2
        LabelVolumeType::PointType point2;
        point2[0] = centroid[0] - eigenVectors(0,sIndex)*axisLength/6;
        point2[1] = centroid[1] - eigenVectors(1,sIndex)*axisLength/6;
        point2[2] = centroid[2] - eigenVectors(2,sIndex)*axisLength/6;
        //std::cout << "\tPoint2: " << point2<< std::endl;

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
              double d1 = eigenVectors(0,lIndex) * centroid[0] + eigenVectors(1,lIndex) * centroid[1] + eigenVectors(2,lIndex)* centroid[2];
              // val1 is positive for left and negative for right.
              double val1 = eigenVectors(0,lIndex) * index[0] + eigenVectors(1,lIndex) * index[1] + eigenVectors(2,lIndex) * index[2] - d1;
              if (flipped[lIndex])
              {
                val1 *= -1;
              }
              double d2 = eigenVectors(0,pIndex) * centroid[0] + eigenVectors(1,pIndex) * centroid[1] + eigenVectors(2,pIndex)* centroid[2];
              // val2 is positive for posterior and negative for anterior.
              double val2 = eigenVectors(0,pIndex) * index[0] + eigenVectors(1,pIndex) * index[1] + eigenVectors(2,pIndex) * index[2] - d2;
              if (flipped[pIndex])
              {
                val2 *= -1;
              }

              double d3 = eigenVectors(0,sIndex) * point1[0] + eigenVectors(1,sIndex) * point1[1] + eigenVectors(2,sIndex)* point1[2];
              double val3 = eigenVectors(0,sIndex) * index[0] + eigenVectors(1,sIndex) * index[1] + eigenVectors(2,sIndex) * index[2] - d3;

              double d4 = eigenVectors(0,sIndex) * point2[0] + eigenVectors(1,sIndex) * point2[1] + eigenVectors(2,sIndex)* point2[2];
              double val4 = eigenVectors(0,sIndex) * index[0] + eigenVectors(1,sIndex) * index[1] + eigenVectors(2,sIndex) * index[2] - d4;
              if (flipped[sIndex])
              {
                val3 *= -1;
                val4 *= -1;
              }
             // label 14 Left Anterior Base
              if (val1 >0 && val2<0 && val4<0)
              {
                  imageIterator.Set(14);
              }

             // label 15 Left Anterior Middle
              else if (val1 >0 && val2<0 && val3<0 && val4>0)
              {
                  imageIterator.Set(15);
              }
             // label 16 Left Anterior Apex
              else if (val1 >0 && val2<0 && val3>0)
              {
                  imageIterator.Set(16);
              }
             // label 17 Left Posterior Base
              else if (val1 >0 && val2>0 && val4<0)
              {
                  imageIterator.Set(17);
              }
             // label 18 Left Posterior Middle
              else if (val1 >0 && val2>0 && val3<0 && val4>0)
              {
                  imageIterator.Set(18);
              }
             // label 19 Left Posterior Apex
              else if (val1 >0 && val2>0 && val3>0)
              {
                  imageIterator.Set(19);
              }
             // label 20 Right Anterior Base
              else if (val1 <0 && val2<0 && val4<0)
              {
                  imageIterator.Set(20);
              }
             // label 21 Right Anterior Middle
              else if (val1 <0 && val2<0 && val3<0 && val4>0)
              {
                  imageIterator.Set(21);
              }
             // label 22 Right Anterior Apex
             else  if (val1 <0 && val2<0 && val3>0)
              {
                  imageIterator.Set(22);
              }
             // label 23 Right Posterior Base
              else if (val1 <0 && val2>0 && val4<0)
              {
                  imageIterator.Set(23);
              }
             // label 24 Right Posterior Middle
              else if (val1 <0 && val2>0 && val3<0 && val4>0)
              {
                  imageIterator.Set(24);
              }
             // label 25 Right Posterior Apex
              else if (val1 <0 && val2>0 && val3>0)
              {
                  imageIterator.Set(25);
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
          writer->SetInput(outputLabelVolume);
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
