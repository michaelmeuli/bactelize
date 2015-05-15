/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/


#include "bactelize.h"

extern int cellchannel;
extern int bacteriachannel;
extern int lysosomechannel;
extern int binaryLowerThresholdBacteria;
extern int cellInsideThreshold;
extern double minNumberOfmm3;  
extern double maxNumberOfmm3;
extern int maxclustersize;
extern std::ofstream fileout;
extern std::string fileoutName;
extern int numberOfStreamDivisions;
extern int numberOfBins;
extern double minSearchRadius;


ImageType2D::Pointer maxintprojection(ImageType3D::Pointer inputImageMIP, unsigned int projectionDirection) {
  unsigned int i, j;
  unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i ) {
    if (i != projectionDirection) {
      direction[j] = i;
      j++;
      }
    }
  ImageType2D::RegionType region;
  ImageType2D::RegionType::SizeType size2DMIP;
  ImageType2D::RegionType::IndexType index2DMIP;
  ImageType3D::RegionType requestedRegion = inputImageMIP -> GetRequestedRegion();
  index2DMIP[ direction[0] ]    = requestedRegion.GetIndex()[ direction[0] ];
  index2DMIP[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size2DMIP[ direction[0] ]     = requestedRegion.GetSize()[  direction[0] ];
  size2DMIP[ 1- direction[0] ]  = requestedRegion.GetSize()[  direction[1] ];
  region.SetSize( size2DMIP );
  region.SetIndex( index2DMIP );
  ImageType2D::Pointer outputImageMIP = ImageType2D::New();
  outputImageMIP->SetRegions(region);
  outputImageMIP->Allocate();
 
  SliceIteratorTypeInput  inputIt(  inputImageMIP, inputImageMIP->GetRequestedRegion() );
  LinearIteratorTypeInput outputIt( outputImageMIP, outputImageMIP->GetRequestedRegion() );
  inputIt.SetFirstDirection(  direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );
  outputIt.GoToBegin();
  while (!outputIt.IsAtEnd()) {
    while (!outputIt.IsAtEndOfLine()) {
      outputIt.Set(itk::NumericTraits<unsigned short>::NonpositiveMin());
      ++outputIt;
      }
    outputIt.NextLine();
    }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while(!inputIt.IsAtEnd()) {
    while (!inputIt.IsAtEndOfSlice()) {
      while (!inputIt.IsAtEndOfLine()) {
        outputIt.Set(vnl_math_max(outputIt.Get(), inputIt.Get()));
        ++inputIt;
        ++outputIt;
        }
      outputIt.NextLine();
      inputIt.NextLine();
      }
    outputIt.GoToBegin();
    inputIt.NextSlice();
    }
  return outputImageMIP;
  }


BinaryImageType2D::Pointer maxintprojection(BinaryImageType3D::Pointer inputImageMIP, unsigned int projectionDirection) {
  unsigned int i, j;
  unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i ) {
    if (i != projectionDirection) {
      direction[j] = i;
      j++;
      }
    }
  BinaryImageType2D::RegionType region;
  BinaryImageType2D::RegionType::SizeType size2DMIP;
  BinaryImageType2D::RegionType::IndexType index2DMIP;
  BinaryImageType3D::RegionType requestedRegion = inputImageMIP -> GetRequestedRegion();
  index2DMIP[direction[0]]    = requestedRegion.GetIndex()[direction[0]];
  index2DMIP[1-direction[0]] = requestedRegion.GetIndex()[direction[1]];
  size2DMIP[direction[0]]     = requestedRegion.GetSize()[direction[0]];
  size2DMIP[1-direction[0]]  = requestedRegion.GetSize()[direction[1]];
  region.SetSize(size2DMIP);
  region.SetIndex(index2DMIP);
  BinaryImageType2D::Pointer outputImageMIP = BinaryImageType2D::New();
  outputImageMIP->SetRegions(region);
  outputImageMIP->Allocate();
 
  SliceIteratorTypeBinary  inputIt(  inputImageMIP, inputImageMIP->GetRequestedRegion() );
  LinearIteratorTypeBinary outputIt( outputImageMIP, outputImageMIP->GetRequestedRegion() );
  inputIt.SetFirstDirection(  direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );
  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() ) {
    while ( ! outputIt.IsAtEndOfLine() ) {
      outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
      ++outputIt;
      }
    outputIt.NextLine();
    }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while(!inputIt.IsAtEnd()) {
    while (!inputIt.IsAtEndOfSlice()) {
      while (!inputIt.IsAtEndOfLine()) {
        outputIt.Set(vnl_math_max(outputIt.Get(), inputIt.Get()));
        ++inputIt;
        ++outputIt;
        }
      outputIt.NextLine();
      inputIt.NextLine();
      }
    outputIt.GoToBegin();
    inputIt.NextSlice();
    }
  return outputImageMIP;
  }


void dumpmetadatadic(ImageType5D::Pointer image5D) {
  std::cout << "--== Metadata from dictionary ==--" << std::endl;
  image5D->Update();
  itk::MetaDataDictionary imgMetaDictionary = image5D->GetMetaDataDictionary();
  std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
  for (std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin(); itKey != imgMetaKeys.end(); ++itKey) {
    std::string tmp;
    itk::ExposeMetaData<std::string>( imgMetaDictionary, *itKey, tmp );
    std::cout << "\t" << *itKey << " ---> " << tmp << std::endl;
    }
  }


void dumpimageio(ReaderType::Pointer reader) {
  const itk::ImageIOBase * imageIO = reader->GetImageIO();
  itk::ImageIORegion region = imageIO->GetIORegion();
  int regionDim = region.GetImageDimension();
  std::cout << "--== Metadata from ImageIOBase ==--" << std::endl;
  for (int i = 0; i < regionDim; i++) {
    std::cout << "\tDimension " << i + 1 << " Size: "
              << region.GetSize(i) << std::endl;
    }
  for (int i = 0; i < regionDim; i++) {
    if (region.GetSize(i) > 1) {
      std::cout << "\tSpacing " << i + 1 << ": "
                << imageIO->GetSpacing(i) << std::endl;
      }
    }
  std::cout << "\tByte Order: "
            << imageIO->GetByteOrderAsString(imageIO->GetByteOrder())
            << std::endl;
  std::cout << "\tPixel Stride: " << imageIO->GetPixelStride() << std::endl;
  std::cout << "\tPixel Type: "
            << imageIO->GetPixelTypeAsString(imageIO->GetPixelType())
            << std::endl;
  std::cout << "\tImage Size (in pixels): "
            << imageIO->GetImageSizeInPixels() << std::endl;
  std::cout << "\tPixel Type: "
            << imageIO->GetComponentTypeAsString(imageIO->GetComponentType())
            << std::endl;
  std::cout << "\tRGB Channel Count: "
            << imageIO->GetNumberOfComponents() << std::endl;
  std::cout << "\tNumber of Dimensions: "
            << imageIO->GetNumberOfDimensions() << std::endl;
  }


void printImInfo(BinaryImageType3D::Pointer image3D) {
  ImageType3D::RegionType region3D = image3D->GetLargestPossibleRegion();
  int regionDim = region3D.GetImageDimension();
  const ImageType3D::SpacingType& sp = image3D->GetSpacing();
  const ImageType3D::PointType& orgn = image3D->GetOrigin();
  for (int i = 0; i < regionDim; i++) {
    std::cout << "Dimension " << i + 1 << " Size: "
              << region3D.GetSize(i) << std::endl;
    }
  for (int i = 0; i < regionDim; i++) {
    if ( region3D.GetSize(i) > 1 ) {
      std::cout << "Spacing " << i + 1 << ": "
                << std::setprecision(3) << sp[i] * 1000 << "um" << std::endl;
      }
    }
  for (int i = 0; i < regionDim; i++) {
    std::cout << "Origin " << i + 1 << ": "
              << orgn[i] << std::endl;
    }
  }


ImageType3D::Pointer extractchannel(ImageType5D::Pointer image5D, int channelnr) {
  // Get the 3D of second channell
  ExtractFilterType5D4D::Pointer extractfilter5D4Dch2 = ExtractFilterType5D4D::New();
  extractfilter5D4Dch2->InPlaceOn();
  extractfilter5D4Dch2->SetDirectionCollapseToSubmatrix();
  ImageType5D::RegionType region5D = image5D->GetLargestPossibleRegion();
  ImageType5D::SizeType size5Dch = region5D.GetSize();
  size5Dch[4] = 0; 
  ImageType5D::IndexType start5Dch = region5D.GetIndex();
  start5Dch[4] = channelnr;
  ImageType5D::RegionType region5Dch2;
  region5Dch2.SetSize(  size5Dch  );
  region5Dch2.SetIndex( start5Dch );
  extractfilter5D4Dch2->SetExtractionRegion( region5Dch2 );
  extractfilter5D4Dch2->SetInput( image5D );

  ExtractFilterType4D3D::Pointer extractfilter4D3Dt1 = ExtractFilterType4D3D::New();
  extractfilter4D3Dt1->InPlaceOn();
  extractfilter4D3Dt1->SetDirectionCollapseToSubmatrix();
  extractfilter5D4Dch2->Update();
  ImageType4D::RegionType region4D = extractfilter5D4Dch2->GetOutput()->GetLargestPossibleRegion();
  ImageType4D::SizeType size4Dt1 = region4D.GetSize();
  size4Dt1[3] = 0;
  ImageType4D::IndexType start4Dt1 = region4D.GetIndex();
  start4Dt1[3] = 0;
  ImageType4D::RegionType region4Dt1;
  region4Dt1.SetSize(  size4Dt1  );
  region4Dt1.SetIndex( start4Dt1 );
  extractfilter4D3Dt1->SetExtractionRegion( region4Dt1 );
  extractfilter4D3Dt1->SetInput( extractfilter5D4Dch2->GetOutput() );
  extractfilter4D3Dt1->Update();
  return extractfilter4D3Dt1->GetOutput();
  }


void printHistogramNormalized(ImageType3D::Pointer inputImage) {
  NormalizeFilterType::Pointer  normalizeFilter = NormalizeFilterType::New();
  normalizeFilter->SetInput(inputImage);
  RescaleFilterTypeNormalized::Pointer rescaleNormalized = RescaleFilterTypeNormalized::New();
  rescaleNormalized->SetInput( normalizeFilter->GetOutput() ); 
  rescaleNormalized->SetOutputMinimum(0);
  rescaleNormalized->SetOutputMaximum(4095);
  HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  typedef HistogramFilterType::HistogramSizeType   SizeType;
  SizeType size(1);
  size[0] =  numberOfBins;     
  histogramFilter->SetHistogramSize( size );
  histogramFilter->SetMarginalScale( 10.0 ); 
  HistogramFilterType::HistogramMeasurementVectorType lowerBound( 1 );
  HistogramFilterType::HistogramMeasurementVectorType upperBound( 1 );
  lowerBound[0] = 0;
  upperBound[0] = 4095;
  histogramFilter->SetHistogramBinMinimum( lowerBound );
  histogramFilter->SetHistogramBinMaximum( upperBound ); 
  histogramFilter->SetInput(rescaleNormalized->GetOutput());
  histogramFilter->Update();
  const HistogramType * histogram = histogramFilter->GetOutput();
  const unsigned int histogramSize = histogram->Size();
  std::cout << "Histogram size " << histogramSize << std::endl;
  for( unsigned int bin=0; bin < histogramSize; bin++ )
    {
    std::cout << "bin = " << std::setw(3) << bin << 
      "        measurement = " << std::setw(10) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, 0) <<
      "        frequency = " << std::setw(10) << histogram->GetFrequency( bin, 0 ) << std::endl;    
    }    
  }


void write2D(ImageType2D::Pointer image2Dbacteria, std::string filenamepath) {
  RescaleFilterTypeWriter::Pointer rescaleFilter = RescaleFilterTypeWriter::New();
  rescaleFilter->SetInput( image2Dbacteria ); 
  rescaleFilter->SetOutputMinimum( 0 );
  rescaleFilter->SetOutputMaximum( 255 );
  WriterType::Pointer writer = WriterType::New();
  TIFFIOType::Pointer tiffIO = TIFFIOType::New();
  writer->SetImageIO(tiffIO);
  writer->SetFileName( filenamepath );             
  writer->SetInput(rescaleFilter->GetOutput());
  writer->Update();
  }


void write2D(BinaryImageType2D::Pointer image2Dbacteria, std::string filenamepath) {
  WriterType::Pointer writer = WriterType::New();
  TIFFIOType::Pointer tiffIO = TIFFIOType::New();
  writer->SetImageIO(tiffIO);
  writer->SetFileName( filenamepath );             
  writer->SetInput(image2Dbacteria);
  writer->Update();
  }


std::string getFilename(std::string inputFileName, int seriesnr, int seriesCount, std::string suffix) {
  std::string inputfilename = itksys::SystemTools::GetFilenameName(inputFileName);   
  std::string filenamebase = inputfilename.substr(0, inputfilename.find_first_of('.'));
  int sigFigs = 0;
  while (seriesCount >= 10) {
    seriesCount /= 10;
    sigFigs++;
    }
  std::string filename = filenamebase;
  std::stringstream ssout;
  ssout << '_';
  int currentSigFigs = 0;
  int currentSeries = seriesnr;
  while (currentSeries >= 10) {
    currentSeries /= 10;
    currentSigFigs++;
    }
  for (int i=0; i<(sigFigs - currentSigFigs); i++) {
    ssout << 0;
    }
  ssout << seriesnr;
  ssout << suffix;
  filename.append(ssout.str());
  return filename;
  }


ImageSizeType getImSize(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter) {
  binaryImageToShapeLabelMapFilter->Update();
  BinaryImageToShapeLabelMapFilterType::OutputImageType::RegionType region = binaryImageToShapeLabelMapFilter->GetOutput()->GetLargestPossibleRegion();
  const BinaryImageToShapeLabelMapFilterType::OutputImageType::SpacingType& spacing = binaryImageToShapeLabelMapFilter->GetOutput()->GetSpacing();
  ImageSizeType imSize;
  for (unsigned int i = 0; i < imSize.Size(); i++) {
    imSize[i] = spacing[i] * region.GetSize(i);
    }
  return imSize;
  }


void printObjectLabelAndCentroid(BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject, ImageSizeType imSize) {
  BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();   
  std::cout << "Object with label " << static_cast<int>(labelObject->GetLabel()) << " and centroid:  ";
  for (unsigned int j = 0; j < imSize.Size(); j++) {
    std::cout << centroid[j] / imSize[j] * 100 << "%  ";
    }
  std::cout << std::endl;
  }


void printObjectInfo(BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject) {
  std::cout << "NumberOfPixels: " << labelObject->GetNumberOfPixels() << std::endl;
  std::cout << "PhysicalSize: " << (labelObject->GetPhysicalSize()) * 1000000000 << "um3" << std::endl;
  std::cout << "Elongation: " << labelObject->GetElongation() << std::endl;
  std::cout << "Roundness: " << labelObject->GetRoundness() << std::endl;
  std::cout << "EquivalentEllipsoidDiameter: " << (labelObject->GetEquivalentEllipsoidDiameter() * 1000) << "um" << std::endl;
  std::cout << "EquivalentSphericalRadius: " << (labelObject->GetEquivalentSphericalRadius() * 1000) << "um" << std::endl;
  }


void printShapeLabelObjects(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter) {
  binaryImageToShapeLabelMapFilter->Update();
  ImageSizeType imSize = getImSize(binaryImageToShapeLabelMapFilter);
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();      
    printObjectLabelAndCentroid(labelObject, imSize);  
    printObjectInfo(labelObject);
    std::cout << std::endl;
    }    
  }


SampleType::Pointer getCentroidsAsSample(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter) {
  binaryImageToShapeLabelMapFilter->Update();
  SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize(3);
  MeasurementVectorType mv;
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();     
    mv[0] = centroid[0];
    mv[1] = centroid[1];
    mv[2] = centroid[2];
    sample->PushBack(mv);
    } 
  return sample;
  }


void excludeIfSet(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, std::set<int> labelsToRemove) {
  std::set<int>::const_iterator itr;
  for (itr = labelsToRemove.begin(); itr != labelsToRemove.end(); ++itr) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetLabelObject(*itr);
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();   
    ImageSizeType imSize = getImSize(binaryImageToShapeLabelMapFilter);
    std::cout << "Removing:  ";
    printObjectLabelAndCentroid(labelObject, imSize);   
    binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel(*itr);
    }
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects remaining in shape label map." << std::endl;
  std::cout << std::endl;
  binaryImageToShapeLabelMapFilter->Update();
}


BinaryImageType3D::Pointer getBinaryIm(ImageType3D::Pointer image3Dbacteria) {
  NormalizeFilterType::Pointer  normalizeFilter = NormalizeFilterType::New();
  normalizeFilter->SetInput( image3Dbacteria );
  RescaleFilterTypeNormalized::Pointer rescaleNormalized = RescaleFilterTypeNormalized::New();
  rescaleNormalized->SetInput( normalizeFilter->GetOutput() ); 
  rescaleNormalized->SetOutputMinimum(0);
  rescaleNormalized->SetOutputMaximum(4095);
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  MedianFilterType::InputSizeType radius; 
  radius.Fill(1);
  medianFilter->SetRadius(radius);
  medianFilter->SetInput( rescaleNormalized->GetOutput() );
  BinaryFilterType::Pointer binaryfilter = BinaryFilterType::New();
  binaryfilter->SetInput( medianFilter->GetOutput() );
  binaryfilter->SetOutsideValue(0);
  binaryfilter->SetInsideValue(255);
  binaryfilter->SetLowerThreshold(binaryLowerThresholdBacteria);
  binaryfilter->Update();
  BinaryImageType3D::Pointer binaryimage3Dbacteria = binaryfilter->GetOutput(); 
  return binaryimage3Dbacteria;
  }


void excludeSmallObjects(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, double minNumberOfmm3) {
  binaryImageToShapeLabelMapFilter->Update();
  std::vector<unsigned long> labelsToRemove;
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    if (labelObject->GetPhysicalSize() < minNumberOfmm3) {       
      labelsToRemove.push_back(labelObject->GetLabel());
      }
    }    
  std::cout << "Removing " << labelsToRemove.size() << " objects from shape label map." << std::endl;
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i) {
    binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
    }
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects remaining in shape label map." << std::endl;
  }
  

void excludeClusters(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, int maxclustersize) {
  binaryImageToShapeLabelMapFilter->Update();
  SampleType::Pointer sample = getCentroidsAsSample(binaryImageToShapeLabelMapFilter);
  TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample(sample);
  treeGenerator->SetBucketSize(16);
  treeGenerator->Update();
  TreeType::Pointer tree = treeGenerator->GetOutput();
  MeasurementVectorType queryPoint;
  DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
  DistanceMetricType::OriginType origin(3);

  std::set<int> labelsToRemove;
  ImageSizeType imSize = getImSize(binaryImageToShapeLabelMapFilter);
  for (unsigned int i = 0; i < sample->Size(); i++) {
    queryPoint = sample->GetMeasurementVector(i);     
    for ( unsigned int j = 0; j < sample->GetMeasurementVectorSize(); ++j ) {
      origin[j] = queryPoint[j];
      }
    distanceMetric->SetOrigin(origin);  
    TreeType::InstanceIdentifierVectorType neighbors;

    double maxDiameter = 0.0;  
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::VectorType ellipsoidVector = labelObject->GetEquivalentEllipsoidDiameter();
    for ( unsigned int vd = 0; vd < ellipsoidVector.GetVectorDimension(); ++vd ) {
      if (maxDiameter < ellipsoidVector[vd]) {
        maxDiameter = ellipsoidVector[vd];
        }
      }
    double radius = minSearchRadius;   
    if (maxDiameter/2 > minSearchRadius) {
      radius = maxDiameter/2;
      } 

    tree->Search(queryPoint, radius, neighbors); 
    if ((neighbors.size() > maxclustersize) || (labelObject->GetPhysicalSize() > maxNumberOfmm3)) {
      std::cout << "Centroid of object with label " << static_cast<int>(labelObject->GetLabel()) << " is set as query point." << std::endl;
      if (radius > minSearchRadius) {
        std::cout << "Search radius was increased to " << radius * 1000 << "um" << std::endl;
        }
      if (labelObject->GetPhysicalSize() > maxNumberOfmm3) {
        std::cout << "Object with label " << static_cast<int>(labelObject->GetLabel()) << " is too big to be one bacteria: " 
                  << (labelObject->GetPhysicalSize()) * 1000000000 << "um3" << std::endl;
        }
      if (neighbors.size() > maxclustersize) {
        std::cout << "Maximum cluster size was exceeded: " << neighbors.size() << " bacteria close to each other." << std::endl;
        }
      std::cout << "Object(s) which are set to be removed in this query:" << std::endl;
      for (unsigned int k = 0; k < neighbors.size(); ++k) {  
        BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObjectNeighbor = 
	    				binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(neighbors[k]);
        printObjectLabelAndCentroid(labelObjectNeighbor, imSize);   
        std::cout << "  Distance to query point: "; 
        std::cout << distanceMetric->Evaluate(tree->GetMeasurementVector(neighbors[k])) * 1000 << "um" << std::endl;
        labelsToRemove.insert(static_cast<int>(labelObjectNeighbor->GetLabel()));         
        }
      std::cout << std::endl;
      }    
    }  
  excludeIfSet(binaryImageToShapeLabelMapFilter, labelsToRemove);
  }





void excludeOutside(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, 
	            ImageType3D::Pointer image3Dcell) {
  binaryImageToShapeLabelMapFilter->Update();
  std::vector<unsigned long> labelsToRemove;
  StatisticsLabelMapFilterType::Pointer statisticsLabelMapFilter = StatisticsLabelMapFilterType::New();
  statisticsLabelMapFilter->SetInput1(binaryImageToShapeLabelMapFilter->GetOutput());
  statisticsLabelMapFilter->SetInput2(image3Dcell);
  statisticsLabelMapFilter->InPlaceOn();
  statisticsLabelMapFilter->Update();

  ImageSizeType imSize = getImSize(binaryImageToShapeLabelMapFilter);
  for(unsigned int i = 0; i < statisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    StatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = 
      statisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    printObjectLabelAndCentroid(labelObject, imSize);   
    std::cout << "  Mean value of object with label " << static_cast<int>(labelObject->GetLabel()) << " in cellchannel: " 
              << labelObject->GetMean() << std::endl; 
    if (labelObject->GetMean() < cellInsideThreshold) { 
      printObjectLabelAndCentroid(labelObject, imSize); 
      std::cout << "  Object with label " << static_cast<int>(labelObject->GetLabel()) << " will be excluded." << std::endl;  
      std::cout << "    Mean value in cellchannel is: " << labelObject->GetMean() << std::endl; 
      std::cout << "    cellInsideThreshold is: " << cellInsideThreshold << std::endl; 
      printObjectLabelAndCentroid(labelObject, imSize); 
      labelsToRemove.push_back(labelObject->GetLabel());
      }
    }
  std::cout << "Removing " << labelsToRemove.size() << " objects from shape label map." << std::endl;
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i) {
    binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
    }
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() 
            << " objects remaining in shape label map." << std::endl;
  }





void writeMeanResults(StatisticsLabelMapFilterType::Pointer statisticsLabelMapFilter, 
                      BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, 
	              ImageType3D::Pointer image3Dred, 
                      std::string fulloutputfilenameResults, std::string seriesName) {
  statisticsLabelMapFilter->SetInput1(binaryImageToShapeLabelMapFilter->GetOutput());
  statisticsLabelMapFilter->SetInput2(image3Dred);
  statisticsLabelMapFilter->InPlaceOn();
  statisticsLabelMapFilter->Update();

  fileout.open(fulloutputfilenameResults.c_str(), std::ofstream::app); 
  fileout << seriesName << "\t";
  unsigned int bacteriacount = statisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
  for(unsigned int i = 0; i < statisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    StatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObjectMe = 
      statisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    double mean = labelObjectMe->GetMean();
    fileout << mean << "\t";
    std::cout << "Mean value of object with label " << static_cast<int>(labelObjectMe->GetLabel()) << " in lysosomechannel: " 
              << mean << std::endl; 
    }
  fileout << "\n";
  fileout.close();
  std::cout << "Total bacteria counted (in statisticsLabelMapFilter): " << bacteriacount << std::endl;
  std::cout << std::endl;  
  }




int processSeries(std::string inputFileName, std::string outputdirectory, bool vflag, bool tflag, int fileNr, int seriesNr) {
  std::cout << "Processing series of file: " << inputFileName << std::endl;
  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  io->DebugOn();
  ReaderType::Pointer reader = ReaderType::New();
  std::cout << "reader->GetUseStreaming(): " << reader->GetUseStreaming() << std::endl;
  std::cout << "done checking streaming usage" << std::endl;
  reader->SetImageIO(io);
  reader->SetFileName(inputFileName);
  StreamingFilter::Pointer streamer = StreamingFilter::New();
  streamer->SetInput(reader->GetOutput());
  streamer->SetNumberOfStreamDivisions(numberOfStreamDivisions);
  ImageType5D::Pointer image5D =  ImageType5D::New();

  typedef itk::ImageFileWriter< ImageType5D > WriterType5D;
  typename WriterType5D::Pointer writer;
  writer = WriterType5D::New();
  writer->SetInput( streamer->GetOutput() );
  itk::SCIFIOImageIO::Pointer ioOut = itk::SCIFIOImageIO::New();
  ioOut->DebugOn();
  writer->SetImageIO( ioOut );

  reader->UpdateOutputInformation();
  io->SetSeries(0);
  reader->Modified(); 
  int seriesEnd = io->GetSeriesCount();
  std::cout << "seriesEnd: " << seriesEnd << std::endl;  
  std::cout << std::endl;

  int seriesnr=0;
  if (tflag) {
    seriesnr = seriesNr;
    seriesEnd = seriesNr+1;
    io->SetSeries(seriesnr);
    reader->Modified();
    }
  while ( seriesnr < seriesEnd ) {
    std::string outputfilenameSeries = getFilename(inputFileName, seriesnr, seriesEnd, ".ome.tiff");
    std::string fulloutputfilenameSeries = outputdirectory + outputfilenameSeries; 
    std::cout << "Writing file: " << fulloutputfilenameSeries << " ..." << std::endl;
    writer->SetFileName( fulloutputfilenameSeries );
    writer->Update();
    std::cout << "Reading seriesnr: " << seriesnr << std::endl;
    image5D = streamer->GetOutput();

    ImageType3D::Pointer image3Dcell = extractchannel(image5D, cellchannel);
    ImageType3D::Pointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
    ImageType3D::Pointer image3Dred = extractchannel(image5D, lysosomechannel);

    ImageType2D::Pointer image2Dcell = maxintprojection(image3Dcell);
    std::string outputfilename2Dcell = getFilename(inputFileName, seriesnr, seriesEnd, "_a_cell.tiff");
    std::string fulloutputfilename2Dcell = outputdirectory + outputfilename2Dcell; 
    std::cout << "Writing file: " << fulloutputfilename2Dcell << " ..." << std::endl;
    write2D(image2Dcell, fulloutputfilename2Dcell);

    ImageType2D::Pointer image2Dbacteria = maxintprojection(image3Dbacteria);
    std::string outputfilename2Dbacteria = getFilename(inputFileName, seriesnr, seriesEnd, "_b_bacteria.tiff");
    std::string fulloutputfilename2Dbacteria = outputdirectory + outputfilename2Dbacteria; 
    std::cout << "Writing file: " << fulloutputfilename2Dbacteria << " ..." << std::endl;
    write2D(image2Dbacteria, fulloutputfilename2Dbacteria);
  
    ImageType2D::Pointer image2Dred = maxintprojection(image3Dred);
    std::string outputfilename2Dred = getFilename(inputFileName,seriesnr, seriesEnd, "_c_lysosome.tiff");
    std::string fulloutputfilename2Dred = outputdirectory + outputfilename2Dred; 
    std::cout << "Writing file: " << fulloutputfilename2Dred << " ..." << std::endl;
    write2D(image2Dred, fulloutputfilename2Dred);
  
    BinaryImageType3D::Pointer binaryimage3Dbacteria = getBinaryIm(image3Dbacteria);

    if (vflag) {
      std::cout << std::endl;
      dumpimageio(reader);
      std::cout << std::endl;
      dumpmetadatadic(image5D);
      std::cout << std::endl;
      printHistogramNormalized(image3Dbacteria);
      std::cout << std::endl;
      printImInfo(binaryimage3Dbacteria);
      }

    BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
    binaryImageToShapeLabelMapFilter->SetInput(binaryimage3Dbacteria);
    binaryImageToShapeLabelMapFilter->Update();
    std::cout << std::endl;
    std::cout << "Number of objects originally found in binaryImageToShapeLabelMapFilter: ";
    std::cout << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects." << std::endl;

    std::cout << "Writing 2D of original binaryImageToShapeLabelMapFilter..." << std::endl;
    LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImageFilter = LabelMapToBinaryImageFilterType::New();
    labelMapToBinaryImageFilter->SetInput(binaryImageToShapeLabelMapFilter->GetOutput());
    labelMapToBinaryImageFilter->Update();
    BinaryImageType2D::Pointer binaryimage2Dbacteria = maxintprojection(labelMapToBinaryImageFilter->GetOutput());
    std::string outputfilenamebinary2Dbacteria = getFilename(inputFileName, seriesnr, seriesEnd, "_d_bacteria_binary.tiff");
    std::string fulloutputfilenamebinary2Dbacteria = outputdirectory + outputfilenamebinary2Dbacteria; 
    std::cout << "Writing file: " << fulloutputfilenamebinary2Dbacteria << " ..." << std::endl;
    write2D(binaryimage2Dbacteria, fulloutputfilenamebinary2Dbacteria);
    std::cout << std::endl;
    std::cout << "Removing objects which are less than " << minNumberOfmm3 * 1000000000 << " um3..." << std::endl;
    excludeSmallObjects(binaryImageToShapeLabelMapFilter, minNumberOfmm3);
    std::cout << std::endl;
    std::cout << "Removing objects outside of macrophages..." << std::endl;
    excludeOutside(binaryImageToShapeLabelMapFilter, image3Dcell);

    std::cout << std::endl;
    if (vflag) {
      std::cout << "Objects remaining after exluding small objects:" << std::endl;
      printShapeLabelObjects(binaryImageToShapeLabelMapFilter);
      }
    std::cout << "Removing clusters and clumps..." << std::endl;
    excludeClusters(binaryImageToShapeLabelMapFilter, maxclustersize);

    std::string fulloutputfilenameResults = outputdirectory + fileoutName;
    std::string seriesName = getFilename(inputFileName, seriesnr, seriesEnd);
    StatisticsLabelMapFilterType::Pointer statisticsLabelMapFilterLysosome = StatisticsLabelMapFilterType::New();
    writeMeanResults(statisticsLabelMapFilterLysosome, binaryImageToShapeLabelMapFilter, image3Dred, 
                     fulloutputfilenameResults, seriesName);
   
    L2ImageType::Pointer labelMapToReadMeanImage = L2ImageType::New();
    labelMapToReadMeanImage->SetInput(statisticsLabelMapFilterLysosome->GetOutput());
    labelMapToReadMeanImage->Update();
    ImageType3D::Pointer image3DbacteriaReadMean = labelMapToReadMeanImage->GetOutput();
    ImageType2D::Pointer image2DReadMean = maxintprojection(image3DbacteriaReadMean);
    RGBFilterType::Pointer colormapImageFilter = RGBFilterType::New();
    colormapImageFilter->SetInput(image2DReadMean);
    colormapImageFilter->SetColormap(RGBFilterType::Jet);
    colormapImageFilter->Update();
    std::string outputfilename2DReadMean = getFilename(inputFileName, seriesnr, seriesEnd, "_e_ReadMean.tiff");
    std::string fulloutputfilename2DReadMean = outputdirectory + outputfilename2DReadMean; 
    WriterTypeRGB::Pointer writerRGB = WriterTypeRGB::New();
    TIFFIOType::Pointer tiffIO = TIFFIOType::New();
    writerRGB->SetImageIO(tiffIO);
    writerRGB->SetFileName(fulloutputfilename2DReadMean);             
    writerRGB->SetInput(colormapImageFilter->GetOutput());
    std::cout << "Writing file: " << fulloutputfilename2DReadMean << " ..." << std::endl;
    writerRGB->Update();

    seriesnr++;
    if (seriesnr < seriesEnd) {
      io->SetSeries(seriesnr);
      reader->Modified();
      }
    std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}






