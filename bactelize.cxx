#include "bactelize.h"

extern int nucleichannel;
extern int bacteriachannel;
extern int lysosomechannel;
extern int binaryLowerThresholdBacteria;
extern double minNumberOfmm3;  
extern double maxNumberOfmm3;
extern int maxclustersize;
extern std::ofstream fileout;
extern std::string fileoutName;
extern int numberOfStreamDivisions;
extern int numberOfBins;
extern double maxSingleObjectDiameter;
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


void printObjectInfo(BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject) {
  std::cout << "Label: " << static_cast<int>(labelObject->GetLabel()) << std::endl; 
  std::cout << "NumberOfPixels: " << labelObject->GetNumberOfPixels() << std::endl;
  std::cout << "PhysicalSize: " << (labelObject->GetPhysicalSize()) * 1000000000 << "um3" << std::endl;
  std::cout << "Elongation: " << labelObject->GetElongation() << std::endl;
  std::cout << "Roundness: " << labelObject->GetRoundness() << std::endl;
  std::cout << "EquivalentEllipsoidDiameter: " << (labelObject->GetEquivalentEllipsoidDiameter() * 1000) << "um" << std::endl;
  std::cout << "EquivalentSphericalRadius: " << (labelObject->GetEquivalentSphericalRadius() * 1000) << "um" << std::endl;
  }


void printShapeLabelObjects(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter) {
  binaryImageToShapeLabelMapFilter->Update();
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    printObjectInfo(labelObject);
    std::cout << std::endl;
    }    
  }

void printShapeLabelObjects(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, ImageSizeType imSize) {
  binaryImageToShapeLabelMapFilter->Update();
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();      
    std::cout << "Centroid in %:  ";
    for (unsigned int i = 0; i < imSize.Size(); i++) {
      std::cout << centroid[i] / imSize[i] * 100 << "  ";
      }
    std::cout << std::endl;
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


void printCentroids(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, ImageSizeType imSize) {
  binaryImageToShapeLabelMapFilter->Update();
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();      
    std::cout << "Centroid of object with label " << static_cast<int>(labelObject->GetLabel()) << " in %:  ";
    for (unsigned int i = 0; i < imSize.Size(); i++) {
      std::cout << centroid[i] / imSize[i] * 100 << "  ";
      }
    std::cout << std::endl;
    } 
  }


void printSampleVectors(SampleType::Pointer sample, ImageSizeType imSize) {       
  for (unsigned int i = 0; i < sample->Size(); i++) {
    MeasurementVectorType queryPoint = sample->GetMeasurementVector(i);
    std::cout << "MeasurementVector of sample i=" << i << " in %:  "; 
    for (unsigned int i = 0; i < imSize.Size(); i++) {
      std::cout << queryPoint[i] / imSize[i] * 100 << "  ";
      }
    std::cout << std::endl; 
    }
  }


void excludeIfSet(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, std::set<int> setToRemove) {
  std::vector<unsigned long> labelsToRemove;
  std::set<int>::const_iterator itr;
  std::cout << "Labels which are set to be removed: ";
  if (setToRemove.empty()) {
    std::cout << "(nothing set to be removed)" << std::endl;
    }
  for (itr = setToRemove.begin(); itr != setToRemove.end(); ++itr) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(*itr);
    std::cout << static_cast<int>(labelObject->GetLabel()) << "  "; 
    labelsToRemove.push_back(labelObject->GetLabel());
    }   
  std::cout << std::endl;
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " shape label map objects." << std::endl;
  std::cout << "Removing " << labelsToRemove.size() << " objects from shape label map." << std::endl;
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i) {
    binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
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
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " shape label map objects." << std::endl;
  std::cout << "Removing " << labelsToRemove.size() << " objects from shape label map." << std::endl;
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i) {
    binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
    }
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects remaining in shape label map." << std::endl << std::endl;
  }
  

void excludeClusters(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, int maxclustersize) {
  binaryImageToShapeLabelMapFilter->Update();
  std::cout << "Adding centroids to sample..." << std::endl;
  SampleType::Pointer sample = getCentroidsAsSample(binaryImageToShapeLabelMapFilter);

  ImageSizeType imSize = getImSize(binaryImageToShapeLabelMapFilter);
  std::cout << "Printing centroids in binaryImageToShapeLabelMapFilter..." << std::endl;
  printCentroids(binaryImageToShapeLabelMapFilter, imSize);
  std::cout << std::endl;
  std::cout << "Printing MeasurementVectors in sample..." << std::endl;
  printSampleVectors(sample, imSize);
  std::cout << std::endl;

  TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
  treeGenerator->SetSample(sample);
  treeGenerator->SetBucketSize(16);
  treeGenerator->Update();
  TreeType::Pointer tree = treeGenerator->GetOutput();
  MeasurementVectorType queryPoint;
  DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
  DistanceMetricType::OriginType origin(3);

  std::set<int> setToRemove;
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
    if (maxDiameter > maxSingleObjectDiameter) {
      std::cout << "Size of object is too big to be one bacteria as " << std::endl; 
      std::cout << " the maximum of the equivalent ellipsoid diameter is: " << maxDiameter * 1000 << "um" << std::endl;
      if (maxDiameter/2 > minSearchRadius) {
        std::cout << "Increasing search radius to " << maxDiameter/2 * 1000 << "um" << std::endl;
        radius = maxDiameter/2;
        }
      }

    tree->Search(queryPoint, radius, neighbors); 
    std::cout << "kd-tree radius search result:" << std::endl;
    std::cout << "query point is object labeled " << static_cast<int>(labelObject->GetLabel()) << std::endl;   
    std::cout << "search radius = " << radius * 1000 << "um" << std::endl;
    if ((neighbors.size() > maxclustersize) || (maxDiameter/2 > minSearchRadius)) {   
      std::cout << "Neighbors which are set to be removed in this query:" << std::endl;
      for (unsigned int k = 0; k < neighbors.size(); ++k) {       
        BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObjectNeighbor = 
	  					binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(neighbors[k]);
        std::cout << "  Neighbor with label " << static_cast<int>(labelObjectNeighbor->GetLabel()) << " and a distance of: "; 
        std::cout << distanceMetric->Evaluate(tree->GetMeasurementVector(neighbors[k])) * 1000 << "um" << std::endl;
        setToRemove.insert(neighbors[k]);         
        }
      }
    else {
      std::cout << "Neighbors which are not set to be removed in this query:" << std::endl;
      for (unsigned int k = 0; k < neighbors.size(); ++k) {           
        BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObjectNeighbor = 
	  					binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(neighbors[k]);
        std::cout << "  Neighbor with label " << static_cast<int>(labelObjectNeighbor->GetLabel()) << " and a distance of: "; 
        std::cout << distanceMetric->Evaluate(tree->GetMeasurementVector(neighbors[k])) * 1000 << "um" << std::endl;
        }
      }
    std::cout << std::endl;    
    }  
  excludeIfSet(binaryImageToShapeLabelMapFilter, setToRemove);
  }


void excludeLargeObjects(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, double maxNumberOfmm3) {
  binaryImageToShapeLabelMapFilter->Update();
  std::vector<unsigned long> labelsToRemove;
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    if (labelObject->GetPhysicalSize() > maxNumberOfmm3) {       
      labelsToRemove.push_back(labelObject->GetLabel());
      }
    }    
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " shape label map objects." << std::endl;
  std::cout << "Removing " << labelsToRemove.size() << " objects from shape label map." << std::endl;
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i) {
    binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
    }
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects remaining in shape label map." << std::endl << std::endl;
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

    ImageType3D::Pointer image3Dnuclei = extractchannel(image5D, nucleichannel);
    ImageType3D::Pointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
    ImageType3D::Pointer image3Dred = extractchannel(image5D, lysosomechannel);

    ImageType2D::Pointer image2Dnuclei = maxintprojection(image3Dnuclei);
    std::string outputfilename2Dnuclei = getFilename(inputFileName, seriesnr, seriesEnd, "_a_nuclei.tiff");
    std::string fulloutputfilename2Dnuclei = outputdirectory + outputfilename2Dnuclei; 
    std::cout << "Writing file: " << fulloutputfilename2Dnuclei << " ..." << std::endl;
    write2D(image2Dnuclei, fulloutputfilename2Dnuclei);

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
      std::cout << std::endl;
      }

    std::string fulloutputfilenameResults = outputdirectory + fileoutName;
    fileout.open(fulloutputfilenameResults.c_str(), std::ofstream::app); 
    std::string seriesName = getFilename(inputFileName, seriesnr, seriesEnd);
    fileout << seriesName << "\t";

    BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
    binaryImageToShapeLabelMapFilter->SetInput(binaryimage3Dbacteria);
    LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImageFilter = LabelMapToBinaryImageFilterType::New();

    std::cout << std::endl;
    std::cout << "Objects originally found in binaryImageToShapeLabelMapFilter..." << std::endl;
    ImageSizeType imSize = getImSize(binaryImageToShapeLabelMapFilter);
    printShapeLabelObjects(binaryImageToShapeLabelMapFilter, imSize);
    std::cout << "Writing 2D of original binaryImageToShapeLabelMapFilter..." << std::endl;
    labelMapToBinaryImageFilter->SetInput(binaryImageToShapeLabelMapFilter->GetOutput());
    labelMapToBinaryImageFilter->Update();
    BinaryImageType2D::Pointer binaryimage2Dbacteria = maxintprojection(labelMapToBinaryImageFilter->GetOutput());
    std::string outputfilenamebinary2Dbacteria = getFilename(inputFileName, seriesnr, seriesEnd, "_d_bacteria_binary.tiff");
    std::string fulloutputfilenamebinary2Dbacteria = outputdirectory + outputfilenamebinary2Dbacteria; 
    std::cout << "Writing file: " << fulloutputfilenamebinary2Dbacteria << " ..." << std::endl;
    write2D(binaryimage2Dbacteria, fulloutputfilenamebinary2Dbacteria);

    std::cout << "Removing objects which are less than " << minNumberOfmm3 * 1000000000 << " um3..." << std::endl;
    excludeSmallObjects(binaryImageToShapeLabelMapFilter, minNumberOfmm3);
    std::cout << "Removing clusters of more than " << maxclustersize << " bacteria..." << std::endl;
    excludeClusters(binaryImageToShapeLabelMapFilter, maxclustersize);
    std::cout << "Removing objects which are bigger than " << maxNumberOfmm3 * 1000000000 << " um3..." << std::endl;
    excludeLargeObjects(binaryImageToShapeLabelMapFilter, maxNumberOfmm3);

    LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
    binaryImageToShapeLabelMapFilter->Update();
    labelMapToLabelImageFilter->SetInput(binaryImageToShapeLabelMapFilter->GetOutput());
    LabelImageToStatisticsLabelMapFilterType::Pointer labelImageToStatisticsLabelMapFilter = LabelImageToStatisticsLabelMapFilterType::New();
    labelImageToStatisticsLabelMapFilter->SetFeatureImage(image3Dred);
    labelImageToStatisticsLabelMapFilter->SetInput(labelMapToLabelImageFilter->GetOutput());
    labelImageToStatisticsLabelMapFilter->Update();


    


    bool er = false;    
    if (binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() == labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()) {
      std::cout << std::endl;
      for(unsigned int i = 0; i < labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
        LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObjectSt = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
        int statLabel = static_cast<int>(labelObjectSt->GetLabel());
        std::cout << std::setw(30) << "Label Statistics LabelMap: " << std::setw(5) << statLabel; 
        LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroidSt = labelObjectSt->GetCentroid ();      
        std::cout << std::setw(20) << "Centroid in %:";
        for (unsigned int j = 0; j < imSize.Size(); j++) {
          std::cout << std::setw(7) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << centroidSt[j] / imSize[j] * 100 << "  ";
          }
        std::cout << std::endl;   
        BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObjectBi = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
        int binLabel = static_cast<int>(labelObjectBi->GetLabel());
        std::cout << std::setw(30) << "Label Binary LabelMap: " << std::setw(5) << binLabel;  
        BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroidBi = labelObjectBi->GetCentroid ();      
        std::cout << std::setw(20) << "Centroid in %:";
        for (unsigned int j = 0; j < imSize.Size(); j++) {
          std::cout << std::setw(7) << std::setprecision(3) << std::setiosflags(std::ios::fixed) << centroidBi[j] / imSize[j] * 100 << "  ";
          }
        std::cout << std::endl << std::endl;
        if (statLabel != binLabel) er = true;        
        }   
      std::cout << std::endl;

      if (er) {
        std::cout << "Error: Not the same label!" << std::endl;
        fileout << "\n";
        fileout.close();
        } 
      else {
        unsigned int bacteriacount = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
        for(unsigned int i = 0; i < labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
          LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
          double mean = labelObject->GetMean();
          fileout << mean << "\t";
          std::cout << "Mean value of object with label " << static_cast<int>(labelObject->GetLabel()) << " in lysosomechannel: " << mean << std::endl; 
          }
        fileout << "\n";
        fileout.close();
        std::cout << "Total bacteria counted (in statistics label map): " << bacteriacount << std::endl;
        std::cout << std::endl;  
        }
      } 
    else {
      std::cout << "Error: Number of LabelObjects not the same in StatisticsLabelMap!" << std::endl;
      std::cout << "NumberOfLabelObjects Shape: " << static_cast<int>(binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()) << std::endl; 
      std::cout << "NumberOfLabelObjects Statistics: " << static_cast<int>(labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()) << std::endl;
      fileout << "\n";
      fileout.close();        
      }



    L2ImageType::Pointer labelMapToReadMeanImage = L2ImageType::New();
    labelMapToReadMeanImage->SetInput(labelImageToStatisticsLabelMapFilter->GetOutput());
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
  std::cout << std::endl;
  return EXIT_SUCCESS;
}






