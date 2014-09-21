#include "bactelize.h"

extern int nucleichannel;
extern int bacteriachannel;
extern int lysosomechannel;
extern int binaryLowerThresholdBacteria;
extern double minNumberOfmm3;  
extern std::ofstream fileout;
extern std::string fileoutName;
extern int numberOfStreamDivisions;
extern int numberOfBins;


ImageType2D::Pointer maxintprojection(ImageType3D::Pointer inputImageMIP, unsigned int projectionDirection) {
  unsigned int i, j;
  unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i ) {
    if (i != projectionDirection) {
      direction[j] = i;
      j++;
      }
    }
  ImageType2D::RegionType region2DMIP;
  ImageType2D::RegionType::SizeType size2DMIP;
  ImageType2D::RegionType::IndexType index2DMIP;
  ImageType3D::RegionType requestedRegion = inputImageMIP -> GetRequestedRegion();
  index2DMIP[ direction[0] ]    = requestedRegion.GetIndex()[ direction[0] ];
  index2DMIP[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size2DMIP[ direction[0] ]     = requestedRegion.GetSize()[  direction[0] ];
  size2DMIP[ 1- direction[0] ]  = requestedRegion.GetSize()[  direction[1] ];
  region2DMIP.SetSize( size2DMIP );
  region2DMIP.SetIndex( index2DMIP );
  ImageType2D::Pointer outputImageMIP = ImageType2D::New();
  outputImageMIP->SetRegions( region2DMIP );
  outputImageMIP->Allocate();
 
  SliceIteratorTypeInput  inputIt(  inputImageMIP, inputImageMIP->GetRequestedRegion() );
  LinearIteratorTypeInput outputIt( outputImageMIP, outputImageMIP->GetRequestedRegion() );
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
  while( !inputIt.IsAtEnd() ) {
    while ( !inputIt.IsAtEndOfSlice() ) {
      while ( !inputIt.IsAtEndOfLine() ) {
        outputIt.Set( vnl_math_max( outputIt.Get(), inputIt.Get() ));
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
  BinaryImageType2D::RegionType region2DMIP;
  BinaryImageType2D::RegionType::SizeType size2DMIP;
  BinaryImageType2D::RegionType::IndexType index2DMIP;
  BinaryImageType3D::RegionType requestedRegion = inputImageMIP -> GetRequestedRegion();
  index2DMIP[ direction[0] ]    = requestedRegion.GetIndex()[ direction[0] ];
  index2DMIP[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size2DMIP[ direction[0] ]     = requestedRegion.GetSize()[  direction[0] ];
  size2DMIP[ 1- direction[0] ]  = requestedRegion.GetSize()[  direction[1] ];
  region2DMIP.SetSize( size2DMIP );
  region2DMIP.SetIndex( index2DMIP );
  BinaryImageType2D::Pointer outputImageMIP = BinaryImageType2D::New();
  outputImageMIP->SetRegions( region2DMIP );
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
  while( !inputIt.IsAtEnd() ) {
    while ( !inputIt.IsAtEndOfSlice() ) {
      while ( !inputIt.IsAtEndOfLine() ) {
        outputIt.Set( vnl_math_max( outputIt.Get(), inputIt.Get() ));
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
  for(std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin(); itKey != imgMetaKeys.end(); ++itKey) {
    std::string tmp;
    itk::ExposeMetaData<std::string>( imgMetaDictionary, *itKey, tmp );
    std::cout << "\t" << *itKey << " ---> " << tmp << std::endl;
    }
  }


void dumpimageio(ReaderType::Pointer reader) {
  const itk::ImageIOBase * imageIO = reader->GetImageIO();
  itk::ImageIORegion regionIO = imageIO->GetIORegion();
  int regionDimIO = regionIO.GetImageDimension();
  std::cout << "--== Metadata from ImageIOBase ==--" << std::endl;
  for(int i = 0; i < regionDimIO; i++) {
    std::cout << "\tDimension " << i + 1 << " Size: "
              << regionIO.GetSize(i) << std::endl;
    }
  for(int i = 0; i < regionDimIO; i++) {
    if ( regionIO.GetSize(i) > 1 ) {
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



void printSpacing(BinaryImageType3D::Pointer image3D) {
  ImageType3D::RegionType region3D = image3D->GetLargestPossibleRegion();
  int regionDimIm = region3D.GetImageDimension();
  const ImageType3D::SpacingType& sp = image3D->GetSpacing();
  const ImageType3D::PointType& orgn = image3D->GetOrigin();
  for(int i = 0; i < regionDimIm; i++) {
    std::cout << "Dimension " << i + 1 << " Size: "
              << region3D.GetSize(i) << std::endl;
    }
  for(int i = 0; i < regionDimIm; i++) {
    if ( region3D.GetSize(i) > 1 ) {
      std::cout << "Spacing " << i + 1 << ": "
                << sp[i] << std::endl;
      }
    }
  for(int i = 0; i < regionDimIm; i++) {
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


void printObjectInfo(LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject) {         
  std::cout << "NumberOfPixels: " << labelObject->GetNumberOfPixels() << std::endl;
  std::cout << "PhysicalSize: " << (labelObject->GetPhysicalSize()) * 1000000000 << "um3" << std::endl;
  std::cout << "Elongation: " << labelObject->GetElongation() << std::endl;
  std::cout << "Roundness: " << labelObject->GetRoundness() << std::endl;
  std::cout << "EquivalentEllipsoidDiameter: " << (labelObject->GetEquivalentEllipsoidDiameter() * 1000) << "um" << std::endl;
  std::cout << "EquivalentSphericalRadius: " << (labelObject->GetEquivalentSphericalRadius() * 1000) << "um" << std::endl;
  }


ImageSizeType getImSize(LabelImageToShapeLabelMapFilterType::Pointer labelImageToShapeLabelMapFilter) {
  LabelImageToShapeLabelMapFilterType::OutputImageType::RegionType region = labelImageToShapeLabelMapFilter->GetOutput()->GetLargestPossibleRegion();
  const ImageType3D::SpacingType& spacing = labelImageToShapeLabelMapFilter->GetOutput()->GetSpacing();
  ImageSizeType imSize;
  for (unsigned int i = 0; i < imSize.Size(); i++) {
    imSize[i] = spacing[i] * region.GetSize(i);
    }
  return imSize;
  }


void printCentroids(LabelImageToShapeLabelMapFilterType::Pointer labelImageToShapeLabelMapFilter, ImageSizeType imSize) {
  for(unsigned int i = 0; i < labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();     
    std::cout << "Centroid of labelObject i=" << i << " in %:  ";
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






int processSeries(std::string inputFileName, std::string outputdirectory, bool vflag, bool tflag, int fileNr, int seriesNr) {
  std::cout << "Processing series of file: " << inputFileName << std::endl;
  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  io->DebugOn();
  typename ReaderType::Pointer reader = ReaderType::New();
  std::cout << "reader->GetUseStreaming(): " << reader->GetUseStreaming() << std::endl;
  std::cout << "done checking streaming usage" << std::endl;
  reader->SetImageIO(io);
  reader->SetFileName(inputFileName);
  typename StreamingFilter::Pointer streamer = StreamingFilter::New();
  streamer->SetInput(reader->GetOutput());
  streamer->SetNumberOfStreamDivisions(numberOfStreamDivisions);
  ImageType5D::Pointer image5D =  ImageType5D::New();
  image5D = streamer->GetOutput();
  reader->UpdateOutputInformation();
  io->SetSeries(0);
  reader->Modified(); 
  int seriesEnd = io->GetSeriesCount();
  std::cout << "seriesEnd: " << seriesEnd << std::endl;  

  int seriesnr=0;
  if (tflag) {
    seriesnr = seriesNr;
    seriesEnd = seriesNr+1;
    io->SetSeries(seriesnr);
    reader->Modified();
    }
  while ( seriesnr < seriesEnd ) {
    std::cout << "Reading seriesnr: " << seriesnr << std::endl;
    image5D->Update();

    ImageType3D::Pointer image3Dnuclei = extractchannel(image5D, nucleichannel);
    ImageType3D::Pointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
    ImageType3D::Pointer image3Dred = extractchannel(image5D, lysosomechannel);
  
    BinaryImageType3D::Pointer binaryimage3Dbacteria = getBinaryIm(image3Dbacteria);

    if (vflag) {
      std::cout << std::endl;
      dumpimageio(reader);
      std::cout << std::endl;
      dumpmetadatadic(image5D);
      std::cout << std::endl;
      printHistogramNormalized(image3Dbacteria);
      std::cout << std::endl;
      printSpacing(binaryimage3Dbacteria);
      std::cout << std::endl;
      }

    std::string fulloutputfilenameResults = outputdirectory + fileoutName;
    fileout.open(fulloutputfilenameResults.c_str(), std::ofstream::app); 
    std::string seriesName = getFilename(inputFileName, seriesnr, seriesEnd);
    fileout << seriesName << "\t";

    BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
    binaryImageToLabelMapFilter->SetInput(binaryimage3Dbacteria);
    LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
    labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
    LabelImageToShapeLabelMapFilterType::Pointer labelImageToShapeLabelMapFilter = LabelImageToShapeLabelMapFilterType::New ();
    labelImageToShapeLabelMapFilter->SetInput( labelMapToLabelImageFilter->GetOutput() );
    labelImageToShapeLabelMapFilter->Update();
    assert (binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() == labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()); 
    std::vector<unsigned long> labelsToRemove1;
    for(unsigned int i = 0; i < labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
      LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
      if (vflag) {
        printObjectInfo(labelObject);
        } 
      if (labelObject->GetPhysicalSize() < minNumberOfmm3) {       
        labelsToRemove1.push_back(labelObject->GetLabel());
        }
      std::cout << std::endl;
      }    
    std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " label map objects." << std::endl;
    std::cout << "Removing " << labelsToRemove1.size() << " objects from label map." << std::endl;
    for(unsigned int i = 0; i < labelsToRemove1.size(); ++i) {
      binaryImageToLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove1[i]);
      }
    std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects remaining in label map." << std::endl << std::endl;







    labelImageToShapeLabelMapFilter->Update();
    assert (binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() == labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()); 
    SampleType::Pointer sample = SampleType::New();
    sample->SetMeasurementVectorSize(3);

    MeasurementVectorType mv;
    std::cout << "Adding centroids to itk::Statistics::ListSample<MeasurementVectorType>::Pointer sample..." << std::endl;
    for(unsigned int i = 0; i < labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
      LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
      LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::CentroidType centroid = labelObject->GetCentroid ();     
      mv[0] = centroid[0];
      mv[1] = centroid[1];
      mv[2] = centroid[2];
      sample->PushBack( mv );
      } 

    ImageSizeType imSize = getImSize(labelImageToShapeLabelMapFilter);
    printCentroids(labelImageToShapeLabelMapFilter, imSize);
    std::cout << std::endl;
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
      LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
      LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType::VectorType ellipsoidVector = labelObject->GetEquivalentEllipsoidDiameter();
      for ( unsigned int vd = 0; vd < ellipsoidVector.GetVectorDimension(); ++vd ) {
        if (maxDiameter < ellipsoidVector[vd]) {
          maxDiameter = ellipsoidVector[vd];
          }
        }
      std::cout << "Maximum of equivalent ellipsoid diameter: " << maxDiameter * 1000 << "um" << std::endl;
      double radius = 0.003;
      if (maxDiameter >= 0.006) {   		//size of object is too big to be one bacteria
        radius = maxDiameter/2 + 0.001;
        setToRemove.insert(i);
        }
      tree->Search(queryPoint, radius, neighbors); 
      std::cout << "kd-tree radius search result:" << std::endl;
      std::cout << "query point = " << i << std::endl;
      std::cout << "search radius = " << radius * 1000 << "um" << std::endl;
      int clustersize = 3;
      if (neighbors.size() >= clustersize) {
        for (unsigned int k = 0; k < neighbors.size(); ++k) {       
          std::cout << "neighbor = " << neighbors[k] << ": " << distanceMetric->Evaluate(tree->GetMeasurementVector(neighbors[k])) * 1000 << "um" << std::endl;
          setToRemove.insert(neighbors[k]);         
          }
        }
      std::cout << std::endl; 
         
      }
    std::cout << std::endl;
 
    std::cout << "setToRemove:" << std::endl;
    std::set<int>::const_iterator itp;
    for (itp = setToRemove.begin(); itp != setToRemove.end(); ++itp) {
      int f = *itp; 
      std::cout << f << "  ";
      }
    std::cout << std::endl;

    std::vector<unsigned long> labelsToRemove2;
    std::set<int>::const_iterator itr;
    for (itr = setToRemove.begin(); itr != setToRemove.end(); ++itr) {
      LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(*itr);
      labelsToRemove2.push_back(labelObject->GetLabel());
      }   
    std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " label map objects." << std::endl;
    std::cout << "Removing " << labelsToRemove2.size() << " objects from label map." << std::endl;
    for(unsigned int i = 0; i < labelsToRemove2.size(); ++i) {
      binaryImageToLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove2[i]);
      }
    std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects remaining in label map." << std::endl;
    std::cout << std::endl;
    labelImageToShapeLabelMapFilter->Update();
    printCentroids(labelImageToShapeLabelMapFilter, imSize);
    std::cout << std::endl;






    LabelImageToStatisticsLabelMapFilterType::Pointer labelImageToStatisticsLabelMapFilter = LabelImageToStatisticsLabelMapFilterType::New();
    labelImageToStatisticsLabelMapFilter->SetFeatureImage(image3Dred);
    labelImageToStatisticsLabelMapFilter->SetInput(labelMapToLabelImageFilter->GetOutput());
    labelImageToStatisticsLabelMapFilter->Update();
    assert (binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() == labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()); 
    unsigned int bacteriacount = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
    for(unsigned int i = 0; i < labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
      LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
      if (vflag) {
        printObjectInfo(labelObject);
        }
      double mean = labelObject->GetMean();
      fileout << mean << "\t";
      std::cout << "Mean value of object with label " << static_cast<int>(labelObject->GetLabel()) << " in lysosomechannel: " << mean << std::endl; 
      std::cout << std::endl;
      }
    fileout << "\n";
    fileout.close();
    std::cout << "Total bacteria counted (in statistics label map): " << bacteriacount << std::endl;
    std::cout << std::endl;

    L2ImageType::Pointer labelMapToReadMeanImage = L2ImageType::New();
    labelMapToReadMeanImage->SetInput( labelImageToStatisticsLabelMapFilter->GetOutput() );
    labelMapToReadMeanImage->Update();
    ImageType3D::Pointer image3DbacteriaReadMean = labelMapToReadMeanImage->GetOutput();

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

    LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImageFilter = LabelMapToBinaryImageFilterType::New();
    labelMapToBinaryImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
    labelMapToBinaryImageFilter->Update();
    BinaryImageType2D::Pointer binaryimage2Dbacteria = maxintprojection(labelMapToBinaryImageFilter->GetOutput());
    std::string outputfilenamebinary2Dbacteria = getFilename(inputFileName, seriesnr, seriesEnd, "_d_bacteria_binary.tiff");
    std::string fulloutputfilenamebinary2Dbacteria = outputdirectory + outputfilenamebinary2Dbacteria; 
    std::cout << "Writing file: " << fulloutputfilenamebinary2Dbacteria << " ..." << std::endl;
    write2D(binaryimage2Dbacteria, fulloutputfilenamebinary2Dbacteria);

    ImageType2D::Pointer image2DReadMean = maxintprojection(image3DbacteriaReadMean);
    RGBFilterType::Pointer colormapImageFilter = RGBFilterType::New();
    colormapImageFilter->SetInput(image2DReadMean);
    colormapImageFilter->SetColormap( RGBFilterType::Jet );
    colormapImageFilter->Update();
    std::string outputfilename2DReadMean = getFilename(inputFileName, seriesnr, seriesEnd, "_e_ReadMean.tiff");
    std::string fulloutputfilename2DReadMean = outputdirectory + outputfilename2DReadMean; 
    WriterTypeRGB::Pointer writerRGB = WriterTypeRGB::New();
    TIFFIOType::Pointer tiffIO = TIFFIOType::New();
    writerRGB->SetImageIO(tiffIO);
    writerRGB->SetFileName( fulloutputfilename2DReadMean );             
    writerRGB->SetInput(colormapImageFilter->GetOutput());
    std::cout << "Writing file: " << fulloutputfilename2DReadMean << " ..." << std::endl;
    writerRGB->Update();

    seriesnr++;
    if ( seriesnr < seriesEnd) {
      io->SetSeries(seriesnr);
      reader->Modified();
      }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  return EXIT_SUCCESS;
}






