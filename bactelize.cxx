#include "bactelize.h"



typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ImageType3D > SliceIteratorType;

ImageType2D::Pointer maxintprojection(ImageType3D::ConstPointer inputImageMIP, unsigned int projectionDirection) {
  unsigned int i, j;
  unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i )
    {
    if (i != projectionDirection)
      {
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
 
  SliceIteratorType  inputIt(  inputImageMIP, inputImageMIP->GetRequestedRegion() );
  LinearIteratorType outputIt( outputImageMIP, outputImageMIP->GetRequestedRegion() );
  inputIt.SetFirstDirection(  direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );
  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() )
    {
    while ( ! outputIt.IsAtEndOfLine() )
      {
      outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
      ++outputIt;
      }
    outputIt.NextLine();
    }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while( !inputIt.IsAtEnd() )
    {
    while ( !inputIt.IsAtEndOfSlice() )
      {
      while ( !inputIt.IsAtEndOfLine() )
        {
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
  // Dump the metadata dictionary
  std::cout << std::endl;
  std::cout << "--== Metadata from dictionary ==--" << std::endl;
  itk::MetaDataDictionary imgMetaDictionary = image5D->GetMetaDataDictionary();
  std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
  for(std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
      itKey != imgMetaKeys.end(); ++itKey)
    {
    std::string tmp;
    itk::ExposeMetaData<std::string>( imgMetaDictionary, *itKey, tmp );
    std::cout << "\t" << *itKey << " ---> " << tmp << std::endl;
    }
  std::cout << std::endl;
  }


void setspacing(ImageType5D::Pointer image5D, float x, float y, float z, float t, float c) {
  // SetSpacing
  std::cout << "--== Correcting spacing and setting origin ==--" << std::endl;
  ImageType5D::SpacingType spacing;
  spacing[0] = x;  
  spacing[1] = y;  
  spacing[2] = z;  
  spacing[3] = t; 
  spacing[4] = c;  
  image5D->SetSpacing( spacing );
  ImageType5D::PointType origin;
  origin.Fill(0.0);
  image5D->SetOrigin( origin );
  ImageType5D::RegionType region5D = image5D->GetLargestPossibleRegion();
  int regionDimIm = region5D.GetImageDimension();
  const ImageType5D::SpacingType& sp = image5D->GetSpacing();
  const ImageType5D::PointType& orgn = image5D->GetOrigin();
  for(int i = 0; i < regionDimIm; i++)
    {
    std::cout << "\tDimension " << i + 1 << " Size: "
              << region5D.GetSize(i) << std::endl;
    }
  for(int i = 0; i < regionDimIm; i++)
  {
    if ( region5D.GetSize(i) > 1 ) {
      std::cout << "\tSpacing " << i + 1 << ": "
                << sp[i] << std::endl;
    }
  }
  for(int i = 0; i < regionDimIm; i++)
    {
    std::cout << "\tOrigin " << i + 1 << ": "
              << orgn[i] << std::endl;
    }
  std::cout << std::endl;
  }



ImageType3D::ConstPointer extractchannel(ImageType5D::Pointer image5D, int channelnr) {
  // Get the 3D of second channell
  ExtractFilterType5D4D::Pointer extractfilter5D4Dch2 = ExtractFilterType5D4D::New();
  extractfilter5D4Dch2->InPlaceOn();
  extractfilter5D4Dch2->SetDirectionCollapseToSubmatrix();
  ImageType5D::RegionType region5D = image5D->GetLargestPossibleRegion();
  ImageType5D::SizeType size5Dch = region5D.GetSize();
  std::cout << "Extract 5D to 4D channell " << channelnr << ": size5Dch= "
	<< size5Dch[0] << ", " << size5Dch[1] << ", " << size5Dch[2] << ", " << size5Dch[3] << ", " << size5Dch[4] << std::endl;  
  size5Dch[4] = 0;
  std::cout << "Extract 5D to 4D channell " << channelnr << ": size5Dch= " 
	<< size5Dch[0] << ", " << size5Dch[1] << ", " << size5Dch[2] << ", " << size5Dch[3] << ", " << size5Dch[4] << std::endl; 
  ImageType5D::IndexType start5Dch = region5D.GetIndex();
  std::cout << "Extract 5D to 4D channell " << channelnr << ": start5Dch= " 
	<< start5Dch[0] << ", " << start5Dch[1] << ", " << start5Dch[2] << ", " << start5Dch[3] << ", " << start5Dch[4] << std::endl; 
  start5Dch[4] = channelnr;
  std::cout << "Extract 5D to 4D channell " << channelnr << ": start5Dch= " 
	<< start5Dch[0] << ", " << start5Dch[1] << ", " << start5Dch[2] << ", " << start5Dch[3] << ", " << start5Dch[4] << std::endl; 
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
  std::cout << "Extract 4D to 3D timepoint 1: size4Dt1= " << size4Dt1[0] << ", " << size4Dt1[1] << ", " << size4Dt1[2] << ", " << size4Dt1[3] << std::endl;  
  size4Dt1[3] = 0;
  std::cout << "Extract 4D to 3D timepoint 1: size4Dt1= " << size4Dt1[0] << ", " << size4Dt1[1] << ", " << size4Dt1[2] << ", " << size4Dt1[3] << std::endl;  
  ImageType4D::IndexType start4Dt1 = region4D.GetIndex();
  std::cout << "Extract 4D to 3D timepoint 1: start4Dt1= " << start4Dt1[0] << ", " << start4Dt1[1] << ", " << start4Dt1[2] << ", " << start4Dt1[3] << std::endl; 
  start4Dt1[3] = 0;
  std::cout << "Extract 4D to 3D timepoint 1: start4Dt1= " << start4Dt1[0] << ", " << start4Dt1[1] << ", " << start4Dt1[2] << ", " << start4Dt1[3] << std::endl; 
  ImageType4D::RegionType region4Dt1;
  region4Dt1.SetSize(  size4Dt1  );
  region4Dt1.SetIndex( start4Dt1 );
  extractfilter4D3Dt1->SetExtractionRegion( region4Dt1 );
  extractfilter4D3Dt1->SetInput( extractfilter5D4Dch2->GetOutput() );
  extractfilter4D3Dt1->Update();
  return extractfilter4D3Dt1->GetOutput();
  }


void printHistogram(ImageType3D::ConstPointer inputImageMIP) {
  HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  typedef HistogramFilterType::HistogramSizeType   SizeType;
  SizeType size( 1 );
  size[0] =  40;        // number of bins for the green channel
  histogramFilter->SetHistogramSize( size );
  histogramFilter->SetMarginalScale( 10.0 ); 
  HistogramFilterType::HistogramMeasurementVectorType lowerBound( 1 );
  HistogramFilterType::HistogramMeasurementVectorType upperBound( 1 );
  lowerBound[0] = 0;
  upperBound[0] = 4095;
  histogramFilter->SetHistogramBinMinimum( lowerBound );
  histogramFilter->SetHistogramBinMaximum( upperBound ); 
  histogramFilter->SetInput(  inputImageMIP  );
  histogramFilter->Update();
  const HistogramType * histogram = histogramFilter->GetOutput();
  const unsigned int histogramSize = histogram->Size();
  std::cout << std::endl << "Histogram size " << histogramSize << std::endl;
  for( unsigned int bin=0; bin < histogramSize; bin++ )
    {
    std::cout << "bin = " << std::setw(3) << bin << 
      "        measurement = " << std::setw(10) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, 0) <<
      "        frequency = " << std::setw(10) << histogram->GetFrequency( bin, 0 ) << std::endl;    
    }    
  }


SeriesReader::SeriesReader(std::string inputFileName)
        : m_io(itk::SCIFIOImageIO::New()), m_reader(ReaderType::New()), m_inputFileName(inputFileName),
          m_streamer(StreamingFilter::New()), m_image5D(ImageType5D::New()), m_seriesStart(0), m_seriesEnd(1) {
    m_io->DebugOn();
    std::cout << "reader->GetUseStreaming(): " << m_reader->GetUseStreaming() << std::endl;
    std::cout << "done checking streaming usage" << std::endl;
    m_reader->SetImageIO( m_io );
    m_reader->SetFileName( m_inputFileName );
    m_streamer->SetInput( m_reader->GetOutput() );
    m_streamer->SetNumberOfStreamDivisions( 4 );
    m_image5D = m_streamer->GetOutput();
    m_reader->UpdateOutputInformation();
    m_seriesEnd = m_io->GetSeriesCount();
    std::cout << "Number of series: " << m_seriesEnd + 1 << std::endl;
    }

ImageType5D::Pointer SeriesReader::get5DImage(int series) {
    m_io->SetSeries(series); 
    m_reader->Modified();  
    m_image5D->Update();
    return m_image5D;
  }

void SeriesReader::dumpimageio() {
  // Dump the metadata naturally contained within ImageIOBase
  const itk::ImageIOBase * imageIO = m_reader->GetImageIO();
  itk::ImageIORegion regionIO = imageIO->GetIORegion();
  int regionDimIO = regionIO.GetImageDimension();
  std::cout << "--== Metadata from ImageIOBase ==--" << std::endl;
  for(int i = 0; i < regionDimIO; i++)
    {
    std::cout << "\tDimension " << i + 1 << " Size: "
              << regionIO.GetSize(i) << std::endl;
    }
  for(int i = 0; i < regionDimIO; i++)
  {
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
  std::cout << std::endl;
  }


