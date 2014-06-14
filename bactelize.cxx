#include "bactelize.h"



typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ImageType3D > SliceIteratorType;

ImageType2D::Pointer maxintprojection(ImageType3D::ConstPointer inputImageMIP) {

  unsigned int projectionDirection = 2;
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
    std::cout << "m_seriesEnd: " << m_seriesEnd << std::endl;
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


