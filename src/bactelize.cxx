#include "itkSCIFIOImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkMetaDataObject.h"
#include "itkStreamingImageFilter.h"

int main( int argc, char * argv [] )
{

  if( argc < 2 )
    {
    std::cerr << "Missing command line arguments" << std::endl;
    std::cerr << "Usage :  Test input-ome-tiff-ImageFileName " << std::endl;
    return -1;
    }

  typedef unsigned int       PixelType;
  const unsigned int          Dimension = 2;
  typedef itk::Image<PixelType, Dimension > ImageType;

  typedef typename itk::ImageFileReader< ImageType > ReaderType;
  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  io->DebugOn();
  typename ReaderType::Pointer reader = ReaderType::New();
  std::cout << "reader->GetUseStreaming(): " << reader->GetUseStreaming() << std::endl;
  std::cout << "done checking streaming usage" << std::endl;
  reader->SetImageIO( io );
  const char * inputFileName  = argv[1];
  reader->SetFileName( inputFileName );

  typedef itk::StreamingImageFilter< ImageType, ImageType > StreamingFilter;
  typename StreamingFilter::Pointer streamer = StreamingFilter::New();
  streamer->SetInput( reader->GetOutput() );
  streamer->SetNumberOfStreamDivisions( 4 );


  try
    {
    streamer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Problem encoutered while reading image file : " << argv[1] << std::endl;
    std::cerr << excp << std::endl;
    return -1;
    }
}
