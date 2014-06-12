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

int main( int argc, char * argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << argv[0] << " image5DFile" << std::endl;
    return -1;
    }

  itk::SCIFIOImageIO::Pointer io = itk::SCIFIOImageIO::New();
  io->DebugOn();
  typename ReaderType::Pointer reader = ReaderType::New();
  std::cout << "reader->GetUseStreaming(): " << reader->GetUseStreaming() << std::endl;
  std::cout << "done checking streaming usage" << std::endl;
  reader->SetImageIO( io );
  const char * inputFileName  = argv[1];
  reader->SetFileName( inputFileName );


  typename StreamingFilter::Pointer streamer = StreamingFilter::New();
  streamer->SetInput( reader->GetOutput() );
  streamer->SetNumberOfStreamDivisions( 4 );

  typename ImageType5D::Pointer image5D = ImageType5D::New();
  image5D = streamer->GetOutput();

  reader->UpdateOutputInformation();
  io->SetSeries(3);    //set series even if different series file is given as command line argument
  reader->Modified();

  try
    {
    image5D->Update();
    }
  catch (itk::ExceptionObject &e)
    {
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }


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


  // Dump the metadata naturally contained within ImageIOBase
  const itk::ImageIOBase * imageIO = reader->GetImageIO();
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


  // SetSpacing
  std::cout << "--== Correcting spacing and setting origin ==--" << std::endl;
  ImageType5D::SpacingType spacing;
  spacing[0] = 0.33; // spacing along X
  spacing[1] = 0.33; // spacing along Y
  spacing[2] = 1.20; // spacing along Z 
  spacing[3] = 1.0;
  spacing[4] = 1.0;
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
  

  // Get the 3D of second channell
  typedef itk::ExtractImageFilter< ImageType5D, ImageType4D > ExtractFilterType5D4D;
  ExtractFilterType5D4D::Pointer extractfilter5D4Dch2 = ExtractFilterType5D4D::New();
  extractfilter5D4Dch2->InPlaceOn();
  extractfilter5D4Dch2->SetDirectionCollapseToSubmatrix();
  ImageType5D::SizeType size5Dch2 = region5D.GetSize();
  std::cout << "Extract 5D to 4D channell 2: size5Dch2= " 
	<< size5Dch2[0] << ", " << size5Dch2[1] << ", " << size5Dch2[2] << ", " << size5Dch2[3] << ", " << size5Dch2[4] << std::endl;  
  size5Dch2[4] = 0;
  std::cout << "Extract 5D to 4D channell 2: size5Dch2= " 
	<< size5Dch2[0] << ", " << size5Dch2[1] << ", " << size5Dch2[2] << ", " << size5Dch2[3] << ", " << size5Dch2[4] << std::endl; 
  ImageType5D::IndexType start5Dch2 = region5D.GetIndex();
  std::cout << "Extract 5D to 4D channell 2: start5Dch2= " 
	<< start5Dch2[0] << ", " << start5Dch2[1] << ", " << start5Dch2[2] << ", " << start5Dch2[3] << ", " << start5Dch2[4] << std::endl; 
  start5Dch2[4] = 1;
  std::cout << "Extract 5D to 4D channell 2: start5Dch2= " 
	<< start5Dch2[0] << ", " << start5Dch2[1] << ", " << start5Dch2[2] << ", " << start5Dch2[3] << ", " << start5Dch2[4] << std::endl; 
  ImageType5D::RegionType region5Dch2;
  region5Dch2.SetSize(  size5Dch2  );
  region5Dch2.SetIndex( start5Dch2 );
  extractfilter5D4Dch2->SetExtractionRegion( region5Dch2 );
  extractfilter5D4Dch2->SetInput( image5D );

  typedef itk::ExtractImageFilter< ImageType4D, ImageType3D > ExtractFilterType4D3D;
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

  ImageType3D::ConstPointer inputImageMIP;
  inputImageMIP = extractfilter4D3Dt1->GetOutput();


  // Histogram of second channell
  typedef itk::Statistics::ImageToHistogramFilter<ImageType3D>   HistogramFilterType;
  HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  typedef HistogramFilterType::HistogramSizeType   SizeType;
  SizeType size( 1 );
  size[0] =  40;        // number of bins for the green channel
  histogramFilter->SetHistogramSize( size );

  histogramFilter->SetMarginalScale( 10.0 ); 
  HistogramFilterType::HistogramMeasurementVectorType lowerBound( 3 );
  HistogramFilterType::HistogramMeasurementVectorType upperBound( 3 );
  lowerBound[0] = 0;
  upperBound[0] = 65536;
  histogramFilter->SetHistogramBinMinimum( lowerBound );
  histogramFilter->SetHistogramBinMaximum( upperBound ); 
  histogramFilter->SetInput(  inputImageMIP  );
  histogramFilter->Update();
  
  typedef HistogramFilterType::HistogramType  HistogramType;
  const HistogramType * histogram = histogramFilter->GetOutput();
  const unsigned int histogramSize = histogram->Size();
  std::cout << std::endl << "Histogram size " << histogramSize << std::endl;
 
  std::cout << std::endl << "Histogram of the green channell" << std::endl;
  for( unsigned int bin=0; bin < histogramSize; bin++ )
    {
    std::cout << "bin = " << std::setw(3) << bin << 
      "        measurement = " << std::setw(10) << std::setprecision(1) << std::setiosflags(std::ios::fixed) <<  histogram->GetMeasurement (bin, 0) <<
      "        frequency = " << std::setw(10) << histogram->GetFrequency( bin, 0 ) << std::endl;	
    }


  ImageType2D::Pointer outputImageMIP = maxintprojection(inputImageMIP);

  typedef itk::Image<unsigned char, 2>  ImageTypeWriter;
  typedef itk::ImageFileWriter< ImageTypeWriter > WriterType;
  typedef itk::RescaleIntensityImageFilter< ImageType2D, ImageTypeWriter >  RescaleFilterType;

  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput( outputImageMIP ); 
  rescaleFilter->SetOutputMinimum( 0 );
  rescaleFilter->SetOutputMaximum( 255 );
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput(rescaleFilter->GetOutput());
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }




  QuickView viewer;
  viewer.AddImage(outputImageMIP.GetPointer(), true, itksys::SystemTools::GetFilenameName(argv[1]));  
  viewer.Visualize();



  return EXIT_SUCCESS;
}
