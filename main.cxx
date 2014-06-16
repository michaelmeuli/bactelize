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


int seriesnr = 0; 
float xspacing = 0.33;
float yspacing = 0.33;
float zspacing = 1.2;
float tspacing = 1.0;
float cspacing = 1.0;
int nucleichannel      = 0;
int bacteriachannel    = 1;
int lysotrackerchannel = 2;


int main( int argc, char * argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << argv[0] << " series.ome.tiff-file" << " outputdirectory" << std::endl;
    return -1;
    }

  const char * inputFileName  = argv[1];
  SeriesReader seriesreader(inputFileName);
  std::cout << "Getting 5D Image of series number: " << seriesnr << std::endl;
  ImageType5D::Pointer image5D = seriesreader.get5DImage(seriesnr);
  seriesreader.dumpimageio();
  dumpmetadatadic(image5D);
  setspacing(image5D, xspacing, yspacing, zspacing, tspacing, cspacing);
  ImageType3D::ConstPointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
  printHistogram(image3Dbacteria);
  ImageType2D::Pointer image2Dbacteria = maxintprojection(image3Dbacteria);

  RescaleFilterTypeWriter::Pointer rescaleFilter = RescaleFilterTypeWriter::New();
  rescaleFilter->SetInput( image2Dbacteria ); 
  rescaleFilter->SetOutputMinimum( 0 );
  rescaleFilter->SetOutputMaximum( 255 );
  WriterType::Pointer writer = WriterType::New();
  std::string outputdirectory = argv[2];
  std::string filename = itksys::SystemTools::GetFilenameName(argv[1]);
  std::string filenamepath = outputdirectory + filename; 
  writer->SetFileName( filenamepath );             //seriesnumber has to be added
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

  NormalizeFilterType::Pointer  normalizeFilter = NormalizeFilterType::New();
  normalizeFilter->SetInput( image3Dbacteria );
  normalizeFilter->Update();


  QuickView viewer;
  viewer.AddImage(image2Dbacteria.GetPointer(), true, itksys::SystemTools::GetFilenameName(argv[1]));  
  viewer.Visualize();

  return EXIT_SUCCESS;
}
