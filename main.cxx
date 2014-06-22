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
    

  std::string fullinputfilename = argv[1];
  std::string outputdirectory = argv[2];
  int seriesnr = 1; 

  SeriesReader seriesreader(fullinputfilename);
  std::cout << "Getting 5D Image of series number: " << seriesnr << std::endl;
  ImageType5D::Pointer image5D = seriesreader.get5DImage(seriesnr);
  seriesreader.dumpimageio();
  dumpmetadatadic(image5D);
  setspacing(image5D, xspacing, yspacing, zspacing, tspacing, cspacing);

  ImageType3D::ConstPointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
  printHistogram(image3Dbacteria);
  ImageType2D::Pointer image2Dbacteria = maxintprojection(image3Dbacteria);
  std::string outputfilename2Dbacteria = seriesreader.getFilename(seriesnr, "_bacteria");
  std::string fulloutputfilename2Dbacteria = outputdirectory + outputfilename2Dbacteria; 
  std::cout << "Writing file: " << fulloutputfilename2Dbacteria << " ..." << std::endl;
  write2D(image2Dbacteria, fulloutputfilename2Dbacteria);

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
  binaryfilter->SetLowerThreshold(200);

  binaryfilter->Update();

  ImageType3D::ConstPointer binaryimage3Dbacteria = binaryfilter->GetOutput();
  ImageType2D::Pointer binaryimage2Dbacteria = maxintprojection(binaryimage3Dbacteria);
  std::string outputfilenamebinary2Dbacteria = seriesreader.getFilename(seriesnr, "_bacteria_binary");
  std::string fulloutputfilenamebinary2Dbacteria = outputdirectory + outputfilenamebinary2Dbacteria; 
  std::cout << "Writing file: " << fulloutputfilenamebinary2Dbacteria << " ..." << std::endl;
  write2D(binaryimage2Dbacteria, fulloutputfilenamebinary2Dbacteria);

  QuickView viewer;
  viewer.AddImage(image2Dbacteria.GetPointer(), true, outputfilename2Dbacteria); 
  viewer.AddImage(binaryimage2Dbacteria.GetPointer(), true, outputfilename2Dbacteria); 
  viewer.Visualize();

  return EXIT_SUCCESS;
}
