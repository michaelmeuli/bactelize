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

int nucleichannel      = 0;
int bacteriachannel    = 1;
int lysosomechannel    = 2;

int main( int argc, char * argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << argv[0] << " series.ome.tiff-file" << " outputdirectory" << std::endl;
    return -1;
    }
    
  float xspacing = 0.33;
  float yspacing = 0.33;
  float zspacing = 1.2;
  float tspacing = 1.0;
  float cspacing = 1.0;
  std::string fullinputfilename = argv[1];
  std::string outputdirectory = argv[2];
  int seriesnr = 1; 
  int binaryLowerThresholdBacteria = 200;

  SeriesReader seriesreader(fullinputfilename);
  std::cout << "Getting 5D Image of series number: " << seriesnr << std::endl;
  ImageType5D::Pointer image5D = seriesreader.get5DImage(seriesnr);
  seriesreader.dumpimageio();
  dumpmetadatadic(image5D);
//  setSpacing(image5D, xspacing, yspacing, zspacing, tspacing, cspacing);

  ImageType3D::Pointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
  ImageType3D::Pointer image3Dred = extractchannel(image5D, lysosomechannel);
  
  printHistogram(image3Dbacteria);
  
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
  
  setSpacing(binaryimage3Dbacteria, xspacing, yspacing, zspacing);
  printSpacing(binaryimage3Dbacteria);
  
  BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
  binaryImageToShapeLabelMapFilter->SetInput(binaryimage3Dbacteria);
  binaryImageToShapeLabelMapFilter->Update();
  // The output of this filter is an itk::ShapeLabelMap, which contains itk::ShapeLabelObject's
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " objects." << std::endl;
  std::vector<unsigned long> labelsToRemove;
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    if (labelObject->GetNumberOfPixels() < 20) {                      // labelObject->GetPhysicalSize() < 10   didn't work as value is always 0
      labelsToRemove.push_back(labelObject->GetLabel());
      }
    }    
  std::cout << "Removing " << labelsToRemove.size() << " objects." << std::endl;
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i) {
    binaryImageToShapeLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
    }
  std::cout << "There are " << binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() 
            << " objects remaining." << std::endl;  
  unsigned int bacteriacount = binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
  unsigned int colcount = 0; 
  for(unsigned int i = 0; i < binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    labelObject->Print(std::cout, 2);
    unsigned long long  pixelValue = 0;  //ImageType3D::PixelType (16bit) was too short!
    for(unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++) {
      pixelValue += image3Dred->GetPixel( labelObject->GetIndex(pixelId) );
      }
    unsigned int mean = pixelValue / (labelObject->Size());
    std::cout << "Mean value in lysosomechannel of object " << i << ": " << mean << std::endl;  
    if ( mean > 100 ) {
      colcount++;
      }
    std::cout << std::endl;    
    }
  std::cout << "Total bacteria counted: " << bacteriacount << std::endl;
  std::cout << "Bacteria colocalizing in lysosomechannel: " << colcount << std::endl;
    
    
  ImageType2D::Pointer image2Dbacteria = maxintprojection(image3Dbacteria);
  std::string outputfilename2Dbacteria = seriesreader.getFilename(seriesnr, "_bacteria");
  std::string fulloutputfilename2Dbacteria = outputdirectory + outputfilename2Dbacteria; 
  std::cout << "Writing file: " << fulloutputfilename2Dbacteria << " ..." << std::endl;
  write2D(image2Dbacteria, fulloutputfilename2Dbacteria);
  
  ImageType2D::Pointer image2Dred = maxintprojection(image3Dred);
  std::string outputfilename2Dred = seriesreader.getFilename(seriesnr, "_lysosome");
  std::string fulloutputfilename2Dred = outputdirectory + outputfilename2Dred; 
  std::cout << "Writing file: " << fulloutputfilename2Dred << " ..." << std::endl;
  write2D(image2Dred, fulloutputfilename2Dred);
  
  BinaryImageType2D::Pointer binaryimage2Dbacteria = maxintprojection(binaryimage3Dbacteria);
  std::string outputfilenamebinary2Dbacteria = seriesreader.getFilename(seriesnr, "_bacteria_binary");
  std::string fulloutputfilenamebinary2Dbacteria = outputdirectory + outputfilenamebinary2Dbacteria; 
  std::cout << "Writing file: " << fulloutputfilenamebinary2Dbacteria << " ..." << std::endl;
  write2D(binaryimage2Dbacteria, fulloutputfilenamebinary2Dbacteria);

  QuickView viewer;
  viewer.AddImage(image2Dbacteria.GetPointer(), true, outputfilename2Dbacteria);
  viewer.AddImage(image2Dred.GetPointer(), true, outputfilename2Dred);
  viewer.AddImage(binaryimage2Dbacteria.GetPointer(), true, outputfilenamebinary2Dbacteria); 
  viewer.Visualize();


  return EXIT_SUCCESS;
}
