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
  int minNumberOfPixels = 100;
  int meanRedThreshold = 50;

  SeriesReader seriesreader(fullinputfilename);
  std::cout << "Getting 5D Image of series number: " << seriesnr << std::endl;
  ImageType5D::Pointer image5D = seriesreader.get5DImage(seriesnr);
  seriesreader.dumpimageio();
  dumpmetadatadic(image5D);

  ImageType3D::Pointer image3Dnuclei = extractchannel(image5D, nucleichannel);
  ImageType3D::Pointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
  ImageType3D::Pointer image3Dred = extractchannel(image5D, lysosomechannel);
  
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

  printHistogram( rescaleNormalized->GetOutput() );
  std::cout << std::endl;
  printSpacing(binaryimage3Dbacteria);
  std::cout << std::endl;
  std::string outputfilenamefileout = outputdirectory + "fileout.txt";
  std::ofstream fileout(outputfilenamefileout.c_str()); 

  BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
  binaryImageToLabelMapFilter->SetInput(binaryimage3Dbacteria);
  LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
  labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
  LabelImageToShapeLabelMapFilterType::Pointer labelImageToShapeLabelMapFilter = LabelImageToShapeLabelMapFilterType::New ();
  labelImageToShapeLabelMapFilter->SetInput( labelMapToLabelImageFilter->GetOutput() );
  labelImageToShapeLabelMapFilter->Update();
  assert (binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() == labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()); 

  std::vector<unsigned long> labelsToRemove;
  for(unsigned int i = 0; i < labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
//  labelObject->Print(std::cout, 5);
    if (labelObject->GetNumberOfPixels() < minNumberOfPixels) {                      // labelObject->GetPhysicalSize() < 10   didn't work as value is always 0
      labelsToRemove.push_back(labelObject->GetLabel());
      }
    }    
  std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " label map objects." << std::endl;
  std::cout << "Removing " << labelsToRemove.size() << " objects from label map." << std::endl;
  for(unsigned int i = 0; i < labelsToRemove.size(); ++i) {
    binaryImageToLabelMapFilter->GetOutput()->RemoveLabel(labelsToRemove[i]);
    }
  std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() 
            << " objects remaining in label map." << std::endl << std::endl;

  LabelImageToStatisticsLabelMapFilterType::Pointer labelImageToStatisticsLabelMapFilter = LabelImageToStatisticsLabelMapFilterType::New();
  labelImageToStatisticsLabelMapFilter->SetFeatureImage(image3Dred);
  labelImageToStatisticsLabelMapFilter->SetInput(labelMapToLabelImageFilter->GetOutput());
  labelImageToStatisticsLabelMapFilter->Update();
  assert (binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() == labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects()); 

  unsigned int bacteriacount = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
  for(unsigned int i = 0; i < labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
//  labelObject->Print(std::cout, 4);
    double mean = labelObject->GetMean();
    std::cout << "Mean value of object with label " << (int)(labelObject->GetLabel()) << " in lysosomechannel: " << mean << std::endl;   //todo c++ style cast
    fileout << mean << "\t";
    }
  std::cout << "Total bacteria counted (in statistics label map): " << bacteriacount << std::endl;
  std::cout << std::endl;

  L2ImageType::Pointer labelMapToReadMeanImage = L2ImageType::New();
  labelMapToReadMeanImage->SetInput( labelImageToStatisticsLabelMapFilter->GetOutput() );
  labelMapToReadMeanImage->Update();
  ImageType3D::Pointer image3DbacteriaReadMean = labelMapToReadMeanImage->GetOutput();
  

  ImageType2D::Pointer image2Dnuclei = maxintprojection(image3Dnuclei);
  std::string outputfilename2Dnuclei = seriesreader.getFilename(seriesnr, "_a_nuclei.tiff");
  std::string fulloutputfilename2Dnuclei = outputdirectory + outputfilename2Dnuclei; 
  std::cout << "Writing file: " << fulloutputfilename2Dnuclei << " ..." << std::endl;
  write2D(image2Dnuclei, fulloutputfilename2Dnuclei);

  ImageType2D::Pointer image2Dbacteria = maxintprojection(image3Dbacteria);
  std::string outputfilename2Dbacteria = seriesreader.getFilename(seriesnr, "_b_bacteria.tiff");
  std::string fulloutputfilename2Dbacteria = outputdirectory + outputfilename2Dbacteria; 
  std::cout << "Writing file: " << fulloutputfilename2Dbacteria << " ..." << std::endl;
  write2D(image2Dbacteria, fulloutputfilename2Dbacteria);
  
  ImageType2D::Pointer image2Dred = maxintprojection(image3Dred);
  std::string outputfilename2Dred = seriesreader.getFilename(seriesnr, "_c_lysosome.tiff");
  std::string fulloutputfilename2Dred = outputdirectory + outputfilename2Dred; 
  std::cout << "Writing file: " << fulloutputfilename2Dred << " ..." << std::endl;
  write2D(image2Dred, fulloutputfilename2Dred);

  LabelMapToBinaryImageFilterType::Pointer labelMapToBinaryImageFilter = LabelMapToBinaryImageFilterType::New();
  labelMapToBinaryImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
  labelMapToBinaryImageFilter->Update();
  BinaryImageType2D::Pointer binaryimage2Dbacteria = maxintprojection(labelMapToBinaryImageFilter->GetOutput());
  std::string outputfilenamebinary2Dbacteria = seriesreader.getFilename(seriesnr, "_d_bacteria_binary.tiff");
  std::string fulloutputfilenamebinary2Dbacteria = outputdirectory + outputfilenamebinary2Dbacteria; 
  std::cout << "Writing file: " << fulloutputfilenamebinary2Dbacteria << " ..." << std::endl;
  write2D(binaryimage2Dbacteria, fulloutputfilenamebinary2Dbacteria);

  ImageType2D::Pointer image2DReadMean = maxintprojection(image3DbacteriaReadMean);
  RGBFilterType::Pointer colormapImageFilter = RGBFilterType::New();
  colormapImageFilter->SetInput(image2DReadMean);
  colormapImageFilter->SetColormap( RGBFilterType::Jet );
  colormapImageFilter->Update();
  std::string outputfilename2DReadMean = seriesreader.getFilename(seriesnr, "_e_ReadMean.tiff");
  std::string fulloutputfilename2DReadMean = outputdirectory + outputfilename2DReadMean; 
  WriterTypeRGB::Pointer writerRGB = WriterTypeRGB::New();
  writerRGB->SetFileName( fulloutputfilename2DReadMean );             
  writerRGB->SetInput(colormapImageFilter->GetOutput());
  std::cout << "Writing file: " << fulloutputfilename2DReadMean << " ..." << std::endl;
  writerRGB->Update();


  return EXIT_SUCCESS;
}



