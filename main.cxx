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
//  setSpacing(image5D, xspacing, yspacing, zspacing, tspacing, cspacing);

  ImageType3D::Pointer image3Dnuclei = extractchannel(image5D, nucleichannel);
  ImageType3D::Pointer image3Dbacteria = extractchannel(image5D, bacteriachannel);
  ImageType3D::Pointer image3Dred = extractchannel(image5D, lysosomechannel);
  
  NormalizeFilterType::Pointer  normalizeFilter = NormalizeFilterType::New();
  normalizeFilter->SetInput( image3Dbacteria );
  RescaleFilterTypeNormalized::Pointer rescaleNormalized = RescaleFilterTypeNormalized::New();
  rescaleNormalized->SetInput( normalizeFilter->GetOutput() ); 
  rescaleNormalized->SetOutputMinimum(0);
  rescaleNormalized->SetOutputMaximum(4095);

  rescaleNormalized->Update();
  printHistogram( rescaleNormalized->GetOutput() );
  std::cout << std::endl;

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
  

  typedef itk::BinaryImageToLabelMapFilter<BinaryImageType3D> BinaryImageToLabelMapFilterType;
  BinaryImageToLabelMapFilterType::Pointer binaryImageToLabelMapFilter = BinaryImageToLabelMapFilterType::New();
  binaryImageToLabelMapFilter->SetInput(binaryimage3Dbacteria);
  binaryImageToLabelMapFilter->Update();
  std::cout << "There are " << binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " label map objects." << std::endl;

  typedef itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, BinaryImageType3D> LabelMapToLabelImageFilterType;
  LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
  labelMapToLabelImageFilter->SetInput(binaryImageToLabelMapFilter->GetOutput());
  labelMapToLabelImageFilter->Update();

  typedef itk::LabelImageToShapeLabelMapFilter <BinaryImageType3D> LabelImageToShapeLabelMapFilterType;
  LabelImageToShapeLabelMapFilterType::Pointer labelImageToShapeLabelMapFilter = LabelImageToShapeLabelMapFilterType::New ();
  labelImageToShapeLabelMapFilter->SetInput( labelMapToLabelImageFilter->GetOutput() );
  labelImageToShapeLabelMapFilter->Update();


  std::cout << "There are " << labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() << " shape label map objects." << std::endl << std::endl;
  std::vector<unsigned long> labelsToRemove;
  for(unsigned int i = 0; i < labelImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    LabelImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    labelObject->Print(std::cout, 5);
    std::cout << std::endl;
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
  binaryImageToLabelMapFilter->Update();
  labelMapToLabelImageFilter->Update();

  typedef itk::LabelImageToStatisticsLabelMapFilter< BinaryImageType3D, ImageType3D > LabelImageToStatisticsLabelMapFilterType;
  LabelImageToStatisticsLabelMapFilterType::Pointer labelImageToStatisticsLabelMapFilter = LabelImageToStatisticsLabelMapFilterType::New();
  labelImageToStatisticsLabelMapFilter->SetFeatureImage(image3Dred);
  labelImageToStatisticsLabelMapFilter->SetInput(labelMapToLabelImageFilter->GetOutput());
  labelImageToStatisticsLabelMapFilter->Update();

// assert binaryImageToLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() == labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects() 
  unsigned int bacteriacount = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();
  unsigned int colcount = 0; 
  for(unsigned int i = 0; i < labelImageToStatisticsLabelMapFilter->GetOutput()->GetNumberOfLabelObjects(); i++) {
    LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType* labelObject = labelImageToStatisticsLabelMapFilter->GetOutput()->GetNthLabelObject(i);
    labelObject->Print(std::cout, 4);
    double mean = labelObject->GetMean();
    std::cout << "Mean value of object with label " << (int)(labelObject->GetLabel()) << " in lysosomechannel: " << mean << std::endl;   //todo c++ style cast
    if ( mean > meanRedThreshold ) {
      colcount++;
      }
    std::cout << std::endl;    
    }
  std::cout << "Total bacteria counted (in statistics label map): " << bacteriacount << std::endl;
  std::cout << "Bacteria colocalizing in lysosomechannel: " << colcount << std::endl;


  typedef itk::LabelMapToAttributeImageFilter< LabelImageToStatisticsLabelMapFilterType::OutputImageType, ImageType3D,
     itk::Functor::MeanLabelObjectAccessor<LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType> > L2ImageType;
  L2ImageType::Pointer labelMapToReadMeanImage = L2ImageType::New();
  labelMapToReadMeanImage->SetInput( labelImageToStatisticsLabelMapFilter->GetOutput() );
  labelMapToReadMeanImage->Update();
  ImageType3D::Pointer image3DbacteriaReadMean = labelMapToReadMeanImage->GetOutput();
  
  typedef itk::ImageFileWriter< ImageType3D > WriterType3D;
  WriterType3D::Pointer writer3D = WriterType3D::New();
  std::string outputfilename3DreadMean = seriesreader.getFilename(seriesnr, "_readMean.vtk");
  std::string fulloutputfilename3DreadMean = outputdirectory + outputfilename3DreadMean; 
  writer3D->SetFileName( fulloutputfilename3DreadMean );             
  writer3D->SetInput(image3DbacteriaReadMean);
  std::cout << "Writing file: " << fulloutputfilename3DreadMean << " ..." << std::endl;
  writer3D->Update();

  typedef itk::VectorImage< InputPixelType, 3 > VectorImageType3D;
  typedef itk::ComposeImageFilter< ImageType3D, VectorImageType3D > ComposeImageFilterType;
  ComposeImageFilterType::Pointer composeImageFilter = ComposeImageFilterType::New();
  composeImageFilter->SetInput( 0, image3Dnuclei );
  composeImageFilter->SetInput( 1, image3Dbacteria );
  composeImageFilter->SetInput( 2, image3Dred );
  composeImageFilter->SetInput( 3, image3DbacteriaReadMean );
  composeImageFilter->Update();
  VectorImageType3D::Pointer vectorImage = composeImageFilter->GetOutput();

  typedef itk::ImageFileWriter< VectorImageType3D > WriterTypeVector;
  WriterTypeVector::Pointer writerVector = WriterTypeVector::New();
  std::string outputfilename3Dvector = seriesreader.getFilename(seriesnr, "_vector.vtk");
  std::string fulloutputfilename3Dvector = outputdirectory + outputfilename3Dvector; 
  writerVector->SetFileName(fulloutputfilename3Dvector);             
  writerVector->SetInput(vectorImage);
  std::cout << "Writing file: " << fulloutputfilename3Dvector << " ..." << std::endl;
  writerVector->Update();


  typedef itk::ImageFileWriter< ImageType3D > WriterTypeSCIFIO;
  typename WriterTypeSCIFIO::Pointer writerSCIFIO = WriterTypeSCIFIO::New();
  std::string outputfilename3DSCIFIO = seriesreader.getFilename(seriesnr, "_readMean.ome.tiff");
  std::string fulloutputfilename3DSCIFIO = outputdirectory + outputfilename3DSCIFIO; 
  writerSCIFIO->SetFileName(fulloutputfilename3DSCIFIO);
  writerSCIFIO->SetInput(image3DbacteriaReadMean);
  itk::SCIFIOImageIO::Pointer ioOut = itk::SCIFIOImageIO::New();
  ioOut->DebugOn();
  writerSCIFIO->SetImageIO(ioOut);
  std::cout << "Writing file: " << fulloutputfilename3DSCIFIO << " ..." << std::endl;
  writerSCIFIO->Update();
  

  typedef itk::ImageFileWriter< VectorImageType3D > WriterTypeVectorSCIFIO;
  typename WriterTypeVectorSCIFIO::Pointer writerVectorSCIFIO = WriterTypeVectorSCIFIO::New();
  std::string outputfilename3DvectorSCIFIO = seriesreader.getFilename(seriesnr, "_vector.ome.tiff");
  std::string fulloutputfilename3DvectorSCIFIO = outputdirectory + outputfilename3DvectorSCIFIO; 
  writerVectorSCIFIO->SetFileName(fulloutputfilename3DvectorSCIFIO);
  writerVectorSCIFIO->SetInput(vectorImage);
  itk::SCIFIOImageIO::Pointer ioOutV = itk::SCIFIOImageIO::New();
  ioOutV->DebugOn();
  writerVectorSCIFIO->SetImageIO(ioOutV);
  std::cout << "Writing file: " << fulloutputfilename3DvectorSCIFIO << " ..." << std::endl;
  writerVectorSCIFIO->Update();



  ImageType2D::Pointer image2Dbacteria = maxintprojection(image3Dbacteria);
  std::string outputfilename2Dbacteria = seriesreader.getFilename(seriesnr, "_bacteria.tiff");
  std::string fulloutputfilename2Dbacteria = outputdirectory + outputfilename2Dbacteria; 
  std::cout << "Writing file: " << fulloutputfilename2Dbacteria << " ..." << std::endl;
  write2D(image2Dbacteria, fulloutputfilename2Dbacteria);
  
  ImageType2D::Pointer image2Dred = maxintprojection(image3Dred);
  std::string outputfilename2Dred = seriesreader.getFilename(seriesnr, "_lysosome.tiff");
  std::string fulloutputfilename2Dred = outputdirectory + outputfilename2Dred; 
  std::cout << "Writing file: " << fulloutputfilename2Dred << " ..." << std::endl;
  write2D(image2Dred, fulloutputfilename2Dred);
  
  BinaryImageType2D::Pointer binaryimage2Dbacteria = maxintprojection(binaryimage3Dbacteria);
  std::string outputfilenamebinary2Dbacteria = seriesreader.getFilename(seriesnr, "_bacteria_binary.tiff");
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
