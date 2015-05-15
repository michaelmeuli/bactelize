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


#ifndef BACTELIZE_H
#define BACTELIZE_H

#include "itkSCIFIOImageIO.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkStreamingImageFilter.h"
#include "itkMetaDataObject.h"
#include "itkMetaDataDictionary.h"
#include "itkImageIOBase.h"
#include "itkExtractImageFilter.h"
#include "itksys/SystemTools.hxx"
#include "vnl/vnl_math.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageToHistogramFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIterator.h"

#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToAttributeImageFilter.h"
#include "itkStatisticsLabelObjectAccessors.h"
#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkTIFFImageIO.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkEuclideanDistanceMetric.h"
#include "itkLabelMap.h"
#include "itkStatisticsLabelObject.h"
#include "itkShapeLabelObject.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <fstream>
#include <set>


typedef unsigned short InputPixelType;
typedef unsigned char  BinaryPixelType;
typedef itk::Image< InputPixelType, 5 > ImageType5D;
typedef itk::Image< InputPixelType, 4 > ImageType4D;
typedef itk::Image< InputPixelType, 3 > ImageType3D;
typedef itk::Image< InputPixelType, 2 > ImageType2D;
typedef itk::Image< double, 3 >   DoubleImageType3D;
typedef itk::Image< BinaryPixelType, 3 >   BinaryImageType3D;
typedef itk::Image< BinaryPixelType, 2 >   BinaryImageType2D;
typedef itk::ImageFileReader< ImageType5D > ReaderType;
typedef itk::StreamingImageFilter< ImageType5D, ImageType5D > StreamingFilter;
typedef itk::ExtractImageFilter< ImageType5D, ImageType4D > ExtractFilterType5D4D;
typedef itk::ExtractImageFilter< ImageType4D, ImageType3D > ExtractFilterType4D3D;
typedef itk::Statistics::ImageToHistogramFilter<ImageType3D>   HistogramFilterType;
typedef HistogramFilterType::HistogramType  HistogramType;
typedef itk::NormalizeImageFilter< ImageType3D, DoubleImageType3D >  NormalizeFilterType;
typedef itk::Image<unsigned char, 2>  ImageTypeWriter;
typedef itk::ImageFileWriter< ImageTypeWriter > WriterType;
typedef itk::RescaleIntensityImageFilter< ImageType2D, ImageTypeWriter >  RescaleFilterTypeWriter;
typedef itk::RescaleIntensityImageFilter< DoubleImageType3D, ImageType3D >  RescaleFilterTypeNormalized;
typedef itk::MedianImageFilter< ImageType3D, ImageType3D > MedianFilterType;
typedef itk::BinaryThresholdImageFilter< ImageType3D, BinaryImageType3D >  BinaryFilterType;
typedef itk::BinaryImageToShapeLabelMapFilter< BinaryImageType3D, itk::LabelMap< itk::StatisticsLabelObject< ImageType3D::PixelType, ImageType3D::ImageDimension > > > BinaryImageToShapeLabelMapFilterType;
typedef itk::StatisticsLabelMapFilter< BinaryImageToShapeLabelMapFilterType::OutputImageType, ImageType3D> StatisticsLabelMapFilterType;
typedef itk::LabelMapToAttributeImageFilter< StatisticsLabelMapFilterType::OutputImageType, ImageType3D,
   itk::Functor::MeanLabelObjectAccessor< StatisticsLabelMapFilterType::OutputImageType::LabelObjectType> > L2ImageType;
typedef itk::LabelMapToBinaryImageFilter< BinaryImageToShapeLabelMapFilterType::OutputImageType, BinaryImageType3D > LabelMapToBinaryImageFilterType;
typedef itk::RGBPixel <unsigned char >   RGBPixelType;
typedef itk::Image< RGBPixelType, 2 >    RGBImageType;
typedef itk::ScalarToRGBColormapImageFilter< ImageType2D, RGBImageType > RGBFilterType;
typedef itk::ImageFileWriter< RGBImageType > WriterTypeRGB;
typedef itk::TIFFImageIO TIFFIOType;
typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorTypeInput;
typedef itk::ImageSliceConstIteratorWithIndex< ImageType3D > SliceIteratorTypeInput;
typedef itk::ImageLinearIteratorWithIndex< BinaryImageType2D > LinearIteratorTypeBinary;
typedef itk::ImageSliceConstIteratorWithIndex< BinaryImageType3D > SliceIteratorTypeBinary;
typedef itk::Vector< double, 3 > MeasurementVectorType;
typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;
typedef itk::FixedArray< double, 3 > ImageSizeType;
typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
typedef TreeGeneratorType::KdTreeType TreeType;
typedef TreeType::KdTreeNodeType      NodeType;
typedef itk::Statistics::EuclideanDistanceMetric< MeasurementVectorType > DistanceMetricType;


ImageType2D::Pointer       maxintprojection(ImageType3D::Pointer, unsigned int projectionDirection = 2);
BinaryImageType2D::Pointer maxintprojection(BinaryImageType3D::Pointer, unsigned int projectionDirection = 2);
void dumpmetadatadic(ImageType5D::Pointer);
void dumpimageio(ReaderType::Pointer);
void printImInfo(BinaryImageType3D::Pointer);
ImageType3D::Pointer extractchannel(ImageType5D::Pointer, int channelnr);
void printHistogramNormalized(ImageType3D::Pointer);
void write2D(ImageType2D::Pointer, std::string filenamepath);
void write2D(BinaryImageType2D::Pointer, std::string filenamepath);
std::string getFilename(std::string inputFileName, int seriesnr, int seriesCount, std::string suffix = "");
ImageSizeType getImSize(BinaryImageToShapeLabelMapFilterType::Pointer, ImageSizeType);
void printObjectLabelAndCentroid(BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType*);
void printObjectInfo(BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType*);
void printShapeLabelObjects(BinaryImageToShapeLabelMapFilterType::Pointer);
SampleType::Pointer getCentroidsAsSample(BinaryImageToShapeLabelMapFilterType::Pointer);
void excludeIfSet(BinaryImageToShapeLabelMapFilterType::Pointer, std::set<int>);
BinaryImageType3D::Pointer getBinaryIm(ImageType3D::Pointer);  
void excludeSmallObjects(BinaryImageToShapeLabelMapFilterType::Pointer, double minNumberOfmm3);
void excludeClusters(BinaryImageToShapeLabelMapFilterType::Pointer, int maxclustersize);
void excludeOutside(BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, 
	            ImageType3D::Pointer image3Dcell);
void writeMeanResults(StatisticsLabelMapFilterType::Pointer statisticsLabelMapFilter, 
                      BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter, 
	              ImageType3D::Pointer image3Dred, 
                      std::string fulloutputfilenameResults, std::string seriesName);
int processSeries(std::string inputFileName, std::string outputdirectory, bool vflag, bool tflag, int fileNr, int seriesNr);


#endif


