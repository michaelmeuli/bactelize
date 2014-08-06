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
#include "itkSCIFIOImageIO.h"
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
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelMapToAttributeImageFilter.h"
#include "itkStatisticsLabelObjectAccessors.h"
#include "itkLabelImageToStatisticsLabelMapFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkLabelMapToBinaryImageFilter.h"

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
typedef itk::BinaryImageToLabelMapFilter<BinaryImageType3D> BinaryImageToLabelMapFilterType;
typedef itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, BinaryImageType3D> LabelMapToLabelImageFilterType;
typedef itk::LabelImageToShapeLabelMapFilter <BinaryImageType3D> LabelImageToShapeLabelMapFilterType;
typedef itk::LabelImageToStatisticsLabelMapFilter< BinaryImageType3D, ImageType3D > LabelImageToStatisticsLabelMapFilterType;
typedef itk::LabelMapToAttributeImageFilter< LabelImageToStatisticsLabelMapFilterType::OutputImageType, ImageType3D,
   itk::Functor::MeanLabelObjectAccessor<LabelImageToStatisticsLabelMapFilterType::OutputImageType::LabelObjectType> > L2ImageType;
typedef itk::LabelMapToBinaryImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, BinaryImageType3D> LabelMapToBinaryImageFilterType;
typedef itk::RGBPixel<unsigned char>   RGBPixelType;
typedef itk::Image<RGBPixelType, 2>    RGBImageType;
typedef itk::ScalarToRGBColormapImageFilter<ImageType2D, RGBImageType> RGBFilterType;
typedef itk::ImageFileWriter<RGBImageType> WriterTypeRGB;


ImageType2D::Pointer       maxintprojection(ImageType3D::Pointer, unsigned int projectionDirection = 2);
BinaryImageType2D::Pointer maxintprojection(BinaryImageType3D::Pointer, unsigned int projectionDirection = 2);
void dumpmetadatadic(ImageType5D::Pointer);
void setSpacing(ImageType5D::Pointer, float x, float y, float z, float t, float c);  //obsolete
void setSpacing(BinaryImageType3D::Pointer image3D, float x, float y, float z);
void printSpacing(BinaryImageType3D::Pointer);
ImageType3D::Pointer extractchannel(ImageType5D::Pointer, int channelnr);
void printHistogram(ImageType3D::Pointer);
void write2D(ImageType2D::Pointer, std::string filenamepath);
void write2D(BinaryImageType2D::Pointer, std::string filenamepath);


class SeriesReader {   
  public:
    SeriesReader(std::string inputFileName);
    ImageType5D::Pointer get5DImage(int series); 
    void dumpimageio();
    int getSeriesStart();
    int getSeriesEnd();
    std::string getFilename(int seriesnr, std::string suffix);

  private:    
    itk::SCIFIOImageIO::Pointer m_io;
    ReaderType::Pointer m_reader;
    std::string m_inputFileName;
    StreamingFilter::Pointer m_streamer;
    ImageType5D::Pointer m_image5D;
    int m_seriesStart;
    int m_seriesEnd;
};


#endif
