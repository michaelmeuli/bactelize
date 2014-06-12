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
#include "QuickView.h"
#include "vnl/vnl_math.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageToHistogramFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

  typedef unsigned short InputPixelType;
  typedef itk::Image< InputPixelType, 5 > ImageType5D;
  typedef itk::Image< InputPixelType, 4 > ImageType4D;
  typedef itk::Image< InputPixelType, 3 > ImageType3D;
  typedef itk::Image< InputPixelType, 2 > ImageType2D;
  typedef itk::ImageFileReader< ImageType5D > ReaderType;
  typedef itk::StreamingImageFilter< ImageType5D, ImageType5D > StreamingFilter;

  typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
  typedef itk::ImageSliceConstIteratorWithIndex< ImageType3D > SliceIteratorType;


ImageType2D::Pointer maxintprojection(ImageType3D::ConstPointer inputImageMIP);


#endif
