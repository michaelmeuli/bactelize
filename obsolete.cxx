void setSpacing(ImageType5D::Pointer image5D, float x, float y, float z, float t, float c) {  //obsolete
  // SetSpacing
  std::cout << "--== Correcting spacing and setting origin ==--" << std::endl;
  ImageType5D::SpacingType spacing;
  spacing[0] = x;  
  spacing[1] = y;  
  spacing[2] = z;  
  spacing[3] = t; 
  spacing[4] = c;  
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
  }


void setSpacing(BinaryImageType3D::Pointer image3D, float x, float y, float z) {  
  std::cout << "--== Setting spacing ==--" << std::endl;
  BinaryImageType3D::SpacingType spacing;
  spacing[0] = x;  
  spacing[1] = y;  
  spacing[2] = z;  
  image3D->SetSpacing( spacing );
  }
