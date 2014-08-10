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
float xspacing = 0.33;
float yspacing = 0.33;
float zspacing = 1.2;
float tspacing = 1.0;
float cspacing = 1.0;
int binaryLowerThresholdBacteria = 200;
int minNumberOfPixels = 100;
std::ofstream fileout;

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

  std::string outputfilenamefileout = outputdirectory + "fileout.txt";
  fileout.open(outputfilenamefileout.c_str()); 
  fileout << "Bactelize!\n";
  fileout.close();

  SeriesReader seriesreader(fullinputfilename, outputdirectory);
  dumpimageio(seriesreader.getReader());
  seriesreader.calculateSeries();
  
 



  return EXIT_SUCCESS;
}



