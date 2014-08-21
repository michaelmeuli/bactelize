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
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <cstdio>


int nucleichannel      = 0;
int bacteriachannel    = 1;
int lysosomechannel    = 2;
int binaryLowerThresholdBacteria = 200;
int minNumberOfPixels = 100;
std::ofstream fileout;
std::string fileoutName = "AA_results.txt";



int getdir(std::string dir, std::vector<std::string> &files) {
  errno = 0;
  DIR *dp;
  struct dirent *dirp;
  if ((dp = opendir(dir.c_str())) == NULL) {
    int errsv = errno;
    std::cout << "Error(" << errsv << ") opening " << dir << std::endl;
    std::cout << "Text version of the error code: " << std::strerror(errsv) << std::endl;
    return EXIT_FAILURE;
    }
  while ((dirp = readdir(dp)) != NULL) {
    struct stat st_buf;
    if (stat((dir+dirp->d_name).c_str(), &st_buf)) {
      int errsv = errno;
      std::cout << "Error(" << errsv << ") getting file status." << std::endl;
      std::cout << "Text version of the error code: " << std::strerror(errsv) << std::endl;
      return EXIT_FAILURE;
      }
    if (S_ISREG(st_buf.st_mode)) {
      files.push_back(dir+dirp->d_name);
      }
    }
  closedir(dp);
  return 0;
}




int main( int argc, char * argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << argv[0] << " inputdirectory" << " outputdirectory" << std::endl;
    return EXIT_FAILURE;
    }    
  std::string inputdirectory = argv[1];
  std::string outputdirectory = argv[2];

  std::string fulloutputfilenameResults = outputdirectory + fileoutName;
  fileout.open(fulloutputfilenameResults.c_str()); 
  fileout << "Bactelize!\n";
  fileout.close();
  
  std::string dir = std::string(inputdirectory);
  std::vector<std::string> files = std::vector<std::string>();
  if (getdir(dir,files)) {
    std::exit(EXIT_FAILURE);   
    }
  std::cout << std::endl;
  std::cout << "Files found in directory:" << std::endl;
  for (unsigned int i = 0;i < files.size();i++) {
    std::cout << "  " << files[i] << std::endl;
    } 
  for (unsigned int i = 0;i < files.size();i++) {
    processSeries(files[i], outputdirectory);
    fileout.open(fulloutputfilenameResults.c_str(), std::ofstream::app); 
    fileout << std::endl << std::endl;
    fileout.close();
    }

  return EXIT_SUCCESS;
}

//    dumpimageio(seriesreader.getReader());





