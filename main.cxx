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
#include <algorithm>


int cellchannel      = 0;
int bacteriachannel    = 1;
int lysosomechannel    = 2;
int binaryLowerThresholdBacteria = 200;
int cellInsideThreshold = 1;
double minNumberOfmm3 = 0.000000001;   // 1um3
double maxNumberOfmm3 = 0.000000004;
int maxclustersize = 2;   // if there are more, they get excluded
std::ofstream fileout;
std::string fileoutName = "AA_results.txt";
int numberOfStreamDivisions = 4;
int numberOfBins = 50;
double maxSingleObjectDiameter = 0.006;
double minSearchRadius = 0.006;



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
  std::sort(files.begin(), files.end());
  closedir(dp);
  return 0;
}



int fail(char *argv[]) {
  std::cerr << "\nUsage: " << argv[0] << " inputdirectory outputdirectory [OPTIONS]\n\n"
  "OPTIONS:\n"
  "-t <n1 n2>\n\tTest file n1 found in inputdirectory and series n2 of this file.\n"
  "-v\n\tVerbose output.\n\n";
  return EXIT_FAILURE;
}



int main(int argc, char *argv[]) {

  if( argc < 3) {
    return fail(argv);
    }

  std::string inputdirectory = argv[1];
  std::string outputdirectory = argv[2];
  std::string inputfile = "";

  bool vflag = false;
  bool tflag = false;
  bool fflag = false;
  int fileNr = 0;
  int seriesNr = 0;

  // parse flags
  for (int i = 3; i < argc; i++) {
    if (strcmp (argv[i], "-v") == 0) {
      vflag = true;
      }
    else if (strcmp (argv[i], "-t") == 0) {
      if (i + 2 >= argc) {
        return fail(argv);
        }
      fileNr = atoi(argv[i+1]);
      seriesNr = atoi(argv[i+2]);
      tflag = true;
      i+=2;
      }
    else if (strcmp (argv[i], "-f") == 0) {
      if (i + 1 >= argc) {
        return fail(argv);
        }
      inputfile = argv[i+1];
      fflag = true;
      i+=1;
      }
    }

  std::stringstream ssout;
  for (int i=0; i<(argc); i++) {
    ssout << argv[i];
    ssout << " ";
  }
  std::string fulloutputfilenameResults = outputdirectory + fileoutName;
  fileout.open(fulloutputfilenameResults.c_str(), std::ofstream::app); 
  fileout << ssout.str() << std::endl;
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
  std::cout << std::endl;
  if (tflag) {
    processSeries(files[fileNr], outputdirectory, vflag, tflag, fileNr, seriesNr);
    }
  else if (fflag) {
    processSeries(inputfile, outputdirectory, vflag, tflag, fileNr, seriesNr);
    }
  else {
    for (unsigned int i = 0; i < files.size();i++) {
      processSeries(files[i], outputdirectory, vflag, tflag, fileNr, seriesNr);
      fileout.open(fulloutputfilenameResults.c_str(), std::ofstream::app); 
      fileout << std::endl << std::endl;
      fileout.close();
      }
    }

  return EXIT_SUCCESS;
}




