#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "error.h"
#include "region.h"

void readRegions(const std::string& input_file, uint32_t max_regions, const std::string& chrom_limit, std::vector<Region>& regions, std::ostream& logger){
  logger << "Reading region file " << input_file << std::endl;
  std::ifstream input(input_file.c_str());
  if (!input.is_open()) 
    printErrorAndDie("Failed to open region file");

  regions.clear();
  std::string line;
  int32_t num_regions = 0;
  while (std::getline(input, line) && regions.size() < max_regions){
    num_regions++;
    std::istringstream iss(line);
    std::string chrom, name;
    int32_t start, stop;
    int period;
    double ref_copy;
    if (!(iss >> chrom >> start >> stop >> period >> ref_copy))
      printErrorAndDie("Improperly formatted region file. \nRequired format is tab-delimited columns CHROM START STOP PERIOD NCOPIES\n Bad line: " + line);
    if (start < 1)      printErrorAndDie("Improperly formatted region file. \n Region has a START < 1, but START must be >= 1\n Bad line: " + line);
    if (stop <= start)  printErrorAndDie("Improperly formatted region file. \n Region has a STOP <= START. Bad line: " + line);
    if (period < 1)     printErrorAndDie("Improperly formatted region file. \n Region has a PERIOD < 1. Bad line: " + line);
    if (period > 9)     printErrorAndDie("Improperly formatted region file. \n Region has a PERIOD > 9. Bad line: " + line);

    if (!chrom_limit.empty() && chrom.compare(chrom_limit) != 0)
      continue;
    if (iss >> name)
      regions.push_back(Region(chrom, start-1, stop, period, name));
    else
      regions.push_back(Region(chrom, start-1, stop, period));
  }
  input.close();
  logger << "Region file contains " << num_regions << " regions";
  if (!chrom_limit.empty())
    logger << ", of which " << regions.size() << " were located on the requested chromosome";
  logger << "\n" << std::endl;

  if (!chrom_limit.empty() && regions.empty())
    printErrorAndDie("Region file " + input_file + " did not contain any regions on the requested chromosome: " + chrom_limit);
}

void orderRegions(std::vector<Region>& regions){
  std::sort(regions.begin(), regions.end());
}
