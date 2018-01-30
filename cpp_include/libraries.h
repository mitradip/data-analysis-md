#ifndef LIBRARIES_H
#define LIBRARIES_H
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <algorithm>
typedef std::string string;
std::string trim_left(const std::string& str)
{
  const std::string pattern = " \f\n\r\t\v";
  return str.substr(str.find_first_not_of(pattern));
}

//
//Right trim
//
std::string trim_right(const std::string& str)
{
  const std::string pattern = " \f\n\r\t\v";
  return str.substr(0,str.find_last_not_of(pattern) + 1);
}



//
//Left and Right trim
//
std::string trim(const std::string& str)
{
  return trim_left(trim_right(str));
}

long int line_skip(string file_name,long int start_pos,long int lines)
{
  std::fstream file;
  file.open(file_name,std::ios::in);
  file.seekg(start_pos,std::ios::beg);
  string xxxx;
  for(long int i=0;i<lines;++i)
    getline(file,xxxx);
  long int stop_pos=file.tellg();
  file.close();
  return stop_pos;
}



#endif
