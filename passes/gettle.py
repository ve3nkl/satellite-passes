"""
Copyright (c) 2016-2020, 2021 Nikolai Ozerov (VE3NKL)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

-----------------------------------------------------------------------------

This script downloads TLE files from the specified sources.

"""

import sys, os, traceback
import argparse
import datetime
import urllib.request as request
import yaml
sys.path.append('../')
import satools.SAConfig as sc


def epoc2dt(epocstr):
  year = int("20" + epocstr[0:2])
  k = epocstr.find(".")
  if k > -1:
    days = int(epocstr[2:k])
    res = datetime.datetime(year, 1, 1, tzinfo=None) + datetime.timedelta(days=days-1)
    return res
  else:
    raise Exception("Invalid epoch value: " + epocstr)

if __name__ == "__main__":

# Parse parameters and store them in the local variables

  parser = argparse.ArgumentParser(description='Update satellite TLE files.')
  parser.add_argument("-v", "--version", action="store_true", default=False, help="display current version")
  parser.add_argument("-c", "--configuration",  action="store", type=str, default="", help="specify path to configuration file")
  parser.add_argument("-o", "--days-old", action="store", type=int, default=2, help="how old the data can become before it has to be updated")
  group = parser.add_mutually_exclusive_group(required=False)
  group.add_argument("-i", "--info",  action="store_true", default=False, help="only display info on existing TLE files, no downloads")
  group.add_argument("-f", "--force",  action="store_true", default=False, help="download TLEs from the web regardless of how recent the existing data is")
  args = parser.parse_args()
  
  if args.version:
    print("gettle.py verion 1.1")
    print("Copyright (c) 2016-2020, 2021 Nikolai Ozerov (VE3NKL)")
    quit(0)

  if args.configuration == "":
    config_file_name = "./data/satellites.conf"
  else:
    config_file_name = args.configuration
  
  try:
    config = sc.SAConfig(config_file_name)
  except Exception as e:
    print( "ERROR: " + str(e) )
    traceback.print_exc()
    print( )
    quit(1)
    
  tle_dir = config.get_tle_directory()
    
  if not os.path.exists(tle_dir):
    print( "ERROR: TLE directory specified in the configuration file does not exist: " + tle_dir + "." )
    print( )
    quit(1) 
    
  tle      = {}                     # New TLEs will be saved here
  xtle     = {}                     # Epoch for existing TLEs will be saved here
  cur_time = datetime.datetime.utcnow()

        
  # Read existing tle files

  tle_c_total = 0
  tle_c_outdated = 0
  oldest = 0
  
  no_tle_data = 0

  for sat in config.get_sat_list():
    infile = os.path.join(tle_dir, sat.nickname + ".tle")

    data_exists = False
    
    if os.path.exists(infile):
    
      lines = None
      with open(infile, "r") as f:
        lines = f.readlines()
        if len(lines)>= 3:
          tle_c_total = tle_c_total + 1
          data_exists = True
      
      if data_exists:
        sat_norad_id = lines[1][2:7]
        sat_epoc = lines[1][18:32]
        if sat.number == sat_norad_id:
          xtle[sat.number] = sat_epoc
          old = (cur_time-epoc2dt(sat_epoc)).days
          if old > oldest:
            oldest = old
        else:
          print( "ERROR: Existing TLE file for satellite " + sat.name + " (file name " + infile + ") has inconsistent NORAD ID (" + sat_norad_id + ")." )
          
        if args.info and data_exists:  
          print( "    Satellite "  + sat.name + " existing TLE data is " + str(old) + " days old." )
      else:
        no_tle_data += 1
        
    else:
      no_tle_data += 1

  if not args.info:           # 'Display info only' was NOT requested

    if not args.force:
      if oldest <= args.days_old and no_tle_data == 0:
        print( "    All TLE files contain up-to-date data. The least recent is " + str(oldest) + " days old." )
        print( "" )
        quit(0)
  
    # TLE source URL

    tle_data = []

    for source in config.get_tle_sources():
      for url in source["urls"]:
  
        print( "Loading TLE source data from URL: " + url + " ..." )
  
        try:
          with request.urlopen(url) as tle_file_response:
            tle_file_data     = tle_file_response.read()
            tle_file_text = tle_file_data.decode('utf-8')
            data_list = tle_file_text.split("\n")
            tle_data = tle_data + data_list
        except Exception as e:
          print( "    ERROR: error reading the specified URL: " + url + "." )
          print( str(e) )
          traceback.print_exc()
          print( "" )
          quit(1)
  
        print( "    successfully loaded." )
      
    print( "Parsing loaded data ..." )
   
  # Parse the loaded file
  
    line_type = 0
    line_title = ""
    line_1 = ""
    line_2 = ""
    cat_num = ""
    n = 0
    for line in tle_data:
      n = n + 1
      line = line.rstrip()
      if len(line) > 1:
        if line[0:2] == "1 ":
          if line_type != 1:
            print( "    ERROR: Parsing error of loaded TLE data. Wrong line sequence (1). Line #" + str(n) + "." )
            print( "    [" + line + "]" )
            line_type = 0            
          else:
            line1 = line
            cat_num = line[2:7]
            epoch   = line[18:32]
            line_type = 2
        elif line[0:2] == "2 ":
          if line_type != 2:
            print( "    ERROR: Parsing error of loaded TLE data. Wrong line sequence (2). Line #" + str(n) + "." )
            print( "    [" + line + "]" )
            line_type = 0            
          else:
            line2 = line
            if cat_num != "":
              if cat_num in tle:
                if tle[cat_num][0] < epoch:
                  tle[cat_num] = [epoch, line_title, line1, line2]
              else:
                tle[cat_num] = [epoch, line_title, line1, line2]
                
            line_type = 0                        
        else:
          if line_type != 0:
            print( "    ERROR: Parsing error of loaded TLE data. Wrong line sequence (3). Line #" + str(n) + "." )
            print( "    [" + line + "]" )
            line_type = 0          
          
          line_title = line
          line_type = 1
  
    print( "    parsing completed." )
    
  # Process the data
  
    for sat in config.get_sat_list():
      outfile = os.path.join(tle_dir, sat.nickname + ".tle")
      
      if sat.number in tle:
        if not (sat.number in xtle) or (tle[sat.number][0] > xtle[sat.number]):
          with open(outfile, "w") as f:
            f.write(tle[sat.number][1] + "\n")
            f.write(tle[sat.number][2] + "\n")
            f.write(tle[sat.number][3] + "\n")          
          
          print( "    TLE data for satellite " + sat.name + " has been replaced with more recent data (" + epoc2dt(tle[sat.number][0]).strftime("%b %d, %Y") + ")." )
        else:
          old = (cur_time-epoc2dt(xtle[sat.number])).days
          if old <= args.days_old:
            print( "    TLE data for satellite " + sat.name + " is up-to-date." )
          else:
            print( "    TLE data for satellite " + sat.name + " is " + str(old) + " days old but more recent data is not available." )
      else:
        print( "ERROR: No TLE data found for satellite " + sat.name + " (catalog number " + str(sat.number) + ")" )
  
