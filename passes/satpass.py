"""
Copyright (c) 2012-2020, 2021 Nikolai Ozerov (VE3NKL)

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

- History -------------------------------------------------------------------
2021/07/28   1.2   Calculate orbit # for elliptical orbits
2021/07/31   1.3   Determine timezone automatically if it is not 
                   specified explicitely
-----------------------------------------------------------------------------

This script calculates and prints details on satellite passes as a yaml document.

"""
import sys, os, traceback
import argparse
import getopt
import yaml
import datetime
import decimal
sys.path.append('../')
import satools.SAConfig as sc
import satools.SAGeo as sg
from psgp4gm import psgp4gm
import pytz
from timezonefinder import TimezoneFinder


"""
  Convert datetime to normalized epoc (see psgp4gm.epoc2norm for details)
"""
def dt2norm(dt):
  year, month, day, hour, minute, second, wday, yday, dst = dt.timetuple()
  if year < 1960 or year > 2059:
    raise Exception("datetime is out of valid range: " + str(dt))
  e = yday + ((( hour * 60 + minute ) * 60 + second) / (24 * 60 * 60))
  return psgp4gm.epoc2norm(year, e)

"""  
  Validate date in the format yyyy-mm-dd
"""
def valid_date(s):
  formats = [
    "%b %d,%Y",
    "%b %d, %Y",
    "%Y/%m/%d",
    "%Y-%m-%d"
  ]
  valid = False
  for pattern in formats:
    try:
      d = datetime.datetime.strptime(s, pattern)
      valid = True
      break
    except ValueError:
      pass
  if valid:
    return d.replace(hour=0, minute=0, second=0, microsecond=0, tzinfo=None)
  msg = "Not a valid date: '{0}'.".format(s)
  raise argparse.ArgumentTypeError(msg)
  
"""  
  Validate hour interval, for example: 1-3, 9-18, 23-05.
"""
def valid_interval(s):
  tokens = s.split("-")
  if len(tokens) == 2:
    frm = tokens[0]
    to  = tokens[1]
    if frm.isdigit() and to.isdigit():
      n_from = int(frm)
      n_to   = int(to)
      if n_from >= 0 and n_from <= 23 and n_to >= 0 and n_to <= 23:
        return (n_from, n_to)
  msg = "Not a valid hour interval: '{0}'.".format(s)
  raise argparse.ArgumentTypeError(msg)
  
"""       
  Validate time zone offset (-12:00 to 12:00 with 15 minute intervals)
"""
def valid_zone_offset(s):
  sign = 1
  if len(s) > 0 and s[0:1] == "-":
    sign = -1
    s = s[1:]
  if len(s) == 4 and s[0:2].isdigit() and ( s[2:4] in ["00", "15", "30", "45"] ):
    hours = int(s[0:2])
    mins  = int(s[2:4])
    total = ( hours * 60 + mins ) * 60
    if sign < 0:
      total = - total
    return datetime.timedelta(seconds=total)
  msg = "Not a valid time zone offset: '{0}'.".format(s)
  raise argparse.ArgumentTypeError(msg)

"""
  Validate maidenhead grid (4, 6 or 8 characters)
"""
def valid_grid(s):
  try:
    g = sg.MHGridSquare(s)
  except Exception as e:
    print( "ERROR: " + str(e) )
    quit(1)
  return g 
  
"""  
  Validate altitude
"""
def valid_altitude(s):
  in_meters = True
  if len(s) >= 2 and s[-1:] == "m":
    s = s[0:-1]
  elif len(s) >= 3 and s[-2:] == "ft":
    in_meters = False
    s = s[0:-2]
  sign = 1
  if len(s) > 1 and s[0:1] == "-":
    sign = -1
    s = s[1:]
  if s.isdigit():
    a = int(s) * sign
    if not in_meters:
      a = a * 0.3048
      in_meters = True
    return a
  msg = "Not a valid altitude: '{0}'.".format(s)
  raise argparse.ArgumentTypeError(msg)

"""
  Main line logic. It is called when the module is called as a command with parameters and it can also be called
  directly when the module is imported. This function returns a tuple containing a return code, a list of lines for
  stdout and a list of lines for stderr.
"""
def main_line(tle_dir, sat_objects, spot_name, spot_object, altitude, time_zone_delta, starting, period, max_elevation, between_hours):

  rc = 0
  out = []
  err = []
    
  # Check all parameters for correctness

  # TLE directory

  if not os.path.exists(tle_dir):
    err.append("Specified TLE directory does not exist: " + tle_dir + ".")
    err.append("")
    rc = 2
    return (rc, out, err)

  # Process satellite TLE files                   

  tle_data = []

  for sat in sat_objects:
    tle_name = sat.nickname + ".tle"
    
    try:
      with open(os.path.join(tle_dir,tle_name), "r") as f:
        t_data = f.readlines()
      if len(t_data) < 3:
        err.append("TLE file " + os.path.join(tle_dir,tle_name) + " contains invalid data.")
        err.append("")
        rc = 2
        return (rc, out, err)
        
      tle_data.append((sat.nickname,t_data))
      
    except:
      err.append("Error opening TLE file " + os.path.join(tle_dir,tle_name) + ".")
      err.append("")
      rc = 2
      return (rc, out, err)                              

  # Calculate starting and ending date/time of the working interval (in UTC)
  
  spot_tz = datetime.timezone(time_zone_delta, name=spot_name)
  utc_tz  = datetime.timezone(datetime.timedelta(seconds=0), name="UTC")

  start_dt = starting.replace(tzinfo=spot_tz)      # Theses dates are in the 'spot' time zone
  end_dt   = starting + datetime.timedelta(days=period+1)
  
  start_utc = start_dt.astimezone(utc_tz).replace(tzinfo=None)
  end_utc   = end_dt.astimezone(utc_tz).replace(tzinfo=None)
  
  current_utc = datetime.datetime.utcnow()      # Current time
  if start_utc < current_utc:                   # Adjust start_utc if necessary
    start_utc = current_utc
  
  # Let's process all satellites ...

  obj = []

  for sat_name, data in tle_data:
    
    try:
    
      psgp4gm.init_wgs72(data[1].encode(), data[2].encode(), spot_object.latitude, spot_object.longitude, altitude)
      
      current_utc = start_utc    
      
      while(current_utc <= end_utc):
        r = psgp4gm.find_rise_set(current_utc.year, current_utc.month, current_utc.day, current_utc.hour, current_utc.minute, current_utc.second, 0, 5)

        rise_datetime = datetime.datetime(r[0][0], r[0][1], r[0][2], r[0][3], r[0][4], r[0][5], 0, None)
        rise_azimuth  = decimal.Decimal(str(r[2]))
        max_datetime  = datetime.datetime(r[3][0], r[3][1], r[3][2], r[3][3], r[3][4], r[3][5], 0,None)
        max_angle     = decimal.Decimal(str(r[4]))
        max_azimuth   = decimal.Decimal(str(r[5]))
        set_datetime  = datetime.datetime(r[6][0], r[6][1], r[6][2], r[6][3], r[6][4], r[6][5], 0, None)
        set_azimuth   = decimal.Decimal(str(r[8]))
      
        if rise_datetime <= end_utc:
      
          rise_datetime2 = rise_datetime.replace(tzinfo=utc_tz).astimezone(spot_tz).replace(tzinfo=None)
          max_datetime2  = max_datetime.replace(tzinfo=utc_tz).astimezone(spot_tz).replace(tzinfo=None)
          set_datetime2  = set_datetime.replace(tzinfo=utc_tz).astimezone(spot_tz).replace(tzinfo=None)
      
          if time_zone_delta.total_seconds == 0:
            tzstr = "+0000"
          else:
            total_seconds = time_zone_delta.seconds + time_zone_delta.days*3600*24
            if total_seconds < 0:
              sign = "-"
              total_seconds = - total_seconds
            else:
              sign = "+"
            tzstr =  sign + str(total_seconds // 3600).rjust(2, "0") + ":"
            tzstr += str(total_seconds % 3600 // 60).rjust(2, "0")
      
          if r[4] >= max_elevation:
            
            if (between_hours[0] == between_hours[1]) or \
               ((between_hours[0] < between_hours[1]) and (set_datetime2.hour >= between_hours[0] and rise_datetime2.hour < between_hours[1])) or \
               ((between_hours[0] > between_hours[1]) and (set_datetime2.hour >= between_hours[0] or rise_datetime2.hour < between_hours[1])):
      
              # Calculate orbit number.
              orbit_s = ""
              orbit_n = psgp4gm.calculate_orbit_number(dt2norm(rise_datetime))
              if orbit_n >= 0:
                orbit_s = str(orbit_n)
              
              o = {'sat-pass':{'QTH':spot_name.strip(),
                   'satellite':data[0].strip(),
                   'satellite-user-name': sat_name,
                   'time-zone': tzstr,
                   'rise-date':rise_datetime2.strftime("%b %d, %Y"),
                   'rise-orbit': orbit_s,
                   'rise-info':{'time':rise_datetime2.strftime("%H:%M:%S"), 'azimuth':str(rise_azimuth.quantize(decimal.Decimal('.1')))},
                   'max-info':{'time':max_datetime2.strftime("%H:%M:%S"), 'azimuth':str(max_azimuth.quantize(decimal.Decimal('.1'))), 'elevation':str(max_angle.quantize(decimal.Decimal('.1')))},
                   'set-info':{'time':set_datetime2.strftime("%H:%M:%S"), 'azimuth':str(set_azimuth.quantize(decimal.Decimal('.1')))} }}
              obj.append(o)
      
        current_utc = set_datetime + datetime.timedelta(minutes=10)
      
    except Exception as e:
      print( "ERROR: " + str(e))
      print( "       during attempt to calculate pass for " + sat_name ) 
      traceback.print_exc()
      print( )
      quit(1)

  out = yaml.dump(obj, default_flow_style=False).split("\n")
  return (rc, out, err)

if __name__ == "__main__":

  # Parse parameters and store them in the local variables

  parser = argparse.ArgumentParser(description='Calculating Satellite Pass Parameters.')
  parser.add_argument('satellite', nargs='+')
  parser.add_argument("-v", "--version", action="store_true", default=False, help="display current version")
  parser.add_argument("-c", "--configuration",  action="store", type=str, default="", help="specify path to configuration file")
  parser.add_argument("-q", "--qth", action="store", type=str, default="", help="name of QTH for pass calculations (defined in configuration file)")
  parser.add_argument("-s", "--starting", action="store", type=valid_date, 
                      default=datetime.datetime.now().replace(hour=0, minute=0, second=0, microsecond=0, tzinfo=None),
                      help="starting date for pass calculations")
  parser.add_argument("-p", "--period", action="store", type=int, default=1, help="number of days to calculate passes for")
  parser.add_argument("-m", "--maximum-elevation", action="store", type=int, default=0, help="ignore passes with maximum elevation below this value")
  parser.add_argument("-b", "--between-hours", action="store", type=valid_interval, default=(0, 0), help="ignore passes outside of the specified interval (hourFrom-hourTo)" )
  parser.add_argument("-z", "--time-zone-offset", action="store", type=valid_zone_offset, help="if specified along with --qth, the sepcified time offset overrides the offset specified in QTH")
  parser.add_argument("-a", "--altitude", action="store", type=valid_altitude, help="altitude, for ex.: 100, 100m or 328ft")
  group = parser.add_mutually_exclusive_group(required=False)
  group.add_argument("-g", "--maidenhead-grid", action="store", type=valid_grid, help="if specified along with --qth, the specifeid grid overrides the location in QTH")
  group.add_argument("-n", "--standard-time", action="store_true", default=False, help="use standard time zone offset defined in QTH definition (see configuration file)")
  group.add_argument("-y", "--daylight-time", action="store_true", default=False, help="use daylight time zone offset defined in QTH definition (see configuration file)")
  
  args = parser.parse_args()
  
  if args.version:
    print("satpass.py verion 1.2")
    print("Copyright (c) 2012-2020, 2021 Nikolai Ozerov (VE3NKL)")
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
  
  sat_objects = []
  if len(args.satellite) == 1 and args.satellite[0].upper() == "ALL":
    for sat in config.get_sat_list():
      sat_objects.append(sat)
  else:
    for sat_name in args.satellite:
      sat = config.get_sat_by_nickname(sat_name)
      if sat is None:
        sat = config.get_sat_by_name(sat_name)
      if sat is not None:
        sat_objects.append(sat)
      else:
        print( "ERROR: Unrecognized satellite name: " + sat_name + "." )
        print( )
        quit(1) 
  
  qth_object = None
  
  if args.qth != "":
    qth_object = config.get_QTH_by_title(args.qth)
    if qth_object is None:
      print( "ERROR: Unrecognized QTH name: " + args.qth + "." )
      print( )
      quit(1) 

  if args.maidenhead_grid is not None:
    spot_name   = args.maidenhead_grid.orig_code
    spot_object = args.maidenhead_grid.getCenterSpot()
  elif qth_object is not None:
    spot_name   = qth_object.title
    spot_object = sg.Spot(qth_object.latitude, qth_object.longitude)
  else:
    print( "ERROR: Either --maidenhead-grid or --qth option needs to be specified." )
    print( )
    quit(1) 
    
  if args.time_zone_offset is not None:
    time_zone_delta = args.time_zone_offset
  elif qth_object is not None:
    if args.standard_time:
      time_zone_delta = qth_object.tz_ostd
    elif args.daylight_time:
      time_zone_delta = qth_object.tz_odst
    else:
      print( "ERROR: When --time_zone-offset is NOT specified and --qth option is specified, either --satndard_time or --daylight-time option is required." )
      print( )
      quit(1) 
  else:
    # Neither time zone nor qth profile is specified, try to detect the time zone based 
    # on the specified location.
    tf = TimezoneFinder()
    timezone_str = tf.timezone_at(lat=spot_object.latitude, lng=spot_object.longitude)
    
    if timezone_str is not None:
      tz_name = timezone_str
      
      tz = None
      try:
        tz = pytz.timezone(timezone_str)
      except pytz.exceptions.UnknownTimeZoneError:
        pass
      
      if tz is not None:
        current_utc       = datetime.datetime.utcnow()
        current_localized = tz.localize(current_utc)
        time_zone_delta = current_localized.utcoffset()
      else:
        print( "ERROR: Time zone could not be determined, specify either --time-zone-offset or --qth option." )
        print( )
        quit(1)    
    else:
      print( "ERROR: Time zone could not be determined, specify either --time-zone-offset or --qth option." )
      print( )
      quit(1)    
    
  if args.altitude is not None:
    altitude = args.altitude
  elif qth_object is not None:
    altitude = qth_object.altitude
  else:
    altitude = 0
  
  starting = args.starting
  period   = args.period
  max_elevation  = args.maximum_elevation
  between_hours  = args.between_hours
  
  """
  # Uncomment for debugging ...
  print( tle_dir )
  print( spot_name )
  print( spot_object )
  print( altitude )
  print( time_zone_delta )
  print( starting )
  print( period )
  print( max_elevation )
  print( between_hours )
  """
  
  rc, out, err = main_line(tle_dir, sat_objects, spot_name, spot_object, altitude, time_zone_delta, starting, period, max_elevation, between_hours)
  if rc == 0:
    for line in out:
      print(line)
  else:
    for line in err:
      print(line, file=sys.stderr)
    quit(rc)
