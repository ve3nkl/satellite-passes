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

Read configuration file and provide access to the information contained in it.

"""

import yaml
import datetime
from satools.SAGeo import MHGridSquare

class Sat:
  
  def __init__(self, name, designator, number, nickname):
    self.name = name
    self.designator = designator
    self.number = number
    self.nickname = nickname

class QTH:
  
  qth_by_title = {}
  
  @classmethod
  def new_qth_latlong(cls, title, latitude, longitude, altitude, tz_offset_dst, tz_offset_std, declination, time_zone_name):
    qth = cls(title, tz_offset_dst, tz_offset_std, declination, time_zone_name)
    qth.set_coordinates_latlong(latitude, longitude, altitude)
    QTH.qth_by_title[title] = qth
    return qth
    
  @classmethod
  def new_qth_maiden(cls, title, maidenhead, altitude, tz_offset_dst, tz_offset_std, declination, time_zone_name):
    qth = cls(title, tz_offset_dst, tz_offset_std, declination, time_zone_name)
    qth.set_coordinates_maiden(maidenhead, altitude)
    QTH.qth_by_title[title] = qth
    return qth
    
  @classmethod
  def zone_offset_to_delta(cls, s):
    sign = 1
    if len(s) > 0 and s[0:1] == "-":
      sign = -1
      s = s[1:]
    if len(s) == 5 and s[0:2].isdigit() and s[2:3] == ":" and ( s[3:5] in ["00", "15", "30", "45"] ):
      hours = int(s[0:2])
      mins  = int(s[3:5])
      total = ( hours * 60 + mins ) * 60
      if sign < 0:
        total = - total
      return datetime.timedelta(seconds=total)
    else:
      return None
  
  def __init__(self, title, tz_offset_dst, tz_offset_std, declination, time_zone_name):
    self.title = title
    self.tz_odst = QTH.zone_offset_to_delta(tz_offset_dst)
    if self.tz_odst is None:
      raise Exception("'time-zone-offset.daylight' value is invalid: " + tz_offset_dst)
    self.tz_ostd = QTH.zone_offset_to_delta(tz_offset_std)
    if self.tz_ostd is None:
      raise Exception("'time-zone-offset.standard' value is invalid: " + tz_offset_std)
    self.declanation = declination
    self.time_zone_name = time_zone_name
    
  def set_coordinates_latlong(self, latitude, longitude, altitude):
    self.latitude  = latitude
    self.longitude = longitude
    self.altitude  = altitude
    self.maidenhead = MHGridSquare.lat_long_to_maiden(latitude, longitude)
    
  def set_coordinates_maiden(self, maidenhead, altitude):
    center = MHGridSquare(maidenhead).getCenterSpot()
    self.latitude  = center.latitude
    self.longitude = center.longitude
    self.altitude  = altitude
    self.maidenhead = maidenhead
    
  def get_declination(self):
    return self.declination

class SAConfig:
  
  def __init__(self, conf_file_name):
    self.conf_file_name = conf_file_name
    
    try:
      f = open(self.conf_file_name, "r")
      conf_data = f.read()
      f.close()
    except Exception:
      raise Exception("Configuration file not found (" + self.conf_file_name + ")")            
      
    try:
      if str(yaml.__version__) >= "5.1":
        self.conf_yaml = yaml.load(conf_data, Loader=yaml.FullLoader)
      else:
        self.conf_yaml = yaml.load(conf_data)
    except yaml.YAMLError as e:
      raise Exception("Configuration file has syntax errors. " + str(e))
            
    if not "awa-sat-config" in self.conf_yaml:
      raise Exception("Configuration file has invalid format")
    
    if "tle-directory" in self.conf_yaml["awa-sat-config"]:
      self.tle_directory = self.conf_yaml["awa-sat-config"]["tle-directory"]
    else:
      raise Exception("'tle-directory' is not defined in configuration file")
      
    for qth_descr in self.conf_yaml["awa-sat-config"]["qth-list"]:
      if "latitude" in qth_descr["location"]:
        QTH.new_qth_latlong(qth_descr["id"],
                            qth_descr["location"]["latitude"],
                            qth_descr["location"]["longitude"],
                            qth_descr["location"]["altitude"],
                            qth_descr["time-zone-offset"]["daylight"],
                            qth_descr["time-zone-offset"]["standard"],
                            qth_descr["magnetic-declination"],
                            qth_descr["time-zone-name"])
      elif "maidenhead" in qth_descr["location"]:
        QTH.new_qth_maiden(qth_descr["id"],
                            qth_descr["location"]["maidenhead"],
                            qth_descr["location"]["altitude"],
                            qth_descr["time-zone-offset"]["daylight"],
                            qth_descr["time-zone-offset"]["standard"],
                            qth_descr["magnetic-declination"],
                            qth_descr["time-zone-name"])
      else:
        pass
    
    self.sat_list = []    
    self.sat_by_name = {}
    self.sat_by_nickname = {}
    for sat_descr in self.conf_yaml["awa-sat-config"]["satellite-list"]:
      self.sat_list.append(Sat(sat_descr["name"], sat_descr["inter-id"], sat_descr["norad-id"], sat_descr["nickname"]))
      self.sat_by_name[sat_descr["name"]] = self.sat_list[-1]
      self.sat_by_nickname[sat_descr["nickname"]] = self.sat_list[-1]
    self.sat_list = sorted(self.sat_list, key=lambda x: x.name)
    
  
  def get_QTH_list(self):
    r = QTH.qth_by_title.values()
    r = sorted(r, key=lambda x: x.title)
    return r      
    
  def get_QTH_by_title(self, title):
    if title in QTH.qth_by_title:
      return QTH.qth_by_title[title]
    else:
      return None
      
  def get_tle_directory(self):
    return self.tle_directory
    
  def get_sat_list(self):
    return self.sat_list
    
  def get_sat_by_name(self, name):
    if name in self.sat_by_name:
      return self.sat_by_name[name]
    else:
      return None
      
  def get_sat_by_nickname(self, nickname):
    if nickname in self.sat_by_nickname:
      return self.sat_by_nickname[nickname]
    else:
      return None
  
  def get_tle_sources(self):
    return self.conf_yaml["awa-sat-config"]["tle-sources"]
