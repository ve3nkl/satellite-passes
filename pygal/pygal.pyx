#cython: language_level=3

"""
Copyright (c) 2022 Nikolai Ozerov (VE3NKL)

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
2022/08/02   1.0   Make changes to increase precision for certain calculations
                   as well as to make python interface more friendly by using
                   classes (and structures)
-----------------------------------------------------------------------------

"""

from cpython.datetime cimport datetime, timedelta

"""
Wrapping up gravity model structure initialization and other C routines.

"""

cdef extern from "string.h":

  void strcpy(char* str1, char* str2)
  
cdef extern from "math.h":

  double sin(double)
  double asin(double)
  double fabs(double)
  double atan2(double, double)
  
cdef extern from "gal_math_macros.h":

  cdef double GAL_R2D, GAL_D2R
  cdef double GAL_PI
  
cdef extern from "gal_pn.h":
  
  void gal_pn(double p[3], double *r, double u[3])
    
cdef extern from "gal_pm.h":
  
  double gal_pm(double p[3])
  
cdef extern from "gal_rxp.h":
  
  void gal_rxp(double m[3][3], double p1[3], double p2[3])  
  
cdef extern from "gal_pn.h":
  
  void gal_pn(double p[3], double* mod, double u[3])    
  
cdef extern from "gal_pmp.h":
  
  void gal_pmp(double p1[3], double p2[3], double p3[3])    
  
cdef extern from "gal_pxp.h":          # Vector cross-product
  
  void gal_pxp(double a[3], double b[3], double axb[3])
  
cdef extern from "gal_pdp.h":          # Vector dot-product
  
  double gal_pdp(double a[3], double b[3])
  
cdef extern from "gal_sepp.h":
  
  double gal_sepp(double a[3], double b[3])
  
cdef extern from "gal_status_t.h":

  ctypedef struct gal_status_t:
    int errc
    pass
    
cdef extern from "gal_stsalloc.h":

  gal_status_t* gal_stsalloc()
  
cdef extern from "gal_stsfree.h":

  void gal_stsfree(gal_status_t* status)  
  
cdef extern from "gal_emparams.h":

  void gal_emparams(int em, double* sma, double* f, gal_status_t* status)
  
cdef extern from "gal_me2bf.h":

  void gal_me2bf(double merf[2][3], double w, double omega, double bfrf[2][3])
  
cdef extern from "gal_gmst82.h":

  double gal_gmst82(double e1, double e2)

cdef extern from "gal_gm.h":

  ctypedef struct gal_gm_t:
    double gm
    pass
    
  ctypedef enum:
    GAL_GMEA_WGS72
    pass    
    
  ctypedef enum:    
    GAL_NORMALIZED
    pass
    
cdef extern from "gal_frame_macros.h":

  ctypedef enum:
    GAL_EMEA_WGS1984
    pass
    
cdef extern from "gal_gmfree.h":

  void gal_gmfree(gal_gm_t* gm)
    
cdef extern from "gal_gmget_wgs72.h":

  gal_gm_t* gal_gmget_wgs72(int maxn, int maxm, int norm, gal_status_t* status)

cdef extern from "gal_sgp4gm.h":

  void gal_sgp4gm(gal_gm_t *gm, double *tumin, double *mu, double *re, double *xke, double *j2, double *j3, double *j4, double *j3oj2, gal_status_t* status)
  
cdef extern from "gal_sgp4.h":

  ctypedef struct gal_sgp4_t:
    int satnum
    double time
    double pv[2][3]
    pass
  
  ctypedef struct gal_tle_t:
    int             satnum         # US Strategic Command Object Number          
    char            classification # Security Classification                     
    char            intldesg[12]   # International Designator
    int             epochyr        # Epoch Year                                  
    double          epochdays      # Epoch Day of Year (plus Fraction)           
    double          date1          # UTC date part 1 - normalized                
    double          date2          # UTC date part 2 - normalized                
    double          ndot           # Mean motion derivative (rev/day /2)         
    double          nddot          # Mean motion second derivative (rev/day2 /6) 
    double          bstar          # Bstar / Drag Term                           
    int             ephtype        # Ephermeris Type                             
    int             setnum         # Element set number                          
    double          inclo          # Inclination                                 
    double          nodeo          # Right Ascension of Ascending Node (deg)     
    double          ecco           # Eccentricity                                
    double          argpo          # Argument of Perigee (deg)                   
    double          mo             # Mean Anamaly (deg)                          
    double          no             # Mean Motion (rev/day)                       
    int             revnum         # Epoch Revolution Number                       

  void gal_sgp4(gal_sgp4_t *sgp4, double epoch1, double epoch2, double pv[2][3], gal_status_t* status)
  
  void gal_sgp4init(gal_gm_t *gm, gal_tle_t *tle, gal_sgp4_t *sgp4, gal_status_t* status )
  
cdef extern from "gal_sgp4rs.h":  

  void gal_sgp4rs(double utc1, double utc2, double latitude, double longitude, double height, double minel,
                 double timefwd, double gm, double re, double inf, gal_sgp4_t *sgp4, double passinfo[3][3],
                 gal_status_t* status)
  
cdef extern from "gal_tledec.h":
  
  ctypedef struct gal_tle_t:
    pass
  
  void gal_tledec(char *card1, char *card2, gal_tle_t *tle, gal_status_t* status)
  
cdef extern from "gal_days2jd.h":
  
  void gal_days2jd(int year, double days, double *djm0, double *djm, gal_status_t *status)  
  
cdef extern from "gal_cal2jd.h":

  void gal_cal2jd(int iy, int im, int id, double *djm0, double *djm, gal_status_t* status)
  
cdef extern from "gal_jd2cal.h":

  void gal_jd2cal(double dj1, double dj2, int *iy, int *im, int *id, double *fd, gal_status_t* status)

cdef extern from "gal_gc2gd.h":
    
  void gal_gc2gd(int n, double xyz[3], double *elong, double *phi, double *height, gal_status_t *status)
    
cdef extern from "gal_emdetails.h":    

  void gal_emdetails(int em, int *body, char *name, double *sma, double *inf, gal_status_t* status) 

cdef extern from "gal_gd2gce.h":

  void gal_gd2gce(double a, double f, double elong, double phi, double height, double xyz[3], gal_status_t* status)
  
cdef extern from "gal_sezm.h":

  void gal_sezm(double latitude, double longitude, double sezmatrix[3][3])
  
cdef extern from "gal_anp.h":

  double gal_anp(double angle)
  
cdef extern from "gal_anpm.h":

  double gal_anpm(double angle)  

cdef extern from "gal_gsupv00.h":
  
  int gal_gsupv00(double epoch1, double epoch2, double pv[2][3], gal_status_t* status)
  
cdef extern from "gal_c2tpv00a.h":
  
  void gal_c2tpv00a(double gcrf[2][3], double utc1, double utc2, double dut1, double lod, double xp, double yp, double itrf[2][3], gal_status_t *status)

cdef extern from "gal_gd2gc.h":
  void gal_gd2gc(int n, double elong, double phi, double height, double xyz[3], gal_status_t *status)
  
cdef extern from "gal_ta2ea.h":
  
  double gal_ta2ea(double ta, double ecc)
  

cdef class SunPosition:
  
   cdef int    visible
   cdef double azimuth
   cdef double elevation
   
   def __init__(self, visible: int, azimuth: double, elevation: double):
     self.visible   = visible
     self.azimuth   = azimuth
     self.elevation = elevation


cdef class RiseSetInfo:

  cdef int      aos_orbit_number
  cdef datetime aos_time     
  cdef double   aos_azimuth  
  cdef datetime max_time     
  cdef double   max_azimuth  
  cdef double   max_elevation
  cdef datetime los_time     
  cdef double   los_azimuth  
  
  def __init__(self, aos_orbit_number: int, aos_time: datetime, aos_azimuth: double, 
               max_time: datetime, max_azimuth: double, max_elevation: double,
               los_time: datetime, los_azimuth: double):
    self.aos_orbit_number = aos_orbit_number
    self.aos_time         = aos_time     
    self.aos_azimuth      = aos_azimuth  
    self.max_time         = max_time     
    self.max_azimuth      = max_azimuth  
    self.max_elevation    = max_elevation
    self.los_time         = los_time     
    self.los_azimuth      = los_azimuth  
  
  def get_orbit_number(self):
    return self.aos_orbit_number
    
  def get_aos_time(self):
    return self.aos_time
    
  def get_aos_azimuth(self):
    return self.aos_azimuth
    
  def get_max_time(self):
    return self.max_time
    
  def get_max_azimuth(self):
    return self.max_azimuth
    
  def get_max_elevation(self):
    return self.max_elevation
    
  def get_los_time(self):
    return self.los_time
    
  def get_los_azimuth(self):
    return self.los_azimuth
    
  def __str__(self):
    return "           " + "         " + "   " + str(round(self.max_elevation,1)).rjust(4) + "\n" + \
           self.aos_time.strftime("%Y/%m/%d %H:%M:%S ")+ \
           self.max_time.strftime("%H:%M:%S ")+self.los_time.strftime("%H:%M:%S\n") + \
           ("#"+str(self.aos_orbit_number)).rjust(10)+"   "+str(round(self.aos_azimuth,1)).rjust(5)+"    "+ \
           str(round(self.max_azimuth,1)).rjust(5)+"    "+str(round(self.los_azimuth,1)).rjust(5)
           

cdef class VisibleSatPositionError:
  
  cdef timedelta time_error       # If None, do not compare time
  cdef double    azimuth_error    # if negative, do not compare azimuth
  cdef double    elevation_error  # if negative, do not compare elevation
  cdef double    distance_error   # if negative, do not compare distance
  
  def __init__(self, time_error, azimuth_error, elevation_error, distance_error):
    self.time_error      = time_error
    self.azimuth_error   = azimuth_error
    self.elevation_error = elevation_error
    self.distance_error  = distance_error
    
  def get_time_error(self):
    return self.time_error
    
  def get_azimuth_error(self):
    return self.azimuth_error
    
  def get_elevation_error(self):
    return self.elevation_error
    
  def get_distance_error(self):
    return self.distance_error
    

cdef class VisibleSatPosition:
  
  cdef datetime time
  cdef double   azimuth
  cdef double   elevation
  cdef double   distance
  
  def __init__(self):
    self.time      = None     # While time remains None the object is considered
    self.azimuth   = 0.0      # to remain NOT initialized.
    self.elevation = 0.0
    self.distance  = 0.0
    
  def is_initialized(self):
    return self.time is not None
    
  def get_time(self):
    return self.time
    
  def get_azimuth(self):
    return self.azimuth
    
  def get_elevation(self):
    return self.elevation
    
  def get_distance(self):
    return self.distance
    
  def copy_from(self, other: VisibleSatPosition):
    self.time      = other.time
    self.azimuth   = other.azimuth
    self.elevation = other.elevation
    self.distance  = other.distance
    
  def close_to(self, other, allowed_error):
    if allowed_error.get_time_error() is not None:
      if self.time is not None and other.get_time() is not None:
        if self.time >= other.get_time():
          if self.time - other.get_time() > allowed_error.get_time_error():
            return False
        else:
          if other.get_time() - self.time > allowed_error.get_time_error():
            return False
      elif self.time is None and other.get_time() is None: # Both objects are NOT
        return True                                        # initialized, they are 
      else:                                                # within the error threshold
        return False
        
    if allowed_error.get_azimuth_error() >= 0:
      if self.azimuth >= other.get_azimuth():
        if self.azimuth - other.get_azimuth() > allowed_error.get_azimuth_error():
          return False
      else:
        if other.get_azimuth() - self.azimuth > allowed_error.get_azimuth_error():
          return False
          
    if allowed_error.get_elevation_error() >= 0:
      if self.elevation >= other.get_elevation():
        if self.elevation - other.get_elevation() > allowed_error.get_elevation_error():
          return False
      else:
        if other.get_elevation() - self.elevation > allowed_error.get_elevation_error():
          return False
          
    if allowed_error.get_distance_error() >= 0:
      if self.distance >= other.get_distance():
        if self.distance - other.get_distance() > allowed_error.get_distance_error():
          return False
      else:
        if other.get_distance() - self.distance > allowed_error.get_distance_error():
          return False
          
    return True
    
    
  def __str__(self):
    return self.time.strftime("%Y/%m/%d %H:%M:%S.%f")[0:22]+" "+str(round(self.azimuth,2))+" "+str(round(self.elevation,2))


cdef class VisibleSatOrbit:
  
  cdef double r_latitude
  cdef double r_longitude
  cdef double height
  cdef double omega, mu, sma, inf, re, gm
  cdef double observer[3]
  cdef double sez[3][3]
  
  cdef gal_tle_t tle
  cdef gal_sgp4_t sgp4
  
  cdef bytes b_card1, b_card2
  
  cdef double tle_norm_epoch


  """
    Initialize object using tle cards and current station position on Earth
  """
  def __init__(self, b_card1: bytes, b_card2: bytes, latitude: double, longitude: double, height: double):
  
    cdef gal_gm_t* gm72  
      
    cdef int rc
    cdef double tumin
    cdef double xke
    cdef double j2
    cdef double j3
    cdef double j4
    cdef double j3oj2 
    
    cdef gal_status_t* status
  
    cdef double epoch2
    cdef double r
    
    cdef char* card1
    cdef char* card2
    
    self.r_latitude  = latitude * GAL_D2R
    self.r_longitude = longitude * GAL_D2R
    self.height    = height
    
    self.b_card1 = b_card1
    self.b_card2 = b_card2
    
    omega = 7.292115146706979e-05
    status = gal_stsalloc()
  
    gm72 = gal_gmget_wgs72(180, 180, GAL_NORMALIZED, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_gmget_wgs72 returned code " + str(rc)
      
    self.gm = gm72.gm
  
    gal_emparams(GAL_EMEA_WGS1984, &self.sma, &self.inf, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_emparams returned code " + str(rc)
    
    gal_sgp4gm(gm72, &tumin, &self.mu, &self.re, &xke, &j2, &j3, &j4, &j3oj2, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_sgp4gm returned code " + str(rc)  
    
    gal_gd2gce(self.sma, self.inf, self.r_longitude, self.r_latitude, self.height, self.observer, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_gd2gce returned code " + str(rc)    
    
    gal_sezm(self.r_latitude, self.r_longitude, self.sez)
    
    card1 = b_card1
    card2 = b_card2
    gal_tledec(card1, card2, &self.tle, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_tledec returned code " + str(rc)      
    
    """
    Convert TLE epoc time to normalized form.
    """
    gal_days2jd(self.tle.epochyr, self.tle.epochdays, &self.tle_norm_epoch, &epoch2, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_days2jd returned code " + str(rc)
    self.tle_norm_epoch = self.tle_norm_epoch + epoch2
    
    gal_sgp4init(gm72, &self.tle, &self.sgp4, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_sgp4init returned code " + str(rc)        
      
    gal_gmfree(gm72)    
    gal_stsfree(status)
    
  """
    Find the next satellite rise information starting at some point in time.
    Ignore passes that do not satisfy the minimum elevation criterion.
  """
  def find_rise_set(self, t: datetime, min_elevation: int, days_fwd: int):
  
    cdef gal_sgp4_t sgp4l
  
    cdef gal_status_t* status
    
    cdef double r_min_elevation
    cdef double days
    
    cdef int iy1, im1, id1, iy2, im2, id2, iy3, im3, id3 
    cdef double time1, time2, time3
    cdef int ihh1, imm1, iss1, ihh2, imm2, iss2, ihh3, imm3, iss3   
    
    cdef double passinfo[3][3]  
    cdef double epoch1, epoch2
    cdef int  rc
    
    status = gal_stsalloc()  
  
    gal_cal2jd(t.year, t.month, t.day, &epoch1, &epoch2, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_cal2jd returned code " + str(rc)              
      
    epoch2 = epoch2 + t.hour / 24.0 + t.minute / 24.0 / 60.0 + t.second / 24.0 /60.0 / 60.0
    
    sgp4l = self.sgp4
    
    r_min_elevation = min_elevation * GAL_D2R
    days            = days_fwd
  
    gal_sgp4rs(epoch1, epoch2, 
               self.r_latitude, self.r_longitude, self.height, 
               r_min_elevation, days, 
               self.gm, self.sma, 1.0 / self.inf, 
               &sgp4l, passinfo, status)
               
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_sgp4rs returned code " + str(rc)                
    
    gal_jd2cal(passinfo[0][0], 0.0, &iy1, &im1, &id1, &time1, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_jd2cal returned code (1) " + str(rc)                  
    ihh1 = (int)(24.0 * time1)
    imm1 = (int)((time1 * 24.0 - ihh1) * 60.0)
    iss1 = (int)((time1 * 24.0 * 60.0  - ihh1 * 60.0 - imm1) * 60.0)
    
    gal_jd2cal(passinfo[1][0], 0.0, &iy2, &im2, &id2, &time2, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_jd2cal returned code (2)" + str(rc)                    
    ihh2 = (int)(24.0 * time2)
    imm2 = (int)((time2 * 24.0 - ihh2) * 60.0)
    iss2 = (int)((time2 * 24.0 * 60.0  - ihh2 * 60.0 - imm2) * 60.0)  
    
    gal_jd2cal(passinfo[2][0], 0.0, &iy3, &im3, &id3, &time3, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_jd2cal returned code (3)" + str(rc)                    
    ihh3 = (int)(24.0 * time3)
    imm3 = (int)((time3 * 24.0 - ihh3) * 60.0)
    iss3 = (int)((time3 * 24.0 * 60.0  - ihh3 * 60.0 - imm3) * 60.0)   
    
    gal_stsfree(status)   
    
    aos_time      = datetime(iy1, im1, id1, ihh1, imm1, iss1)
    aos_azimuth   = float(passinfo[0][1] * GAL_R2D)
    los_time      = datetime(iy3, im3, id3, ihh3, imm3, iss3)
    los_azimuth   = float(passinfo[2][2] * GAL_R2D)
    max_time      = datetime(iy2, im2, id2, ihh2, imm2, iss2)
    max_azimuth   = float(passinfo[1][2] * GAL_R2D)
    max_elevation = float(passinfo[1][1] * GAL_R2D)
  
    return RiseSetInfo(
             self.calculate_orbit_number(passinfo[0][0]),
             datetime(iy1, im1, id1, ihh1, imm1, iss1),
             passinfo[0][2] * GAL_R2D,
             datetime(iy2, im2, id2, ihh2, imm2, iss2),
             passinfo[1][2] * GAL_R2D,
             passinfo[1][1] * GAL_R2D,
             datetime(iy3, im3, id3, ihh3, imm3, iss3),
             passinfo[2][2] * GAL_R2D
          )
    

  """
    Calculate satellite position as observed from a location on Earth at the
    time moment 'time'
  """
  def get_position(self, time: datetime, vsp: VisibleSatPosition):
  
    cdef double e1, e2  
    cdef gal_status_t* status  
    cdef double pv[2][3]
    cdef double s
    
    cdef double efg[2][3]
    
    cdef double sr[3]
    cdef double srunit[3]
    cdef double srmod
    cdef double rho[3]
    cdef double az, sel, elev, dist
    
    cdef int rc
    
    status = gal_stsalloc()  
    
    gal_cal2jd(time.year, time.month, time.day, &e1, &e2, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_cal2jd returned code " + str(rc)      
      
    e1 = e1 + e2
    e2 = (time.hour + time.minute / 60.0 + time.second / 3600.0 + time.microsecond / 3600000000.0) / 24.0        
    
    gal_sgp4(&self.sgp4, e1, e2, pv, status)
    if (status[0].errc <> 0):
      rc = status[0].errc
      gal_stsfree(status)
      raise Exception, "gal_sgp4 returned code " + str(rc)            
    
    gal_me2bf(pv, gal_gmst82(e1, e2), self.omega, efg)  
    gal_pmp(efg[0], self.observer, sr)
    
    gal_pn(sr, &srmod, srunit)
    gal_rxp(self.sez, sr, rho)
    
    vsp.azimuth  = atan2(rho[1], -rho[0])
    vsp.distance = gal_pm(sr)
    sel = rho[2] / vsp.distance
    
    if (fabs(sel) > 1.0):
      sel = sel / fabs(sel)
    
    vsp.elevation = asin(sel)
    
    vsp.azimuth   = gal_anp(vsp.azimuth) * GAL_R2D
    vsp.elevation = gal_anpm(vsp.elevation) * GAL_R2D
    
    gal_stsfree(status)
    vsp.time = time
    
  """
    Calculate orbit number for the given date/time.
  """
  def calculate_orbit_number(self, norm_epoch):
    
    cdef double ecco_anomaly
    cdef double mean_anomaly_rad
    cdef double mean_anomaly_difference
    cdef double ascending_node_true_anomaly
    cdef double norm_epoch_at_ascending_node
    cdef int orbit_number
  
    if self.tle.ecco < 1:
      ascending_node_true_anomaly = 360.0 - self.tle.argpo
      ecco_anomaly = gal_ta2ea(ascending_node_true_anomaly/180.0*GAL_PI,self.tle.ecco)
      mean_anomaly_rad = ecco_anomaly - self.tle.ecco*sin(ecco_anomaly)
      if mean_anomaly_rad < 0:
        mean_anomaly_rad += 2.0 * GAL_PI
      ascending_node_mean_anomaly = mean_anomaly_rad*180.0/GAL_PI
      if ascending_node_mean_anomaly <= self.tle.mo:
        mean_anomaly_difference = self.tle.mo - ascending_node_mean_anomaly
      else:
        mean_anomaly_difference = self.tle.mo + 360.0 - ascending_node_mean_anomaly
      
      norm_epoch_at_ascending_node = self.tle_norm_epoch - (mean_anomaly_difference / 360.0 / self.tle.no)
      orbit_number = self.tle.revnum + (norm_epoch - norm_epoch_at_ascending_node) * self.tle.no
    else:
      orbit_number = -1
      
    return orbit_number


"""

  Calculate Sun position in the sky as viewed by an observer on Earth at the
  specified time and location.

"""

def sun_position(time: datetime, latitude: double, longitude: double, height: double):
  
  cdef gal_status_t* status
  cdef double epoch1, epoch2
  cdef double pv[2][3]
  cdef double itrf[2][3]
  cdef double xyz[3]
  cdef double sun[3]
  cdef double elong, phi, h
  
  cdef double a, b, c      # Ellipsoide coordinates
  cdef double t
  cdef double p[3]
  cdef double n[3]
  cdef double u[3]
  cdef double w[3]
  cdef double r
  
  cdef double tv[3]        # Normal vector to the tangent plane
  
  cdef double alpha, gamma
  cdef int visible
  cdef int rc
  
  elong = longitude * GAL_PI / 180.0
  phi   = latitude * GAL_PI / 180.0
  h     = height
  
  status = gal_stsalloc()
  
  gal_cal2jd(time.year, time.month, time.day, &epoch1, &epoch2, status)
  if status[0].errc <> 0:
    rc = status[0].errc
    gal_stsfree(status)
    raise Exception, "gal_cal2jd returned code " + str(rc)
  epoch2 = epoch2 + (time.hour + time.minute / 60.0 + time.second / 3600.0 + time.microsecond / 1000000.0) / 24.0 
   
  gal_gsupv00(epoch1, epoch2, pv, status)
  
  if status[0].errc <> 0:
    rc = status[0].errc
    gal_stsfree(status)
    raise Exception, "Date/Time in invalid range for gal_gsupv00."
    
  gal_c2tpv00a(pv, epoch1, epoch2, 0.5, 0.0, 0.0, 0.0, itrf, status)
  if status[0].errc <> 0:
    rc = status[0].errc
    gal_stsfree(status)
    raise Exception, "c2tpv00a return code " + str(rc)

# Let's get constants for the ellipsoide formula in the Geocentric reference frame

  gal_gd2gc(GAL_EMEA_WGS1984, 0.0, 0.0, h, xyz, status)
  if status[0].errc <> 0:
    rc = status[0].errc
    gal_stsfree(status)
    raise Exception, "gd2gc return code " + str(rc)  
  a = xyz[0]
  
  gal_gd2gc(GAL_EMEA_WGS1984, GAL_PI/2.0, 0.0, h, xyz, status)
  if status[0].errc <> 0:
    rc = status[0].errc
    gal_stsfree(status)
    raise Exception, "gd2gc return code " + str(rc)  
  b = xyz[1]
  
  gal_gd2gc(GAL_EMEA_WGS1984, 0.0, GAL_PI/2.0, h, xyz, status)
  if status[0].errc <> 0:
    rc = status[0].errc
    gal_stsfree(status)
    raise Exception, "gd2gc return code " + str(rc)  
  c = xyz[2]
   
# Convert the point on the Earth to the Geocentric frame reference
  
  gal_gd2gc(GAL_EMEA_WGS1984, elong, phi, h, xyz, status)
  if status[0].errc <> 0:
    rc = status[0].errc
    gal_stsfree(status)
    raise Exception, "gd2gc return code " + str(rc)
     
# Obtain Sun position vector
  
  sun[0] = itrf[0][0] - xyz[0]
  sun[1] = itrf[0][1] - xyz[1]
  sun[2] = itrf[0][2] - xyz[2]
  
# Calculate the parameter for a point on the tangent plane to the ellipsoide

  t = -(xyz[0]*sun[0]/a**2 + xyz[1]*sun[1]/b**2 + xyz[2]*sun[2]/c**2)/((xyz[0]/a**2)**2 + (xyz[1]/b**2)**2 + (xyz[2]/c**2)**2)
  
# Calculate a projection of the 'sun' vestor on the tangent plane  
  
  p[0] = sun[0] + t*xyz[0]/a**2
  p[1] = sun[1] + t*xyz[1]/b**2
  p[2] = sun[2] + t*xyz[2]/c**2
  
# Determine sun elevation  

  alpha = gal_sepp(p, sun)
  
# Calculate another parameter t to determine what side of the tangent plane the sun is

  t = 1.0/( (xyz[0] + sun[0])*xyz[0]/a**2 + (xyz[1] + sun[1])*xyz[1]/b**2 + (xyz[2] + sun[2])*xyz[2]/c**2 )
  if t >= 0:
    visible = 1
  else:
    visible = 0
    
# Calculate parameter t for a projection of the nothern pole (0, 0, c) to the tangent plane

  t =  ( xyz[0]**2/a**2 + xyz[1]**2/b**2 + xyz[2]**2/c**2 )/( (xyz[0]/a**2)**2 + (xyz[1]/b**2)**2 + (xyz[2]/c**2)**2 )
  
# Calculate a projection of the vector going from the given point on the surface to the projection of the northern pole
# on the tangent plane. This vector gives the direction to North on the plane and can be used to determine azimuth
# of the sun.
  
  n[0] = t*xyz[0]/a**2 - xyz[0]
  n[1] = t*xyz[1]/b**2 - xyz[1]
  n[2] = t*xyz[2]/c**2 - xyz[2]
  
# Determine sun azimuth 

  tv[0] = -xyz[0]/a**2
  tv[1] = -xyz[1]/b**2
  tv[2] = -xyz[2]/c**2

  gal_pn(tv, &r, u)
  
  gal_pxp(n, p, w)

  gamma = atan2(gal_pdp(u,w), gal_pdp(n,p))
  if gamma < 0:
    gamma = GAL_PI * 2.0 + gamma
    
  gal_stsfree(status)  
  
  return SunPosition(visible, alpha * 180.0 / GAL_PI, gamma * 180.0 / GAL_PI)
  
