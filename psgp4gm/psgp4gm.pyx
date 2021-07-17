#cython: language_level=3

"""
Copyright (c) 2014-2020, 2021 Nikolai Ozerov (VE3NKL)

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

Wrapping up gravity model structure initialization.

"""

cdef extern from "string.h":

  void strcpy(char* str1, char* str2)
  
cdef extern from "math.h":

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
  
  ctypedef struct gal_tle_t  

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
  
cdef extern from "gal_cal2jd.h":

  void gal_cal2jd(int iy, int im, int id, double *djm0, double *djm, gal_status_t* status)
  
cdef extern from "gal_jd2cal.h":

  void gal_jd2cal(double dj1, double dj2, int *iy, int *im, int *id, double *fd, gal_status_t* status)
  
cdef extern from "gal_frame_macros.h":

  ctypedef enum:
    GAL_EMEA_WGS1984
    pass
    
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
  
# Define Global variables  

cdef gal_sgp4_t sgp4orig

cdef gal_sgp4_t sgp4

cdef double omega
cdef double mu, sma, inf, re, gm

cdef double observer[3]

cdef double sez[3][3]

cdef double rlat, rlon, h
    
def init_wgs72(card1, card2, latitude, longitude, height):

  global sgp4, sgp4orig
  global omega
  global mu, sma, inf, re, gm
  global observer
  global sez
  global rlat, rlon, h

  cdef int rc 
  cdef gal_status_t* status
  cdef gal_tle_t tle
  cdef gal_gm_t* gm72
  
  cdef double tumin
  cdef double xke
  cdef double j2
  cdef double j3
  cdef double j4
  cdef double j3oj2  
  
  omega = 7.292115146706979e-05
  rlat = latitude * GAL_D2R
  rlon = longitude * GAL_D2R
  h    = height
  
  status = gal_stsalloc()
    
  gm72 = gal_gmget_wgs72(180, 180, GAL_NORMALIZED, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_gmget_wgs72 returned code " + str(status[0].errc)
  
  gm = gm72.gm
  
  gal_emparams(GAL_EMEA_WGS1984, &sma, &inf, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_emparams returned code " + str(status[0].errc)
  
  gal_sgp4gm(gm72, &tumin, &mu, &re, &xke, &j2, &j3, &j4, &j3oj2, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_sgp4gm returned code " + str(status[0].errc)  
  
  gal_gd2gce(sma, inf, rlon, rlat, h, observer, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_gd2gce returned code " + str(status[0].errc)    
  
  gal_sezm(rlat, rlon, sez)
  
  gal_tledec(card1, card2, &tle, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_tledec returned code " + str(status[0].errc)      

  gal_sgp4init(gm72, &tle, &sgp4orig, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_sgp4init returned code " + str(status[0].errc)        
  
  gal_gmfree(gm72)
  sgp4 = sgp4orig
  
  gal_stsfree(status)

  return
  
def get_position(y, m, d, hh, mm, ss):

  global sgp4
  global omega
  global mu, sma, inf, re, gm
  global observer
  global sez
  global rlat, rlon, h  

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
  
  cdef int year, month, day
  cdef int hr, min, sec
  
  status = gal_stsalloc()  
    
  year = y
  month = m
  day = d
  hr = hh
  min = mm
  sec = ss  
  
  gal_cal2jd(year, month, day, &e1, &e2, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_cal2jd returned code " + str(status[0].errc)      
    
  e2 = e2 + hr / 24.0 + min / 24.0 / 60.0 + sec / 24.0 /60.0 / 60.0        
  
  gal_sgp4(&sgp4, e1, e2, pv, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_sgp4 returned code " + str(status[0].errc)            
  
  gal_me2bf(pv, gal_gmst82(e1, e2), omega, efg)  
  gal_pmp(efg[0], observer, sr)
  
  gal_pn(sr, &srmod, srunit)
  gal_rxp(sez, sr, rho)
  
  az  = atan2(rho[1], -rho[0])
  dist = gal_pm(sr)
  sel = rho[2] / dist
  
  if (fabs(sel) > 1.0):
    sel = sel / fabs(sel)
  
  elev = asin(sel)
  
  az = gal_anp(az)
  elev = gal_anpm(elev)
  
  gal_stsfree(status)
  return (e1, e2, az * GAL_R2D, elev * GAL_R2D, sr[0], sr[1], sr[2], dist)


def find_rise_set(y, m, d, hh, mm, ss, minel, daysfwd):

  global sgp4orig
  global omega
  global mu, sma, inf, re, gm
  global observer
  global sez
  global rlat, rlon, h

  cdef gal_sgp4_t sgp4l

  cdef gal_status_t* status
  
  cdef double rminel
  cdef double days
  cdef int year, month, day
  cdef int hr, min, sec
  
  cdef int iy1, im1, id1, iy2, im2, id2, iy3, im3, id3 
  cdef double time1, time2, time3
  cdef int ihh1, imm1, iss1, ihh2, imm2, iss2, ihh3, imm3, iss3   
  
  cdef double passinfo[3][3]  
  cdef double epoch1, epoch2
  
  year = y
  month = m
  day = d
  hr = hh
  min = mm
  sec = ss
  
  status = gal_stsalloc()  

  gal_cal2jd(year, month, day, &epoch1, &epoch2, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_cal2jd returned code " + str(status[0].errc)              
    
  epoch2 = epoch2 + hr / 24.0 + min / 24.0 / 60.0 + sec / 24.0 /60.0 / 60.0
  
  sgp4l = sgp4orig
  
  rminel     = minel * GAL_D2R
  days       = daysfwd

  gal_sgp4rs(epoch1, epoch2, rlat, rlon, h, rminel, days, gm, sma, 1.0 / inf, &sgp4l, passinfo, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_sgp4rs returned code " + str(status[0].errc)                
  
  gal_jd2cal(passinfo[0][0], 0.0, &iy1, &im1, &id1, &time1, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_jd2cal returned code " + str(status[0].errc)                  
  ihh1 = (int)(24.0 * time1)
  imm1 = (int)((time1 * 24.0 - ihh1) * 60.0)
  iss1 = (int)((time1 * 24.0 * 60.0  - ihh1 * 60.0 - imm1) * 60.0)
  
  gal_jd2cal(passinfo[1][0], 0.0, &iy2, &im2, &id2, &time2, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_jd2cal returned code " + str(status[0].errc)                    
  ihh2 = (int)(24.0 * time2)
  imm2 = (int)((time2 * 24.0 - ihh2) * 60.0)
  iss2 = (int)((time2 * 24.0 * 60.0  - ihh2 * 60.0 - imm2) * 60.0)  
  
  gal_jd2cal(passinfo[2][0], 0.0, &iy3, &im3, &id3, &time3, status)
  if (status[0].errc <> 0):
    raise Exception, "gal_jd2cal returned code " + str(status[0].errc)                    
  ihh3 = (int)(24.0 * time3)
  imm3 = (int)((time3 * 24.0 - ihh3) * 60.0)
  iss3 = (int)((time3 * 24.0 * 60.0  - ihh3 * 60.0 - imm3) * 60.0)   
  
  gal_stsfree(status)   

  return ((iy1, im1, id1, ihh1, imm1, iss1, passinfo[0][0]), passinfo[0][1] * GAL_R2D, passinfo[0][2] * GAL_R2D, (iy2, im2, id2, ihh2, imm2, iss2,  passinfo[1][0]), passinfo[1][1] * GAL_R2D, passinfo[1][2] * GAL_R2D, (iy3, im3, id3, ihh3, imm3, iss3,  passinfo[2][0]), passinfo[2][1] * GAL_R2D, passinfo[2][2] * GAL_R2D)


def sun_position(y, m, d, hh, mm, ss, lng, lat, height):
  
  cdef int year, month, day
  cdef int hr, min, sec
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
  
  year = y
  month = m
  day = d
  hr = hh
  min = mm
  sec = ss
  
  elong = lng * GAL_PI / 180.0
  phi   = lat * GAL_PI / 180.0
  h     = height
  
  status = gal_stsalloc()
  
  gal_cal2jd(year, month, day, &epoch1, &epoch2, status)
  if status[0].errc <> 0:
    raise Exception, "gal_cal2jd returned code " + str(status[0].errc)
  epoch2 = epoch2 + hr / 24.0 + min / 24.0 / 60.0 + sec / 24.0 /60.0 / 60.0
   
  gal_gsupv00(epoch1, epoch2, pv, status)
  
  if status[0].errc <> 0:
    raise Exception, "Date/Time in invalid range for gal_gsupv00."
    
  gal_c2tpv00a(pv, epoch1, epoch2, 0.5, 0.0, 0.0, 0.0, itrf, status)
  
  if status[0].errc <> 0:
    raise Exception, "c2tpv00a return code " + str(status[0].errc)

# Let's get constants for the ellipsoide formula in the Geocentric reference frame

  gal_gd2gc(GAL_EMEA_WGS1984, 0.0, 0.0, h, xyz, status)
  if status[0].errc <> 0:
    raise Exception, "gd2gc return code " + str(status[0].errc)  
  a = xyz[0]
  
  gal_gd2gc(GAL_EMEA_WGS1984, GAL_PI/2.0, 0.0, h, xyz, status)
  if status[0].errc <> 0:
    raise Exception, "gd2gc return code " + str(status[0].errc)  
  b = xyz[1]
  
  gal_gd2gc(GAL_EMEA_WGS1984, 0.0, GAL_PI/2.0, h, xyz, status)
  if status[0].errc <> 0:
    raise Exception, "gd2gc return code " + str(status[0].errc)  
  c = xyz[2]
   
# Convert the point on the Earth to the Geocentric frame reference
  
  gal_gd2gc(GAL_EMEA_WGS1984, elong, phi, h, xyz, status)

  if status[0].errc <> 0:
    raise Exception, "gd2gc return code " + str(status[0].errc)
     
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
  
# Calculate a projection of the vector going from the given point om the surface to the projection of the northern pole
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
  
  return (visible, alpha * 180.0 / GAL_PI, gamma * 180.0 / GAL_PI)
  
