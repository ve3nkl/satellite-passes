"""
Copyright (c) 2016-2020, 2021, 2022 Nikolai Ozerov (VE3NKL)

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

Perform calculations using maidenhead grid squares as well as Lattitude,Longitude
coordinates.

"""

from geographiclib.geodesic import Geodesic
import math

"""

   The Spot object represents a single spot on the Earth surface
   defines by its Latitude and Longitude coordinates.
   The Latitude is a float number from -90 (South pole) to 90
   (North pole). The Longitude is a float number from -180 to
   180. Negative numbers represent the Western hemisphere.

"""

class Spot:
  
  def __init__(self, latitude, longitude):
    self.latitude = latitude
    self.longitude = longitude
    
  @classmethod
  def middle(cls, spot1, spot2):
    return cls((spot1.latitude + spot2.latitude)/2.0, (spot1.longitude + spot2.longitude)/2.0)
    
  def distanceToAnotherSpot(self, anotherSpot):
    geod = Geodesic.WGS84
    g = geod.Inverse(self.latitude, self.longitude, anotherSpot.latitude, anotherSpot.longitude)
    return g["s12"]
    
  def distanceAzimuthToAnotherSpot(self, anotherSpot):
    geod = Geodesic.WGS84
    g = geod.Inverse(self.latitude, self.longitude, anotherSpot.latitude, anotherSpot.longitude)
    return (g["s12"],g["azi2"])
    
  def get_MHGrid(self):
    return MHGridSquare.lat_long_to_maiden(self.latitude, self.longitude)
    
  def __str__(self):
    return str(self.latitude) + ", " + str(self.longitude)
  
"""

   The MHGridSquare object represents a Mainheaden Grid Square. It
   is defined by a 4, 6 or 8-character string.

"""  
  
class MHGridSquare:
  
  L1 = "ABCDEFGHIJKLMNOPQR"
  L2 = "0123456789"
  L3 = "ABCDEFGHIJKLMNOPQRSTUVWX"
  
  def __init__(self, code):
    self.orig_code = code
    self.code = code.upper()
    
    B = -90                        # Bottom Latitude
    V = 180                        # Vertical size of the square
    L = -180                       # Left Longitude
    H = 360                        # Horizontal size of the square
    
    if len(self.code) == 4 or len(self.code) == 6 or len(self.code) == 8:
      
      n1 = MHGridSquare.L1.find(self.code[0:1])
      if n1 > -1:
        H = H * 1.0 / len(MHGridSquare.L1)
        L = L + n1 * H

        n2 = MHGridSquare.L1.find(self.code[1:2])
        if n2 > -1:
          V = V * 1.0 / len(MHGridSquare.L1)
          B = B + n2 * V          
          
          n1 = MHGridSquare.L2.find(self.code[2:3])
          if n1 > -1:
            H = H * 1.0 / len(MHGridSquare.L2)
            L = L + n1 * H

            n2 = MHGridSquare.L2.find(self.code[3:4])
            if n2 > -1:
              V = V * 1.0 / len(MHGridSquare.L2)
              B = B + n2 * V              
              
              if len(self.code) == 6 or len(self.code) == 8:              
                
                n1 = MHGridSquare.L3.find(self.code[4:5])
                if n1 > -1:
                  H = H * 1.0 / len(MHGridSquare.L3)
                  L = L + n1 * H

                  n2 = MHGridSquare.L3.find(self.code[5:6])
                  if n2 > -1:
                    V = V * 1.0 / len(MHGridSquare.L3)
                    B = B + n2 * V          
                    
                    if len(self.code) == 8:
                    
                      n1 = MHGridSquare.L2.find(self.code[6:7])
                      if n1 > -1:
                        H = H * 1.0 / len(MHGridSquare.L2)
                        L = L + n1 * H

                        n2 = MHGridSquare.L2.find(self.code[7:8])
                        if n2 > -1:
                          V = V * 1.0 / len(MHGridSquare.L2)
                          B = B + n2 * V             
                          
                        else:
                          raise Exception("Incorrect Maidenhead string (character 8): " + code)          
                      else:
                        raise Exception("Incorrect Maidenhead string (character 7): " + code)                        
                  else:
                    raise Exception("Incorrect Maidenhead string (character 6): " + code)          
                else:
                  raise Exception("Incorrect Maidenhead string (character 5): " + code)                
                                
            else:
              raise Exception("Incorrect Maidenhead string (character 4): " + code)          
          else:
            raise Exception("Incorrect Maidenhead string (character 3): " + code)                        
        else:
          raise Exception("Incorrect Maidenhead string (character 2): " + code)          
      else:
        raise Exception("Incorrect Maidenhead string (character 1): " + code)
    else:
      raise Exception("Incorrect Maidenhead string length: " + code)

    self.BottomLeft = Spot(B, L)
    self.BottomRight = Spot(B, L + H)
    self.TopLeft = Spot(B + V, L)
    self.TopRight = Spot(B + V, L + H)

  def getCornerCoordinates(self):
    return (                                  
      self.BottomLeft,
      self.BottomRight,
      self.TopLeft,    
      self.TopRight
    )
    
  def getCenterSpot(self):
    return Spot.middle(self.BottomLeft, self.TopRight)
    
  def distanceBetweenCenters(self, anotherSquare):
    return Spot.middle(self.BottomLeft, self.TopRight).distanceToAnotherSpot(Spot.middle(anotherSquare.BottomLeft, anotherSquare.TopRight))
    
  def distanceToSpot(self, spot):
    
    """
       Start with determining the relative positions of the two 
       squares.
    """
    
    RV = 0
    
    if self.TopLeft.latitude < spot.latitude:
      RV = -1                        # Our square is below the spot
    elif self.BottomLeft.latitude > spot.latitude:
      RV = 1                         # Our square is above the spot
    else:
      RV = 0                         # The Parallel line coming through
                                     # the spot crosses the square
    RH = 0
    if self.BottomLeft.longitude <= spot.longitude and spot.longitude <= self.BottomRight.longitude:
      RH = 0                         # The Meridian line coming through
                                     # the spot crosses the square
    else:
      w = self.TopRight.longitude - self.TopLeft.longitude
      
      if spot.longitude > self.TopRight.longitude:
        d = spot.longitude - self.TopRight.longitude
        if d <= 360 - d - w:
          RH = 1                     # The spot is to the right of the 
        else:                        # square
          RH = -1                    # The spot is to the left of the
      else:                          # square
        d = self.TopLeft.longitude - spot.longitude
        if d <= 360 - d - w:
          RH = -1                    # The spot is to the left of the
        else:                        # square
          RH = 1                     # The spot is to the right of the
                                     # square
    if RV == 0 and RH == 0:
      return (0, 0)                  # The spot is inside the square
    elif RV == 0:
      if RH == 1:
        dmax1 = spot.distanceToAnotherSpot(self.TopLeft)
        dmax2 = spot.distanceToAnotherSpot(self.BottomLeft)
        dmin = spot.distanceToAnotherSpot(Spot(spot.latitude, self.TopRigh.longitude))
      else:
        dmax1 = spot.distanceToAnotherSpot(self.TopRight)
        dmax2 = spot.distanceToAnotherSpot(self.BottomRight)
        dmin = spot.distanceToAnotherSpot(Spot(spot.latitude, self.TopLeft.longitude))
        
      if dmax1 >= dmax2:
        return (dmin, dmax1)
      else:
        return (dmin, dmax2)        
    elif RH == 0:
      if RV == 1:
        dmax1 = spot.distanceToAnotherSpot(self.BottomLeft)
        dmax2 = spot.distanceToAnotherSpot(self.BottomRight)
        dmin = spot.distanceToAnotherSpot(self.TopLeft.latitude, spot.longitude)
      else:
        dmax1 = spot.distanceToAnotherSpot(self.TopLeft)
        dmax2 = spot.distanceToAnotherSpot(self.TopRight)
        dmin = spot.distanceToAnotherSpot(self.BottomLeft.latitude, spot.longitude)        
        
      if dmax1 >= dmax2:
        return (dmin, dmax1)
      else:
        return (dmin, dmax2)            
        
    elif RV == 1 and RH == 1:
      return (spot.distanceToAnotherSpot(self.TopRight), spot.distanceToAnotherSpot(self.BottomLeft))
    elif RV == 1 and RH == -1:
      return (spot.distanceToAnotherSpot(self.TopLeft), spot.distanceToAnotherSpot(self.BottomRight))
    elif RV == -1 and RH == 1:
      return (spot.distanceToAnotherSpot(self.BottomRight), spot.distanceToAnotherSpot(self.TopLeft))
    else:
      return (spot.distanceToAnotherSpot(self.BottomLeft), spot.distanceToAnotherSpot(self.TopRight))
      
  def distanceToAnotherSquare(self, anotherSquare):
    dmin1, dmax1 = self.distanceToSpot(anotherSquare.BottomLeft)  
    dmin2, dmax2 = self.distanceToSpot(anotherSquare.BottomRight)  
    dmin3, dmax3 = self.distanceToSpot(anotherSquare.TopLeft)  
    dmin4, dmax4 = self.distanceToSpot(anotherSquare.TopRight)  
    
    dmin = dmin1
    if dmin > dmin2:
      dmin = dmin2
    if dmin > dmin3:
      dmin = dmin3
    if dmin > dmin4:
      dmin = dmin4            
      
    dmax = dmax1
    if dmax < dmax2:
      dmax = dmax2
    if dmax < dmax3:
      dmax = dmax3
    if dmax < dmax4:
      dmax = dmax4                  
      
    return (dmin, dmax)
      
  @classmethod
  def lat_long_to_maiden(cls, latitude, longitude):
    lo1 = int((longitude + 180.0) / 20.)
    if lo1 < 0:
      lo1 = 0
    elif lo1 > 17:
      lo1 = 17
    
    lo1rem = (longitude + 180.0) - lo1 * 20.0
    
    lo2 = int(lo1rem / 2.0)
    if lo2 < 0:
      lo2 = 0
    elif lo2 > 9:
      lo2 = 9
      
    lo2rem = lo1rem - lo2 * 2.0
    
    lo3 = int(lo2rem * 12.0)
    if lo3 < 0:
      lo3 = 0
    elif lo3 > 23:
      lo3 = 23
    
    lo3rem = lo2rem * 12.0 - lo3
    
    lo4 = int(lo3rem * 10)
    if lo4 < 0:
      lo4 = 0
    elif lo4 > 9:
      lo4 = 9
    
    la1 = int((latitude + 90.0) / 10.)
    if la1 < 0:
      la1 = 0
    elif la1 > 17:
      la1 = 17
    
    la1rem = (latitude + 90.0) - la1 * 10.0
    
    la2 = int(la1rem)
    if la2 < 0:
      la2 = 0
    elif la2 > 9:
      la2 = 9
      
    la2rem = la1rem - la2
    
    la3 = int(la2rem * 24.0)
    if la3 < 0:
      la3 = 0
    elif la3 > 23:
      la3 = 23
    
    la3rem = la2rem * 24.0 - la3
    la4 = int(la3rem * 10)
    if la4 < 0:
      la4 = 0
    elif la4 > 9:
      la4 = 9
      
    return MHGridSquare.L1[lo1] + MHGridSquare.L1[la1] + str(lo2) + str(la2) + MHGridSquare.L3.lower()[lo3] + MHGridSquare.L3.lower()[la3] + str(lo4) + str(la4)

  """
    Determine a direction of an arrow than is drawn from grid square 1 to grid square 2. Both grid squares are represented
    by strings (of the same length).
  """
  @classmethod
  def arrow_direction(cls, grid_str1, grid_str2):
    if len(grid_str1) != len(grid_str2):
      raise Exception("Grid squares are represented by strings of different length.")
    pos = 0
    same_horizontal_position = True
    same_vertical_position   = True
    rH = ""
    rV = ""
    for c1, c2 in zip(grid_str1.upper(), grid_str2.upper()):
      pos += 1
      if c1 == c2:
        continue
      if pos % 2 == 1: # Horizontal direction
        if pos == 1:   #   First character for this direction
          if c1 > c2 and ord(c1) - ord(c2) < 9 or c1 < c2 and ord(c1)-ord('A')+ord('R')-ord(c2)+1 < 9:
            rH = "W"
            same_horizontal_position = False
          elif c1 < c2 and ord(c2) - ord(c1) < 9 or c2 < c1 and ord(c2)-ord('A')+ord('R')-ord(c1)+1 < 9:
            rH = "E"
            same_horizontal_position = False
          else:
            raise Exception("Horizontal distance is too high.")
        else:          #   Second or higher character for this direction
          if same_horizontal_position:
            if c1 > c2:
              rH = "W"
              same_horizontal_position = False
            else:
              rH = "E"
              same_horizontal_position = False
          
      else:            # Vertical direction
        if pos == 2:   #   First character for this direction
          if c1 > c2 and ord(c1) - ord(c2) < 9 or c1 < c2 and ord(c1)-ord('A')+ord('R')-ord(c2)+1 < 9:
            rV = "S"
            same_vertical_position = False
          elif c1 < c2 and ord(c2) - ord(c1) < 9 or c2 < c1 and ord(c2)-ord('A')+ord('R')-ord(c1)+1 < 9:
            rV = "N"
            same_vertical_position = False
          else:
            raise Exception("Vertical distance is too high.")
        else:          #   Second character for this direction
          if same_vertical_position:
            if c1 > c2:
              rV = "S"
              same_vertical_position = False
            else:
              rV = "N"
              same_vertical_position = False
              
      if not same_horizontal_position and not same_vertical_position:
        break    

    return (rH, rV)
