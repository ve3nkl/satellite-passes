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

- History -------------------------------------------------------------------
2021/07/28   1.2   Show orbit # if available
-----------------------------------------------------------------------------

This script reads satellite pass info from stdin and generates PDF document
with visual representation of the passes.


"""

import os, sys
import argparse
import math
import datetime
import decimal
import yaml
from pyx import *

def draw_handle(scl, cx, cy, w, h):

  xbl = cx - w / 2.
  ybl = cy
  xbr = cx + w / 2.
  ybr = cy
  xtl = cx - w / 2.
  ytl = cy + h - w / 2.
  xtr = cx + w / 2.
  ytr = cy + h  - w / 2.
  xtcr = cx + w / 2.
  ytcr = cy + h + w / 2.
  xtcl = cx - w / 2.
  ytcl = cy + h + w / 2.  
  
  handle = path.path(path.moveto(xbl / scl, ybl / scl), path.lineto(xbr / scl, ybr / scl),
    path.lineto(xtr / scl, ytr / scl), path.curveto(xtcr / scl, ytcr / scl, xtcl / scl, ytcl / scl, xtl / scl, ytl / scl), path.closepath())

  return handle
  
def genPass(cnvs, scl, cx, cy, title, tz_offset_str, satellite, rise_dt, set_dt, rise_azimuth, set_azimuth, max_elev, max_azimuth, max_dt, orbit_n):
  
  title_bump = 2.0

  main_circle_r = 85

  sat_circle_r = 30.
  sat_circle_x = cx
  sat_circle_y = cy

  pfx = "-" + str(cx) + "-" + str(cy)

  total_size = cy + 100

  watch_circle_r = 70.
  watch_scale1_len = 0.1
  watch_scale2_len = 0.05

  hour_handle_width  = 5.0
  hour_handle_height = (watch_circle_r - sat_circle_r) * 0.4 + sat_circle_r

  minute_handle_width  = 3.0
  minute_handle_height = (watch_circle_r - sat_circle_r) * 0.8 + sat_circle_r

  pilot_circle_20 = (70.0 / 90.0) * sat_circle_r
  pilot_circle_40 = (50.0 / 90.0) * sat_circle_r
  pilot_circle_60 = (30.0 / 90.0) * sat_circle_r
  pilot_circle_80 = (10.0 / 90.0) * sat_circle_r

  max_elev_str = str(decimal.Decimal(str(max_elev)).quantize(decimal.Decimal('.1')))
  max_elev_str = max_elev_str.split(".")[0]

  date_circle_r = 10
  date_circle_x = sat_circle_x + watch_circle_r - date_circle_r - 2
  date_circle_y = sat_circle_y
  date_circle_text = str(rise_dt.day)
  if len(date_circle_text) == 1:
    date_circle_text = "0" + date_circle_text

  title_lower = rise_dt.strftime("%b %d, %Y") + "   TZ Offset " + tz_offset_str

  if rise_dt.hour >= 12:
    rise_hour_angle = (rise_dt.hour - 12.0 + rise_dt.minute / 60. + rise_dt.second / 3600.) / 12. * 360.
    rise_ampm = "PM"
  else:
    rise_hour_angle = (rise_dt.hour - 12.0 + rise_dt.minute / 60. + rise_dt.second / 3600.) / 12. * 360.
    rise_ampm = "AM"

  rise_minute_angle = (rise_dt.minute + rise_dt.second / 60.) / 60. * 360.
  set_minute_angle  = (set_dt.minute + set_dt.second / 60.) / 60. * 360.

  rise_x = sat_circle_r * math.sin(rise_azimuth / 360. * 2.0 * math.pi) + sat_circle_x
  rise_y = sat_circle_r * math.cos(rise_azimuth / 360. * 2.0 * math.pi) + sat_circle_y

  set_x = sat_circle_r * math.sin(set_azimuth / 360. * 2.0 * math.pi) + sat_circle_x
  set_y = sat_circle_r * math.cos(set_azimuth / 360. * 2.0 * math.pi) + sat_circle_y

  middle_x = (set_x + rise_x) / 2.0
  middle_y = (set_y + rise_y) / 2.0

  # Calculate middle angle ...

  a = rise_azimuth
  b = set_azimuth
  if a > b:  # Make sure that a <= b
    c = a
    a = b
    b = c

  middle_angle = (a + b) / 2.0

  if b - a > 180.0: # If the difference is greater than half-circle then we are on the wrong side
    middle_angle = middle_angle + 180.0

  if middle_angle > 360.0:
    middle_angle = middle_angle - 360.0

  middle_arc_x = sat_circle_r * math.sin(middle_angle / 360. * 2.0 * math.pi) + sat_circle_x
  middle_arc_y = sat_circle_r * math.cos(middle_angle / 360. * 2.0 * math.pi) + sat_circle_y

  t_x = (sat_circle_x - middle_arc_x) * (max_elev / 90.0) + middle_arc_x
  t_y = (sat_circle_y - middle_arc_y) * (max_elev / 90.0) + middle_arc_y

  q_x = (t_x - middle_x) * 2.0 + middle_x
  q_y = (t_y - middle_y) * 2.0 + middle_y

# Watch dial circle.

  circle = path.circle(sat_circle_x / scl, sat_circle_y / scl , watch_circle_r / scl)
  cnvs.stroke(circle)
  
# Let's mark hours and minutes.

  ref_points = []

  for i in range(0, 60):

    px1 = math.sin(2. * math.pi / 60. * i) * watch_circle_r + sat_circle_x
    py1 = math.cos(2. * math.pi / 60. * i) * watch_circle_r + sat_circle_y

    if i % 5 == 0:
      px2 = (sat_circle_x - px1) * watch_scale1_len + px1
      py2 = (sat_circle_y - py1) * watch_scale1_len + py1
      w = 1

      ref_points.append((px2, py2))
      
    else:
      px2 = (sat_circle_x - px1) * watch_scale2_len + px1
      py2 = (sat_circle_y - py1) * watch_scale2_len + py1
      w = 0.5

    line = path.line(px1 / scl, py1 / scl, px2 / scl, py2 / scl)
    if w == 1:
      cnvs.stroke(line)
    else:
      cnvs.stroke(line)      
    
# Set numbers for hours. Satellite name instead of 1 and 11, AM or PM instead of 6 and W(est) and E(ast) instead of 9 and 3 correspondingly.

  cnvs.text(ref_points[0][0] / scl, (ref_points[0][1]-24) / scl, satellite, [text.size.Large, text.valign.middle, text.halign.center])
  
  cnvs.text((ref_points[0][0]) / scl, (ref_points[0][1]-5) / scl, "12", [text.size.normalsize, text.valign.middle, text.halign.center])  
  cnvs.text((ref_points[1][0]-5) / scl, (ref_points[1][1]-4) / scl, "1", [text.size.normalsize, text.valign.middle, text.halign.center])
  cnvs.text((ref_points[2][0]-5) / scl, (ref_points[2][1]-3) / scl, "2", [text.size.normalsize, text.valign.middle, text.halign.center])

  cnvs.text((ref_points[3][0]-25) / scl, (ref_points[3][1]) / scl, "E", [text.size.large, text.valign.middle, text.halign.center])

  cnvs.text((ref_points[4][0]-5) / scl, (ref_points[4][1]+3) / scl, "4", [text.size.normalsize, text.valign.middle, text.halign.center])
  cnvs.text((ref_points[5][0]-5) / scl, (ref_points[5][1]+4) / scl, "5", [text.size.normalsize, text.valign.middle, text.halign.center])

  cnvs.text(ref_points[6][0] / scl, (ref_points[6][1]+8) / scl, str(rise_ampm), [text.size.large, text.valign.middle, text.halign.center])

  cnvs.text((ref_points[7][0]+5) / scl, (ref_points[7][1]+3) / scl, "7", [text.size.normalsize, text.valign.middle, text.halign.center])
  cnvs.text((ref_points[8][0]+5) / scl, (ref_points[8][1]+4) / scl, "8", [text.size.normalsize, text.valign.middle, text.halign.center])

  cnvs.text((ref_points[9][0]+22) / scl, (ref_points[9][1]) / scl, "W", [text.size.large, text.valign.middle, text.halign.center])

  cnvs.text((ref_points[10][0]+5) / scl, (ref_points[10][1]-4) / scl, "10", [text.size.normalsize, text.valign.middle, text.halign.center])
  cnvs.text((ref_points[11][0]+5) / scl, (ref_points[11][1]-3) / scl, "11", [text.size.normalsize, text.valign.middle, text.halign.center])

# Calendar window with the date.

  circle = path.circle(date_circle_x / scl, date_circle_y / scl , date_circle_r / scl)
  cnvs.stroke(circle, [deco.filled([color.rgb(1.0, 1.0, 1.0)])])
  cnvs.text((date_circle_x) / scl, (date_circle_y) / scl, date_circle_text, [text.size.normalsize, text.valign.middle, text.halign.center])  

# Hour Handle.

  cnvs.stroke(draw_handle(scl, sat_circle_x, sat_circle_y, hour_handle_width, hour_handle_height),
    [deco.filled([color.grey(0.5), color.transparency(0.5)]), trafo.rotate(- rise_hour_angle, sat_circle_x / scl, sat_circle_y / scl)])

# Minute Handle.

  cnvs.stroke(draw_handle(scl, sat_circle_x, sat_circle_y, minute_handle_width, minute_handle_height),
    [deco.filled([color.grey(0.5), color.transparency(0.5)]), trafo.rotate(- rise_minute_angle, sat_circle_x / scl, sat_circle_y / scl)])

# The second minute handle to mark the end of the pass time.

  cnvs.stroke(draw_handle(scl, sat_circle_x, sat_circle_y, minute_handle_width, minute_handle_height),    
    [deco.filled([color.rgb(1.0, 1.0, 1.0), color.transparency(0.5)]), trafo.rotate(- set_minute_angle, sat_circle_x / scl, sat_circle_y / scl)])

# Draw the main outter circle.

  main_circle = path.circle(sat_circle_x / scl, sat_circle_y / scl, main_circle_r / scl)
  cnvs.stroke(main_circle)

# Insert text in the upper part of the ring.

  p1 = path.path(path.arcn(sat_circle_x / scl, sat_circle_y / scl, (watch_circle_r + (main_circle_r - watch_circle_r) / 2.) / scl, 170, 10))
  cnvs.draw(p1, [ deco.curvedtext(title, textattrs=[text.vshift.mathaxis, text.size.footnotesize]) ])
  
# Insert text in the lower part of the ring.

  p2 = path.path(path.arc(sat_circle_x / scl, sat_circle_y / scl, (watch_circle_r + (main_circle_r - watch_circle_r) / 2.) / scl, 190, 350))
  cnvs.draw(p2, [ deco.curvedtext(title_lower, textattrs=[text.vshift.mathaxis, text.size.footnotesize]) ])

# Add navigational circles inside the satellite circle.

  nav_circle = path.circle(sat_circle_x / scl, sat_circle_y / scl, sat_circle_r / scl)
  cnvs.stroke(nav_circle, [deco.filled([color.rgb(1.0, 1.0, 1.0)])])

  pi20_circle = path.circle(sat_circle_x / scl, sat_circle_y / scl, pilot_circle_20 / scl)
  cnvs.stroke(pi20_circle, [style.linewidth(0.01), style.linestyle.dashed])

  pi40_circle = path.circle(sat_circle_x / scl, sat_circle_y / scl, pilot_circle_40 / scl)
  cnvs.stroke(pi40_circle, [style.linewidth(0.01), style.linestyle.dashed])

  pi60_circle = path.circle(sat_circle_x / scl, sat_circle_y / scl, pilot_circle_60 / scl)
  cnvs.stroke(pi60_circle, [style.linewidth(0.01), style.linestyle.dashed])      

  cnvs.text((sat_circle_x) / scl, (sat_circle_y) / scl, max_elev_str, [text.size.normalsize, text.valign.middle, text.halign.center, color.rgb(0., 0., 0.)])

# Draw the satellite pass.

  vhx = (set_x - rise_x) / 4.
  vhy = (set_y - rise_y) / 4.
  vvx = (t_x - middle_x) / 3.
  vvy = (t_y - middle_y) / 3.

  p_x = rise_x + vhx + vvx * 4.
  p_y = rise_y + vhy + vvy * 4.

  r_x = rise_x + vhx * 3 + vvx * 4.
  r_y = rise_y + vhy * 3 + vvy * 4.   

#  p_x = rise_x * 0.25 + q_x * 0.75
#  p_y = rise_y * 0.25 + q_y * 0.75
#  r_x = set_x * 0.25 + q_x * 0.75
#  r_y = set_y * 0.25 + q_y * 0.75  

  traj = path.path(path.moveto(rise_x / scl, rise_y / scl), path.curveto(p_x / scl, p_y / scl, r_x / scl, r_y / scl, set_x / scl, set_y / scl))
  cnvs.stroke(traj, [deco.earrow([deco.stroked([style.linejoin.round])], size=0.2)])

# Insert text reflecting exact AOS info

  aos_txt_x = sat_circle_x / scl - main_circle_r / scl
  aos_txt_y = sat_circle_y / scl - main_circle_r / scl  
  cnvs.text(aos_txt_x, aos_txt_y, rise_dt.isoformat()[11:19], [text.size.footnotesize, text.valign.bottom, text.halign.left])
  cnvs.text(aos_txt_x, aos_txt_y + 10 / scl, str(rise_azimuth), [text.size.footnotesize, text.valign.bottom, text.halign.left])
  cnvs.text(aos_txt_x, aos_txt_y + 20 / scl, "A", [text.size.footnotesize, text.valign.bottom, text.halign.left])  

# Insert text reflecting exact LOS info

  los_txt_x = sat_circle_x / scl + main_circle_r / scl
  los_txt_y = sat_circle_y / scl - main_circle_r / scl  
  cnvs.text(los_txt_x, los_txt_y, set_dt.isoformat()[11:19], [text.size.footnotesize, text.valign.bottom, text.halign.right])
  cnvs.text(los_txt_x, los_txt_y + 10 / scl, str(set_azimuth), [text.size.footnotesize, text.valign.bottom, text.halign.right])
  cnvs.text(los_txt_x, los_txt_y + 20 / scl, "L", [text.size.footnotesize, text.valign.bottom, text.halign.right])  

# Insert text reflecting exact MAX info

  max_txt_x = sat_circle_x / scl - main_circle_r / scl
  max_txt_y = sat_circle_y / scl + main_circle_r / scl
  max_elev_str2 = str(decimal.Decimal(str(max_elev)).quantize(decimal.Decimal('.01')))[0:-1]
  cnvs.text(max_txt_x, max_txt_y, max_dt.isoformat()[11:19], [text.size.footnotesize, text.valign.top, text.halign.left])
  cnvs.text(max_txt_x, max_txt_y - 10 / scl, str(max_azimuth), [text.size.footnotesize, text.valign.top, text.halign.left])  
  cnvs.text(max_txt_x, max_txt_y - 20 / scl, max_elev_str2, [text.size.footnotesize, text.valign.top, text.halign.left])

# Insert orbit #

  if orbit_n != "":
    orb_txt_x = sat_circle_x / scl + main_circle_r / scl
    orb_txt_y = sat_circle_y / scl + main_circle_r / scl
    cnvs.text(orb_txt_x, orb_txt_y, "Orbit", [text.size.footnotesize, text.valign.top, text.halign.right])
    cnvs.text(orb_txt_x, orb_txt_y - 10 / scl, orbit_n, [text.size.footnotesize, text.valign.top, text.halign.right])

  return 

if __name__ == "__main__":
  
  parser = argparse.ArgumentParser(description='Represent satellite pass info visually and save it as a PDF file.')
  parser.add_argument("-v", "--version", action="store_true", default=False, help="display current version")
  parser.add_argument("-o", "--output",  action="store", type=str, help="specify name for the output pdf file")
  args = parser.parse_args()
  
  if args.version:
    print("passpdf.py verion 1.2")
    print("Copyright (c) 2014-2020, 2021 Nikolai Ozerov (VE3NKL)")
    quit(0)

  try:
    data = sys.stdin.read()
  except:
    print( "Error reading stdin." )
    print( )
    quit(1)               

  y = yaml.safe_load(data)
  z = []
  
  rise_dt_first  = None
  qth_name_first = None

  for p in y:

    rise_dt = None
    try:
      rise_dt = datetime.datetime.strptime(str(p["sat-pass"]["rise-date"]) + " " + str(p["sat-pass"]["rise-info"]["time"]), "%b %d, %Y %H:%M:%S")
      rise_dt.replace(tzinfo=None)
    except:
      sys.stderr.write('Error retrieving rise date/time.: [' + str(p["sat-pass"]["rise-date"]) + '] / [' + str(p["sat-pass"]["rise-info"]["time"]) + ']\n')
      quit(1)

    set_dt = None
    try:
      set_dt = datetime.datetime.strptime(str(p["sat-pass"]["rise-date"]) + " " + str(p["sat-pass"]["set-info"]["time"]), "%b %d, %Y %H:%M:%S")
      set_dt.replace(tzinfo=None)
    except:
      sys.stderr.write('Error retrieving set time: [' + str(p["sat-pass"]["set-info"]["time"]) + ']\n')
      quit(1)

    max_dt = None
    try:
      max_dt = datetime.datetime.strptime(str(p["sat-pass"]["rise-date"]) + " " + str(p["sat-pass"]["max-info"]["time"]), "%b %d, %Y %H:%M:%S")
      max_dt.replace(tzinfo=None)
    except:
      sys.stderr.write('Error retrieving max time: [' + str(p["sat-pass"]["max-info"]["time"]) + ']\n')
      quit(1)

    # We used the same date for the set time. If the pass crossed midnight set_dt will be < rise_dt. To correct
    # this we add 24 hours to set_dt.

    if set_dt < rise_dt:
      set_dt = set_dt + datetime.timedelta(hours=24)

    if max_dt < rise_dt:
      max_dt = max_dt + datetime.timedelta(hours=24)

    z.append((rise_dt, set_dt, max_dt, p))
    if rise_dt_first is None:
      rise_dt_first  = rise_dt
      qth_name_first = p["sat-pass"]["QTH"]

  z.sort(key=lambda tuple: tuple[0])

  i = 0
  j = 0
  row_len = 4
  rows_per_page = 5
  outline_side = 5
  outline_side_size = 180
  
  scale = 40.

  rows = (len(y) + row_len - 1) // row_len

  d = document.document()

  c = canvas.canvas()
  
  print( "Number of passes: " + str(len(z)) )

  for (rise_dt, set_dt, max_dt, p) in z:

    cx = j * outline_side_size + (outline_side_size // 2)
    cy = outline_side_size * row_len - i * outline_side_size + (outline_side_size // 2)
       
    genPass(c, scale, 
            cx, cy,
            p["sat-pass"]["QTH"],
            p["sat-pass"]["time-zone"],
            p["sat-pass"]["satellite-user-name"].upper(),
            rise_dt,
            set_dt,
            float(p["sat-pass"]["rise-info"]["azimuth"]),
            float(p["sat-pass"]["set-info"]["azimuth"]),
            float(p["sat-pass"]["max-info"]["elevation"]),
            float(p["sat-pass"]["max-info"]["azimuth"]),
            max_dt,
            p["sat-pass"]["rise-orbit"])

    j = j + 1
    if j >= row_len:
      i = i + 1
      j = 0
      if i >= rows_per_page:
        d.append(document.page(c, paperformat=document.paperformat.Letter, centered=0, fittosize=1, margin=2.5))
        i = 0
        j = 0
        c = canvas.canvas()

  if i > 0 or j > 0:
    x = outline_side_size*row_len-1                    # Make sure that the entire canvas
    line = path.line(x / scale, 0, (x - 1) / scale, 0) # is used for drawing by inserting
    c.stroke(line)                                     # a tiny speck in the bottom right
                                                       # corner
    d.append(document.page(c, paperformat=document.paperformat.Letter, centered=0, fittosize=1, margin=2.5))
  
  # We construct the output file name automatically based on the AOS time of the earliest pass and a name
  # of the location we calculated the pass parameters for.
  
  if args.output:
    output_pdf_file_name = args.output
  else:
    output_pdf_file_name = rise_dt_first.strftime("%Y%m%d-%H%M%S-") + qth_name_first.replace(" ", "-") + ".pdf"
    
  d.writePDFfile(output_pdf_file_name)
  
  
