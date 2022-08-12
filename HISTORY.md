### 2022/08/12 - v. 1.4 - Use new PyGal interface to LibGAL, new *--sort-by* option for satpass.py
The new PyGal interface to LibGAL uses the Object Oriented approach and dramatically simplifies the calling code. satpass.py now has sorting functionality. All generated passes are sorted by time (by default) but can be sorted by satellite name with the new *--sort-by* option.

### 2021/07/28 - v. 1.3 - Automatic Time Zone
When time zone is not specified (nether explicitly nor using a QTH name), the code tries to derive the time zone from the specified location.

### 2021/07/28 - v. 1.2 - Add orbit number
Calculate orbit number for elliptical orbits. Make it available in the pdf file. 