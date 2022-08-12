# Calculating satellite passes - SatPass
### What is it for?
Calculate parameters of the coming satellite passes: 

- the exact time and azimuth of a satellite rising above the horizon (AOS)
- the time, azimuth and elevation of the satellite pass at the moment when it reaches its highest point
- the time and azimuth of the satellite when it goes below the horizon (LOS)
- satellite orbit number

These parameters can be represented in a visual form and a series of these visual representations for a number of coming satellite passes can be saved as a pdf file. Here is an example:

![](sample.jpg)

The number in the middle of the small circle is the maximum satellite elevation for the pass. The "filled" clock hands represent the time of AOS, the third hand represents the time of LOS. The other elements are self-explanatory.

### How can I Try It?
Installation of the code, especially the General Astrodynamics Library the code depends on, can be a time consuming task. The easiest way to try the code is to download and run a Docker image. It does not require any installation (except Docker, if it is not yet installed). Also, it should work on any operating system supporting Docker.

```
  docker pull awasys/satpass:v14
```
Run the image in a terminal:

```
  docker run -it --entrypoint /bin/bash awasys/satpass:v14
```

Now skip to the *How to use* section and play with the code, but return back here to find out how to copy the pdf files you produced to your host file system (because they are created in the docker container filesystem).

After your pdf files are generated, go to your host computer shell (the one where you can issue the docker command) and issue the following command:

```
docker ps -a

CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
44b70c6a3635        satpass:v14         "/bin/bash"         7 minutes ago       Up 7 minutes                            stupefied_rhodes
```

The output will look close to the one shown above, copy the container name where you played with the code. It will be displayed under the NAMES title. In this example this name is *stupefied_rhodes*. Now to copy your pdf file from the docker container to your host filesystem, use the following command:

```
docker cp stupefied_rhodes:/opt/satpass/passes/20220812-155434-FN03gu.pdf ./
```

Replace the container name and pdf file name with the names in your environment. This command will copy the pdf file into your current local directory. Now you can view or print it.

### What does it consist of?
There are 3 utilities written in Python:

**gettle.py** is used to download satellite TLE files. You have to keep your locally saved TLE files up-to-date to produce reliable satellite pass information.

**satpass.py** performs calculations for one or more satellites that are going to be observed from a single location. 

**passpdf.py** creates a pdf file with visual representations of the passes based on the data calculated by satpass.py.

In addition to Python 3 (at least version 3.9 recommended) and a number of Python modules the **satpass.py** utility requires a shared library to be installed: [General Astrodynamics Library - libgal](http://www.amsat-bda.org/GAL_Home.html). A copy of the distribution of the library version 0.6.0 is included in directory LibGAL. One of the required Python modules (PyX) needs TeX (LaTeX) to be installed in order to be able to generate *pdf* files.

### Where can it run?
Since the required library is written in C++ and the utilities are written in Python 3 (well, there is also a Cython component), in theory, it can run anywhere. However, it was tested under Linux on some Raspberry Pi models (including Pi Zero) and under Mac OS with Python, Cython and various Python packages installed using macports. 

### How to install?
There are 3 major steps:

1. LibGAL

  It is recommended that you start with installing LibGAL first. Installation instructions are included within its distribution package. The link mentioned above contains some documentation on the library.

2. Python Module *pygal*

  After you generate the **libgal** library, go to the **pygal** folder and update the **build.sh** script so that the variables point to the correct locations in your environment. There is no script for Windows (sorry), but, perhaps, a similar install could be done under Cygwin. You will need Cython to compile the **pygal** package. After successful run of the build you should have the **pygal.so** library that will be recognized by Python as an external module. It serves as a wrapper around **libgal**, so that the Python code can make use of this library.

3. Python Code

  Since the utilities are written in Python 3, they can simply be copied to a directory of your choice and run from there. Note, that there is a number of Python modules that have to be installed to satisfy the dependencies. One of them, PyX, will require TeX (LaTeX) product to be installed.

### How to configure?
There is a configuration file that is used to keep a list of satellites, an optional list of geographic locations and some other running options. The file is in the **yaml** format and can be created and/or edited by any text editor. A sample configuration file is provided: **satellites.conf** in the **data** sub-directory. When TLE data is downloaded from the internet, a small file containing TLE information is created for each satellite. These files will be created on the first successful run of the **gettle.py** utility. The location of the files is specified in the configuration. With provided configuration file this location is the same **data** sub-directory where the configuration file is located.

### How to use?

Start with refreshing (or populating, if this is your first time) the TLE data. To do it simply switch to the **passes** directory and run the command:

```
python gettle.py
```
When it gets its job done successfully, run the following command to create a *pdf* file with satellite passes:

```
python satpass.py -g FN03gu -a 270m -z -0400 so50 | python passpdf.py
``` 
The pass information will be calculated for satellite SO-50 observed from geographic location in the centre of the FN03gu grid square at the hight of 270m (you can also specify 886ft). All dates and times will be local for the timezone with the offset -4 hours. If you don't specify time zone offset, the code will attempt to determine it automatically. The name of the pdf file will be chosen automatically and will include the date and time of the beginning of first pass and the location name (grid square in this case). 

There are additional options for all utilities, use **--help** to display them (for example, **python gettle.py --help**).

Please, keep in mind that generating **pdf** files may take some time depending on how modern your computer is. 

### Docker File
The **docker** directory contains file **Dockerfile** describing how to build a Docker image. The file can be viewed to understand code dependencies and the steps to install it. The file can also be used to build a Docker image. To do this, open a terminal, switch to the directory containing this README file and run the following command:

```
  docker build -t satpass -f ./Docker/Dockerfile .
```
To run the built image:

```
  docker run -it --entrypoint /bin/bash satpass
```

### What is the License?
The Python code is distributed under the "MIT" (some also call it "Expat") license. However, because of the fact that the **libgal** library is distributed under GNU GPL 2, the entire package is governed by the **GNU GPL 2** license.
