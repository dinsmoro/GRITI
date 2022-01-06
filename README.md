# GRITI: General Resource for Ionospheric Transient Investigations

## About GRITI
**GRITI** is an open-source, Python 3-based analysis tool for ionospheric activity. For several data sources it supports data download, loading, and analysis while for other data sources it only covers data loading and analysis. It is designed for flexibility and expandability in adding new data sources and analysis methods.

## Using GRITI
Everything in **GRITI** is launched from **GRITI_main.py**. The settings are from the beginning from the file to the line **`#!!!END OF SETTINGS!!!`**. The settings are grouped under data sources or analysis methods, if the analysis method spans multiple data sources. Data sources are loaded on run based on the analysis methods chosen in **GRITI_main.py**'s settings.

In **GRITI_main_config.ini you must set paths and relevant login info** (if you don't plan on using a data source that needs login info, then you don't need to fill that login info in).

## Supported Data Sources
### Automatic Download
- delta-vTEC (vertical total electron content) based on [Madrigal's TEC dataset](http://cedar.openmadrigal.org/list)
- [Kp Index](http://www.gfz-potsdam.de/en/kp-index/)
- [NASA OMNI Data](https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/hroformat.txt)

### User Must Download from Source
- [Haystack ISR](http://cedar.openmadrigal.org/list) (incoherent scatter radar)
- [Canadian magnetometer network](https://www.geomag.nrcan.gc.ca/data-donnee/sd-en.php)
- [AMPERE](http://ampere.jhuapl.edu/dataget/index.html)-derived ionospheric model estimates (model estimates not yet live)
- [LISN's TEC dataset](http://lisn.igp.gob.pe/) supports LISN's modified post-processed RINEX format

### Not Fully Implemented
- [Pokerflat ISR](http://cedar.openmadrigal.org/list) (needs more methods for data processing to achieve higher output quality)

## Selection of Supported Analysis Methods
- Keograms \[delta-vTEC, AMPERE-derived ionospheric model estimates, and Canadian magnetometer network]
- Radius-around-a-point averaging \[delta-vTEC]
- Movies & snapshots of movies \[delta-vTEC, AMPERE-derived ionospheric model estimates]
- RTI (range-time-intensity) \[ISR]
- FFT, Lomb-Scargle, & CPSD (cross-power spectral density) \[all]
- Walking FFT \[delta-vTEC, AMPERE-derived ionospheric model estimates, NASA OMNI data]

## Installation
To get the code from this repository, clone it using Git or click the green Code button at the top and "Download ZIP". All of the functions need to be in the same directory as **GRITI_main.py** but the data sources folder and other output folders can be set elsewhere in **GRITI_main_settings.ini**.

**GRITI** has dependencies on the following Python 3 packages: NumPy, Matplotlib, Scipy, h5py, Numba, Cartopy, Basemap, Astropy, timezonefinder, pytz, and html2text.

[Anaconda](https://www.anaconda.com/products/individual), a Python 3 distribution that includes many useful scientific packages, comes with many of those needed packages automatically. It doesn't come with:
- [Cartopy](https://scitools.org.uk/cartopy/docs/latest/) install with **conda** using `conda install -c anaconda cartopy`.
- [Basemap](https://matplotlib.org/basemap/) !deprecated, I am in the process of converting to Cartopy! install with **conda** using `conda install -c anaconda basemap`, if needed install higher-res maps with `conda install -c conda-forge basemap-data-hires`, and see [this link](https://stackoverflow.com/questions/52295117/basemap-import-error-in-pycharm-keyerror-proj-lib/53751941) for `KeyError: 'PROJ_LIB'` issues that may arise.
- [timezonefinder](https://timezonefinder.readthedocs.io/en/latest/) install with **conda** using `conda install -c conda-forge timezonefinder`.
- [html2text](https://github.com/Alir3z4/html2text/) install with **conda** using `conda install -c conda-forge html2text`.

_Note: On Windows, use Anaconda Prompt to input **conda** commands._

Anaconda comes with [Spyder](https://www.spyder-ide.org/), a Python IDE, that is a very easy place to run **GRITI** from. Open **GRITI_main.py** in Spyder, set the relevant settings, and hit the green "Run File" arrow to use **GRITI** quickly.

**GRITI** has been tested with Python 3.7, so later versions should work and earlier Python 3 versions may work. Python 2 probably won't work.

**GRITI** has only been tested with Windows 10, and paths may not be properly be setup to work in all instances (among other things, I'm sure). [Raise an issue](https://github.com/dinsmoro/GRITI/issues/new) on this repo to let me know if any issues arise from other OSes.

**GRITI** is still under some construction (like the Basemap to Cartopy conversion) and some functions may not work as-is. [Raise an issue](https://github.com/dinsmoro/GRITI/issues/new) on this repo to let me know about it.

**Find any issue?** [Raise that issue!](https://github.com/dinsmoro/GRITI/issues/new)

## Contributing
Follow this [code of conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/) and raise an issue or make a pull request with the new addition/fix. Make sure the new feature/fix is commented!

**GRITI** is simple as far as a software suite goes, and **GRITI_main.py** ideally should only have function calls. (I develop in **GRITI_main.py** before packaging as a function, so **GRITI_main.py** isn't only functions at the moment - but I wish it was!)

## GRITI Literature
### Cite this code with
Dinsmore, R., Mathews, J.D., Urbina, J., 2021. General resource for ionospheric transient investigations (GRITI): An open-source code developed in support of the Dinsmore et al. (2021) results. MethodsX 8, 101456. https://doi.org/10.1016/j.mex.2021.101456

### Papers that have used GRITI
Dinsmore, R., Mathews, J.D., Coster, A., Robinson, R.M., Sarkhel, S., Erickson, P.J., Urbina, J., 2021. Multi-instrument observations of SCIPS: 1. ISR and GPS TEC results. Journal of Atmospheric and Solar-Terrestrial Physics 213, 105515. https://doi.org/10.1016/j.jastp.2020.105515

---
**GRITI** was made by [Ross Dinsmore](https://github.com/dinsmoro), who was supported by the US National Science Foundation under Grant No. AGS-1241407 to The Pennsylvania State University. 
