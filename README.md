# AMULET
<!-- ABOUT THE PROJECT -->
## About The Project
AMULET (Analysis MUon LifETime) consists in a series of C++ and python scripts written to analyze data relative to the life-times of cosmic muons.
The project has two main parts:
* Reading data taken with a caen digitizer (CAEN DT5751). This is the most important part and it is actually the main focus of this project.
* Analyze data to produce plots, perform fit, and do everything is needed to measure in the most accurate way the lifetime of cosmic muons decaying both in carbon (scintillator material) and other materials. 

<!-- GETTING STARTED -->
## Getting Started
Using thes script is simple!

First of all ensure to have the proper prerequisites, otherwise install them
Download all the project and compile the script by executing the ./compile_scripts.sh script that is present in the folder called digitizer_scripts/

There are really only two things that is important that you do properly
* Using the correct prerequisites (in particular ROOT)
* **Work using the correct directories and the correct naming schemes for your files!**


### Prerequisites
* [ROOT CERN toolkit](https://root.cern/): this project is based on ROOT and, in particular, it uses RDataFrames to generate .root files containing TTree with all the relevant quantities that can be extracted from xml RAW data. Since when I developed the project (2021) RDataFrames and VecOps::RVec where under construction it is particularly important that, if you encounter problems, you install the same version I used: **ROOT v6.24/00**.
* python3
* lxml module of python3
  ```
  pip3 install lxml
  ```


<!-- USAGE EXAMPLES -->
## Usage
_Contact me_

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated** :smile: .

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


<!-- CONTACT -->
## Contact

Your Name - massimogirola1@gmail.com - m.girola2@campus.unimib.it

Project Link: [https://github.com/mgirola/amulet](https://github.com/mgirola/amulet)


<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [ROOT CERN](https://root.cern/)
* [ROOT forum](https://root-forum.cern.ch/)

