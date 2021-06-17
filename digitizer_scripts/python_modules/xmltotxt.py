import sys
#sys.path.append('../python_modules') #if you change position of the file tou need to uncomment and modify this line
from readxmlclass.DT5751readPersonalized import DT5751reader as digiread

digidata = digiread(sys.argv[1]) #read data from file at directory specified in sys.argv[1]
digidata.generate_datafile(sys.argv[2]) #write data from file at directory specified in sys.argv[2]
