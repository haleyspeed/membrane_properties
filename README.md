# membrane_properties.py
Python program to analyze Input Resistance and Capacitance from seal tests in whole-cell patch clamp

### About
Before every whole-cell experiment, the seal teast is recorded (-10mV step for 200ms, 10 sweeps at least) to measure the input resistance and capacitance offline. Using clampfit, the last 10 sweeps of the seal test are averaged and the area between the step onset and offset is taken as the charge (Q). The membrane capacitance is then roughly estimated by the equation:

  C = Q/V where Q = area in pA * ms and V = -10mV 
  
  
  The conversions to A * s and V are written into the program. The input resistance is calculated by the equation:
  
  V = IR where V = -10mV and I = difference in current from before the step onset and the steady-state current 75-100ms after onset.
  
  
  
### Assumptions
  CSV to be read into the file should contain the following columns of data:
  
  - mouse	
  - cell	
  - genotype	
  - treatment	
  - current(pA)	
  - tau1(1/s)	
  - charge(pA*s)	
  
  
  For my own records, I also include:
  - file  (name of the seal test file and channel if doing multiple simultaneous recordings)
  - cursors (positions of cursors for setting charge/area measurements or input current measurements)

  
  
### Things You Will Need to Change
  
  - file_name
  - os.chdir (working directory) 
  - column names if you don't use the ones I've listed above
  
  
### Output
  
  Four csv files:
  
  * file_name_per_cell.csv (Data per cell  with input resistance and capacitance calculations)
  * file_name_per_mouse.csv (Data aggregated per mouse  with input resistance and capacitance calculations [with desc stats per mouse])
  * file_name_per_cell_desc.csv (Descriptive statistics with each cell as an n of 1)
  * file_name_per_mouse_desc.csv (Descriptive statistics with average of cells per mouse as an n of 1 [more accurate])
  
  
  
  
### Important Notes
  
  - Save directory is the same as the source directory. To save in a different directory, place another os.chir() statement with the new output directory above the "Save Data" Section of the script.
  
  
### Reference Information

Author: Haley E. Speed, PhD

Department of Neurology and Neurotherapeutics

University of Texas Southwestern Medical Center, Dallas, TX

Copyright 2017
