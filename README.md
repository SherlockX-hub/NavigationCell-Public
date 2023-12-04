# NavigationCell-Public
Some methods to analyze some cells in navigation.

## BNT_modified:
Modified and corrected functions in the BNT package based on personal requests. Replace the original ones before using.  
BNT package refers to Behavioural Neurology Toolbox, (c) Vadim Frolov 2018 (https://bitbucket.org/cnc-ntnu/bnt/src/master).

## cellDetect:
Integrate the detection of place cell, grid cell, border cell, head direction cell and speed cell into one code. Users can execute the code with personalized inputs.  
Note: the threshold is 95th and stability is not used for detecting cells. If one wants to have more criteria, one can write custom scripts to detect cells by using output values.

## decode:
Use Naive Bayes Decoder to decode the position or head direction of mice. Modified the code from Kai Gao.
