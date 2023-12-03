# NavigationCell-Public
Some methods to analysis some cells in navigation.

## all:
Integrate the detection of place cell, grid cell, border cell, head direction cell and speed cell into one code. Users can excute the code by personalized inputs.
Note: the threshold is 95th and stability is not used for detecting cells. If one wants to have more criteria, they can write custom scripts to detect cells by using ouput values.
BNT package refers to Behavioural Neurology Toolbox, (c) Vadim Frolov 2018 (https://bitbucket.org/cnc-ntnu/bnt/src/master).

## decode:
Use Naive Bayes Decoder to decode the position or head direction of mice. Modified the code from Kai Gao.