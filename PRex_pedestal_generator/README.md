# UVa decoder
A modified version of the UVa MPD4 decoder. 

## main change

* changed the bank number and the tag 
* change the data search structure of the 'inputhandler.cpp' so as to decode the PRex data 
* changed the grephic out style


## pedestal generator instructions

### 1. check the whether the mapping match the prex mapping
### 2. change the config/gem.cfg file 
    # runType
    RUNTYPE: CalPedestal
    # pedestal
    SAVEPED: ./Pedestal/     // this is the file name of the generated pedestal
    # # Input File for physics analysis; NOT FOR OTHER TYPE ANALYSIS

    INPUTFILE:   //  the pedestal raw run file name


### 3. generate the pedestal
    ./mpd4_decoder

### 4. save the file to the prex database
* it will generate a file named PRex_Pedestal.txt
* this file is conpatible with the standard prex pedestal, need to replace the pedestal in the prex database
