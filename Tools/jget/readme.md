# file name            save path [save in current path if not specified]
prexRHRS_20532.dat.0
prexRHRS_20532.dat.1
prexRHRS_20532.dat.2 


# Automaticly get the data file from the mss

### input format
#### 1. run number 

---
    ./jget 20110
    ./jget 2033
---

#### 2. run name 

---
    ./jget prexRHRS_20532.dat      will read all the files that match prexRHRS_20532.dat.*
    ./jget prexRHRS_20532.dat.0
---

#### 3. run list
 Read file list name must end with .txt, otherwise it will not take it
---
    ./jget runlist.txt      
---

### format of the run list ffile 

        run nmae (or run number) start run ID (optional) end run ID(optional)

       * prexRHRS_1111.dat,     0 , 1
        
       * prexRHRS_1111.dat
         
       * 1111，     0 , 1
    
       * 1111
       * 1111， 1 （will only read the one with file run ID and run ID） 

