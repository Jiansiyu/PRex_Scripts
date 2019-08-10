'''
* jget scipt used for get the file from the mss 
*     
*      Siyu Jian (sj9va@virginia.edu)
*
* By defaul, the looking path is located at  /mss/halla/happexsp/raw
* full name of the file need to be specified
* it also can take input / multi input 
*
* input(it can take any of those):
*       run list file [filename, save path[optinal]]
*       run number
*       run name 
*       run name with path
*       save path [optinal, if not specified, it will save to the current path]
*
* othertings:
*             It will check the existance of the file in the MSS, if not exist, skip it and write an log
*
'''
import os
import glob
class mssget(object):
    '''
    '''
    def __init__(self):
        self.mssPath="/mss/halla/happexsp/raw/" 
        pass
    
    def ReadListFile(self,RunListFileName="runList.txt"):
        '''
        read the file list from the file 
        add support for multi format support

        formate can be 

        run nmae (or run number) start run ID (optional) end run ID(optional)

        prexRHRS_1111.dat,     0 , 1
        or
        prexRHRS_1111.dat
        or 
        1111     0 , 1
        '''

        # parser the file 
        # only take txt file 
        if ".txt" in RunListFileName:
            pass
    
    def ReadList(self, avg=""):
        '''
        need multi format support
        input number? input filename? input filename with path? or combined
        '''
        fileList=[]
        
        if os.path.isfile() and ".txt" in arg:
            # this is a readlist txt file contains the files that want to readout
            pass
        
        elif avg.isnumeric():
            # this is the number that want to convert into an file name
            # and get all the files in the mss
            filename=self.ReadListWNumber(runNumber=int(avg))
            for item in self._GetAllRunsWithName(filename=filename):
                fileList.append(item)
        
        else:
            for item in self.ReadListWname(filename=avg)
                fileList.append(item)
        return fileList


    def ReadListWNumber(self,runNumber=123):
        '''
        read the file list from the input number
        '''
        filename=""
        if runNumber < 20000:
            filename="${1}prexRHRS_${2}.dat".format(self.mssPath,runNumber)
            
        else:
            filename="${1}prexLHRS_${2}.dat".format(self.mssPath,runNumber)
            
        return self.ReadListWname(filename=filename)

    def ReadListWname(self, filename="aa.dat.0"):
        '''
        Get all file list
        if the file end with run number, it will only read out this one 
        if the file end with '.dat', it will return the list all the runs ith different runID
        
        '''
        # form full filename

        if self._isFileWpath(filename=filename):
            filename="${1}${2}".format(self.mssPath,filename)

        filenamewithoutNumber=""
        if filename[-1].isdigit():
            filenamewithoutNumber=os.path.splitext(filename)[0]
        else:
            filenamewithoutNumber=filename
        
        # get all the files in the mss
        return glob.glob("${}.*".format(filenamewithoutNumber))
    
    def ReadListName(self,ListFilename=""):
        '''
        read the file list from the input number
        need to check whether this is a filename with path or not
        '''
        pass
    
    def _isFileWpath(self, filename=""):
        '''
        Check the input number wether have path
        '''
        return '/' in filename

    def IsFileExist(self,filename="text.txt"):
        '''
        Check file exist or not
        '''
        return os.path.isfile(filename)
    
    def GetFiles(self, target="", SavePath=""):
        jgetCMD="jget ${1} ${2}".format(target,SavePath)
        print("Get file with command ${}".format(jgetCMD))
        // start get the file and save the file into the folder 
    
    def _WriteToLog(self,log="", logLevel="1",LogType=""):
        pass
    
if __name__ == "__main__":
    '''
    if no input 
    if input numbe (or array of number )
    if input is a filename(with or with out path)
    '''
    pass