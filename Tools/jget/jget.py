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

class mssget(object):
    '''
    '''
    def __init__(self):
        pass
    
    def ReadListFile(self):
        '''
        read the file list from the file 
        '''
        pass
    
    def ReadList(self):
        '''
        need multi format support
        input number? input filename? input filename with path? or combined
        '''
    def ReadListNumber(self):
        '''
        read the file list from the input number
        '''
        pass
    
    def ReadListName(self):
        '''
        read the file list from the input number
        need to check whether this is a filename with path or not
        '''
        pass
    def _isFileWpath(self):
        '''
        Check the input number wether have path
        '''
        pass

    def IsFileExist(self,filename="text.txt"):
        pass
    
    def GetFiles(self):
        pass
    
    def _WriteToLog(self):
        pass
    
if __name__ == "__main__":
    '''
    if no input 
    if input numbe (or array of number )
    if input is a filename(with or with out path)
    '''
    pass