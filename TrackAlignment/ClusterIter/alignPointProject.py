import numpy

if __name__ == "__main__":
    # create UID, 
    DataPointCounter = 0
    FilenameCounter  = 0

    aRangeMin=-0.02
    aRangeMax=0.02
    aBin=0.005

    bRangeMin=-0.02
    bRangeMax=0.02
    bBin=0.005
    
    cRangeMin=-0.02
    cRangeMax=0.02
    cBin=0.005

    PointList=[]
    for dx0 in numpy.arange(aRangeMin,aRangeMax+aBin,aBin):
        x0=round(dx0,3)
        for dy0 in numpy.arange(bRangeMin,bRangeMax+bBin,bBin):
            y0=round(dy0,3)
            for dz0 in numpy.arange(cRangeMin,cRangeMax+cBin,bBin):
                z0=round(dz0,3)
                line="{}, {}, {}".format(x0,y0,z0)
                PointList.append(line)    # generate point condidate
                DataPointCounter=DataPointCounter+1
                print(DataPointCounter)
    
'''  
    PointInFile=100000
    f=open("{}.txt".format(FilenameCounter),"w")
    for item0 in PointList:
        for item1 in PointList:
            for item2 in PointList:
                for item3 in PointList:
                    for item4 in PointList:
                        dataPoint='{}, {}, {}, {}, {}, {}'.format(DataPointCounter, item0, item1,item2, item3, item4)
                        DataPointCounter=DataPointCounter+1
                        print(DataPointCounter)
'''