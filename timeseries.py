'''
Created on 16 Feb 2021

@author: thomasgumbricht
'''

# Standard library imports

from tqdm import tqdm

# Third party imports

import numpy as np
 
# Package application imports

import geoimagine.support.karttur_dt as mj_dt

from geoimagine.gis import kt_gis as ktgis

from geoimagine.timeseries.numbautil import TimeSeriesNumba

class TScoreFuncs:
    ''' Time Series Core functions
    '''
    def __init__(self):
        '''Empty call to access the functions
        '''
        pass
    
    def _SetMinMax(self,a):
        '''Numpy array min, max and range extraction
        '''
        self.min = np.amin(a)
        self.max = np.amax(a)
        self.range = self.max-self.min
           
class TSdata(TScoreFuncs):
    '''Class for managing time series data read and write etc
    '''    
    def __init__(self, pp, session):
        TScoreFuncs.__init__(self)
        
        self.session = session
                
        self.pp = pp  
        
        self.verbose = self.pp.process.verbose 
        
        self.session._SetVerbosity(self.verbose)

    def _ResetSrcDstDatesLater(self,locus):
        ''' Adjust soruce and destination data based on what already exists
        '''
        
        # Adjust the output related to the input            
        if 'resample' in self.pp.process.processid.lower():
            
            firstSrcDate = self.process.srcperiod.datumD[self.srcDateL[0]]['acqdate']
            
            lastSrcDate = self.process.srcperiod.datumD[self.srcDateL[-1]]['acqdate']
            
            if self.process.srcperiod.timestep in ['M','monthlyday']:
            
                lastSrcDate = mj_dt.AddMonth(lastSrcDate, 1)

            datumL = list(self.pp.dstLayerD[locus].keys())

            acqdateL = [self.process.dstperiod.datumD[item]['acqdate'] for item in datumL]
            #self._maskLayer = True
            #self._maskLayerD = {'datum':self.srcDateL[0],'comp':srccompL[0]}

            #Loop all datums in output and make sure all the required input is available 
            okDateL = []
            notOkDateL = []
            if self.process.dstperiod.timestep[len(self.process.dstperiod.timestep)-1] =='D':
                for x,datum in enumerate(acqdateL):
                    startdate = datum
                    enddate = mj_dt.DeltaTime(startdate, self.process.dstperiod.periodstep-1)
                    if startdate >= firstSrcDate and enddate <= lastSrcDate:
                        okDateL.append(datumL[x])
                    else:
                        notOkDateL.append(datumL[x])
            elif self.process.dstperiod.timestep == 'M':
                for x,datum in enumerate(acqdateL):
                    startdate = datum
                    enddate = mj_dt.AddMonth(startdate, 1)
                    enddate = mj_dt.DeltaTime(enddate , -1)
                    if startdate >= firstSrcDate and enddate <= lastSrcDate:
                        okDateL.append(datumL[x])       
                    else:
                        print ('adding not ok', datumL[x],startdate, firstSrcDate, enddate, lastSrcDate)
                        notOkDateL.append(datumL[x])
                        if self.process.proc.acceptmissing:
                            firstSrcDate = startdate
                            lastSrcDate = enddate
                            
            elif self.process.dstperiod.timestep == 'A':
                for x,datum in enumerate(acqdateL):
                    startdate = datum
                    enddate = mj_dt.AddYear(startdate, 1)
                    enddate = mj_dt.DeltaTime(enddate , -1)
                    if startdate >= firstSrcDate and enddate <= lastSrcDate:
                        okDateL.append(datumL[x])       
                    else:
                        if self.process.proc.acceptmissing:
                            firstSrcDate = startdate
                            lastSrcDate = enddate
                        else:
                            print (startdate, firstSrcDate)
                            print (enddate, lastSrcDate)
                            print (datumL[x])
                            notOkDateL.append(datumL[x])
                            if self.process.proc.acceptmissing:
                                firstSrcDate = startdate
                                lastSrcDate = enddate

            else:
                exit('Unknown timestep in PrcessTimeSeries',self.process.dstperiod.timestep)

            if self.verbose:
                print ('notOkDateL',notOkDateL)
            print ('dsttadums', self.pp.dstLayerD[locus])
            #Remove the destination dates that are not OK unless acceptmissing is True
            if not self.process.proc.acceptmissing:
                for item in notOkDateL:
                    self.pp.dstLayerD[locus].pop(item)   
        
    def _OpenSrcLayers(self,locus):
        '''Open timeseries source datalayers for reading
        '''
        
        self.OpenSrcLayers = []
        
        self.AcqDateStrSrcLayerD = {}
        
        for comp in self.pp.srcCompL:
        
            firstLayer = True
            
            n = 0
            
            '''
            self.LayerForAcqDateStrD = {}
            
            self.srcDateD[comp] = []
            '''
            
            for datum in self.pp.srcLayerD[locus]:
                
                # if this data is non existent, continue to next
                if datum in self.pp.srcLayerDateNonExistD[locus][comp]:
                    
                    # src layer does not exist, set to none
                    self.pp.srcLayerD[locus][datum][comp] = False
                    
                    continue
                
                if self.pp.srcLayerD[locus][datum][comp].FPN in self.OpenSrcLayers:
                    
                    #Duplicate layer reading, for seasonal or static data combined with timesereis
    
                    self.pp.srcLayerD[locus][datum][comp].copy = True
                    
                    #self.AcqDateStrSrcLayerD[datum] = datum
                
                else:
                    n += 1
                    self.pp.srcLayerD[locus][datum][comp].copy = False
                    
                    self.pp.srcLayerD[locus][datum][comp].srcDS, self.pp.srcLayerD[locus][datum][comp].layer = ktgis.RasterOpenGetFirstLayer(self.pp.srcLayerD[locus][datum][comp].FPN,{'mode':'read'})
                    
                    # Append the openend layer to the OpenLayers list
                    self.OpenSrcLayers.append(self.pp.srcLayerD[locus][datum][comp].FPN)
                    
                    #self.LayerForAcqDateStrD[datum] = datum
                    
                    self.pp.srcLayerD[locus][datum][comp].layer.GetGeometry()
                    
                    lins = self.pp.srcLayerD[locus][datum][comp].layer.lins
                    
                    cols = self.pp.srcLayerD[locus][datum][comp].layer.cols
                    
                    if comp in self.srcLayerReaderD:
                        
                        if lins != self.srcLayerReaderD[comp]['lins']:
                            
                            exitstr = ('ERROR in number of lins in timeseries.ProcessTimeSeries._OpenSourceLayers')
                            
                            exit(exitstr)
                            
                        if cols != self.srcLayerReaderD[comp]['cols']:
                            
                            exitstr = ('ERROR in number of cols in timeseries.ProcessTimeSeries._OpenSourceLayers')
                            
                            exit(exitstr)
                            
                    else:
                        self.srcLayerReaderD[comp] = {'lins':lins,'cols':cols}
                    
                    l = self.pp.srcLayerD[locus][datum][comp].layer
                    
                    self.srcCellNullD[comp] = self.pp.srcLayerD[locus][datum][comp].layer.cellnull
                    
                    #if hasattr(self.pp.srcLayerD[locus][datum][comp].comp,'id'):
                    
                    #    self.idD[self.pp.srcLayerD[locus][datum][comp].comp.id] = comp
                    
                    if firstLayer: 
                        
                        self.geoFormatD[comp] = {'lins':l.lins,'cols':l.cols,'projection':l.projection,'geotrans':l.geotrans,'cellsize':l.cellsize}
                        
                        self.geoFormatD[locus] = {'lins':l.lins,'cols':l.cols,'projection':l.projection,'geotrans':l.geotrans,'cellsize':l.cellsize}
                        
                        self.firstFPN = self.pp.srcLayerD[locus][datum][comp].FPN
                        
                        firstLayer = False
                        
                    else:
                        
                        gfD = {'lins':l.lins,'cols':l.cols,'projection':l.projection,'geotrans':l.geotrans,'cellsize':l.cellsize}
                        
                        for item in self.geoFormatD[locus]:
                        
                            if self.geoFormatD[locus][item] != gfD[item]:
                            
                                if item == 'cellsize':
                                
                                    if round(self.geoFormatD[locus][item],4) != round(gfD[item],4):
                                    
                                        print ('layers can not be processed together (%s = %s and %s = %s' %(item, self.geoFormatD[locus][item], item, gfD[item]))
                                        
                                        exit()
                                
                                elif item == 'geotrans':
                                
                                    pass
                                
                                else:
                                
                                    print ('layers can not be processed together (%s = %s and %s = %s' %(item, self.geoFormatD[locus][item], item, gfD[item]))
                                    
                                    print ('layer 0 ', self.firstFPN)
                                   
                                    print ('layer 1', self.pp.srcLayerD[locus][datum][comp].FPN)
                                    
                                    print (self.geoFormatD[locus][item])
                                    
                                    print (gfD[item])
                                    
                                    exit()
                                                                
    def _CreateOpenDstRasterFiles(self):
        '''Create and open the destination layers
        '''
        
        for datum in self.pp.dstLayerD[self.locus]:
            
            for comp in self.pp.dstLayerD[self.locus][datum]:
                
                #Transfer the geoformat for this locus
                for item in self.geoDstFormatD:
                    
                    setattr(self.pp.dstLayerD[self.locus][datum][comp], item, self.geoDstFormatD[item])  
                
                if not self.pp.dstLayerD[self.locus][datum][comp]._Exists() or self.pp.process.overwrite:
                    
                    self.pp.dstLayerD[self.locus][datum][comp].dstDS = ktgis.RasterCreateWithFirstLayer(self.pp.dstLayerD[self.locus][datum][comp].FPN, self.pp.dstLayerD[self.locus][datum][comp])
                       
                    self.dstDateD[comp].append(datum)   
                
                else:
 
                    #self.session._InsertLayer(self.pp.dstLayerD[locus][datum][comp],self.process.overwrite,self.process.delete)
                    
                    self.pp.dstLayerD[self.locus][datum][comp] = False
                    
                    self.dstDateD[comp].append(False)
    
    def _SetReadBlocks(self,lins,cols,blocklines):
        imgSize =  lins*cols
        blockSize = cols*blocklines
        
        readitems = int(lins/blocklines)
        
        if readitems < lins/blocklines:
            lastFullReadItem = readitems
            lastBlockSize = imgSize - readitems*blockSize
            lastBlockLines = int(lastBlockSize/cols)
            
            readitems += 1
        else:
            lastBlockLines = 0
            lastFullReadItem  = readitems + blocklines #never reached
        return (imgSize, blockSize, readitems, lastBlockLines, lastFullReadItem)
    
    def _SetReadWriteDim(self):
        '''
        '''

        # Get min and max lins and cols
        maxlins = 0
        maxcols = 0
        minlins = 999999
        mincols = 999999
        
        # Loop over the soruce compostions 
        for comp in self.srcLayerReaderD:
            
            if self.srcLayerReaderD[comp]['lins'] > maxlins:
                
                maxlins = self.srcLayerReaderD[comp]['lins']
            
            if self.srcLayerReaderD[comp]['cols'] > maxcols:
            
                maxcols = self.srcLayerReaderD[comp]['cols']
            
            if self.srcLayerReaderD[comp]['lins'] < minlins:
            
                minlins = self.srcLayerReaderD[comp]['lins']
            
            if self.srcLayerReaderD[comp]['cols'] < mincols:
            
                mincols = self.srcLayerReaderD[comp]['cols']
            
        if minlins < maxlins:
            
            #the ratio must be an exact integer
            linsratio = float(maxlins)/float(minlins)
            
            if not linsratio.is_integer():
                
                exit('The ratio of layers at different resolutions must be an integer')
            
            colsratio = int(maxcols/mincols) 
            
            linsratio = int(linsratio)
            
            if not colsratio == linsratio:
            
                exit('The ratio of layers at different resolutions must be an integer')  
            
            blocklines = linsratio
            
            self.allSameSize = False
        
        elif mincols < maxcols:
        
            exit('The ratio of layers at different resolutions must be an integer') 
        
        else:
        
            blocklines = min(10,maxlins)
            
            self.allSameSize = True
            
        self.blocklines = blocklines

        #Set the read rules for min, only if min is less than max, otherwise the blockreading collapses
        if not self.allSameSize:
        
            imgSize, blockSize, readitems, lastBlockLines, lastFullReadItem = self._SetReadBlocks(minlins,mincols,1)
            
            blockratio = int(maxlins/minlins)
            
            for comp in self.srcLayerReaderD:
                
                if self.srcLayerReaderD[comp]['lins'] == minlins:
                    
                    self.srcLayerReaderD[comp]['imgsize'] = imgSize
                    
                    self.srcLayerReaderD[comp]['blocklines'] = 1
                    
                    self.srcLayerReaderD[comp]['blocksize'] = blockSize
                    
                    self.srcLayerReaderD[comp]['blockratio'] = blockratio
                    
                    self.srcLayerReaderD[comp]['readitems'] = readitems
                    
                    self.srcLayerReaderD[comp]['lastBlockLines'] = lastBlockLines
                    
                    lastBlockSize = imgSize - (readitems-1)*blockSize
                    
                    self.srcLayerReaderD[comp]['lastBlockSize'] = lastBlockSize

        #Set the read rules for max
        imgSize, blockSize, readitems, lastBlockLines, lastFullReadItem = self._SetReadBlocks(maxlins,maxcols,blocklines)
        
        self.readitems = readitems
        
        self.lastFullReadItem = lastFullReadItem
        
        blockratio = 1
        
        for comp in self.srcLayerReaderD:
        
            if self.srcLayerReaderD[comp]['lins'] == maxlins:
            
                self.srcLayerReaderD[comp]['imgsize'] = imgSize
                
                self.srcLayerReaderD[comp]['blocklines'] = blocklines
                
                self.srcLayerReaderD[comp]['blocksize'] = blockSize
                
                self.srcLayerReaderD[comp]['blockratio'] = blockratio
                
                self.srcLayerReaderD[comp]['readitems'] = readitems
                
                self.srcLayerReaderD[comp]['lastBlockLines'] = lastBlockLines
                
                self.srcLayerReaderD[comp]['lastFullReadItem'] = lastFullReadItem
                
                lastBlockSize = imgSize - (readitems-1)*blockSize

                self.srcLayerReaderD[comp]['lastBlockSize'] = lastBlockSize
                
                #Set the destination geoformat after the largest input
                
                self.geoDstFormatD = self.geoFormatD[comp]
         
    def _CreateDstNullArrays(self,n=0):
        ''' Create null arrays for each destination layer
        '''
        self.dstNullDarr = {}
        
        for comp in self.pp.dstCompL:
            
            cellnull = self.pp.dstCompD[comp].cellnull
            
            if not n:
                
                n = len(self.dstDateD[comp])
                
            self.dstNullDarr[comp] = np.ones(n, np.float32)
            
            self.dstNullDarr[comp] *= cellnull
               
    def _IniBlockRead(self,l):
        ''' Initate the block reading by setting the dimensions, then read the blocks
        '''

        self.wl = l*self.blocklines #self.wl for read and writeline in order for gdal to get info on where in the file to qrite

        #Read the data
        self.srcRD = {}
        
        for comp in self.pp.srcCompL:
                    
            blockSize = self.srcLayerReaderD[comp]['blocksize']
            
            blockLines = self.srcLayerReaderD[comp]['blocklines']
            
            #Check the reading to not go beyond the image 
            if l == self.lastFullReadItem:

                blockLines = self.srcLayerReaderD[comp]['lastBlockLines']
                
                blockSize = self.srcLayerReaderD[comp]['lastBlockSize']

            #Also reset the gobal bolcklines for writing
            self.blocklines = blockLines
            
            # Read the data from file
            self._BlockRead(self.locus,comp,blockSize,blockLines)

            #Set null to nan for the processing
            cellnull = self.srcCellNullD[comp]
            
            self.srcRD[comp][self.srcRD[comp]==cellnull]=np.nan

            '''
            if 'spl3smp' in self.process.proc.srccompD[comp]['product'].lower():
                if self.process.proc.processid.lower() not in ['seasonfilltsmodissingletile','seasonfilltssmap']:
                    #soil moisture = 0.02 is a default fallback in SMAP that is useless, better set to null in the resample
                    self.srcRD[comp][self.srcRD[comp] <= 0.02]=np.nan
            '''
           
    def _BlockRead(self,locus,comp,blockSize,blockLines):
        '''Blockread the source data
        '''
        
        # all time series processing is done with data as Float 32
        self.srcRD[comp] = np.ones( ( len(self.pp.srcPeriod.datumL), blockSize), np.float32)
        
        self.srcRD[comp] *= self.srcCellNullD[comp]
        
        x = 0

        for datum in self.pp.srcLayerD[locus]:
            
            if self.pp.srcLayerD[locus][datum][comp]:
                
                if self.pp.srcLayerD[locus][datum][comp].copy:
                
                    copydatum = self.AcqDateStrSrcLayerD[datum]
                    
                    srcDatum = self.LayerForAcqDateStrD[copydatum]
                    
                    self.srcRD[comp][x] = self.pp.srcLayerD[locus][srcDatum][comp].layer.NPBLOCK
                
                else:

                    startlin = int(self.wl/self.srcLayerReaderD[comp]['blockratio'])
                    
                    readcols = self.srcLayerReaderD[comp]['cols']
                    
                    self.pp.srcLayerD[locus][datum][comp].layer.ReadBlock(0,startlin,readcols,blockLines)

                    self.srcRD[comp][x] = self.pp.srcLayerD[locus][datum][comp].layer.NPBLOCK

            else:
                
                pass
            
            x += 1
            
    def _WriteBlock(self,row,datum,dstcomp):
        ''' Write data blocks to destination layers
        '''
        if self.blocklines == 1:
            
            R2D = np.atleast_2d(row)
            
        else:
            
            R2D = np.reshape(row, (self.blocklines,-1))

        self.pp.dstLayerD[self.locus][datum][dstcomp].dstDS.WriteBlock(0,self.wl,R2D) 
            
    def _CloseDstRasterFiles(self):
        ''' Close the created destination layers
        '''
        for datum in self.pp.dstLayerD[self.locus]:
            
            for comp in self.pp.dstLayerD[self.locus][datum]:
                
                if self.pp.dstLayerD[self.locus][datum][comp]: #Non created files set to False
                    
                    self.pp.dstLayerD[self.locus][datum][comp].dstDS._SetStats()
                                        
                    self.pp.dstLayerD[self.locus][datum][comp].dstDS.CloseDS()
                                       
                    self.pp.dstLayerD[self.locus][datum][comp] = None
                    
                    #self.session._InsertLayer(self.process.dstLayerD[self.locus][datum][comp],self.process.overwrite,self.process.delete)
    
class ProcessTimeSeries(TSdata,TimeSeriesNumba):
    ''' Main class for Time Series Processing
    ''' 
     
    def __init__(self, pp, session):
        '''
        '''
        
        # Initiate TScommon
        TSdata.__init__(self, pp, session)
        
        # Initate Time Series Numba
        TimeSeriesNumba.__init__(self)
        
        if self.verbose > 0:

            print ('        ProcessTimeSeries',self.pp.process.processid) 

        # timeseries processing is run location by location
        for locus in self.pp.dstLayerD:
            
            cont = True
            
            for comp in self.pp.dstCompL:
              
                if len(self.pp.dstLayerCreateD[locus][comp]) > 0:
                        
                    cont = False
                        
            if cont:
                
                continue
            
            self.locus = locus
            
            # Reset source and destination dates from the availbale data 
            # Skipped for now as uncldear  
            # self._ResetSrcDstDates(locus)
                                    
            # Layer Reader dict for storing metadata on how to read each composition
            self.srcLayerReaderD = {}
                
            # Dict for source data null
            self.srcCellNullD = {}
            
            # Dict for georformats of comp and locus
            self.geoFormatD = {}
            
            # Dict for all destination dates (both existing and non-existing)
            self.dstDateD = {}
            
            for comp in self.pp.dstCompL:
                
                self.dstDateD[comp] = []
             
            self._OpenSrcLayers(locus)
                   
            # Get the dimensions and reading/writing style for this locus
            self._SetReadWriteDim()
            
            # Create the output data
            self._CreateOpenDstRasterFiles()
            
            # Set destination nullarrays
            self._CreateDstNullArrays()
            
            # Loop the processes using data blocks
            self._ProcessLoop()
            
            # Close the destination raster files
            self._CloseDstRasterFiles()
            
            #self._CloseSrcRasterFiles()
                                                 
    def _ProcessLoop(self):
        ''' Process the timeseries block by block
        '''
        
        # Read the source data in chunks (blocks of lines)
        for l in tqdm(range(self.readitems)):
                        
            # Initiate and the read each block chunk
            self._IniBlockRead(l)
            
            # Direct to defined process
            if self.pp.process.processid.lower() == 'timeseriesimagemend':
                
                # Set the active composition
                self.activeComp = self.pp.dstCompL[0] 
        
                # Get the date list from the "local" class list (not the pp.process list(
                self.dstDateL = self.dstDateD[self.activeComp]
                
                # Start the processing
                self._TSFillNoData()
                
            else:
                
                exitstr = '%s process not available in Timeseries.ProcessTimesereis._ProcessLoop' %(self.pp.process.processid)
                
                exit(exitstr)
        
    def _TSFillNoData(self):
        '''Fill timeseries nodata by interpolation
        '''
        
        # Retrieve the array for the acitve compoent
        R = self.srcRD[self.activeComp]
                       
        # Use the numpy function apply_along_axis together with the numpy prucedure _FillAlongAxis
        Rr = np.apply_along_axis( self._FillAlongAxis, 0, R )
        
        # Rearrange the data and loop over all dates, write to file
        for x,row in enumerate(Rr):
            
            datum = self.dstDateL[x]
            
            if datum:
                
                self._WriteBlock(row,datum,self.activeComp)