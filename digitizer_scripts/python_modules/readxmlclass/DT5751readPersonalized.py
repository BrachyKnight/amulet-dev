import os, errno 
import gc
from lxml import etree as ET

class DT5751reader(object):
    def __init__(self, filename):
        if not os.path.isfile(filename):
            raise IOError(errno.ENOENT, os.strerror(errno.ENOENT), filename)

        #self.__context = ET.iterparse(filename,tag='trace')
        
        self.__context = ET.iterparse(filename,
                                      #events=('start','end'),
                                      tag=('event'))
        
        return


    def get(self):

        try:
            elem = next(self.__context)[1]
        except Exception as ex:
            if type(ex).__name__ in ['StopIteration', 'XMLSyntaxError']:
                return None
            else: 
                #print (type(ex).__name__, ex.args)
                raise ex
            
        ret=dict(elem.attrib)
        ret.setdefault('channels', {})
        
        for e in elem:
            if e.tag == 'trace':
                ret['channels'].update({int(e.attrib['channel']): list(map(float, e.text.rstrip('\n ').split())) })
        
        return ret
        
    def generate_datafile(self, txtOutDir):
        g = open(txtOutDir,"w+")
        Nlines = 0
        while True:
            evt = self.get()
            if evt is None: break
            chs = evt.get("channels")
            chskeys = chs.keys()
            evtid = evt.get("id")
            if len(set(len(chs.get(chskey)) for chskey in chskeys)) != 1: #check all channels have same numbers of samples
                raise IndexError('cannels have different numbers of total samples')
            for i,_ in enumerate(chs.get(0)):
                line = ""
                for chskey in chskeys:
                    line += (str(round(chs.get(chskey)[i])) + "\t")
                g.write(line+"\n")
                Nlines += 1
        Nevts = evtid
        if Nlines != int(Nevts)*len(chs.get(0)):
            print( "\nNevents =", Nevts, "\nNlines =", Nlines, "\n(samples in one event) =",len(chs.get(0)), "\nNevents*(samples in one event) =",  int(Nevts)*len(chs.get(0)))
            raise IndexError('N lines in file is not equal to Nevents*(samples in one event)')
        g.close()
        f = open(txtOutDir.replace(".txt","_infos.txt"),"w+") #TODO find a way to save also other infos in the file (in particular sample frequency)
        f.write(str(Nevts)+"\n"+str(len(chs.get(0)))+"\n"+str(len(chskeys))+"\n") #save infos in another file: Nevts, SamplesInOneEvent, Nchannels
        f.close()
