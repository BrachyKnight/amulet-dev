import os, errno 
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

    
