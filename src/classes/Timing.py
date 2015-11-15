'''
Created on 07/02/2013

@author: http://effbot.org/librarybook/timing-example-2.py
'''
import time
from sys import platform

class Timing:
    def __init__(self):
        if platform == "win32":
            self.timer = time.clock
        else:
            self.timer = time.time
        
        self.t0 = self.t1 = 0

    def start(self):
        self.t0 = self.timer()
    
    def finish(self):
        self.t1 = self.timer()
    
    def seconds(self):
        return int(self.t1 - self.t0)
    
    def milli(self):
        return int((self.t1 - self.t0) * 1000)
    
    def micro(self):
        return int((self.t1 - self.t0) * 1000000)
