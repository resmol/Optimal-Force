#!/usr/bin/env python3.6
"""
Created on Wed Feb  2 09:36:24 2022

@author: Alejandro Jodra
"""

class MsgError(Exception):
    def __init__(self, errmsg):
        self.errmsg = errmsg
        
    def __str__(self):
        return "\n\n  ERROR: %s"%self.errmsg

class ElseError(MsgError):
    def __init__(self, option, desc):
        self.errmsg = "Option %s not implemented for %s!"%(option, desc)
    
class NIError(MsgError):
    def __init__(self):
        self.errmsg = "Functionality no implemented yet!"
        
class PureVirtualError(MsgError):
    """
    Use this to mark pure virtual functions in an abstract base class
    These have to be redefined by the inherited class.
    """
    def __init__(self):
        self.errmsg = "Virtual function is pure!"

class NegativeRootError(MsgError):
    """
    Use this to notice an error in the root calculation
    """
    def __init__(self,value):
        self.errmsg = "Yoy can not do a square root of a negative number! %f"%value
    
    def __call__(self):
        return True
