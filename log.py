"""
log.py - a logging system for the ASR pipeline
"""
import os
from tools import *

def init_log(ap, overwrite=True):
    """
    Initialize the log in the directory 'dir'
    If the log already exists, then it will be overwritten.
    """
    
    action = "a"
    if overwrite:
        action = "w"        
    ap.params["logpath"]  = "log.txt"
    fout = open(ap.params["logpath"], action)
    fout.close()
    
    ap.params["errpath"] = "logerr.txt"
    fout = open(ap.params["errpath"], action)
    fout.close()
    return ap

def write_log(ap, message):
    """
    Writes to the log file
    """
    fout = open(ap.params["logpath"], "a")
    fout.write(ap.params["checkpoint"].__str__() + "\t" + message + "\n")
    fout.close()


def write_error(ap, message):
    """
    Writes to the log file
    """
    fout = open(ap.params["errpath"], "a")
    fout.write(ap.params["pending_checkpoint"].__str__() + "\t" + message + "\n")
    fout.close()