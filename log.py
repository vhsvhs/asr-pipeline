"""
log.py - a logging system for the ASR pipeline
"""
import os

LOGPATH = None
ERRPATH = None

def init_log(dir):
    """
    Initialize the log in the directory 'dir'
    If the log already exists, then it will be overwritten.
    """
    if os.path.exists(dir):
        LOGPATH = dir + "/log.txt"
        fout = open(LOGPATH, "w")
        fout.close()
        ERRPATH = dir + "/logerr.txt"
        fout = open(ERRPATH, "w")
        fout.close()
    else:
        print "ERROR: The output directory", dir, "doesn't exist!"

def write_log(checkpoint, message):
    """
    Writes to the log file
    """
    fout = open(LOGPATH, "a")
    fout.write(checkpoint.__str__() + "\t" + message)
    fout.close()


def write_error(message):
    """
    Writes to the log file
    """
    fout = open(ERRPATH, "a")
    fout.write(message)
    fout.close()