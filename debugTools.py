#
# Include these methods in your Python project.
#
#

import datetime

from errorMessages import *

debugOn = True
#debugOn = False
def printDebug(message):
	if debugOn:
		f = open("m2p2.log", "a")
		f.write( getTime() + " " + message + "\n")
		f.close()

def printError(message):
	print "(ERROR) " + message
	f = open("m2p2.log", "a")
	f.write( getTime() + " (ERROR)" + message + "\n")
	f.close()
	
def printTip(message):
	print "(TIP) " + message
	f = open("m2p2.log", "a")
	f.write( getTime() + " (TIP)" + message + "\n")
	f.close()

def getTime():
	return datetime.datetime.now().__str__()

