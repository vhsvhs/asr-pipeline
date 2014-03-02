#
# Builds an HTML report based on the results of the asr pipeline
#
#

import os, re, sys

from configuration import *
from tools import *
from html_helper import *



#########################
#
# main
#
write_css()
write_index()
write_alignments()
write_treesancs()
write_ancestors_indi() # write individual ancestor pages
write_anccomp()
















