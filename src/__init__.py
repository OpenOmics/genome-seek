import os, sys
# Makes relative imports to work in Python 3.6
# without the need of '.' before the name of the
# package or py file.
# Allows for consistent syntax of relative imports 
# across python2 and python3.
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

# Ground source of truth for version information
try:
    # Import from root of project directory
    version = open('VERSION', 'r').readlines()[0].strip()
except IOError:
    # When namespace is __main__ 
    version = open('../VERSION', 'r').readlines()[0].strip()