from SCons.Builder import Builder
from SCons.Action import Action
from SCons.Subst import scons_subst
import re
from glob import glob
from functools import partial
import logging
import os.path
import os
import cPickle as pickle
import math
import xml.etree.ElementTree as et
import gzip
import subprocess
import shlex
import time
import shutil
import tempfile
import codecs
import locale
import bisect
import arpabo
from arpabo import ProbabilityList, Arpabo, Pronunciations, Vocabulary
from os.path import join as pjoin


def run_experiment(target, source, env):
    return None

def TOOLS_ADD(env):
    env.Append(BUILDERS = {"RunExperiment" : Builder(action=run_experiment),
                           })
