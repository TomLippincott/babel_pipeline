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
from morfessor import BaselineModel, AnnotatedCorpusEncoding, AnnotationsModelUpdate, LexiconEncoding, CorpusEncoding, Encoding, MorfessorException, MorfessorIO

def transcript_to_morfessor(target, source, env):
    return None

def TOOLS_ADD(env):
    env.Append(BUILDERS = {"TranscriptToMorfessor" : Builder(action=transcript_to_morfessor),
                           })
