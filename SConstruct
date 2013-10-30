import sys
import os.path
import logging
from glob import glob
import asr_tools
import kws_tools
import torque_tools

#
# load variable definitions from custom.py, and define them for SCons (seems like it should logically
# happen in the reverse order, but anyways...)
#
vars = Variables("custom.py")
vars.AddVariables(
    ("OUTPUT_WIDTH", "", 130),
    ("ATILLA_PATH", "", ""),
    ("SEQUITUR_PATH", "", ""),
    ("ATILLA_INTERPRETER", "", "${ATILLA_PATH}/tools/attila/attila"),
    ("OUTPUT_PATH", "", ""),
    ("ADD_WORDS", "", "/usr/bin/add_words"),
    ("BABEL_REPO", "", None),
    ("BABEL_RESOURCES", "", None),
    ("F4DE", "", None),
    ("INDUS_DB", "", None),
    ("JAVA_NORM", "", "${BABEL_REPO}/KWS/examples/babel-dryrun/javabin"),
    ("OVERLAY", "", None),
    ("LIBRARY_OVERLAY", "", "${OVERLAY}/lib:${OVERLAY}/lib64"),
    ("EXPERIMENTS", "", {}),
    ("LANGUAGES", "", {}),
    ("LOG_LEVEL", "", logging.INFO),
    ("LOG_DESTINATION", "", sys.stdout),
    )

#
# create the actual build environment we'll be using
#
env = Environment(variables=vars, ENV={}, TARFLAGS="-c -z", TARSUFFIX=".tgz",
                  tools=["default", "textfile"] + [x.TOOLS_ADD for x in [asr_tools, kws_tools, torque_tools]],
                  BUILDERS={"CopyFile" : Builder(action="cp ${SOURCE} ${TARGET}")}
                  )

#
# initialize the Python logging system (though we don't really use it in this build, could be useful later)
#
if isinstance(env["LOG_DESTINATION"], basestring):
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=env["LOG_LEVEL"], filename=env["LOG_DESTINATION"])
else:
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=env["LOG_LEVEL"])

#
# each Builder emits a string describing what it's doing (target, source, etc), but with thousands of
# input files, this can be huge.  If the string is larger than OUTPUT_WIDTH, replace some of its
# characters with an ellipsis.  You can get less-truncated output by running e.g. "scons -Q OUTPUT_WIDTH=300"
#
def print_cmd_line(s, target, source, env):
    if len(s) > int(env["OUTPUT_WIDTH"]):
        print s[:int(env["OUTPUT_WIDTH"]) - 10] + "..." + s[-7:]
    else:
        print s

env['PRINT_CMD_LINE_FUNC'] = print_cmd_line




#
# check whether custom.py has defined all the variables we need
#
# undefined = [x for x in vars.keys() if x not in env]
# if len(undefined) != 0:
#     print ""
# One or more parameters (%s) are undefined.  Please edit custom.py appropriately.
# To get started, 'cp custom.py.example custom.py'
# """ % (",".join(undefined))
#     env.Exit()



def experiment(substitutions={}):
    files = {"LANGUAGE_MODEL_FILE" : None, #os.path.join(models, "models", "lm.3gm.arpabo.gz"),
             "PRONUNCIATIONS_FILE" : None, #os.path.join(models, "models", "dict.test"),
             "VOCABULARY_FILE" : None, #os.path.join(models, "models", "vocab"),
             "DATABASE_FILE" : os.path.join(models, "segment", "babel106.dev.LimitedLP.seg.v1.db"),
             "MEL_FILE" : os.path.join(models, "models", "mel"),
             "PHONE_FILE" : os.path.join(models, "models", "pnsp"),
             "PHONE_SET_FILE" : os.path.join(models, "models", "phonesset"),
             "TAGS_FILE" : os.path.join(models, "models", "tags"),
             "PRIORS_FILE" : os.path.join(models, "models", "priors"),
             "TREE_FILE" : os.path.join(models, "models", "tree"),
             "TOPO_FILE" : os.path.join(models, "models", "topo.tied"),
             "TOPO_TREE_FILE" : os.path.join(models, "models", "topotree"),
             "WARP_FILE" : os.path.join(models, "adapt", "warp.lst"),
             "LDA_FILE" : os.path.join(models, "models", "30.mat"),
             }
    directories = {"PCM_PATH" : data,
                   "OUTPUT_PATH" : None, #os.path.join(output_path, "baseline"),
                   "CMS_PATH" : os.path.join(models, "adapt", "cms"),
                   "FMLLR_PATH" : os.path.join(models, "adapt", "fmllr"),
                   "MODEL_PATH" : os.path.join(models, "models"),
                   }
    parameters = {"SAMPLING_RATE" : 8000,                    
                  "FEATURE_TYPE" : "plp",
                  "AC_WEIGHT" : 0.053,
                  "MAX_ERROR" : 15000,
                  "USE_DISPATCHER" : False,
                  }

    for k in files.keys():
        files[k] = substitutions.get(k, files[k])

    for k in directories.keys():
        directories[k] = substitutions.get(k, directories[k])

    for k in parameters.keys():
        parameters[k] = parameters.get(k, parameters[k])

    return [env.Value(x) for x in [files, directories, parameters]]


def submit(path, setup):
    if not env["HAS_TORQUE"]:
        return
    construct = env.SubmitJob(os.path.join(path, "construct.timestamp"), 
                              [setup, env.Value({"name" : "construct",
                                                               "commands" : ["${ATTILA_PATH}/tools/attila/attila construct.py"],
                                                               "path" : path,
                                                               "other" : ["#PBS -W group_list=yeticcls"],
                                                               "interval" : 10,
                                                               })])

    test = env.SubmitJob(os.path.join(path, "test.timestamp"), 
                         [construct, env.Value({"name" : "test",
                                                "commands" : ["${ATTILA_PATH}/tools/attila/attila test.py -w %f -n %s -j $${PBS_ARRAYID} -l 1" % (.13, 50)],
                                                "path" : path,
                                                "array" : 20,
                                                "other" : ["#PBS -W group_list=yeticcls"],
                                                "interval" : 120,
                                                })])

    score = env.ScoreResults(os.path.join(path, "ctm", "scoring", "babel106.dev.sys"), 
                             [env.Value("106.dev"), env.Value(os.path.abspath(os.path.join(path, "ctm")))])
    env.Depends(score, test)
    return score





#
#
#
for (language, language_id), packs in env["LANGUAGES"].iteritems():
    for pack, locations in packs.iteritems():
        data = locations["data"]
        models = locations["models"]
        locale = locations["locale"]
        base_path = os.path.join("work", language, pack)
        #output_path = os.path.join(env["OUTPUT_PATH"], language, pack)
        #print output_path

        # baseline experiment
        baseline_path = os.path.join("work", "experiments", language, pack, "baseline", "baseline", "baseline")
        baseline_vocabulary = env.File(os.path.join(models, "models", "vocab"))
        baseline_pronunciations = env.File(os.path.join(models, "models", "dict.test"))
        baseline_languagemodel = env.File(os.path.join(models, "models", "lm.3gm.arpabo.gz"))

        baseline_experiment = env.CreateSmallASRDirectory(Dir(baseline_path), 
                                                          experiment({"VOCABULARY_FILE" : baseline_vocabulary.rstr(),
                                                                      "PRONUNCIATIONS_FILE" : baseline_pronunciations.rstr(),
                                                                      "LANGUAGE_MODEL_FILE" : baseline_languagemodel.rstr(),
                                                                      "OUTPUT_PATH" : baseline_path,
                                                                      })
                                                          )
        #scores.append(submit(baseline_path, baseline_experiment))
