import sys
import os.path
from os.path import join as pjoin
import logging
from functools import partial
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
    ("ATTILA_PATH", "", ""),
    ("SEQUITUR_PATH", "", ""),
    ("ATTILA_INTERPRETER", "", "${ATTILA_PATH}/tools/attila/attila"),
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
    ("HAS_TORQUE", "", False),
    ("PYTHON_INTERPRETER", "", None),
    ("SCORE_SCRIPT", "", None),
    ("INDUS_DB", "", None),
    ("SCLITE_BINARY", "", None),
    ("JOBS", "", 20),
    ("IBM_MODELS", "", None),
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


def run_experiment(jobs=20, default_files={}, default_directories={}, default_parameters={}, **args):
    files = {k : v for k, v in default_files.iteritems()}
    directories = {k : v for k, v in default_directories.iteritems()}
    parameters = {k : v for k, v in default_parameters.iteritems()}
    for k, v in args.iteritems():
        if k.endswith("FILE"):
            files[k] = v
        elif k.endswith("PATH"):
            directories[k] = v
        else:
            parameters[k] = v

    # files = {k : f for k, f in ["MEL", "PHONE", "PHONE_SET", "TAGS", "PRIORS", "TREE", "TOPO", "TOPO_TREE", "WARP", "LDA"]}
    # files = {"VOCABULARY_FILE" : vocabulary,
    #          "PRONUNCIATIONS_FILE" : pronunciations,
    #          "LANGUAGE_MODEL_FILE" : language_model,
    #          "DATABASE_FILE" : database,
    #          "MEL_FILE" : os.path.join(model_path, "mel"),
    #          "PHONE_FILE" : os.path.join(model_path, "pnsp"),
    #          "PHONE_SET_FILE" : os.path.join(model_path, "phonesset"),
    #          "TAGS_FILE" : os.path.join(model_path, "tags"),
    #          "PRIORS_FILE" : os.path.join(model_path, "priors"),
    #          "TREE_FILE" : os.path.join(model_path, "tree"),
    #          "TOPO_FILE" : os.path.join(model_path, "topo.tied"),
    #          "TOPO_TREE_FILE" : os.path.join(model_path, "topotree"),
    #          "WARP_FILE" : os.path.join(model_path, "warp.lst"),
    #          "LDA_FILE" : os.path.join(model_path, "30.mat"),
    #          }
    # directories = {"PCM_PATH" : data,
    #                "OUTPUT_PATH" : path,
    #                "CMS_PATH" : os.path.join(adapt_path, "adapt", "cms"),
    #                "FMLLR_PATH" : os.path.join(adapt_path, "adapt", "fmllr"),
    #                "MODEL_PATH" : os.path.join(model_path),
    #                }
    # parameters = {"SAMPLING_RATE" : 8000,                    
    #               "FEATURE_TYPE" : "plp",
    #               "AC_WEIGHT" : 0.053,
    #               "MAX_ERROR" : 15000,
    #               "USE_DISPATCHER" : False,
    #               }

    experiment = env.CreateSmallASRDirectory(Dir(args["OUTPUT_PATH"]), [env.Value(x) for x in [files, directories, parameters]])

    if not env["HAS_TORQUE"]:
        return experiment

    construct = env.SubmitJob(os.path.join(path, "construct.timestamp"), 
                              [env.Dir(path), env.Value({"name" : "construct",
                                                         "commands" : ["${ATTILA_PATH}/tools/attila/attila construct.py"],
                                                         "path" : args["OUTPUT_PATH"],
                                                         "other" : ["#PBS -W group_list=yeticcls"],
                                                         "interval" : 10,
                                                         })])

    test = env.SubmitJob(os.path.join(path, "test.timestamp"), 
                         [construct, env.Value({"name" : "test",
                                                "commands" : ["${ATTILA_PATH}/tools/attila/attila test.py -w %f -n %s -j $${PBS_ARRAYID} -l 1" % (.13, jobs)],
                                                "path" : args["OUTPUT_PATH"],
                                                "array" : jobs,
                                                "other" : ["#PBS -W group_list=yeticcls"],
                                                "interval" : 120,
                                                })])

    #score = env.ScoreResults(os.path.join(path, "ctm", "scoring", "babel106.dev.sys"), 
    #                         [env.Value("106.dev"), env.Value(os.path.abspath(os.path.join(path, "ctm")))])
    #env.Depends(score, test)
    return test


    # files = {"VOCABULARY_FILE" : vocabulary,
    #          "PRONUNCIATIONS_FILE" : pronunciations,
    #          "LANGUAGE_MODEL_FILE" : language_model,
    #          "DATABASE_FILE" : database,
    #          "MEL_FILE" : os.path.join(model_path, "mel"),
    #          "PHONE_FILE" : os.path.join(model_path, "pnsp"),
    #          "PHONE_SET_FILE" : os.path.join(model_path, "phonesset"),
    #          "TAGS_FILE" : os.path.join(model_path, "tags"),
    #          "PRIORS_FILE" : os.path.join(model_path, "priors"),
    #          "TREE_FILE" : os.path.join(model_path, "tree"),
    #          "TOPO_FILE" : os.path.join(model_path, "topo.tied"),
    #          "TOPO_TREE_FILE" : os.path.join(model_path, "topotree"),
    #          "WARP_FILE" : os.path.join(model_path, "warp.lst"),
    #          "LDA_FILE" : os.path.join(model_path, "30.mat"),
    #          }
    # directories = {"PCM_PATH" : data,
    #                "OUTPUT_PATH" : path,
    #                "CMS_PATH" : os.path.join(adapt_path, "adapt", "cms"),
    #                "FMLLR_PATH" : os.path.join(adapt_path, "adapt", "fmllr"),
    #                "MODEL_PATH" : os.path.join(model_path),
    #                }
    # parameters = {"SAMPLING_RATE" : 8000,                    
    #               "FEATURE_TYPE" : "plp",
    #               "AC_WEIGHT" : 0.053,
    #               "MAX_ERROR" : 15000,
    #               "USE_DISPATCHER" : False,
    #               }



#
#
#
experiments = []
general_run = partial(run_experiment, SAMPLING_RATE=8000, FEATURE_TYPE="plp", AC_WEIGHT=.13, MAX_ERROR=15000, USE_DISPATCHER=False)
for (language, language_id), packs in env["LANGUAGES"].iteritems():
    for pack, locations in packs.iteritems():

        model_path = pjoin(env["IBM_MODELS"], str(language_id), pack, "models")
        adapt_path = pjoin(env["IBM_MODELS"], str(language_id), pack, "adapt")
        language_pack_run = partial(general_run,
                                    PCM_PATH=locations["data"],
                                    CMS_PATH=pjoin(env["IBM_MODELS"], str(language_id), pack, "adapt", "cms"),
                                    FMLLR_PATH=pjoin(env["IBM_MODELS"], str(language_id), pack, "adapt", "fmllr"),
                                    MODEL_PATH=model_path,
                                    MEL_FILE=pjoin(model_path, "mel"),
                                    PHONE_FILE=pjoin(model_path, "pnsp"),
                                    PHONE_SET_FILE=pjoin(model_path, "phonesset"),
                                    TAGS_FILE=pjoin(model_path, "tags"),
                                    PRIORS_FILE=pjoin(model_path, "priors"),
                                    TREE_FILE=pjoin(model_path, "tree"),
                                    TOPO_FILE=pjoin(model_path, "topo.tied"),
                                    TOPO_TREE_FILE=pjoin(model_path, "topotree"),
                                    WARP_FILE=pjoin(adapt_path, "warp.lst"),
                                    LDA_FILE=pjoin(model_path, "30.mat"),                                             
                                    )
                                       
                                       #MODEL_PATH=pjoin(env["IBM_MODELS"], str(language_id), pack, "models"),

        #ibm_path = os.path.join(env["IBM_MODELS"], str(language_id), pack)
        #model_path = os.path.join(ibm_path, "models")
        #adapt_path = os.path.join(ibm_path, "adapt")
        #data = locations["data"]
        #locale = locations["locale"]


        #output_path = os.path.join("work", language, pack)
        # baseline experiment
        #baseline_path = os.path.join("work", "experiments", language, pack, "baseline", "baseline", "baseline")
        #baseline_vocabulary = env.File(os.path.join(models, "models", "vocab"))
        #baseline_pronunciations = env.File(os.path.join(models, "models", "dict.test"))
        #baseline_languagemodel = env.File(os.path.join(models, "models", "lm.3gm.arpabo.gz"))

        experiments.append(language_pack_run(OUTPUT_PATH=os.path.join("work", language, pack, "baseline"),
                                             VOCABULARY_FILE=locations["vocabulary"], 
                                             PRONUNCIATIONS_FILE=locations["pronunciations"], 
                                             LANGUAGE_MODEL_FILE=locations["language_model"], 
                                             DATABASE_FILE=locations["database"],
                                             ))

#         baseline_experiment = env.CreateSmallASRDirectory(Dir(baseline_path), 
#                                                           experiment({"VOCABULARY_FILE" : baseline_vocabulary.rstr(),
#                                                                       "PRONUNCIATIONS_FILE" : baseline_pronunciations.rstr(),
#                                                                       "LANGUAGE_MODEL_FILE" : baseline_languagemodel.rstr(),
#                                                                       "OUTPUT_PATH" : baseline_path,
#                                                                       })
#                                                           )
#         scores.append(submit(baseline_path, baseline_experiment, jobs=env["JOBS"]))


#         # triple-oracle experiment
#         oracle_pronunciations, oracle_pnsp, oracle_tags = env.AppenToAttila([os.path.join(base_path, x) for x in ["oracle_pronunciations.txt", "oracle_pnsp.txt", "oracle_tags.txt"]],
#                                                                         [os.path.join(data, "conversational", "reference_materials", "lexicon.txt"),
#                                                                          os.path.join(data, "scripted", "reference_materials", "lexicon.txt"),
#                                                                          env.Value(locale)])
#         oracle_text, oracle_text_words = env.CollectText([os.path.join(base_path, x) for x in ["oracle_text.txt", "oracle_text_words.txt"]], 
#                                                          [env.Dir(x) for x in glob(os.path.join(data, "*/*/transcription"))])
#         oracle_vocabulary = env.PronunciationsToVocabulary(os.path.join(base_path, "oracle_vocabulary.txt"), oracle_pronunciations)
#         oracle_languagemodel = env.IBMTrainLanguageModel(os.path.join(base_path, "oracle_lm.3gm.arpabo.gz"), [oracle_text, oracle_text_words, env.Value(3)])


#         oracle_path = os.path.join("work", "experiments", language, pack, "oracle", "oracle", "oracle")
#         triple_oracle_experiment = env.CreateSmallASRDirectory(Dir(os.path.join("work", "experiments", language, pack, "oracle", "oracle", "oracle")),
#                                                                experiment({"VOCABULARY_FILE" : oracle_vocabulary[0].rstr(),
#                                                                            "PRONUNCIATIONS_FILE" : oracle_pronunciations.rstr(),
#                                                                            "LANGUAGE_MODEL_FILE" : oracle_languagemodel[0].rstr(),
#                                                                            "OUTPUT_PATH" : oracle_path,
#                                                                            })
#                                                                )

#         scores.append(submit(oracle_path, triple_oracle_experiment, jobs=env["JOBS"]))

#         # babelgum experiments
#         for model, (probs, prons) in locations.get("babelgum", {}).iteritems():
#             for size in [50000]: #[5000, 10000, 50000]:
#                 babelgum_probabilities, babelgum_pronunciations = env.BabelGumLexicon([os.path.join(base_path, x) for x in ["babelgum_%s_%d_probabilities.txt" % (model, size), 
#                                                                                                                             "babelgum_%s_%d_pronunciations.txt" % (model, size)]], 
#                                                                                       [probs, prons, env.Value(size)])
#                 env.NoClean([babelgum_probabilities, babelgum_pronunciations])


#                 babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities = \
#                     env.FilterBabelGum([os.path.join(base_path, "babelgum_rightwords_%d_%s.txt" % (size, x)) for x in ["pronunciations", "probabilities"]],
#                                        [babelgum_pronunciations, babelgum_probabilities, oracle_vocabulary])


#                 for weight in [.1]: #[.01, .05, .1, .4]:

#                     # adding all babelgum vocabulary (uniform)
#                     babelgum_uniform_vocabulary, babelgum_uniform_pronunciations, babelgum_uniform_languagemodel = env.AugmentLanguageModel(
#                         [os.path.join(base_path, "babelgum_uniform_%s_%d_%f_%s" % (model, size, weight, x)) for x in ["vocabulary.txt", "pronunciations.txt", "lm.3gm.arpabo.gz"]],
#                         [baseline_pronunciations, baseline_languagemodel, babelgum_pronunciations, env.Value(weight)]
#                         )


#                     babelgum_uniform_path = os.path.join("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum", "babelgum", "uniform")
#                     babelgum_uniform_experiment = env.CreateSmallASRDirectory(Dir(babelgum_uniform_path),
#                                                                               experiment({"VOCABULARY_FILE" : babelgum_uniform_vocabulary.rstr(),
#                                                                                           "PRONUNCIATIONS_FILE" : babelgum_uniform_pronunciations.rstr(),
#                                                                                           "LANGUAGE_MODEL_FILE" : babelgum_uniform_languagemodel.rstr(),
#                                                                                           "OUTPUT_PATH" : babelgum_uniform_path,
#                                                                                           })
#                                                                               )
#                     scores.append(submit(babelgum_uniform_path, babelgum_uniform_experiment, jobs=env["JOBS"]))

#                     # adding all babelgum vocabulary (weighted)
#                     babelgum_weighted_vocabulary, babelgum_weighted_pronunciations, babelgum_weighted_languagemodel = env.AugmentLanguageModel(
#                         [os.path.join(base_path, "babelgum_weighted_%s_%d_%f_%s" % (model, size, weight, x)) for x in ["vocabulary.txt", "pronunciations.txt", "lm.3gm.arpabo.gz"]],
#                         [baseline_pronunciations, baseline_languagemodel, babelgum_pronunciations, babelgum_probabilities, env.Value(weight)]
#                         )

#                     babelgum_weighted_path = os.path.join("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum", "babelgum", "babelgum")
#                     babelgum_weighted_experiment = env.CreateSmallASRDirectory(Dir(babelgum_weighted_path),
#                                                                               experiment({"VOCABULARY_FILE" : babelgum_weighted_vocabulary.rstr(),
#                                                                                           "PRONUNCIATIONS_FILE" : babelgum_weighted_pronunciations.rstr(),
#                                                                                           "LANGUAGE_MODEL_FILE" : babelgum_weighted_languagemodel.rstr(),
#                                                                                           "OUTPUT_PATH" : babelgum_weighted_path,
#                                                                                           })
#                                                                               )
#                     scores.append(submit(babelgum_weighted_path, babelgum_weighted_experiment, jobs=env["JOBS"]))
                    

#                     # adding just correct babelgum terms (uniform)
#                     babelgum_rightwords_uniform_vocabulary, babelgum_rightwords_uniform_pronunciations, babelgum_rightwords_uniform_languagemodel = env.AugmentLanguageModel(
#                         [os.path.join(base_path, "babelgum_rightwords_uniform_%s_%d_%f_%s" % (model, size, weight, x)) for x in ["vocabulary.txt", "pronunciations.txt", "lm.3gm.arpabo.gz"]],
#                         [baseline_pronunciations, baseline_languagemodel, babelgum_rightwords_pronunciations, env.Value(weight)]
#                         )

#                     babelgum_rightwords_uniform_path = os.path.join("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum_corrected", "babelgum", "uniform")
#                     babelgum_rightwords_uniform_experiment = env.CreateSmallASRDirectory(Dir(babelgum_rightwords_uniform_path),
#                                                                                          experiment({"VOCABULARY_FILE" : babelgum_rightwords_uniform_vocabulary.rstr(),
#                                                                                                      "PRONUNCIATIONS_FILE" : babelgum_rightwords_uniform_pronunciations.rstr(),
#                                                                                                      "LANGUAGE_MODEL_FILE" : babelgum_rightwords_uniform_languagemodel.rstr(),
#                                                                                                      "OUTPUT_PATH" : babelgum_rightwords_uniform_path,
#                                                                                                      })
#                                                                                          )
#                     scores.append(submit(babelgum_rightwords_uniform_path, babelgum_rightwords_uniform_experiment, jobs=env["JOBS"]))


#                     # adding just correct babelgum terms (weighted)
#                     babelgum_rightwords_weighted_vocabulary, babelgum_rightwords_weighted_pronunciations, babelgum_rightwords_weighted_languagemodel = env.AugmentLanguageModel(
#                         [os.path.join(base_path, "babelgum_rightwords_weighted_%s_%d_%f_%s" % (model, size, weight, x)) for x in ["vocabulary.txt", "pronunciations.txt", "lm.3gm.arpabo.gz"]],
#                         [baseline_pronunciations, baseline_languagemodel, babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities, env.Value(weight)]
#                         )

#                     babelgum_rightwords_weighted_path = os.path.join("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum_corrected", "babelgum", "weighted")
#                     babelgum_rightwords_weighted_experiment = env.CreateSmallASRDirectory(Dir(babelgum_rightwords_weighted_path),
#                                                                                          experiment({"VOCABULARY_FILE" : babelgum_rightwords_weighted_vocabulary.rstr(),
#                                                                                                      "PRONUNCIATIONS_FILE" : babelgum_rightwords_weighted_pronunciations.rstr(),
#                                                                                                      "LANGUAGE_MODEL_FILE" : babelgum_rightwords_weighted_languagemodel.rstr(),
#                                                                                                      "OUTPUT_PATH" : babelgum_rightwords_weighted_path,
#                                                                                                      })
#                                                                                          )
#                     scores.append(submit(babelgum_rightwords_weighted_path, babelgum_rightwords_weighted_experiment, jobs=env["JOBS"]))

# env.CollateResults("work/results.txt", scores)
