import sys
import os.path
from os.path import join as pjoin
import logging
from functools import partial
from glob import glob
import vocabulary_tools
import asr_tools
import kws_tools
import torque_tools
import morfessor_tools
import babel_tools
import trmorph_tools
import re
from common_tools import meta_basename, meta_splitext

#
# load variable definitions from custom.py, and define them for SCons (seems like it should logically
# happen in the reverse order, but anyways...)
#
vars = Variables("custom.py")
vars.AddVariables(
    # just a limit to make SCons' output more readable
    ("OUTPUT_WIDTH", "", 130),

    # these variables determine what experiments are performed
    ("LANGUAGES", "", {}),
    ("RUN_ASR", "", False),
    ("RUN_KWS", "", False),
    ("EVALUATE_PRONUNCIATIONS", "", False),
    ("EXPANSION_SIZES", "", []),
    ("EXPANSION_WEIGHTS", "", []),

    # these variables determine how parallelism is exploited
    ("HAS_TORQUE", "", False),
    ("MAXIMUM_LOCAL_JOBS", "", 4),
    ("MAXIMUM_TORQUE_JOBS", "", 200),

    # these variables define the locations of various tools and data
    ("BASE_PATH", "", None),
    ("OVERLAY", "", "${BASE_PATH}/local"),
    ("LANGUAGE_PACKS", "", "${BASE_PATH}/language_transcripts"),
    ("IBM_MODELS", "", "${BASE_PATH}/ibm_models"),
    ("LORELEI_SVN", "", "${BASE_PATH}/lorelei_svn"),
    ("ATTILA_PATH", "", "${BASE_PATH}/VT-2-5-babel"),
    ("VOCABULARY_EXPANSION_PATH", "", "${BASE_PATH}/vocabulary_expansions"),
    ("INDUSDB_PATH", "", "${BASE_PATH}/lorelei_resources/IndusDB"),
    ("SEQUITUR_PATH", "", ""),
    ("ATTILA_INTERPRETER", "", "${ATTILA_PATH}/tools/attila/attila"),
    ("F4DE", "", None),
    ("JAVA_NORM", "", "${BABEL_REPO}/KWS/examples/babel-dryrun/javabin"),
    ("OVERLAY", "", None),
    ("LIBRARY_OVERLAY", "", "${OVERLAY}/lib:${OVERLAY}/lib64:${LORELEI_TOOLS}/boost_1_49_0/stage/lib/"),
    ("LOG_LEVEL", "", logging.INFO),
    ("LOG_DESTINATION", "", sys.stdout),    
    ("PYTHON_INTERPRETER", "", None),
    ("SCORE_SCRIPT", "", None),
    ("SCLITE_BINARY", "", "${BASE_PATH}/sctk-2.4.5/bin/sclite"),
    ("LORELEI_TOOLS", "", "${BASE_PATH}/lorelei_tools"),
    
    # these variables all have default definitions in terms of the previous, but may be overridden as needed
    ("PYTHON", "", "/usr/bin/python"),
    ("PERL", "", "/usr/bin/perl"),
    ("PERL_LIBRARIES", "", os.environ.get("PERL5LIB", "")),
    ("BABEL_BIN_PATH", "", "${LORELEI_SVN}/tools/kws/bin64"),
    ("BABEL_SCRIPT_PATH", "", "${LORELEI_SVN}/tools/kws/bin64"),
    ("F4DE", "", "${BABEL_RESOURCES}/F4DE"),
    ("INDUS_DB", "", "${BABEL_RESOURCES}/IndusDB"),
    ("WRD2PHLATTICE", "", "${BABEL_BIN_PATH}/wrd2phlattice"),
    ("BUILDINDEX", "", "${BABEL_BIN_PATH}/buildindex"),
    ("BUILDPADFST", "", "${BABEL_BIN_PATH}/buildpadfst"),
    ("FSTCOMPILE", "", "${OVERLAY}/bin/fstcompile"),
    ("QUERY2PHONEFST", "", "${BABEL_BIN_PATH}/query2phonefst"),
    ("STDSEARCH", "", "${BABEL_BIN_PATH}/stdsearch"),
    ("MERGESEARCHFROMPARINDEXPRL", "", "${BABEL_SCRIPT_PATH}/mergeSearchFromParIndex.prl"),
    ("MERGESCORESSUMPOSTNORMPL", "", "${BABEL_SCRIPT_PATH}/merge.scores.sumpost.norm.pl"),
    ("PRINTQUERYTERMLISTPRL", "", "${BABEL_SCRIPT_PATH}/printQueryTermList.prl"),
    ("F4DENORMALIZATIONPY", "", "${BABEL_SCRIPT_PATH}/F4DENormalization.py"),
    ("JAVA_NORM", "", "${BABEL_REPO}/examples/babel-dryrun/javabin"),    
    ("KWSEVALPL", "", "${F4DE}/KWSEval/tools/KWSEval/KWSEval.pl"),    
    )

#
# create the actual build environment we'll be using
#
env = Environment(variables=vars, ENV={}, TARFLAGS="-c -z", TARSUFFIX=".tgz",
                  tools=["default", "textfile"] + [x.TOOLS_ADD for x in [asr_tools, kws_tools, torque_tools, babel_tools, trmorph_tools, vocabulary_tools]],
                  BUILDERS={"CopyFile" : Builder(action="cp ${SOURCE} ${TARGET}")}
                  )

#
# calculate how much parallelism can exist for each SCons instance
#
num_instances = env.GetOption("num_jobs")
env.Replace(LOCAL_JOBS_PER_SCONS_INSTANCE=int(max(env["MAXIMUM_LOCAL_JOBS"] / num_instances, 1)))
env.Replace(TORQUE_JOBS_PER_SCONS_INSTANCE=int(max(env["MAXIMUM_TORQUE_JOBS"] / num_instances, 1)))

Help(vars.GenerateHelpText(env))

# throw an error if a file substitution cannot be made
AllowSubstExceptions()

#
# initialize the Python logging system (though we don't really use it in this build, could be useful later)
#
if isinstance(env["LOG_DESTINATION"], basestring):
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=int(env["LOG_LEVEL"]), filename=env["LOG_DESTINATION"])
else:
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=int(env["LOG_LEVEL"]))

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


"""
file conventions: text.gz, vocabulary.gz, probabilities.gz, language_model.gz, pronunciations.gz, frequencies.gz
"""


#
# morph, lm rerank, lm rerank w averaging, lm reranking for morpheme boundaries
#
experiments = []
properties = {}
figures = {}
results = {}
g2p_evals = []
for language, config in env["LANGUAGES"].iteritems():
    language_id = config["LANGUAGE_ID"]
    locale = config["LOCALE"]
    limited_basic_expansions_file = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "expansions", "simple.ex"))
    limited_bigram_expansions_file = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "expansions", "bigram.ex"))

    limited_basic_expansions = env.SplitExpansion([limited_basic_expansions_file, env.Value(50000)],
                                                  BASE_PATH=pjoin("work", "expansions", language, "limited", "basic"))

    limited_bigram_expansions = env.SplitExpansion([limited_bigram_expansions_file, env.Value(50000)],
                                                   BASE_PATH=pjoin("work", "expansions", language, "limited", "bigram"))

    oracle_text = env.CollectText(pjoin("work", "texts", language, "oracle_text.txt.gz"),
                                  [env.subst("${LANGUAGE_PACKS}/%d.tgz" % language_id), env.Value(".*(sub-train|dev)/transcription.*txt")],
                                  )

    training_text = env.CollectText(pjoin("work", "texts", language, "training_text.txt.gz"),
                                    [env.subst("${LANGUAGE_PACKS}/%d.tgz" % language_id), env.Value(".*sub-train/transcription.*txt")],
                                    )

    dev_text = env.CollectText(pjoin("work", "texts", language, "dev_text.txt.gz"),
                               [env.subst("${LANGUAGE_PACKS}/%d.tgz" % language_id), env.Value(".*dev/transcription.*txt")],
                               )

    oracle_vocabulary_file = env.TextToVocabulary(pjoin("work", "vocabularies", language, "oracle_vocabulary.txt.gz"),
                                                  oracle_text)

    training_vocabulary_file = env.TextToVocabulary(pjoin("work", "vocabularies", language, "training_vocabulary.txt.gz"),
                                                    training_text)

    dev_vocabulary_file = env.TextToVocabulary(pjoin("work", "vocabularies", language, "dev_vocabulary.txt.gz"),
                                               dev_text)

    #env.FilterBy(pjoin("work", "vocabularies", language, "oov_vocabulary.txt.gz"),
    #             [oracle_vocaulary_file, 

    for x in limited_basic_expansions + limited_bigram_expansions:
        (base, ext) = meta_splitext(x.rstr())
        exp_vocab = env.ProbabilityListToVocabulary("%s_vocabulary.gz" % (base),
                                                    x)
        env.FilterBy("%s_in_train_or_dev_vocabulary.gz" % (base),
                     [exp_vocab, oracle_vocabulary_file])

        if language == "turkish":
            env.TurkishFilter("%s_in_language_vocabulary.gz" % (base),
                              x)

    continue
    if os.path.exists(pjoin(env.subst(env["IBM_MODELS"]), str(language_id))):
        limited_pronunciations_file = env.File(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "dict.test")) #config["pronunciations"]
        limited_vocabulary_file = env.File(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "vocab"))
        limited_language_model_file = env.Glob(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "lm.*"))[0]
        

        for pack in ["LLP"]:

            # baseline experiment
            (asr_output, asr_score) = env.RunASR("baseline", LANGUAGE_ID=language_id, PACK=pack, ACOUSTIC_WEIGHT=config["ACOUSTIC_WEIGHT"])
            kws_score = env.RunKWS("baseline", asr_output, LANGUAGE_ID=language_id, PACK=pack)

            continue
            # triple-oracle experiment
            oracle_pronunciations, oracle_pnsp, oracle_tags = env.AppenToAttila([pjoin("work", "oracles", language, pack, x) for x in 
                                                                                 ["oracle_pronunciations.txt", "oracle_pnsp.txt", "oracle_tags.txt"]],
                                                                                ["${LANGUAGE_PACKS}/${LANGUAGE_ID}/conversational/reference_materials/lexicon.txt",
                                                                                 "${LANGUAGE_PACKS}/${LANGUAGE_ID}/scripted/reference_materials/lexicon.txt",
                                                                                 env.Value({"LOCALE" : locale, "SKIP_ROMAN" : config["SKIP_ROMAN"]})],
                                                                                LANGUAGE_ID=language_id)
            
            oracle_text, oracle_text_words = env.CollectText([pjoin("work", "oracles", language, pack, x) for x in ["oracle_text.txt", "oracle_text_words.txt"]], 
                                                             [env.Dir(x) for x in glob(env.subst("${LANGUAGE_PACKS}/%d/*/*/transcription" % language_id))],
                                                             )
            oracle_vocabulary = env.PronunciationsToVocabulary(pjoin("work", "oracles", language, pack, "oracle_vocabulary.txt"), oracle_pronunciations)

            markov = 3
            oracle_language_model = env.IBMTrainLanguageModel(pjoin("work", "oracles", language, pack, "oracle_lm.%dgm.arpabo.gz" % (markov)), 
                                                              [oracle_text, oracle_text_words, env.Value(markov)])

            (asr_output, asr_score) = env.RunASR("triple_oracle", LANGUAGE_ID=language_id, PACK=pack, ACOUSTIC_WEIGHT=config["ACOUSTIC_WEIGHT"],
                                                 VOCABULARY_FILE=oracle_vocabulary[0],
                                                 PRONUNCIATIONS_FILE=oracle_pronunciations,
                                                 LANGUAGE_MODEL_FILE=oracle_language_model[0],
                                                 PHONE_FILE=oracle_pnsp)
            kws_score = env.RunKWS("triple_oracle", asr_output, LANGUAGE_ID=language_id, PACK=pack)
            
            continue
            for expansion in [limited_basic_expansions[1]]:

                g2p = env.RunG2P(pjoin("work", "expansions", language, pack, "%s_g2p.gz" % (os.path.splitext(os.path.basename(expansion.rstr()))[0])), 
                                [expansion, pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "tools/g2p/model-6")])
                all_pronunciations = env.G2PToBabel(pjoin("work", "expansions", language, pack, "%s_pronunciations.gz" % (os.path.splitext(os.path.basename(expansion.rstr()))[0])), 
                                                   [g2p, env.Value(config["PHONEME_SWAP"])])

                for size in env["EXPANSION_SIZES"]:
                    name = pjoin("work", "asr", language, pack, "babelgum_%d" % size)                
                    words, pronunciations = env.TopWords([pjoin("work", "expansions", language, pack, name, "probability_list.gz"), 
                                                          pjoin("work", "expansions", language, pack, name, "pronunciations.gz")], 
                                                         [expansion, all_pronunciations, env.Value({"COUNT" : size})])

                    for weight in env["EXPANSION_WEIGHTS"]:
                        name = "babelgum_%d_%f" % (size, weight)
                        babelgum_vocabulary, babelgum_pronunciations, babelgum_language_model = env.AugmentLanguageModel(
                            [pjoin("work", "expansions", language, pack, name, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(limited_language_model_file.rstr())]],
                            [limited_pronunciations_file, limited_language_model_file, pronunciations, env.Value(weight)]
                            )

                        (asr_output, asr_score) = env.RunASR(name, LANGUAGE_ID=language_id, PACK=pack, ACOUSTIC_WEIGHT=config["ACOUSTIC_WEIGHT"],
                                                             VOCABULARY_FILE=babelgum_vocabulary,
                                                             PRONUNCIATIONS_FILE=babelgum_pronunciations,
                                                             LANGUAGE_MODEL_FILE=babelgum_language_model,
                                                             )

                        kws_score = env.RunKWS(name, asr_output, LANGUAGE_ID=language_id, PACK=pack)
