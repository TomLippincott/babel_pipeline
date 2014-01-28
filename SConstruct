import sys
import os.path
from os.path import join as pjoin
import logging
from functools import partial
from glob import glob
import asr_tools
import kws_tools
import torque_tools
import morfessor_tools
import babel_tools
import re


#
# load variable definitions from custom.py, and define them for SCons (seems like it should logically
# happen in the reverse order, but anyways...)
#
vars = Variables("custom.py")
vars.AddVariables(
    ("OUTPUT_WIDTH", "", 130),

    ("RUN_ASR", "", False),
    ("RUN_KWS", "", False),
    ("EVALUATE_PRONUNCIATIONS", "", False),

    ("EXPANSION_SIZES", "", []),
    ("EXPANSION_WEIGHTS", "", []),

    ("LANGUAGE_PACKS", "", None),
    ("IBM_MODELS", "", None),
    ("LORELEI_SVN", "", None),
    ("ATTILA_PATH", "", None),
    ("VOCABULARY_EXPANSION_PATH", "", None),
    ("INDUSDB_PATH", "", None),

    ("SEQUITUR_PATH", "", ""),
    ("ATTILA_INTERPRETER", "", "${ATTILA_PATH}/tools/attila/attila"),
    ("OUTPUT_PATH", "", ""),
    ("ADD_WORDS", "", "/usr/bin/add_words"),
    ("BABEL_REPO", "", None),
    ("BABEL_RESOURCES", "", None),
    ("F4DE", "", None),
    ("JAVA_NORM", "", "${BABEL_REPO}/KWS/examples/babel-dryrun/javabin"),
    ("OVERLAY", "", None),
    ("LIBRARY_OVERLAY", "", "${OVERLAY}/lib:${OVERLAY}/lib64:${LORELEI_TOOLS}/boost_1_49_0/stage/lib/"),
    ("EXPERIMENTS", "", {}),
    ("LANGUAGES", "", {}),
    ("LOG_LEVEL", "", logging.INFO),
    ("LOG_DESTINATION", "", sys.stdout),
    ("HAS_TORQUE", "", False),
    ("PYTHON_INTERPRETER", "", None),
    ("SCORE_SCRIPT", "", None),
    ("INDUS_DB", "", None),
    ("SCLITE_BINARY", "", None),
    ("JOBS", "", 4),
    ("LOCAL_JOBS", "", 4),
    ("BABEL_REPO", "", None),
    ("BABEL_RESOURCES", "", None),
    ("LORELEI_TOOLS", "", None),
    ("G2P", "", None),
 
    ("OVERLAY", "", None),

    # these variables all have default definitions in terms of the previous, but may be overridden as needed
    ("PYTHON", "", "/usr/bin/python"),
    ("PERL", "", "/usr/bin/perl"),
    ("PERL_LIBRARIES", "", os.environ.get("PERL5LIB", "")),
    ("BABEL_BIN_PATH", "", "${BABEL_REPO}/tools/kws/bin64"),
    ("BABEL_SCRIPT_PATH", "", "${BABEL_REPO}/tools/kws/bin64"),
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
                  tools=["default", "textfile"] + [x.TOOLS_ADD for x in [asr_tools, kws_tools, torque_tools, morfessor_tools, babel_tools]],
                  BUILDERS={"CopyFile" : Builder(action="cp ${SOURCE} ${TARGET}")}
                  )

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
    exp_id = config["BIG_ID"]

    # full_transcripts, limited_transcripts, dev_transcripts = env.SplitTrainDev(
    #     [pjoin("work", "transcripts", language, "%s_transcripts.txt.gz" % (name)) for name in ["full", "limited", "development"]],
    #     env.Dir(pjoin(env["LANGUAGE_PACKS"], str(language_id))))

    # full_vocabulary, limited_vocabulary, dev_vocabulary = [env.TranscriptsToVocabulary(pjoin("work", "vocabularies", language, "%s.txt.gz" % (name)), transcripts) for name, transcripts in zip(["full", "limited", "development"], [full_transcripts, limited_transcripts, dev_transcripts])]
    
    limited_basic_expansions_file = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "expansions", "simple.ex"))
    # limited_bigram_expansions_file = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "expansions", "bigram.ex"))

    #limited_basic_expansions = env.SplitExpansion([limited_basic_expansions_file, env.Value(100000)],
    #                                              BASE_PATH=pjoin("work", "expansions", language, "limited", "basic"))

    # limited_bigram_expansions = env.SplitExpansion([limited_bigram_expansions_file, env.Value(100000)],
    #                                                BASE_PATH=pjoin("work", "expansions", language, "limited", "bigram"))
    
    # morfessor_input = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "train.sorted"))
    # morfessor_output = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "morf.out"))
    # prefixes = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "prefix.lex"))
    # stems = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "stem.lex"))
    # suffixes = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "suffix.lex"))
    
    # figures[(language, "Limited")] = env.PlotReduction(limited_basic_expansions + limited_bigram_expansions + [limited_vocabulary, dev_vocabulary, env.Value({"bins" : 1000})])

    # evaluate G2P performance
    # if os.path.exists(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "tools/g2p/model-6")) and env["EVALUATE_PRONUNCIATIONS"]:
    #     limited_pronunciations_file = env.File(pjoin(env["LORELEI_SVN"], str(language_id), "LimitedLP", "models", "dict.test"))
    #     g2p = env.RunG2P(pjoin("work", "pronunciations", language, "g2p_candidate_pronunciations.gz"),
    #                      [limited_pronunciations_file, pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "tools/g2p/model-6")])

    #     all_pronunciations = env.G2PToBabel(pjoin("work", "pronunciations", language, "candidate_pronunciations.gz"), # % (os.path.splitext(os.path.basename(expansion.rstr()))[0])), 
    #                                         [g2p, env.Value(config["PHONEME_SWAP"])])

    #     g2p_evals += [language, env.PronunciationPerformance(pjoin("work", "pronunciations", language, "performance.txt"), [all_pronunciations, limited_pronunciations_file])]                      



    #if os.path.exists(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "tools/g2p/model-6")):    


    # properties[(language, "Limited")] = {"prefixes" : prefixes,
    #                                      "stems" : stems,
    #                                      "suffixes" : suffixes,
    #                                      "limited_vocabulary" : limited_vocabulary,
    #                                      "dev_vocabulary" : dev_vocabulary,
    #                                      "morfessor_output" : morfessor_output,
    #                                      "morfessor_input" : morfessor_input,                                         
    #                                      #"babelgum_pronunciations" : None,
    #                                      }

    # #continue
    results[(language, "Limited")] = {}

    if os.path.exists(pjoin(env["IBM_MODELS"], str(language_id))):

        #full_pronunciations_file = env.File(pjoin(env["LORELEI_SVN"], str(language_id), "FullLP", "models", "dict.test"))
        #limited_pronunciations_file = env.File(pjoin(env["LORELEI_SVN"], str(language_id), "LimitedLP", "models", "dict.test"))

        exp_path = pjoin("work", "pronunciations", language, "limited")
        #g2p = env.RunG2P(pjoin(exp_path, "basic.gz"), # % (os.path.splitext(os.path.basename(expansion.rstr()))[0])), 
        #                 [limited_basic_expansions_file, pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "tools/g2p/model-6")])
        #env.NoClean(g2p)
        #all_pronunciations = env.G2PToBabel(pjoin(exp_path, "basic_prons.gz"), # % (os.path.splitext(os.path.basename(expansion.rstr()))[0])), 
        #                                    [g2p, env.Value(config["PHONEME_SWAP"])])



        limited_pronunciations_file = env.File(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "dict.test")) #config["pronunciations"]
        limited_vocabulary_file = env.File(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "vocab"))
        limited_language_model_file = env.Glob(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "lm.*"))[0]



        for pack in ["LLP"]:
            (asr_output, asr_score) = env.RunASR("baseline", LANGUAGE_ID=language_id, PACK=pack, ACOUSTIC_WEIGHT=config["ACOUSTIC_WEIGHT"])
            kws_score = env.RunKWS("baseline", asr_output, LANGUAGE_ID=language_id, PACK=pack)


            #PRONUNCIATIONS_FILE=limited_pronunciations_file.rstr(),
            #                             VOCABULARY_FILE=limited_vocabulary_file.rstr(),
            #                             LANGUAGE_MODEL_FILE=limited_language_model_file.rstr())
            #baseline_output = env.RunKWS(pjoin(experiment_path, "kws"), [])

        # baseline experiment
        #(baseline_asr_score, baseline_asr_output) = language_pack_asr_run(OUTPUT_PATH=pjoin("work", language, pack, "baseline"),
        #                                                                  VOCABULARY_FILE=limited_vocabulary_file, 
        #                                                                  PRONUNCIATIONS_FILE=limited_pronunciations_file,
        #                                                                  LANGUAGE_MODEL_FILE=limited_language_model_file, 
        #                                                                  )
        
        # baseline experiment
        #baseline_kws_output = language_pack_kws_run(baseline_asr_output, OUTPUT_PATH=pjoin("work", language, pack, "baseline", "kws"),
        #                                            VOCABULARY_FILE=limited_vocabulary_file, 
        #                                            PRONUNCIATIONS_FILE=limited_pronunciations_file,
        #                                            LANGUAGE_MODEL_FILE=limited_language_model_file, 
        #                                            )

        #baseline_oov_kws_output = language_pack_kws_run(baseline_asr_output, OUTPUT_PATH=pjoin("work", language, pack, "baseline", "oov_kws"),
        #                                                VOCABULARY_FILE=limited_vocabulary_file, 
        #                                                PRONUNCIATIONS_FILE=limited_pronunciations_file,
        #                                                LANGUAGE_MODEL_FILE=limited_language_model_file, 
        #                                                OOV_ONLY=True,
        #                                                )

#         #results[(language, "Limited")]["Baseline"] = baseline_results
#         #{"ASR" : env.File(pjoin("work", language, pack, "baseline", "asr", "scoring", "babel.sys")), 
#         #                                              "KWS" : env.File(pjoin("work", language, pack, "baseline", "kws", "output", "Ensemble.AllOccur.results.txt")),
#         #                                              }

#         # triple-oracle experiment (weighted)
#         exp_path = pjoin("work", language, pack, "triple_oracle")
#         oracle_pronunciations, oracle_pnsp, oracle_tags = env.AppenToAttila([pjoin(exp_path, x) for x in 
#                                                                              ["oracle_pronunciations.txt", "oracle_pnsp.txt", "oracle_tags.txt"]],
#                                                                             [pjoin(data_path.rstr(), "conversational", "reference_materials", "lexicon.txt"),
#                                                                              pjoin(data_path.rstr(), "scripted", "reference_materials", "lexicon.txt"),
#                                                                              env.Value({"LOCALE" : locale, "SKIP_ROMAN" : config["SKIP_ROMAN"]})])

#         oracle_text, oracle_text_words = env.CollectText([pjoin(exp_path, x) for x in ["oracle_text.txt", "oracle_text_words.txt"]], 
#                                                          [env.Dir(x) for x in glob(pjoin(data_path.rstr(), "*/*/transcription"))])
#         oracle_vocabulary = env.PronunciationsToVocabulary(pjoin(exp_path, "oracle_vocabulary.txt"), oracle_pronunciations)
#         oracle_language_model = env.IBMTrainLanguageModel(pjoin(exp_path, "oracle_lm.%dgm.arpabo.gz" % (markov)), 
#                                                           [oracle_text, oracle_text_words, env.Value(markov)])

#         # oracle_results = language_pack_run(OUTPUT_PATH=exp_path,
#         #                                    VOCABULARY_FILE=oracle_vocabulary[0],
#         #                                    PRONUNCIATIONS_FILE=oracle_pronunciations,
#         #                                    LANGUAGE_MODEL_FILE=oracle_language_model[0],
#         #                                    PHONE_FILE=oracle_pnsp,
#         #                                    )
        
#         #results[(language, "Limited")]["Oracle"] = oracle_results

#         # babelgum experiments
        continue
        for expansion in [limited_basic_expansions[1]]:
            continue
            exp_path = pjoin("work", language, pack)                
            g2p = env.RunG2P(pjoin(exp_path, "%s_g2p.gz" % (os.path.splitext(os.path.basename(expansion.rstr()))[0])), 
                            [expansion, pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "tools/g2p/model-6")])
            all_pronunciations = env.G2PToBabel(pjoin(exp_path, "%s_pronunciations.gz" % (os.path.splitext(os.path.basename(expansion.rstr()))[0])), 
                                               [g2p, env.Value(config["PHONEME_SWAP"])])

            for size in env["EXPANSION_SIZES"]:
                exp_path = pjoin("work", language, pack, "babelgum_%d" % size)                
                words, pronunciations = env.TopWords([pjoin(exp_path, "probability_list.gz"), pjoin(exp_path, "pronunciations.gz")], 
                                                     [expansion, all_pronunciations, env.Value({"COUNT" : size})])

                for weight in env["EXPANSION_WEIGHTS"]:
                    exp_path = pjoin("work", language, pack, "babelgum_%d_%f" % (size, weight))
                    babelgum_vocabulary, babelgum_pronunciations, babelgum_language_model = env.AugmentLanguageModel(
                        [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(limited_language_model_file.rstr())]],
                        [limited_pronunciations_file, limited_language_model_file, pronunciations, env.Value(weight)]
                        )
                
                    (babelgum_asr_score, babelgum_asr_output) = language_pack_asr_run(OUTPUT_PATH=exp_path,
                                                                                      VOCABULARY_FILE=babelgum_vocabulary,
                                                                                      PRONUNCIATIONS_FILE=babelgum_pronunciations,
                                                                                      LANGUAGE_MODEL_FILE=babelgum_language_model,
                                                                                      )

                    babelgum_kws_output = language_pack_kws_run(babelgum_asr_output, OUTPUT_PATH=pjoin(exp_path, "kws"),
                                                                VOCABULARY_FILE=babelgum_vocabulary,
                                                                PRONUNCIATIONS_FILE=babelgum_pronunciations,
                                                                LANGUAGE_MODEL_FILE=babelgum_language_model, 
                                                                )

                    # babelgum_kws_output = language_pack_kws_run(babelgum_asr_output, OUTPUT_PATH=pjoin("work", language, pack, "babelgum", "kws"),
                    #                                             VOCABULARY_FILE=babelgum_vocabulary,
                    #                                             PRONUNCIATIONS_FILE=babelgum_pronunciations,
                    #                                             LANGUAGE_MODEL_FILE=babelgum_language_model,
                    #                                             REMOVE_VOCABULARY_FILE=limited_vocabulary_file,
                    #                                             )

                    #babelgum_oov_kws_output = language_pack_kws_run(babelgum_asr_output, OUTPUT_PATH=pjoin("work", language, pack, "babelgum", "oov_kws"),
                    #                                                VOCABULARY_FILE=limited_vocabulary_file, 
                    #                                                PRONUNCIATIONS_FILE=limited_pronunciations_file,
                    #                                                LANGUAGE_MODEL_FILE=limited_language_model_file, 
                    #                                                OOV_ONLY=True,
                    #                                                )


#                     ###
#                     # exp_path = pjoin("work", language, pack, "babelgum_%d_%f_correctprons" % (size, weight))
#                     # vocabulary, pronunciations = env.ReplacePronunciations([pjoin(exp_path, x) for x in ["vocabulary.gz", "corrected_pronunciations.gz"]],
#                     #                                            [pronunciations, oracle_pronunciations])

#                     # babelgum_vocabulary, babelgum_pronunciations, babelgum_language_model = env.AugmentLanguageModel(
#                     #     [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(limited_language_model_file.rstr())]],
#                     #     [limited_pronunciations_file, limited_language_model_file, pronunciations, env.Value(weight)]
#                     #     )

#                     # babelgum_results = language_pack_run(OUTPUT_PATH=exp_path,
#                     #                                      VOCABULARY_FILE=babelgum_vocabulary,
#                     #                                      PRONUNCIATIONS_FILE=babelgum_pronunciations,
#                     #                                      LANGUAGE_MODEL_FILE=babelgum_language_model,
#                     #                                      )

#                     # exp_path = pjoin("work", language, pack, "weighted_babelgum_%d_%f" % (size, weight))
#                     # weighted_babelgum_vocabulary, weighted_babelgum_pronunciations, weighted_babelgum_language_model = env.AugmentLanguageModel(
#                     #     [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(limited_language_model_file.rstr())]],
#                     #     [limited_pronunciations_file, limited_language_model_file, pronunciations, words, env.Value(.1)]
#                     #     )

#                     # weighted_babelgum_results = language_pack_run(OUTPUT_PATH=exp_path,
#                     #                                               VOCABULARY_FILE=weighted_babelgum_vocabulary,
#                     #                                               PRONUNCIATIONS_FILE=weighted_babelgum_pronunciations,
#                     #                                               LANGUAGE_MODEL_FILE=weighted_babelgum_language_model,
#                     #                                               )


#         # for method, (probs, prons) in config["expansions"].iteritems():
#         #     for size in [50000]:
#         #         base_path = pjoin("work", language, pack, method, str(size))
#         #         babelgum_probabilities, babelgum_pronunciations = env.BabelGumLexicon([pjoin(base_path, x) for x in ["expanded_%s_%d_probabilities.txt" % (method, size), 
#         #                                                                                                              "expanded_%s_%d_pronunciations.txt" % (method, size)]], 
#         #                                                                               [probs, prons, env.Value(size)])
#         #         #env.PlotProbabilities(pjoin(base_path, "probabilities.png"), probs)

#         #         for weight in [.1]:

#         #             # double-oracle experiment (uniform)
#         #             exp_path = pjoin("work", language, pack, "double_oracle")
#         #             babelgum_vocabulary, babelgum_pronunciations, babelgum_language_model = env.AugmentLanguageModel(
#         #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(config["language_model"])]],
#         #                 [oracle_pronunciations, config["language_model"], oracle_pronunciations, env.Value(weight)]
#         #                 )



#         #             # adding correct babelgum vocabulary (uniform)
#         #             exp_path = pjoin(base_path, "correct_uniform=%f" % weight)
#         #             babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities = \
#         #                 env.FilterBabelGum([pjoin(exp_path, "babelgum_rightwords_%d_%s.txt" % (size, x)) for x in ["pronunciations", "probabilities"]],
#         #                                    [babelgum_pronunciations, babelgum_probabilities, oracle_vocabulary])

#         #             babelgum_rightwords_uniform_vocabulary, babelgum_rightwords_uniform_pronunciations, babelgum_rightwords_uniform_language_model = env.AugmentLanguageModel(
#         #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(baseline_language_model)]],
#         #                 [baseline_pronunciations, baseline_language_model, babelgum_rightwords_pronunciations, env.Value(weight)]
#         #                 )



#         #             experiments.append(language_pack_run(OUTPUT_PATH=exp_path,
#         #                                                  VOCABULARY_FILE=babelgum_rightwords_uniform_vocabulary,
#         #                                                  PRONUNCIATIONS_FILE=babelgum_rightwords_uniform_pronunciations,
#         #                                                  LANGUAGE_MODEL_FILE=babelgum_rightwords_uniform_language_model,
#         #                                                  ))


#         #             # adding correct babelgum vocabulary (weighted)
#         #             exp_path = pjoin(base_path, "correct_weighted=%f" % weight)

#         #             babelgum_rightwords_weighted_vocabulary, babelgum_rightwords_weighted_pronunciations, babelgum_rightwords_weighted_language_model = env.AugmentLanguageModel(
#         #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(baseline_language_model)]],
#         #                 [baseline_pronunciations, baseline_language_model, babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities, env.Value(weight)]
#         #                 )

#         #             experiments.append(language_pack_run(OUTPUT_PATH=exp_path,
#         #                                                  VOCABULARY_FILE=babelgum_rightwords_weighted_vocabulary,
#         #                                                  PRONUNCIATIONS_FILE=babelgum_rightwords_weighted_pronunciations,
#         #                                                  LANGUAGE_MODEL_FILE=babelgum_rightwords_weighted_language_model,
#         #                                                  ))

# env.Textfile(pjoin("work", "pronunciations", "evaluations.txt"), g2p_evals, SUBST_DICT={})
# #env.BuildSite(target=[], 
# #              source=[env.Value(properties), env.Value(figures), env.Value({k : {kk : {"ASR" : vv[0], "KWS" : vv[1]} for kk, vv in v.iteritems()} for k, v in results.iteritems()})], BASE_PATH="work/babel_site")

# #env.BuildPropertyTables(target=["work/%s_property_table.tex" % (x) for x in ["language", "morfessor", "babelgum"]], 
# #                        source=[env.Value(properties)]) #, env.Value({k : {kk : {"ASR" : vv[0], "KWS" : vv[1]} for kk, vv in v.iteritems()} for k, v in results.iteritems()})])

#env.BuildExtrinsicTables("work/extrinsic_table.tex", source=[env.Value(results)])
# #env.Value({k : {kk : {"ASR" : vv[0], "KWS" : vv[1]} for kk, vv in v.iteritems()} for k, v in results.iteritems()}))

# #                     babelgum_rightwords_uniform_path = pjoin("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum_corrected", "babelgum", "uniform")
# #                     babelgum_rightwords_uniform_experiment = env.CreateSmallASRDirectory(Dir(babelgum_rightwords_uniform_path),
# #                                                                                          experiment({"VOCABULARY_FILE" : babelgum_rightwords_uniform_vocabulary.rstr(),
# #                                                                                                      "PRONUNCIATIONS_FILE" : babelgum_rightwords_uniform_pronunciations.rstr(),
# #                                                                                                      "LANGUAGE_MODEL_FILE" : babelgum_rightwords_uniform_languagemodel.rstr(),
# #                                                                                                      "OUTPUT_PATH" : babelgum_rightwords_uniform_path,
# #                                                                                                      })
# #                                                                                          )
# #                     scores.append(submit(babelgum_rightwords_uniform_path, babelgum_rightwords_uniform_experiment, jobs=env["JOBS"]))


# #                     # adding just correct babelgum terms (weighted)
# #                     babelgum_rightwords_weighted_vocabulary, babelgum_rightwords_weighted_pronunciations, babelgum_rightwords_weighted_languagemodel = env.AugmentLanguageModel(
# #                         [pjoin(base_path, "babelgum_rightwords_weighted_%s_%d_%f_%s" % (model, size, weight, x)) for x in ["vocabulary.txt", "pronunciations.txt", "lm.3gm.arpabo.gz"]],
# #                         [baseline_pronunciations, baseline_languagemodel, babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities, env.Value(weight)]
# #                         )

# #                     babelgum_rightwords_weighted_path = pjoin("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum_corrected", "babelgum", "weighted")
# #                     babelgum_rightwords_weighted_experiment = env.CreateSmallASRDirectory(Dir(babelgum_rightwords_weighted_path),
# #                                                                                          experiment({"VOCABULARY_FILE" : babelgum_rightwords_weighted_vocabulary.rstr(),
# #                                                                                                      "PRONUNCIATIONS_FILE" : babelgum_rightwords_weighted_pronunciations.rstr(),
# #                                                                                                      "LANGUAGE_MODEL_FILE" : babelgum_rightwords_weighted_languagemodel.rstr(),
# #                                                                                                      "OUTPUT_PATH" : babelgum_rightwords_weighted_path,
# #                                                                                                      })
# #                                                                                          )
# #                     scores.append(submit(babelgum_rightwords_weighted_path, babelgum_rightwords_weighted_experiment, jobs=env["JOBS"]))

# #env.CollateResults("work/results.txt", experiments)
