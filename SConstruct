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
    ("JOBS", "", 20),
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
    args["ASR_OUTPUT_PATH"] = pjoin(args["OUTPUT_PATH"], "asr")
    args["KWS_OUTPUT_PATH"] = pjoin(args["OUTPUT_PATH"], "kws")
    for k, v in args.iteritems():
        if k.endswith("FILE"):
            files[k] = v
        elif k.endswith("PATH"):
            directories[k] = v
        else:
            parameters[k] = v

    #
    # ASR 
    #
    experiment = env.CreateASRExperiment(Dir(args["ASR_OUTPUT_PATH"]), [env.Value(x) for x in [files, directories, parameters]])
    
    #construct = env.Construct()
    
    #test = env.Test()
    
    #if env["HAS_TORQUE"]:
    # construct = env.SubmitJob(pjoin(args["ASR_OUTPUT_PATH"], "construct.timestamp"), 
    #                           [env.Value({"name" : "construct",
    #                                       "commands" : ["${ATTILA_PATH}/tools/attila/attila construct.py"],
    #                                       "path" : args["ASR_OUTPUT_PATH"],
    #                                       "other" : ["#PBS -W group_list=yeticcls"],
    #                                       "interval" : 10,
    #                                       })])

    # test = env.SubmitJob(pjoin(args["ASR_OUTPUT_PATH"], "test.timestamp"), 
    #                      [construct, env.Value({"name" : "test",
    #                                             "commands" : ["${ATTILA_PATH}/tools/attila/attila test.py -w %f -n %s -j $${PBS_ARRAYID} -l 1" % (args["ACOUSTIC_WEIGHT"], 
    #                                                                                                                                               jobs)],
    #                                             "path" : args["ASR_OUTPUT_PATH"],
    #                                             "array" : jobs,
    #                                             "other" : ["#PBS -W group_list=yeticcls"],
    #                                             "interval" : 120,
    #                                             })])
    if not env["HAS_TORQUE"]:
        return (env.File(pjoin(args["ASR_OUTPUT_PATH"], "scoring", "babel.sys")),
                env.File(pjoin(args["KWS_OUTPUT_PATH"], "output", "Full-Occur-MITLLFA3-AppenWordSeg.sum.txt")),
                )


    asr_output = env.RunASRExperiment(source=experiment, ACOUSTIC_WEIGHT=args["ACOUSTIC_WEIGHT"])
    #return asr_output

    asr_score = env.ScoreResults(env.Dir(pjoin(args["ASR_OUTPUT_PATH"], "scoring")),
                                 #asr_output)
                                 [env.Dir(os.path.abspath(pjoin(args["ASR_OUTPUT_PATH"], "ctm"))), files["STM_FILE"]])

    #env.Depends(asr_score, test)
    #return asr_score
    #return asr_score
    #else:
    #asr_score = env.File(pjoin(args["ASR_OUTPUT_PATH"], "scoring", "babel.sys"))
    #return experiment

    #
    # KEYWORD SEARCH
    #

    # just make some local variables from the experiment definition (for convenience)
    iv_dict = args["VOCABULARY_FILE"]
    oov_dict = args["OOV_DICTIONARY_FILE"]
    dbfile = args["DATABASE_FILE"]
    kw_file = args["KEYWORDS_FILE"]
    language_id = args["LANGUAGE_ID"]

    iv_query_terms, oov_query_terms, term_map, word_to_word_fst, kw_file = env.QueryFiles([pjoin(args["KWS_OUTPUT_PATH"], x) for x in ["iv_queries.txt", 
                                                                                                                                       "oov_queries.txt",
                                                                                                                                       "term_map.txt",
                                                                                                                                       "word_to_word.fst",
                                                                                                                                       "kwfile.xml"]], 
                                                                                          [kw_file, iv_dict, env.Value(language_id)])

    # JOBS LATTICE_DIRECTORY KW_FILE RTTM_FILE
    base_path = args["OUTPUT_PATH"]
    args["LATTICE_DIRECTORY"] = pjoin(args["ASR_OUTPUT_PATH"], "lat")
    args["JOBS"] = 4
    args["KW_FILE"] = args["KEYWORDS_FILE"]

    
    full_lattice_list = env.LatticeList(pjoin(args["KWS_OUTPUT_PATH"], "lattice_list.txt"),
                                        [dbfile, env.Value(args["LATTICE_DIRECTORY"])])

    lattice_lists = env.SplitList([pjoin(args["KWS_OUTPUT_PATH"], "lattice_list_%d.txt" % (n + 1)) for n in range(args["JOBS"])], full_lattice_list)

    wordpron = env.WordPronounceSymTable(pjoin(args["KWS_OUTPUT_PATH"], "in_vocabulary_symbol_table.txt"),
                                         iv_dict)

    isym = env.CleanPronounceSymTable(pjoin(args["KWS_OUTPUT_PATH"], "cleaned_in_vocabulary_symbol_table.txt"),
                                      wordpron)

    mdb = env.MungeDatabase(pjoin(args["KWS_OUTPUT_PATH"], "munged_database.txt"),
                            [dbfile, full_lattice_list])

    padfst = env.BuildPadFST(pjoin(args["KWS_OUTPUT_PATH"], "pad_fst.txt"),
                             wordpron)

    env.Depends([wordpron, iv_query_terms, full_lattice_list], asr_score)

    full_data_list = env.CreateDataList(pjoin(args["KWS_OUTPUT_PATH"], "full_data_list.txt"),
                                        [mdb] + [env.Value({"oldext" : "fsm.gz", 
                                                            "ext" : "fst",
                                                            "subdir_style" : "hub4",
                                                            "LATTICE_DIR" : args["LATTICE_DIRECTORY"],
                                                            })], BASE_PATH=args["KWS_OUTPUT_PATH"])

    ecf_file = env.ECFFile(pjoin(args["KWS_OUTPUT_PATH"], "ecf.xml"), mdb)

    data_lists = env.SplitList([pjoin(args["KWS_OUTPUT_PATH"], "data_list_%d.txt" % (n + 1)) for n in range(args["JOBS"])], full_data_list)
    
    p2p_fst = env.FSTCompile(pjoin(args["KWS_OUTPUT_PATH"], "p2p_fst.txt"),
                             [isym, word_to_word_fst])

    wtp_lattices = []


    for i, (data_list, lattice_list) in enumerate(zip(data_lists, lattice_lists)):
        wp = env.WordToPhoneLattice(pjoin(args["KWS_OUTPUT_PATH"], "lattices", "lattice_generation-%d.stamp" % (i + 1)), 
                                    [data_list, lattice_list, wordpron, iv_dict, env.Value({"PRUNE_THRESHOLD" : -1,
                                                                                            "EPSILON_SYMBOLS" : "'<s>,</s>,~SIL,<HES>'",
                                                                                            })])

        fl = env.GetFileList(pjoin(args["KWS_OUTPUT_PATH"], "file_list-%d.txt" % (i + 1)), 
                             [data_list, wp])
        idx = env.BuildIndex(pjoin(args["KWS_OUTPUT_PATH"], "index-%d.fst" % (i + 1)),
                             fl)

        wtp_lattices.append((wp, data_list, lattice_list, fl, idx))

    merged = {}
    for query_type, query_file in zip(["in_vocabulary", "out_of_vocabulary"], [iv_query_terms, oov_query_terms]):
        queries = env.QueryToPhoneFST(pjoin(args["KWS_OUTPUT_PATH"], query_type, "query.fst"), 
                                      [p2p_fst, isym, iv_dict, query_file, env.Value({"n" : 1, "I" : 1, "OUTDIR" : pjoin(args["KWS_OUTPUT_PATH"], query_type, "queries")})])
        searches = []
        for i, (wtp_lattice, data_list, lattice_list, fl, idx) in enumerate(wtp_lattices):
            searches.append(env.StandardSearch(pjoin(args["KWS_OUTPUT_PATH"], query_type, "search_output-%d.txt" % (i + 1)),
                                               [data_list, isym, idx, padfst, queries, env.Value({"PRECISION" : "'%.4d'", "TITLE" : "std.xml", "LANGUAGE_ID" : language_id})]))



        qtl, res_list, res, ures = env.Merge([pjoin(args["KWS_OUTPUT_PATH"], query_type, x) for x in ["ids_to_query_terms.txt", "result_file_list.txt", "search_results.xml", "unique_search_results.xml"]], 
                                             [query_file] + searches + [env.Value({"MODE" : "merge-default",
                                                                                   "PADLENGTH" : 4,                                    
                                                                                   "LANGUAGE_ID" : language_id})])

        merged[query_type] = ures
        om = env.MergeScores(pjoin(args["KWS_OUTPUT_PATH"], query_type, "results.xml"), 
                             res)

    iv_oov = env.MergeIVOOV(pjoin(args["KWS_OUTPUT_PATH"], "iv_oov_results.xml"), 
                            [merged["in_vocabulary"], merged["out_of_vocabulary"], term_map, args["KW_FILE"]])

    norm = env.Normalize(pjoin(args["KWS_OUTPUT_PATH"], "norm.kwslist.xml"), 
                         [iv_oov, kw_file])

    normSTO = env.NormalizeSTO(pjoin(args["KWS_OUTPUT_PATH"], "normSTO.kwslist.xml"), 
                               norm)

    kws_score = env.Score(pjoin(args["KWS_OUTPUT_PATH"], "output", "Ensemble.AllOccur.results.txt"), 
                          [normSTO, kw_file, env.Value({"RTTM_FILE" : args["RTTM_FILE"], "ECF_FILE" : ecf_file[0].rstr(), "EXPID" : args["EXPID"]})])

    return (asr_score, kws_score)


#
# morph, lm rerank, lm rerank w averaging, lm reranking for morpheme boundaries
#
experiments = []
general_run = partial(run_experiment, SAMPLING_RATE=8000, FEATURE_TYPE="plp", AC_WEIGHT=.13, MAX_ERROR=15000, USE_DISPATCHER=False)
properties = {}
figures = {}
results = {}
for (language, language_id, expid), config in env["LANGUAGES"].iteritems():
    # prefix, stem, suffix, vocabulary, avg. character length, avg. morph length

    language_id = config["LANGUAGE_ID"]
    locale = config["LOCALE"]

    full_transcripts, limited_transcripts, dev_transcripts = env.SplitTrainDev(
        [pjoin("work", "transcripts", language, "%s_transcripts.txt.gz" % (name)) for name in ["full", "limited", "development"]],
        env.Dir(pjoin(env["LANGUAGE_PACKS"], str(language_id))))

    full_vocabulary, limited_vocabulary, dev_vocabulary = [env.TranscriptsToVocabulary(pjoin("work", "vocabularies", language, "%s_vocabularies.txt.gz" % (name)), transcripts) for name, transcripts in zip(["full", "limited", "development"], [full_transcripts, limited_transcripts, dev_transcripts])]
    
    limited_basic_expansions = env.SplitExpansion([pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "expansions", "simple.ex"), env.Value(100000)],
                                                  BASE_PATH=pjoin("work", "expansions", language, "limited", "basic"))

    limited_bigram_expansions = env.SplitExpansion([pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "expansions", "bigram.ex"), env.Value(100000)],
                                                   BASE_PATH=pjoin("work", "expansions", language, "limited", "bigram"))
    
    morfessor_input = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "train.sorted"))
    morfessor_output = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "morf.out"))
    prefixes = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "prefix.lex"))
    stems = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "stem.lex"))
    suffixes = env.File(pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "data", "suffix.lex"))


    figures[(language, "Limited")] = env.PlotReduction(limited_basic_expansions + limited_bigram_expansions + [limited_vocabulary, dev_vocabulary, env.Value({"bins" : 1000})])
    properties[(language, "Limited")] = {"prefixes" : prefixes,
                                         "stems" : stems,
                                         "suffixes" : suffixes,
                                         "limited_vocabulary" : limited_vocabulary,
                                         "dev_vocabulary" : dev_vocabulary,
                                         "morfessor_output" : morfessor_output,
                                         "morfessor_input" : morfessor_input,
                                         }

    results[(language, "Limited")] = {}
    #continue
    if os.path.exists(pjoin(env["IBM_MODELS"], str(language_id))):

        #full_pronunciations_file = env.File(pjoin(env["LORELEI_SVN"], str(language_id), "FullLP", "models", "dict.test"))
        limited_pronunciations_file = env.File(pjoin(env["LORELEI_SVN"], str(language_id), "LimitedLP", "models", "dict.test"))

        limited_vocabulary_file = env.File(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "vocab"))
        limited_language_model_file = env.Glob(pjoin(env["IBM_MODELS"], str(language_id), "LLP", "models", "lm.*"))[0]

        markov = int(re.match(r".*lm\.(\d+)gm.*", limited_language_model_file.rstr()).group(1))

        for pack in ["LLP"]:
            baseline_pronunciations = env.File(pjoin(env["IBM_MODELS"], str(language_id), "models", "dict.test")) #config["pronunciations"]
            baseline_language_model = env.Glob(pjoin(env["IBM_MODELS"], str(language_id), "models", "lm*")) #config["language_model"]

            data_path = env.Dir(pjoin(env["LANGUAGE_PACKS"], str(language_id)))        
            segmentation_file = env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "segment", "*dev.*"))[0]

            model_path = env.Dir(pjoin(env["IBM_MODELS"], str(language_id), pack, "models"))
            mel_file, phone_file, phone_set_file, tags_file, priors_file, tree_file, topo_file, topo_tree_file, lda_file = \
                [env.File(pjoin(model_path.rstr(), x)) for x in ["mel", "pnsp", "phonesset", "tags", "priors", "tree", "topo.tied", "topotree", "30.mat"]]

            adapt_path = env.Dir(pjoin(env["IBM_MODELS"], str(language_id), pack, "adapt"))
            warp_file = env.File(pjoin(adapt_path.rstr(), "warp.lst"))
            stm_file = env.Glob(pjoin(env["INDUSDB_PATH"], "babel%d*" % (language_id), "*stm"))[0]
            rttm_file = env.Glob(pjoin(env["INDUSDB_PATH"], "babel%d*" % (language_id), "*dev.rttm"))[0]
            keywords = env.Glob(pjoin(env["INDUSDB_PATH"], "babel%d*conv-dev.kwlist.xml" % (language_id)))[0]
            oov_dict = env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "kws-resources", "kws-resources*", "dict.OOV.v2p"))[0]
            language_pack_run = partial(general_run,
                                        jobs=env["JOBS"],
                                        PCM_PATH=data_path,
                                        CMS_PATH=pjoin(env["IBM_MODELS"], str(language_id), pack, "adapt", "cms"),
                                        FMLLR_PATH=pjoin(env["IBM_MODELS"], str(language_id), pack, "adapt", "fmllr"),
                                        MODEL_PATH=model_path,
                                        MEL_FILE=mel_file,
                                        PHONE_FILE=phone_file,
                                        PHONE_SET_FILE=phone_set_file,
                                        TAGS_FILE=tags_file,
                                        PRIORS_FILE=priors_file,
                                        TREE_FILE=tree_file,
                                        TOPO_FILE=topo_file,
                                        TOPO_TREE_FILE=topo_tree_file,
                                        LDA_FILE=lda_file,
                                        WARP_FILE=warp_file,
                                        LANGUAGE_ID=str(language_id),
                                        EXPID=expid,
                                        DATABASE_FILE=segmentation_file,
                                        ACOUSTIC_WEIGHT=config["ACOUSTIC_WEIGHT"],
                                        RTTM_FILE=rttm_file,
                                        STM_FILE=stm_file,
                                        KEYWORDS_FILE=keywords,
                                        OOV_DICTIONARY_FILE=oov_dict,
                                        TRFS_FILE=env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "models", "*.trfs")),
                                        TR_FILE=env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "models", "*.tr")),
                                        CTX_FILE=env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "models", "*.ctx")),
                                        GS_FILE=env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "models", "*.gs")),
                                        MS_FILE=env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "models", "*.ms")),
                                        FS_FILE=env.Glob(pjoin(env["IBM_MODELS"], str(language_id), pack, "models", "*.fs")),
                                        TXT_PATH=pjoin(env["IBM_MODELS"], str(language_id), pack, "SI", "cons"),
                                        )
        
        # baseline experiment
        baseline_results = language_pack_run(OUTPUT_PATH=pjoin("work", language, pack, "baseline"),
                                             VOCABULARY_FILE=limited_vocabulary_file, 
                                             PRONUNCIATIONS_FILE=limited_pronunciations_file,
                                             LANGUAGE_MODEL_FILE=limited_language_model_file, 
                                             )

        results[(language, "Limited")]["Baseline"] = baseline_results
        #{"ASR" : env.File(pjoin("work", language, pack, "baseline", "asr", "scoring", "babel.sys")), 
        #                                              "KWS" : env.File(pjoin("work", language, pack, "baseline", "kws", "output", "Ensemble.AllOccur.results.txt")),
        #                                              }

        # triple-oracle experiment (weighted)
        exp_path = pjoin("work", language, pack, "triple_oracle")
        oracle_pronunciations, oracle_pnsp, oracle_tags = env.AppenToAttila([pjoin(exp_path, x) for x in 
                                                                             ["oracle_pronunciations.txt", "oracle_pnsp.txt", "oracle_tags.txt"]],
                                                                            [pjoin(data_path.rstr(), "conversational", "reference_materials", "lexicon.txt"),
                                                                             pjoin(data_path.rstr(), "scripted", "reference_materials", "lexicon.txt"),
                                                                             env.Value({"LOCALE" : locale, "SKIP_ROMAN" : config["SKIP_ROMAN"]})])

        oracle_text, oracle_text_words = env.CollectText([pjoin(exp_path, x) for x in ["oracle_text.txt", "oracle_text_words.txt"]], 
                                                         [env.Dir(x) for x in glob(pjoin(data_path.rstr(), "*/*/transcription"))])
        oracle_vocabulary = env.PronunciationsToVocabulary(pjoin(exp_path, "oracle_vocabulary.txt"), oracle_pronunciations)
        oracle_language_model = env.IBMTrainLanguageModel(pjoin(exp_path, "oracle_lm.%dgm.arpabo.gz" % (markov)), 
                                                          [oracle_text, oracle_text_words, env.Value(markov)])

        oracle_results = language_pack_run(OUTPUT_PATH=exp_path,
                                           VOCABULARY_FILE=oracle_vocabulary[0],
                                           PRONUNCIATIONS_FILE=oracle_pronunciations,
                                           LANGUAGE_MODEL_FILE=oracle_language_model[0],
                                           PHONE_FILE=oracle_pnsp,
                                           )
        
        results[(language, "Limited")]["Oracle"] = oracle_results

        # babelgum experiments
        for size in [50000]:
            for expansion in [limited_basic_expansions[1]]:
                exp_path = pjoin("work", language, pack, "babelgum_%d" % size)                
                words = env.TopWords(pjoin(exp_path, "probability_list.gz"), [expansion, env.Value({"COUNT" : size})])
                g2p = env.RunG2P(pjoin(exp_path, "g2p.gz"), [words, pjoin(env["VOCABULARY_EXPANSION_PATH"], "%s-subtrain" % language, "tools/g2p/model-6")])
                pronunciations = env.G2PToBabel(pjoin(exp_path, "pronunciations.gz"), g2p)
                babelgum_vocabulary, babelgum_pronunciations, babelgum_language_model = env.AugmentLanguageModel(
                    [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(limited_language_model_file.rstr())]],
                    [limited_pronunciations_file, limited_language_model_file, pronunciations, env.Value(.1)]
                    )
                
                babelgum_results = language_pack_run(OUTPUT_PATH=exp_path,
                                                     VOCABULARY_FILE=babelgum_vocabulary,
                                                     PRONUNCIATIONS_FILE=babelgum_pronunciations,
                                                     LANGUAGE_MODEL_FILE=babelgum_language_model,
                                                     )

        # for method, (probs, prons) in config["expansions"].iteritems():
        #     for size in [50000]:
        #         base_path = pjoin("work", language, pack, method, str(size))
        #         babelgum_probabilities, babelgum_pronunciations = env.BabelGumLexicon([pjoin(base_path, x) for x in ["expanded_%s_%d_probabilities.txt" % (method, size), 
        #                                                                                                              "expanded_%s_%d_pronunciations.txt" % (method, size)]], 
        #                                                                               [probs, prons, env.Value(size)])
        #         #env.PlotProbabilities(pjoin(base_path, "probabilities.png"), probs)

        #         for weight in [.1]:

        #             # double-oracle experiment (uniform)
        #             exp_path = pjoin("work", language, pack, "double_oracle")
        #             babelgum_vocabulary, babelgum_pronunciations, babelgum_language_model = env.AugmentLanguageModel(
        #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(config["language_model"])]],
        #                 [oracle_pronunciations, config["language_model"], oracle_pronunciations, env.Value(weight)]
        #                 )



        #             # adding all babelgum vocabulary (uniform)
        #             exp_path = pjoin(base_path, "uniform=%f" % weight)
        #             babelgum_vocabulary, babelgum_pronunciations, babelgum_language_model = env.AugmentLanguageModel(
        #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(config["language_model"])]],
        #                 [config["pronunciations"], config["language_model"], babelgum_pronunciations, env.Value(weight)]
        #                 )

        #             experiments.append(language_pack_run(OUTPUT_PATH=exp_path,
        #                                                  VOCABULARY_FILE=babelgum_vocabulary,
        #                                                  PRONUNCIATIONS_FILE=babelgum_pronunciations,
        #                                                  LANGUAGE_MODEL_FILE=babelgum_language_model,
        #                                                  ))



        #             # adding all babelgum vocabulary (weighted)
        #             exp_path = pjoin(base_path, "weighted=%f" % weight)
        #             babelgum_weighted_vocabulary, babelgum_weighted_pronunciations, babelgum_weighted_language_model = env.AugmentLanguageModel(
        #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(baseline_language_model)]],
        #                 [baseline_pronunciations, baseline_language_model, babelgum_pronunciations, babelgum_probabilities, env.Value(weight)]
        #                 )

        #             experiments.append(language_pack_run(OUTPUT_PATH=exp_path,
        #                                                  VOCABULARY_FILE=babelgum_weighted_vocabulary,
        #                                                  PRONUNCIATIONS_FILE=babelgum_weighted_pronunciations,
        #                                                  LANGUAGE_MODEL_FILE=babelgum_weighted_language_model,
        #                                                  ))

        #             # adding correct babelgum vocabulary (uniform)
        #             exp_path = pjoin(base_path, "correct_uniform=%f" % weight)
        #             babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities = \
        #                 env.FilterBabelGum([pjoin(exp_path, "babelgum_rightwords_%d_%s.txt" % (size, x)) for x in ["pronunciations", "probabilities"]],
        #                                    [babelgum_pronunciations, babelgum_probabilities, oracle_vocabulary])

        #             babelgum_rightwords_uniform_vocabulary, babelgum_rightwords_uniform_pronunciations, babelgum_rightwords_uniform_language_model = env.AugmentLanguageModel(
        #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(baseline_language_model)]],
        #                 [baseline_pronunciations, baseline_language_model, babelgum_rightwords_pronunciations, env.Value(weight)]
        #                 )



        #             experiments.append(language_pack_run(OUTPUT_PATH=exp_path,
        #                                                  VOCABULARY_FILE=babelgum_rightwords_uniform_vocabulary,
        #                                                  PRONUNCIATIONS_FILE=babelgum_rightwords_uniform_pronunciations,
        #                                                  LANGUAGE_MODEL_FILE=babelgum_rightwords_uniform_language_model,
        #                                                  ))


        #             # adding correct babelgum vocabulary (weighted)
        #             exp_path = pjoin(base_path, "correct_weighted=%f" % weight)

        #             babelgum_rightwords_weighted_vocabulary, babelgum_rightwords_weighted_pronunciations, babelgum_rightwords_weighted_language_model = env.AugmentLanguageModel(
        #                 [pjoin(exp_path, x) for x in ["vocabulary.txt", "pronunciations.txt", os.path.basename(baseline_language_model)]],
        #                 [baseline_pronunciations, baseline_language_model, babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities, env.Value(weight)]
        #                 )

        #             experiments.append(language_pack_run(OUTPUT_PATH=exp_path,
        #                                                  VOCABULARY_FILE=babelgum_rightwords_weighted_vocabulary,
        #                                                  PRONUNCIATIONS_FILE=babelgum_rightwords_weighted_pronunciations,
        #                                                  LANGUAGE_MODEL_FILE=babelgum_rightwords_weighted_language_model,
        #                                                  ))


env.BuildSite(target=[], source=[env.Value(properties), env.Value(figures), env.Value({k : {kk : {"ASR" : vv[0], "KWS" : vv[1]} for kk, vv in v.iteritems()} for k, v in results.iteritems()})], BASE_PATH="work/babel_site")

#                     babelgum_rightwords_uniform_path = pjoin("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum_corrected", "babelgum", "uniform")
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
#                         [pjoin(base_path, "babelgum_rightwords_weighted_%s_%d_%f_%s" % (model, size, weight, x)) for x in ["vocabulary.txt", "pronunciations.txt", "lm.3gm.arpabo.gz"]],
#                         [baseline_pronunciations, baseline_languagemodel, babelgum_rightwords_pronunciations, babelgum_rightwords_probabilities, env.Value(weight)]
#                         )

#                     babelgum_rightwords_weighted_path = pjoin("work", "experiments_%d_%f" % (size, weight), language, pack, "babelgum_corrected", "babelgum", "weighted")
#                     babelgum_rightwords_weighted_experiment = env.CreateSmallASRDirectory(Dir(babelgum_rightwords_weighted_path),
#                                                                                          experiment({"VOCABULARY_FILE" : babelgum_rightwords_weighted_vocabulary.rstr(),
#                                                                                                      "PRONUNCIATIONS_FILE" : babelgum_rightwords_weighted_pronunciations.rstr(),
#                                                                                                      "LANGUAGE_MODEL_FILE" : babelgum_rightwords_weighted_languagemodel.rstr(),
#                                                                                                      "OUTPUT_PATH" : babelgum_rightwords_weighted_path,
#                                                                                                      })
#                                                                                          )
#                     scores.append(submit(babelgum_rightwords_weighted_path, babelgum_rightwords_weighted_experiment, jobs=env["JOBS"]))

#env.CollateResults("work/results.txt", experiments)
