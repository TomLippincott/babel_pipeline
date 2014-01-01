from SCons.Builder import Builder
from SCons.Action import Action
from SCons.Subst import scons_subst
from SCons.Node import Node, NodeList
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
from babel import ProbabilityList, Arpabo, Pronunciations, Vocabulary, FrequencyList, MorfessorOutput, ASRResults, KWSResults
from os.path import join as pjoin
from common_tools import meta_open

def run_asr(target, source, env):
    
    return None

def run_asr_emitter(target, source, env):
    new_sources = []
    new_targets = pjoin(env["BASE_PATH"], "out.txt")
    return new_targets, new_sources

def run_kws(target, source, env):
    return None

def run_kws_emitter(target, source, env):
    return target, source

def build_extrinsic_tables(target, source, env):
    return None

def build_extrinsic_tables_emitter(target, source, env):
    new_sources = [files_to_strings(source[0].read())] + leaves(source[0].read())
    return target, new_sources

def build_latex(target, source, env):
    properties, results = [x.read() for x in source[0:2]]
    languages = set([x[0] for x in properties.keys()])
    lookup = {"PRE" : "Prefixes",
              "STM" : "Stems",
              "SUF" : "Suffixes",
              }
    packs = ["Limited"]
    language_table, morfessor_table, babelgum_table = {}, {}, {}


    for language in languages:
        language_properties = properties[(language, "Limited")]
        with meta_open(language_properties["prefixes"]) as prefix_fd, meta_open(language_properties["stems"]) as stem_fd, meta_open(language_properties["suffixes"]) as suffix_fd:
            pre, stm, suf = [
                [l.strip().split()[0] for l in prefix_fd if "<epsilon>" not in l],
                [l.strip().split()[0] for l in stem_fd if "<epsilon>" not in l],
                [l.strip().split()[0] for l in suffix_fd if "<epsilon>" not in l],
                ]
            babelgum_table[language] = [len(pre), "%.2f" % (sum(map(len, pre)) / max(1.0, float(len(pre)))), 
                                         len(stm), "%.2f" % (sum(map(len, stm)) / float(len(stm))),
                                         len(suf), "%.2f" % (sum(map(len, suf)) / float(len(suf))),
                                         ]            

        with meta_open(language_properties["limited_vocabulary"]) as lim_fd, meta_open(language_properties["dev_vocabulary"]) as dev_fd:
            lim_vocab = set(FrequencyList(lim_fd).make_conservative().keys())
            dev_vocab = set(FrequencyList(dev_fd).make_conservative().keys())
            lim_vocab_size = len(lim_vocab)
            dev_vocab_size = len(dev_vocab)
            both_vocabs = len(lim_vocab.union(dev_vocab))
            avg_len_lim_vocab = sum(map(len, lim_vocab)) / float(len(lim_vocab))
            avg_len_dev_vocab = sum(map(len, dev_vocab)) / float(len(dev_vocab))
            language_table[language] = [lim_vocab_size, "%.2f" % (avg_len_lim_vocab), 
                                        dev_vocab_size, "%.2f" % (avg_len_dev_vocab), 
                                        len([x for x in dev_vocab if x not in lim_vocab])]
        with meta_open(language_properties["morfessor_input"]) as input_fd, meta_open(language_properties["morfessor_output"]) as output_fd:
            input_vocab = FrequencyList({w : int(c) for c, w in [x.strip().split() for x in input_fd]})
            morf_output = MorfessorOutput(output_fd)
            pre = morf_output.morphs["PRE"]
            stm = morf_output.morphs["STM"]
            suf = morf_output.morphs["SUF"]
            morfessor_table[language] = [len(pre), "%.2f" % (sum(map(len, pre)) / max(1.0, float(len(pre)))), 
                                         len(stm), "%.2f" % (sum(map(len, stm)) / float(len(stm))),
                                         len(suf), "%.2f" % (sum(map(len, suf)) / float(len(suf))),
                                         ]
    with meta_open(target[0].rstr(), "w") as ofd:
        body = "\n".join([r"  %s & %s \\" % (l.title(), " & ".join(map(str, v))) for l, v in sorted(language_table.iteritems())])
        ofd.write(r"""
%%language properties
\begin{tabular}{|l|r|r|r|r|r|}
  \hline
  Language & \multicolumn{2}{|c|}{Training} & \multicolumn{2}{|c|}{Development} & OOV \\
  & Count & Avg. length & Count & Avg. length & \\
  \hline
%s
  \hline
\end{tabular}
""" % body)

        body = "\n".join([r"  %s & %s \\" % (l.title(), " & ".join(map(str, v))) for l, v in sorted(morfessor_table.iteritems())])
        ofd.write(r"""
%%morfessor properties
\begin{tabular}{|l|r|r|r|r|r|r|}
  \hline
  Language & \multicolumn{2}{|c|}{Prefixes} & \multicolumn{2}{|c|}{Stems} & \multicolumn{2}{|c|}{Suffixes} \\
  & Count & Avg. length & Count & Avg. length & Count & Avg. length \\
  \hline
%s
  \hline
\end{tabular}
""" % body)

        body = "\n".join([r"  %s & %s \\" % (l.title(), " & ".join(map(str, v))) for l, v in sorted(babelgum_table.iteritems())])
        ofd.write(r"""
%%babelgum properties
\begin{tabular}{|l|r|r|r|r|r|r|}
  \hline
  Language & \multicolumn{2}{|c|}{Prefixes} & \multicolumn{2}{|c|}{Stems} & \multicolumn{2}{|c|}{Suffixes} \\
  & Count & Avg. length & Count & Avg. length & Count & Avg. length \\
  \hline
%s
  \hline
\end{tabular}
""" % body)
        pass

    return None

def build_latex_emitter(target, source, env):
    properties,results = [x.read() for x in source[0:3]]
    new_sources = [env.Value(files_to_strings(properties)),
                   env.Value(files_to_strings(results)),
                   ] + sum(map(leaves, [properties, results]), [])
    return target, new_sources


def build_site(target, source, env):
    properties, figures, results = [x.read() for x in source[0:3]]
    languages = set([x[0] for x in figures.keys()])
    lookup = {"PRE" : "Prefixes",
              "STM" : "Stems",
              "SUF" : "Suffixes",
              }
    packs = ["Limited"]
    base_path = os.path.dirname(target[0].rstr())
    try:
        os.makedirs(pjoin(base_path, "images"))
    except:
        pass
    with meta_open(pjoin(base_path, "theme.css"), "w") as ofd:
        ofd.write("body {text-align : center; vertical-align : top;}\n")
        ofd.write("table {text-align : center; vertical-align : top;}\n")
        ofd.write("tr {text-align : center; vertical-align : top;}\n")
        ofd.write("td {text-align : center; vertical-align : top;}\n")
    with meta_open(target[0].rstr(), "w") as ofd:
        xml = et.TreeBuilder()
        xml.start("html", {})
        xml.start("head", {}), xml.start("link", {"rel" : "stylesheet", "type" : "text/css", "href" : "theme.css"}), xml.end("link"), xml.end("head")
        xml.start("body", {}), xml.start("table", {})
        for language in languages:
            language_properties = properties[(language, "Limited")]
            with meta_open(language_properties["prefixes"]) as prefix_fd, meta_open(language_properties["stems"]) as stem_fd, meta_open(language_properties["suffixes"]) as suffix_fd:
                babel_output = {"PRE" : [l.strip().split()[0] for l in prefix_fd if "<epsilon>" not in l],
                                "STM" : [l.strip().split()[0] for l in stem_fd if "<epsilon>" not in l],
                                "SUF" : [l.strip().split()[0] for l in suffix_fd if "<epsilon>" not in l],
                                }
            with meta_open(language_properties["limited_vocabulary"]) as lim_fd, meta_open(language_properties["dev_vocabulary"]) as dev_fd:
                lim_vocab = set(FrequencyList(lim_fd).make_conservative().keys())
                dev_vocab = set(FrequencyList(dev_fd).make_conservative().keys())
                lim_vocab_size = len(lim_vocab)
                dev_vocab_size = len(dev_vocab)
                both_vocabs = len(lim_vocab.union(dev_vocab))
                avg_len_lim_vocab = sum(map(len, lim_vocab)) / float(len(lim_vocab))
                avg_len_dev_vocab = sum(map(len, dev_vocab)) / float(len(dev_vocab))
            with meta_open(language_properties["morfessor_input"]) as input_fd, meta_open(language_properties["morfessor_output"]) as output_fd:
                input_vocab = FrequencyList({w : int(c) for c, w in [x.strip().split() for x in input_fd]})
                morf_output = MorfessorOutput(output_fd)
                morf_analysis_counts = len(morf_output)
                morf_morph_counts = {k : len(v) for k, v in morf_output.morphs.iteritems()}
                morf_morph_lengths = {k : sum(map(len, v)) / float(len(v)) for k, v in morf_output.morphs.iteritems() if len(v) > 0}
                

            xml.start("tr", {}), xml.start("td", {}), xml.start("h1", {}), xml.data(language.title()), xml.end("h1"), xml.end("td"), xml.end("tr")

            xml.start("table", {})
            xml.start("tr", {}), [(xml.start("td", {}), xml.start("h3", {}), xml.data("%s information" % (x)), xml.end("h3"), xml.end("td")) for x in ["Language", "Morfessor", "BabelGUM"]], xml.end("tr")
            
            xml.start("tr", {})

            # Language information
            xml.start("td", {})
            xml.start("table", {"border" : "1"})
            #xml.start("tr", {}), [(xml.start("td", {}), xml.data(x), xml.end("td")) for x in ["Pack", "Vocabulary size"]], xml.end("tr")
            xml.start("tr", {}), xml.start("td", {}), xml.end("td"), xml.start("td", {}), xml.data("Count"), xml.end("td"), xml.start("td", {}), xml.data("Average length"), xml.end("td"), xml.end("tr")
            xml.start("tr", {}), xml.start("td", {}), xml.data("Limited vocab"), xml.end("td"), xml.start("td", {}), xml.data("%d" % (lim_vocab_size)), xml.end("td"), xml.start("td", {}), xml.data("%.2f" % (avg_len_lim_vocab)), xml.end("td"), xml.end("tr")
            xml.start("tr", {}), xml.start("td", {}), xml.data("Dev vocab"), xml.end("td"), xml.start("td", {}), xml.data("%d" % (dev_vocab_size)), xml.end("td"), xml.start("td", {}), xml.data("%.2f" % (avg_len_dev_vocab)), xml.end("td"), xml.end("tr")
            
            #for name, values in [(lookup[x[0]], x[1]) for x in sorted(morf_output.morphs.iteritems())]:
            #    xml.start("tr", {})
                #xml.start("td", {}), xml.data(name), xml.end("td")
                #xml.start("td", {}), xml.data(str(len(values))), xml.end("td")
                #xml.start("td", {}), xml.data(avg_len), xml.end("td")
            #    xml.end("tr")
            xml.end("table")
            xml.end("td")
            
            # Morfessor information
            xml.start("td", {})
            xml.start("table", {"border" : "1"})
            xml.start("tr", {}), [(xml.start("td", {}), xml.data(x), xml.end("td")) for x in ["Type", "Count", "Average length"]], xml.end("tr")
            for name, values in [(lookup[x[0]], x[1]) for x in sorted(morf_output.morphs.iteritems())]:
                if len(values) > 0:
                    avg_len = "%.2f" % (sum(map(len, values)) / float(len(values)))
                else:
                    avg_len = ""
                xml.start("tr", {})
                xml.start("td", {}), xml.data(name), xml.end("td")
                xml.start("td", {}), xml.data(str(len(values))), xml.end("td")
                xml.start("td", {}), xml.data(avg_len), xml.end("td")
                xml.end("tr")
            xml.end("table")
            xml.end("td")
            
            # BabelGUM information
            xml.start("td", {})
            xml.start("table", {"border" : "1"})
            xml.start("tr", {}), [(xml.start("td", {}), xml.data(x), xml.end("td")) for x in ["Type", "Count", "Average length"]], xml.end("tr")
            for name, values in [(lookup[x[0]], x[1]) for x in sorted(babel_output.iteritems())]:
                if len(values) > 0:
                    avg_len = "%.2f" % (sum(map(len, values)) / float(len(values)))
                else:
                    avg_len = ""
                xml.start("tr", {})
                xml.start("td", {}), xml.data(name), xml.end("td")
                xml.start("td", {}), xml.data(str(len(values))), xml.end("td")
                xml.start("td", {}), xml.data(avg_len), xml.end("td")
                xml.end("tr")
            xml.end("table")
            xml.end("td")
            xml.end("tr")
            
            xml.end("table")
            
            # graphs of IV increase and OOV reduction, type-based and token-based
            xml.start("tr", {}), xml.start("td", {}), xml.start("h3", {}), xml.data("Intrinsic performance evaluation"), xml.end("h3"), xml.end("td"), xml.end("tr")
            xml.start("tr", {}), xml.start("td", {})
            xml.start("table", {})
            for pack in packs:
                image_file = "%s_%s.png" % (language, pack)
                shutil.copy(figures[(language, pack)], pjoin(base_path, "images", image_file))
                xml.start("tr", {}), xml.start("td", {}), xml.start("img", {"src" : pjoin("images", image_file)}), xml.end("img"), xml.end("td"), xml.end("tr")
            xml.end("table")
            xml.end("td"), xml.end("tr")
            
            # word error rate for ASR and maximum term-weighted value for KWS
            xml.start("tr", {}), xml.start("td", {}), xml.start("h3", {}), xml.data("Extrinsic performance evaluation"), xml.end("h3"), xml.end("td"), xml.end("tr")
            xml.start("tr", {}), xml.start("td", {}), xml.start("table", {"border" : "1"})
            xml.start("tr", {}), [(xml.start("td", {}), xml.data(x), xml.end("td")) for x in ["Augmentation", "Error", "Substitutions", "Deletions", "Insertions", "PMiss", "MTWV"]], xml.end("tr")
            for name, values in sorted(results[(language, "Limited")].iteritems()):
                with meta_open(values["ASR"]) as asr_fd, meta_open(values["KWS"]) as kws_fd:
                    asr = ASRResults(asr_fd)
                    kws = KWSResults(kws_fd)
                    xml.start("tr", {})
                    xml.start("td", {}), xml.data(name), xml.end("td")
                    [(xml.start("td", {}), xml.data("%.3f" % (asr.get(x))), xml.end("td")) for x in ["error", "substitutions", "deletions", "insertions"]]
                    [(xml.start("td", {}), xml.data("%.3f" % (kws.get(x))), xml.end("td")) for x in ["pmiss", "mtwv"]]
                    xml.end("tr")
            xml.end("table")
            xml.end("td"), xml.end("tr")
            
        xml.end("table"), xml.end("body")
        xml.end("html")
        ofd.write(et.tostring(xml.close()))
    return None

def files_to_strings(data):
    if isinstance(data, dict):
        return {k : files_to_strings(v) for k, v in data.iteritems()}
    elif isinstance(data, NodeList):
        return data[0].rstr()
    elif isinstance(data, Node):
        return data.rstr()
    else:
        raise Exception(type(data))

def strings_to_files(data):
    if isinstance(data, dict):
        return {k : files_to_strings(v) for k, v in data.iteritems()}
    elif not isinstance(data, basestring):
        return File(data)

def leaves(data):
    if isinstance(data, dict):
        return sum(map(leaves, data.values()), [])
    else:
        return [data]

def build_site_emitter(target, source, env):
    properties, figures, results = [x.read() for x in source[0:3]]
    new_targets = pjoin(env["BASE_PATH"], "index.html")
    new_sources = [env.Value(files_to_strings(properties)),
                   env.Value(files_to_strings(figures)),
                   env.Value(files_to_strings(results)),
                   ] + sum(map(leaves, [properties, figures, results]), [])    #[v[0].rstr() for v in sum([x.values() for x in properties.values()], []) + figures.values() + results.values()]
    return new_targets, new_sources

def TOOLS_ADD(env):
    env.Append(BUILDERS = {"RunASR" : Builder(action=run_asr, emitter=run_asr_emitter),
                           "RunKWS" : Builder(action=run_kws, emitter=run_kws_emitter),
                           "BuildSite" : Builder(action=build_site, emitter=build_site_emitter),
                           "BuildLatex" : Builder(action=build_latex, emitter=build_latex_emitter),
                           "BuildExtrinsicTables" : Builder(action=build_extrinsic_tables, emitter=build_extrinsic_tables_emitter),
                           })
