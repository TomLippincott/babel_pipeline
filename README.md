babel_pipeline
==============

Language packs

It is assumed that each language is in the directory ${LANGUAGE_PACKS}/${LANGUAGE_ID}, with basic subdirectory structure:

  conversational/
    dev/
      audio/
      transcription/
    training/
      audio/
      transcription/
    reference_materials/
  scripted/
    dev/
      audio/
      transcription/
    training/
      audio/
      transcription/
    reference_materials/

There may be additional directories, such as "sub-train", "transcript_roman", etc, but these are not relevant at the moment.