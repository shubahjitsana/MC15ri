#!/usr/bin/env python3
import sys
import basf2_mva

#input and output path
input_file = sys.argv[1]
output_weightfile = sys.argv[2]

general_options = basf2_mva.GeneralOptions()
# general_options.m_datafiles = basf2_mva.vector("scaled/train.root")
general_options.m_datafiles = basf2_mva.vector(f"{input_file}")
general_options.m_treename = "tree"
general_options.m_identifier = f"{output_weightfile}"  # outputted weightfile
general_options.m_variables = basf2_mva.vector(
    "abs_qr",       # calculated using ft.flavorTagger fn
    # "sphericity",   # calculated using ma.buildEventShape fn
    "DeltaZ",       # calculated using vx.TagV fn
    "R2",           # rest of the variables are calculated using ma.buildContinuumSuppression fn
    "thrustBm",  #must    #significance
    "thrustOm",  #must    #significance
    "cosTBTO",  #must    #significance
    "cosTBz",  #must    #significance
    # "KSFWVariables(et)",    #significance
    # "KSFWVariables(mm2)",
    # "KSFWVariables(hso00)",
    # "KSFWVariables(hso01)",
    "KSFWVariables(hso02)",  #must    #significance
    # "KSFWVariables(hso03)",
    # "KSFWVariables(hso04)",    #significance
    # "KSFWVariables(hso10)",    #significance
    "KSFWVariables(hso12)",  #must    #significance
    # "KSFWVariables(hso14)",    #significance
    # "KSFWVariables(hso20)",    #significance
    # "KSFWVariables(hso22)",
    # "KSFWVariables(hso24)",
    "KSFWVariables(hoo0)",  #must    #significance
    # "KSFWVariables(hoo1)",
    # "KSFWVariables(hoo2)",    #significance
    # "KSFWVariables(hoo3)",
    # "KSFWVariables(hoo4)",    #significance
    # "CleoConeCS(1)",
    # "CleoConeCS(2)",
    # "CleoConeCS(3)",
    # "CleoConeCS(4)",
    # "CleoConeCS(5)",
    # "CleoConeCS(6)",
    # "CleoConeCS(7)",
    # "CleoConeCS(8)",
    # "CleoConeCS(9)",
)
general_options.m_target_variable = "isSignal"
fastbdt_options = basf2_mva.FastBDTOptions()

basf2_mva.teacher(general_options, fastbdt_options)