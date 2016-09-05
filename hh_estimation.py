import uncertainpy
import matplotlib.pyplot as plt

parameterlist = [["Cm", 1, None],
                 ["gbar_Na", 120, None],
                 ["gbar_K", 36, None],
                 ["gbar_l", 0.3, None],
                 ["E_Na", 115, None],
                 ["E_K", -12, None],
                 ["E_l", 10.613, None]]


parameters = uncertainpy.Parameters(parameterlist)
model = uncertainpy.HodkinHuxleyModel(parameters=parameters)


model.setAllDistributions(uncertainpy.Distribution(0.15).uniform)


exploration = uncertainpy.UncertaintyEstimation(model,
                                                feature_list="all",
                                                CPUs=7,
                                                output_dir_data="figures/",
                                                output_dir_figures="figures/",
                                                save_figures=True,
                                                figureformat=".pdf")

exploration.allParameters()
