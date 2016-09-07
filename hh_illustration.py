import uncertainpy
import matplotlib.pyplot as plt

plt.xkcd()

parameterlist = [["gbar_Na", 120, None],
                 ["gbar_K", 36, None],
                 ["gbar_l", 0.3, None]]


parameters = uncertainpy.Parameters(parameterlist)
model = uncertainpy.HodkinHuxleyModel(parameters=parameters)


model.run()
uncertainpy.prettyPlot(model.t, model.U, nr_hues=3, sns_style="white", linewidth=2)

parameterlist1 = {"gbar_Na": 0.3*120,
                  "gbar_K": 0.3*36,
                  "gbar_l": 0.3*0.3}

model.setParameterValues(parameterlist1)
model.run()
uncertainpy.prettyPlot(model.t, model.U, new_figure=False, nr_hues=3, sns_style="white", linewidth=2)



parameterlist2 = {"gbar_Na": 1.6*120,
                  "gbar_K": 1.6*36,
                  "gbar_l": 1.6*0.3}

model.setParameterValues(parameterlist2)
model.run()
uncertainpy.prettyPlot(model.t, model.U, new_figure=False, nr_hues=3, sns_style="white", linewidth=2)

plt.savefig("figures/hh_thin.pdf")
# plt.show()
