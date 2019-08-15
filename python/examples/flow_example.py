from simulators import simulators as opm
p = opm.FlowSimulator()
p.setDeckFilenameTo("../../../opm-data/norne/NORNE_ATW2013.DATA")
p.run()
