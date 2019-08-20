import sunbeam
from simulators import simulators

p = simulators.BlackOilSimulator()
p.setDeckFilename('../../../opm-data/norne/NORNE_ATW2013.DATA')

norne = sunbeam.parse('../../../opm-data/norne/NORNE_ATW2013.DATA', recovery=[('PARSE_RANDOM_SLASH', sunbeam.action.ignore)])
p.setEclipseState(norne._state())
p.setDeck(norne._deck())
p.setSchedule(norne._schedule())
p.setSummaryConfig(norne._summary_config())
p.run()