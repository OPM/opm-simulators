import opm.io
from opm.io.parser import Parser, ParseContext
from opm.io.ecl_state import EclipseState
from opm.io.schedule import Schedule
from opm.io.summary import SummaryConfig
from simulators import simulators

parse_context = ParseContext( [("PARSE_RANDOM_SLASH", opm.io.action.ignore)] )
deck = Parser().parse( '../../../opm-data/norne/NORNE_ATW2013.DATA' , parse_context)
state = EclipseState(deck)
schedule = Schedule(deck, state)
summary_config = SummaryConfig(deck, state, schedule)

p = simulators.BlackOilSimulator()
p.setDeck(deck)
p.setEclipseState(state)
p.setSchedule(schedule)
p.setSummaryConfig(summary_config)
p.run()
