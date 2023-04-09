import schemdraw
from schemdraw import Segment, SegmentCircle
import schemdraw.elements as elm

def figure_1():
    with schemdraw.Drawing() as d:
        source = d.add(NoiseSource().label('$V_s$', 'left'))
        d.add(elm.Resistor().at(source.output).label('$R_s$').up())
        d.add(elm.Coax().right().label('$Z_0$'))
        d.add(elm.Resistor().down().label('$R_l$'))
        d.add(elm.Wire('|-').to(source.input))
    return d


class NoiseSource(elm.Element):
    def __init__(self, *d, **kwargs):
        super().__init__(*d, **kwargs)
        self.segments.append(SegmentCircle((0, 0), 0.5))
        self.segments.append(
                Segment(
                    [
                        (-0.25, 0),
                        (-0.2, 0.24),
                        (-0.14, 0.1),
                        (-0.12, -0.2),
                        (-0.08, 0.05),
                        (-0.04, -0.2),
                        (0, 0.22),
                        (0.05, 0.18),
                        (0.09, -0.22),
                        (0.15, -0.04),
                        (0.2, -0.21),
                        (0.25, 0.1)
                    ],
                ),
        )
        self.segments.append(Segment([(0, 0.5), (0, 1)]))
        self.segments.append(Segment([(0, -0.5), (0, -1)]))
        self.anchors = {
                'input': (0, -1),
                'output': (0, 1),
                }
