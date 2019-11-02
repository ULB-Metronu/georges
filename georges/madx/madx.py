from typing import Optional
import cpymad.madx
from georges_core.sequences import Sequence as _Sequence
from georges_core import Kinematics as _Kinematics


class MadX(cpymad.madx.Madx):
    def __init__(self,
                 sequence: Optional[_Sequence] = None,
                 kinematics: Optional[_Kinematics] = None,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.send_sequence(sequence, kinematics)

    def send_sequence(self,
                      sequence: Optional[_Sequence] = None,
                      kinematics: Optional[_Kinematics] = None):
        if sequence is None:
            return
        if kinematics is None:
            kinematics = sequence.kinematics
        df = sequence.to_df(strip_units=True)
        if kinematics is not None:
            self.input(f"""
        BEAM, PARTICLE={kinematics.particule.__name__}, ENERGY={kinematics.etot.m_as('GeV')};
    """.strip())
        self.input(f"{sequence.name or 'SEQ'}: SEQUENCE, L={df.iloc[-1]['AT_EXIT']}, REFER=ENTRY;")
        for name, element in df.iterrows() or []:
            parameters = dict(element[set(list(element.index.values)).intersection(
                set(list(map(lambda _: _.upper(),
                             self._libmadx.get_defined_command(element['CLASS'].lower())['data'].keys())))
            )])
            self.input((f"{name}: {element['CLASS'].lower()}, AT={element['AT_ENTRY']}, " + ', '.join([f"{k}={str(v).strip('([])')}" for k, v in parameters.items()]) + ';').strip())
        self.input("ENDSEQUENCE;")
        self.input(f"USE, SEQUENCE={sequence.name or 'SEQ'};")

    def track(self):
        pass
