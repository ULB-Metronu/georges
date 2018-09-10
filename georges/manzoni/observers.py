import numpy as np

_ = lambda x: x.copy()


class Observer:

    def __init__(self, turn_by_turn_active=False, element_by_element_active=False, **kwargs):
        self._start_data = None
        self._end_data = None
        self._data = []
        self._element_by_element_is_active = element_by_element_active
        self._turn_by_turn_is_active = turn_by_turn_active
        self._track_start = kwargs.get('track_start', _)
        self._track_end = kwargs.get('track_end', _)
        self._turn_by_turn = kwargs.get('turn_by_turn', _)
        self._element_by_element = kwargs.get('element_by_element', _)

    @property
    def start_data(self):
        return self._start_data

    @property
    def end_data(self):
        return self._end_data

    @property
    def data(self):
        return self._data

    @property
    def element_by_element_is_active(self):
        return self._element_by_element_is_active

    @property
    def turn_by_turn_is_active(self):
        return self._turn_by_turn_is_active

    def track_start(self, beam):
        if self._track_start is not None:
            self._start_data = self._track_start(beam)

    def element_by_element(self, turn, element, beam):
        pass

    def turn_by_turn(self, turn, element, beam):
        pass

    def track_end(self, turn, element, beam):
        self._end_data = self._track_end(beam)
        return self


class TurnByTurnObserver(Observer):

    def __init__(self, **kwargs):
        super().__init__(turn_by_turn_active=True, element_by_element_active=False, **kwargs)

    def track_start(self, beam):
        super().track_start(beam)

    def track_end(self, turn, element, beam):
        super().track_end(turn, element, beam)
        self._data = np.array(self.data)
        return self

    def turn_by_turn(self, turn, element, beam):
        self._data.append(self._turn_by_turn(beam))


class ElementByElementObserver(Observer):

    def __init__(self, elements=None, **kwargs):
        super().__init__(turn_by_turn_active=False, element_by_element_active=True, **kwargs)
        self._elements = elements

    @property
    def elements(self):
        return self._elements

    def track_start(self, beam):
        super().track_start(beam)

    def track_end(self, turn, element, beam):
        if self.elements is not None:
            n_elements = len(self.elements)-1
        else:
            n_elements = element
        super().track_end(turn, n_elements, beam)
        self._data = np.array(self.data).reshape(turn+1, n_elements+1, -1, np.shape(beam)[1])
        return self

    def element_by_element(self, turn, element, beam):
        self._data.append(self._element_by_element(beam))
