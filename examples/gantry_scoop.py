import numpy as np
import georges
import georges.manzoni
import pyDOE
import random
from deap import base, creator, tools, algorithms
from scoop import futures


class Helpers:
    def __init__(self):
        self.i = 0
        # C://Users//AGYXV//Documents//Work//cgtr//beamlines//
        self.bl = georges.BeamlineBuilder(path="/Users/chernals/reps/IBA-Optics/beamlines", prefix='SBF112') \
            .add_from_file("ra") \
            .build() \
            .add_drifts(using_collimators=True, pipe_aperture=0.033)

        self.kinematics = georges.physics.kinematics(energy=230)
        self.b = georges.Beam(energy=self.kinematics['energy']).from_twiss_parameters(5e4,
                                                                            BETAX=1.14,
                                                                            ALPHAX=0.0159,
                                                                            EMITX=28.50e-6,
                                                                            BETAY=5.13,
                                                                            ALPHAY=2.59,
                                                                            EMITY=4.49e-6,
                                                                            DPPRMS=0.0,
                                                                            )

        # Predefined settings
        self.settings_100mm = {
            'Q1C': 93.5,
            'Q2C': 47.7,
            'Q3C': 25.0,
            'Q4C': 51.3,
            'Q5C': 46.4,
            'SLITS': 100 / 2 / 1000,
        }

        self.lower_bounds = [
            -20,
            4,
            1,
            -14,
            3,
        ]

        self.upper_bounds = [
            -12,
            12,
            9,
            -6,
            12,
        ]

        # For optimization
        self.variables = [
            ('Q1C', 'K1'),
            ('Q2C', 'K1'),
            ('Q3C', 'K1'),
            ('Q4C', 'K1'),
            ('Q5C', 'K1'),
        ]
        self.context = {
            'PARTICLE': 'PROTON',
            'PC': self.kinematics['momentum'],
            'BRHO': self.kinematics['brho'],
            'IPMQ': 17.5,
            'SLITS': 0.05,
        }
        self.manzoni_line = georges.manzoni.convert_line(self.bl.line, context=self.context, to_numpy=True)
        self.l = georges.manzoni.convert_line(self.bl.line, context=self.context, to_numpy=False)

        self.manzoni_variables = georges.manzoni.transform_variables(self.l, self.variables)

        self.manzoni_beam = self.b.distribution.values
        self.cube = self.scale(pyDOE.lhs(5, 4000, None), [-20, 4, 1, -14, 3], [-12, 12, 9, -6, 12])

    def scale(self, cube, xmin, xmax):
        for i in range(0, cube.shape[1]):
            cube[:, i] *= xmax[i] - xmin[i]
            cube[:, i] += xmin[i]
        return cube

    def initialize(self, cls):
        element = self.cube[self.i, :]
        self.i += 1
        return cls(element)

    def evaluate(self, vs):
        georges.manzoni.adjust_line(self.manzoni_line, self.manzoni_variables, np.array(vs))
        obs = georges.manzoni.ElementByElementObserver(
            element_by_element=lambda b: [1e3 * np.std(b[:, 0]),
                                          1e3 * np.std(b[:, 2]),
                                          b.shape[0],
                                          np.cov(b[:, 0:2].T)[0, 1] / np.sqrt(np.linalg.det(np.cov(b[:, 0], b[:, 1]))),
                                          np.cov(b[:, 2:4].T)[0, 1] / np.sqrt(np.linalg.det(np.cov(b[:, 2], b[:, 3]))),

                                          ]
        )
        o = georges.manzoni.manzoni.track(self.manzoni_line, self.manzoni_beam, observer=obs)
        f = [
            np.abs(o.data[:, -5][0][0]),
            np.abs(o.data[:, -5][0][1]),
            np.abs(o.data[:, -5][0][0] - o.data[:, -5][0][1]) / (o.data[:, -5][0][0] + o.data[:, -5][0][0]),
            1 - o.data[:, -5][0][2] / self.manzoni_beam.shape[0],
            np.abs(o.data[:, -5][0][3]),
            np.abs(o.data[:, -5][0][4])

        ]
        if 4 < f[0] < 5 and 4 < f[1] < 5:
            return f
        else:
            return [
                (f[0] - 4) ** 2 + 6,
                (f[1] - 4) ** 2 + 6,
                f[2] + 6,
                f[3] + 6,
                f[4] + 6,
                f[5] + 6,
            ]

    def algo(self, seed=None, NGEN=20, POP=400, CXPB=0.8):
        random.seed(seed)
        pfrt = tools.ParetoFront()

        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean, axis=0)
        stats.register("min", np.min, axis=0)
        logbook = tools.Logbook()
        logbook.header = "gen", "evals", "min", "avg"

        pop = toolbox.population(n=POP)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in pop if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # This is just to assign the crowding distance to the individuals
        # no actual selection is done
        pop = toolbox.select(pop, len(pop))

        record = stats.compile(pop)
        logbook.record(gen=0, evals=len(invalid_ind), **record)
        print(logbook.stream)
        param = [tuple(ind) for ind in pop if ind.fitness.values]
        fit = [ind.fitness.values for ind in pop if ind.fitness.values]
        fobj = [np.array(i[:]).sum() for i in fit]
        d_index = [(gen, i) for i in range(0, len(param))]
        index = pd.MultiIndex.from_tuples(d_index, names=['gen', 'pop'])
        a = pd.DataFrame(param, columns=['IQ1C', 'IQ2C', 'IQ3C', 'IQ4C', 'IQ5C'])
        b = pd.DataFrame(fit, columns=['SigmaX', 'SigmaY', 'Sym', 'Losses', 'alfaX', 'alfaY'])
        fobj = pd.DataFrame(fobj, columns=['f_obj'])
        c = a.join(b)
        cc = c.join(fobj)
        cc = cc.set_index(index)
        try:
            a = pd.read_csv('Test.csv', index_col=[0, 1])
            b = pd.concat([a, cc])
            b.to_csv('Test.csv')
        except FileNotFoundError:
            cc.to_csv('Test.csv')

        # Begin the generational process
        for gen in range(1, NGEN):
            # Vary the population
            offspring = tools.selTournamentDCD(pop, len(pop))
            offspring = [toolbox.clone(ind) for ind in offspring]

            for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
                if random.random() <= CXPB:
                    toolbox.mate(ind1, ind2)

                toolbox.mutate(ind1)
                toolbox.mutate(ind2)
                del ind1.fitness.values, ind2.fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            # Select the next generation population
            pop = toolbox.select(pop + offspring, POP)
            pfrt.update(pop)
            record = stats.compile(pop)
            logbook.record(gen=gen, evals=len(invalid_ind), **record)
            print(logbook.stream)
            param = [tuple(ind) for ind in pop if ind.fitness.values]
            fit = [ind.fitness.values for ind in pop if ind.fitness.values]
            fobj = [np.array(i[:]).sum() for i in fit]
            d_index = [(gen, i) for i in range(0, len(param))]
            index = pd.MultiIndex.from_tuples(d_index, names=['gen', 'pop'])
            a = pd.DataFrame(param, columns=['IQ1C', 'IQ2C', 'IQ3C', 'IQ4C', 'IQ5C'])
            b = pd.DataFrame(fit, columns=['SigmaX', 'SigmaY', 'Sym', 'Losses', 'alfaX', 'alfaY'])
            fobj = pd.DataFrame(fobj, columns=['f_obj'])
            c = a.join(b)
            cc = c.join(fobj)
            cc = cc.set_index(index)
            try:
                a = pd.read_csv('Test.csv', index_col=[0, 1])
                b = pd.concat([a, cc])
                b.to_csv('Test.csv')
            except FileNotFoundError:
                cc.to_csv('Test.csv')
        return pfrt, pop, logbook


h = Helpers()
creator.create("FitnessMin", base.Fitness, weights=(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0))
creator.create("Individual", list, fitness=creator.FitnessMin)
toolbox = base.Toolbox()
toolbox.register('map', futures.map)
toolbox.register("individual", h.initialize, creator.Individual)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", h.evaluate)
toolbox.register("mate", tools.cxSimulatedBinaryBounded, eta=2, low=h.lower_bounds, up=h.upper_bounds)
toolbox.register("mutate", tools.mutPolynomialBounded, eta=2, low=h.lower_bounds, up=h.upper_bounds,
                 indpb=1 / len(h.variables))
toolbox.register("select", tools.selNSGA2)


if __name__ == '__main__':
    h.algo(NGEN=150, POP=4000, CXPB=0.9)
