import unittest
import matplotlib.pyplot as plt
import numpy as np
import os

from posydon.popsyn.binarypopulation import BinaryPopulation
from posydon.binary_evol import SimulationProperties
from posydon.binary_evol.flow_chart import flow_chart
from posydon.binary_evol.step_end import step_end


class TestBinaryPopulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        np.random.seed(12345)
        cls.POP_KWARGS = {
            "number_of_binaries": int(500),
            "primary_mass_min": 15
        }

        class MyEndStep(step_end):
            def __call__(self, binary):
                step_end.__call__(binary)
                binary.star_1.mass = np.sqrt(binary.star_1.mass)
                # change star_1 mass

        cls.SIM_PROP = SimulationProperties(
                flow = (flow_chart, {}),
                step_HMS_HMS = (MyEndStep, {}),
                step_CO_HeMS = (MyEndStep, {}),
                step_CO_HMS_RLO = (MyEndStep, {}),
                step_detached = (MyEndStep, {}),
                step_CE = (MyEndStep, {}),
                step_SN = (MyEndStep, {}),
                step_end = (MyEndStep, {}),
                )

    def test_init_0(self):
        bin_pop = BinaryPopulation()
        self.assertTrue(isinstance(bin_pop, BinaryPopulation))
        self.assertTrue(hasattr(bin_pop, 'population_properties'))
        self.assertTrue(hasattr(bin_pop, 'entropy'))


    def test_generate(self):
        bin_pop = BinaryPopulation(**self.POP_KWARGS)
        for i in range(bin_pop.number_of_binaries):
            bin_pop.manager.generate(**self.POP_KWARGS)
        self.assertTrue(len(bin_pop) == self.POP_KWARGS["number_of_binaries"])

        self.assertTrue( [b.star_1.mass > self.POP_KWARGS["primary_mass_min"]
                        for b in bin_pop] )

    # def test_generate_initial_binaries(self):
    #     bin_pop = BinaryPopulation(generate_initial_population=False, **self.POP_KWARGS)
    #     bin_pop.generate_initial_binaries()
    #     first_bin = bin_pop[1]
    #     bin_pop.generate_initial_binaries(overwrite=True)
    #     second_bin = bin_pop[1]
    #     self.assertFalse(
    #         first_bin is second_bin, msg="Binaries should not be the same object."
    #     )
    #
    # def test_gen_init_bin_err(self):
    #     bin_pop = BinaryPopulation()
    #     with self.assertRaisesRegex(
    #         ValueError, "set overwrite=True to overwrite existing population"
    #     ):
    #         bin_pop.generate_initial_binaries()
    #
    def test_sim_properties(self):
        bin_pop = BinaryPopulation(population_properties=self.SIM_PROP)
        self.assertTrue(isinstance(bin_pop, BinaryPopulation))
        self.assertTrue(bin_pop.population_properties is self.SIM_PROP)
        bin_pop.population_properties.load_steps()
        return bin_pop

    # def test_evolve(self):
    #     bin_pop = self.test_sim_properties()
    #     test_ids = np.arange(0, 15, 1)
    #     original_bins = bin_pop.copy(ids=test_ids)
    #     bin_pop.evolve()
    #     for b in bin_pop:
    #         with self.subTest("Check event END", binary_ind=b.index):
    #             self.assertTrue(b.event == "END")
    #
    #     for j, b in enumerate(bin_pop[test_ids]):
    #         with self.subTest(
    #             "Check mass changed",
    #             binary_ind=b.index,
    #             test_ind=original_bins[j].index,
    #         ):
    #             self.assertAlmostEqual(
    #                 b.star_1.mass, np.sqrt(original_bins[j].star_1.mass), places=8
    #             )

    # def test_evolve_binary_population(self):
    #     # It is unclear how to test multiprocessing at the moment
    #     POP_KWARGS = {"number_of_binaries": int(500), "primary_mass_min": 15}
    #
    #     def end(binary):
    #         binary.star_1.mass = np.sqrt(binary.star_1.mass)
    #         binary.event = "END"
    #
    #     def get_sim_prop():
    #         SIM_PROP = SimulationProperties(
    #             flow={("H-rich_Core_H_burning", "H-rich_Core_H_burning",
    #             "detached", "ZAMS"): "step_end"}, step_end=end, max_simulation_time=13.7e9)
    #         return SIM_PROP
    #
    #     bin_pop = BinaryPopulation(population_properties=get_sim_prop, **POP_KWARGS)
    #     bin_pop.evolve_binary_population(num_batches=4, verbose=True, use_df=True)

    # def test_evolve_each_binary(self):
    #     bin_pop = self.test_sim_properties()
    #     for num, evolved_bin in enumerate(bin_pop.evolve_each_binary()):
    #         with self.subTest("Evolve generator", num=num):
    #             self.assertTrue(num == evolved_bin.index)
    #             self.assertTrue(evolved_bin.event == "END")

    # def test_copy(self):
    #     bin_pop = self.test_sim_properties()
    #     binary_copy = bin_pop.copy(ids=0)
    #     self.assertFalse(binary_copy is bin_pop[0])
    #     all_binaries_copy = bin_pop.copy()
    #     self.assertFalse(
    #         any([copy_b is b for copy_b, b in zip(all_binaries_copy, bin_pop)])
    #     )

    # TODO: step_times is breaking to_df with only initialized binary / pop
    # def test_to_df(self):
    #     bin_pop = self.test_sim_properties()
    #     self.assertTrue( isinstance(bin_pop.to_df(), pd.DataFrame) )

    # def test_get_bin_by_index(self):
    #     bin_pop = self.test_sim_properties()
    #     test_indicies = [1, 6, 9, 8, 2]
    #     out_bins = bin_pop.get_binaries_by_index(test_indicies)
    #     self.assertTrue([b.index for b in out_bins] == test_indicies)

    # def test_bool_and_len(self):
    #     bin_pop = BinaryPopulation(population_properties=self.SIM_PROP)
    #     self.assertTrue(bool(bin_pop), msg="True if len self > 0")
    #     self.assertTrue(
    #         len(bin_pop) == bin_pop.number_of_binaries, msg="Should be len __binaries"
    #     )

    # def test_get_subpopulation(self):
    #     bin_pop = BinaryPopulation(population_properties=self.SIM_PROP, **self.POP_KWARGS)
    #     for i in range(200, 300):
    #         bin_pop[i].star_2.state = "BH"
    #     subpop = bin_pop.get_subpopulation(star_1_states=None, star_2_states="BH")
    #     self.assertTrue(all([bin.index == 200 + j for j, bin in enumerate(subpop)]))
    #     self.assertTrue(all([bin.star_2.state == "BH" for bin in subpop]))

    # def test_pickle_and_load(self):
    #     bin_pop = BinaryPopulation()
    #     bin_pop.pickle("saved_population.pkl")
    #     self.assertTrue(os.path.isfile("saved_population.pkl"))
    #
    #     loaded_pop = BinaryPopulation.load("saved_population.pkl")
    #     self.assertTrue(isinstance(loaded_pop, BinaryPopulation))
    #
    # def test_unique_sim_prop(self):
    #     bin_pop = BinaryPopulation()
    #     prop = bin_pop.population_properties
    #     self.assertTrue(
    #         all([prop is b.properties for b in bin_pop]),
    #         msg="All binary properties should map to the same object.",
    #     )

    def tearDown(self):
        # remove pickled files
        if os.path.isfile("saved_population.pkl"):
            os.remove("saved_population.pkl")


if __name__ == "__main__":
    unittest.main()
