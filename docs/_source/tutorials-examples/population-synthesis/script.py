from posydon.popsyn.synthetic_population import SyntheticPopulation

if __name__ == "__main__":
    synth_pop = SyntheticPopulation(f"./population_params.ini")
    synth_pop.evolve()
