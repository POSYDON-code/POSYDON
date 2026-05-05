__authors__ = [
    "Eirini Kasdagli <kasdaglie@ufl.edu>",
]

import os
import glob
import argparse
import datetime
import numpy as np
import pandas as pd
from posydon.spectral_synthesis.population_spectra import population_spectra

class SpectralSynthesisRunner():

    def __init__(self, population_file, temp_directory='./spectra_batches',
                num_batches=100, verbose=False,
                JOB_ID=None, RANK=None, size=None,  # mirrors BinaryPopulation [2]
                **kwargs):

        if not os.path.isfile(population_file):
            raise FileNotFoundError(
                f'Population file not found: {population_file}'
            )

        self.population_file = population_file
        self.temp_directory = os.path.abspath(temp_directory)
        self.num_batches = num_batches
        self.verbose = verbose
        self.kwargs = kwargs

        # mirrors: JOB_ID and RANK handling in BinaryPopulation.__init__ [2]
        self.JOB_ID = JOB_ID
        self.RANK = RANK
        self.size = size

        self.batch_files = []

    def create_population_batches(self):
        """
        Create batches of POSYDON h5 populations

        Args:
            population_file : str — path to the population h5 file
            temp_directory  : str — directory to save batch files
            batch_size      : int — number of binaries per batch

        Returns:
            batch_files: list of str — paths to created batch files
        """
        if not os.path.isfile(self.population_file):
            raise FileNotFoundError(f'Population file not found: {self.population_file}')

        temp_directory = os.path.abspath(self.temp_directory)
        if not os.path.exists(temp_directory):
            os.makedirs(temp_directory)


        pop = pd.read_hdf(self.population_file, key='history')
        indices = pop.index.unique()
        indices_split = np.array_split(indices, self.num_batches)
        self.batch_files = []
        
        for i in range(self.num_batches):
            batch_file = os.path.join(temp_directory, f'batch_{i}.h5')
            batch = pop.loc[indices_split[i]]
            batch.to_hdf(batch_file, key='history', format='table')
            self.batch_files.append(batch_file)
        return self.batch_files



    def make_spectra(self,batch_id=None):
        """
        Process a single batch file and save its spectral output.
        Intended to be called by one SLURM job per batch.

        Args:
            batch_id       : int — index of the batch to process
            temp_directory : str — directory containing batch files
        """
        
        if batch_id is None:
            if self.JOB_ID is None:
                raise ValueError(
                    'batch_id is required if not running as a SLURM job array'
                )
            batch_id = self.JOB_ID

        batch_file = os.path.join(self.temp_directory, f'batch_{batch_id}.h5')
        if not os.path.isfile(batch_file):
            raise FileNotFoundError(f'Batch file not found: {batch_file}')

        # mirrors: f"evolution.combined.{rank}.h5" naming in BinaryPopulation [2]
        output_file = f'batch_{batch_id}_spectra.h5'
        

        if self.JOB_ID is None:
            output_file = f'batch_{batch_id}_spectra.h5'
        else:
            # BinaryPopulation._safe_evolve [2]
            output_file = (
                f'batch_{batch_id}_spectra.{self.RANK}.h5'
                if self.RANK is not None
                else f'batch_{batch_id:04d}_spectra.h5'
            )
        print(self.temp_directory)
        output_file_path = os.path.join(self.temp_directory, output_file)
        print(f'[Batch {batch_id}] Starting at {datetime.datetime.now()}')
        start_time = datetime.datetime.now()

        ps = population_spectra(
            population_file=batch_file,
            output_file=output_file,
            output_path = self.temp_directory,
            **self.kwargs
        )
        ps.create_spectrum()

        print(f'[Batch {batch_id}] Finished in {datetime.datetime.now() - start_time}')
        return output_file_path

    def combine_batches(self, keep_batches=False,
                        complib='zlib', complevel=9):
        """Combines all the batches into one and cleans up the batch files
        
                Parameters
            ----------
            keep_batches : bool
                If True, do not delete batch files after combining.
                Mirrors the overwrite flag in PopulationRunner [1].
            complib : str
                Compression library.
                Mirrors combine_saved_files in BinaryPopulation [2].
            complevel : int
                Compression level.
                Mirrors combine_saved_files in BinaryPopulation [2].

            Raises
            ------
            FileNotFoundError
                If no spectra batch files are found in temp_directory.
        """
        pop_batch_pattern = os.path.join(
            self.temp_directory, 'batch_*.h5'
        )
        spectra_batch_pattern = os.path.join(
            self.temp_directory, 'batch_*_spectra*.h5'
        )

        pop_batch_files = sorted(glob.glob(pop_batch_pattern))
        spectra_batch_files = sorted(glob.glob(spectra_batch_pattern))

        if not spectra_batch_files:
            raise FileNotFoundError(
                f'No batch spectra files found in {self.temp_directory}'
            )
        
        output_dir = os.path.dirname(self.population_file)
        output_name = os.path.basename(
            self.population_file
        ).replace('.h5', '_spectra.h5')
        h5file = os.path.join(output_dir, output_name)
        
        if os.path.exists(h5file):
            if self.verbose:
                print(f'Removing pre-existing {h5file}...')
            os.remove(h5file)

        
        combined_flux = None
        all_pop_data = []
        
        for i, f in enumerate(spectra_batch_files):
            if self.verbose:
                print(
                    f'  Reading file {i+1}/{len(spectra_batch_files)}: '
                    f'{os.path.basename(f)}'
                )
            try:
                pop_data = pd.read_hdf(f, key='data')
                all_pop_data.append(pop_data)

                flux = pd.read_hdf(f, key='flux')
                #Combining the flux from all the batches
                if combined_flux is None:
                    # first batch — store wavelength and initialise
                    # combined_flux with zeros of the right shape
                    wavelength = flux['wavelength'].values
                    combined_flux = pd.DataFrame(
                        np.zeros((len(wavelength), len(flux.columns) - 1)),
                        columns=[c for c in flux.columns if c != 'wavelength']
                    )
                for col in combined_flux.columns:
                    combined_flux[col] += flux[col].values

            except Exception as e:
                print(f'  Warning: failed to combine batch {f} — {e}')
                continue

        all_pop_data = pd.concat(all_pop_data, ignore_index=True)

        #Saving into the new file
        with pd.HDFStore(
            h5file, mode='w', complevel=complevel, complib=complib
        ) as store:
            store.append('data', all_pop_data, format='table')
            store.append('flux', combined_flux, format='table')


        #Removing all batches
        all_batch_files = pop_batch_files + spectra_batch_files
        if self.verbose:
            print(
                f'Removing {len(all_batch_files)} batch files '
                f'({len(pop_batch_files)} population + '
                f'{len(spectra_batch_files)} spectra)...'
            )
        for f in all_batch_files:
            try:
                os.remove(f)
            except Exception as e:
                print(f'  Warning: could not remove {f} — {e}')
                              
        # only remove i batches directory is empty 
        if len(os.listdir(self.temp_directory)) == 0:
            os.rmdir(self.temp_directory)
            if self.verbose:
                print(
                    f'Removed empty temp directory: '
                    f'{self.temp_directory}'
                )
        else:
            remaining = os.listdir(self.temp_directory)
            print(
                f'temp_directory not removed — '
                f'{len(remaining)} file(s) remaining: {remaining}'
            )




    def run(self, keep_batches=False):
        """Run the full spectral synthesis pipeline locally.

        Mirrors PopulationRunner.evolve [1] — creates batches,
        processes them sequentially, then combines and cleans up.

        Parameters
        ----------
        keep_batches : bool
            If True, do not delete batch files after combining.
        """
        if self.verbose:
            print('=== Creating batches ===')
        self.create_population_batches()


        # mirrors: for pop in self.binary_populations loop in
        # PopulationRunner.evolve [1]
        for i in range(len(self.batch_files)):
            self.make_spectra(batch_id=i)

        if self.verbose:
            print('\n=== Combining and cleaning up ===')
        self.combine_batches(keep_batches=keep_batches)

        if self.verbose:
            print('\nDone!')