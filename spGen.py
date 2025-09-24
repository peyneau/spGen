"""Synthetic time scan generator.

This program generates a realistic synthetic time scan with both dissolved and
particulate components, simulating nanoparticle events as spherical particles.
The simulation is configured via an input YAML file and produces output files
with the synthetic scan and event information.

Notes
-----

Written by: Pierre-Emmanuel Peyneau
This version: 2025-06-05
Licence: CC-BY 4.0

Parameters
----------

The input parameters are read from an input YAML file and include. An example
YAML file is available in the [README](README.md) file.

id : str
    Tag identifier of the simulation (e.g., 'test1').
n_reads : int
    Number of readings of the synthetic time scan.
dwell_time : float
    Dwell time per reading (seconds).
f_transmission : float
    Probability of ion transmission from plasma to detector.
flux_dissolved : float
    Number of analyte ions produced each second from dissolved species
    (1/second).
flux_particles : float
    Number of nanoparticles entering the plasma each second (1/second).
detector : {'ideal', 'non-ideal'}
    Type of detector.
sir_mean : float
    Mean of the logarithm of the detector response.
sir_std : float
    Standard deviation of the logarithm of the detector response.
size_distribution : {'dirac', 'truncated-normal', 'lognormal'}
    Probability distribution used to draw nanoparticle diameters.
size_particle_mean : float
    Mean of the nanoparticle diameter distribution (cm).
size_particle_std : float
    Standard deviation of the nanoparticle diameter distribution (cm).
shape_distribution : {'inverse-gaussian', 'uniform'}
    Probability distribution for the time dispersion of analyte ions per
    nanoparticle event.
mean : float
    Mean transport time of an ion from plasma to detector, relative to particle
    event start (seconds).
scale : float
    Scale parameter for 'inverse-gaussian' shape_distribution (seconds).
seed : int
    Seed for random number generators.
molar_mass : float
    Molar mass of the monitored isotope (g/mol).
isotopic_abundance : float
    Natural isotope abundance of the monitored isotope.
mass_fraction : float
    Mass fraction of the monitored element in the nanoparticles.
mass_density : float
    Mass density of the nanoparticle (g/cm^3).
timescan_csv : bool
    Whether to save the synthetic time scan as a CSV file.
timescan_csv_header : bool
    Whether to include a header in the CSV file.
timescan_csv_start : float
    Start time of the time scan (seconds).

Written by Pierre-Emmanuel Peyneau This version: 2025-06-05 Licence: CC-BY 4.0


Returns
-------

The program generates the following output files in the specified folder:

output.yaml
    Main characteristics of the simulation
timescan+*.csv
    Synthetic time scan.
events+*.csv:
    list of synthetic particle events contained in the time scan.
"""

__date__ = "2025-06-05"
__codename__ = "spGen"
__version__ = "1-0"


import math
import numpy as np
from scipy import constants
import csv
import yaml
import os
import argparse
import datetime


TIME_PREC = 10
PI = math.pi
AVOGADRO = constants.Avogadro


class InputParams:
    """Class to handle input parameters from a YAML file."""

    def __init__(self, input_dict: dict):
        for key in input_dict:
            setattr(self, key, input_dict[key])


def transportToInverseGaussian(
    length: float,
    velocity: float,
    dispersion: float,
) -> tuple:
    """Transport parameters to inverse Gaussian distribution parameters.

    Parameters
    ----------
    length : float
        Length of the transport path from the plasma to the detector (cm).
    velocity : float
        Mean velocity of the analyte ions in the transport path (cm/s).
    dispersion : float
        Dispersion of the analyte ions in the transport path (cm^2/s).

    Returns
    -------
    tuple
        Mean and scale parameters of the inverse Gaussian distribution.
    """
    return length / velocity, length * length / (2.0 * dispersion)


def inverseGaussianMoments(mu: float, lmbd: float) -> tuple:
    """Moments of the inverse Gaussian distribution.

    Parameters
    ----------
    mu : float
        Theoretical mean of the inverse Gaussian distribution
    lmbd : float
        Theoretical scale parameter of the inverse Gaussian distribution (also
        called shape parameter).

    Returns
    -------
    tuple
        Mean and standard deviation of the inverse Gaussian distribution.
    """
    return mu, (mu**3 / lmbd) ** (0.5)


def poisson_process(rate, maxtime, seed=np.random.default_rng()):
    """Generate a Poisson process with a given rate and maximum time.

    Parameters
    ----------
    rate : float
        Rate of the Poisson process (events per second).
    maxtime : float
        Maximum time for the process (seconds).
    Returns
    -------
    list
        List of event times generated by the Poisson process.
    """
    pptimes = []
    tmax = 0
    if rate > 0:  # otherwise, pptimes remains = []
        while tmax <= maxtime:
            u = seed.uniform(0, 1)
            # Exponential distribution with parameter <rate>
            dt = -(1.0 / rate) * math.log(1 - u)
            tmax += dt
            pptimes.append(tmax)
        # erate the last element which is necessarily > maxtime
        pptimes.pop()
    return pptimes


def parser():
    """Parse command line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="spGen: Synthetic time scan generator"
    )
    parser.add_argument(
        "--folder",
        nargs="+",
        type=str,
        required=True,
        help="List of folder names containing input.yaml files.",
    )
    return parser.parse_args()


def main():

    args = parser()
    folder_names_list = args.folder

    for folder_name in folder_names_list:

        # Folder
        path_input = os.path.join(folder_name, "input.yaml")

        # Open the yaml input file
        with open(path_input, "r") as inputfile:
            yaml_input = yaml.full_load(inputfile)

        # Input parameters: dictionary to class
        params = InputParams(yaml_input)

        # Seeding three independent number generators so that the position of
        # the particle events, the value of <flux_dissolved> and the type of
        # detector are all independent from each other:
        # Dissolved component:
        rngd = np.random.default_rng(params.seed)
        # Particulate component:
        rngp = np.random.default_rng(params.seed)
        # Detector component:
        rngo = np.random.default_rng(params.seed)

        # Dissolved component of the time scan
        n_atoms_dissolved_introduced = rngd.poisson(
            params.flux_dissolved * params.dwell_time, params.n_reads
        )
        n_ions_detector = rngd.binomial(
            n_atoms_dissolved_introduced,
            params.f_transmission,
            params.n_reads,
        )
        timescan_d = n_ions_detector

        # Particulate component of the time scan
        # Numbers of :
        #    | introduced particles
        #    | particle events (nonzero cumulated number of counts)
        timescan_p = np.zeros(params.n_reads)
        ppjumps = poisson_process(
            params.flux_particles * params.dwell_time,
            params.n_reads,
            seed=rngp,
        )
        n_particles = len(ppjumps)
        if params.size_distribution == "dirac":
            diameters = np.ones(n_particles) * params.size_particle_mean
        elif params.size_distribution == "truncated-normal":
            diameters = rngp.normal(
                params.size_particle_mean,
                params.size_particle_std,
                n_particles,
            )
            diameters[diameters < 0.0] = 0.0  # get rid of < 0 sizes
        elif params.size_distribution == "lognormal":
            mu_ln = 2 * math.log(params.size_particle_mean) - 0.5 * math.log(
                params.size_particle_mean**2 + params.size_particle_std**2
            )
            sigma_ln = math.sqrt(
                math.log(
                    1
                    + params.size_particle_std**2
                    / params.size_particle_mean**2
                )
            )
            diameters = rngp.lognormal(mu_ln, sigma_ln, n_particles)
        atoms_per_particle = (
            diameters**3
            * (PI / 6)
            * params.mass_density
            * params.mass_fraction
            * params.isotopic_abundance
            * AVOGADRO
            / params.molar_mass
        )
        # Turn the array of floats into an array of ints
        atoms_per_particle = (atoms_per_particle).astype(int)
        n_atoms_particles_detected = rngp.binomial(
            atoms_per_particle, params.f_transmission, n_particles
        )
        n_atoms_particles_detected_nonzero = n_atoms_particles_detected[
            n_atoms_particles_detected != 0
        ]
        n_particle_events = len(n_atoms_particles_detected_nonzero)

        # Particle event duration
        particle_event_duration = np.zeros(n_particle_events)
        particle_event_start = np.zeros(n_particle_events)
        particle_event_end = np.zeros(n_particle_events)
        i_particle_event = 0  # counter
        for i in range(n_particles):
            # inverse gaussian distribution of the spike duration
            if params.shape_distribution == "inverse-gaussian":
                time = ppjumps[i] * np.ones(
                    n_atoms_particles_detected[i]
                ) + rngp.wald(
                    params.mean / params.dwell_time,
                    params.scale / params.dwell_time,
                    n_atoms_particles_detected[i],
                )
            # uniform distribution of the spike duration
            elif params.shape_distribution == "uniform":
                time = ppjumps[i] * np.ones(
                    n_atoms_particles_detected[i]
                ) + rngp.uniform(
                    0.0,
                    2 * params.mean / params.dwell_time,
                    n_atoms_particles_detected[i],
                )
            time = np.rint(time)
            if time.size > 0:
                particle_event_duration[i_particle_event] = (
                    int(max(time)) - int(min(time)) + 1
                )
                particle_event_start[i_particle_event] = int(min(time))
                particle_event_end[i_particle_event] = int(max(time))
                i_particle_event += 1
            for j in range(len(time)):
                if time[j] < params.n_reads:
                    timescan_p[int(time[j])] += 1

        # Addition of baseline and particle event components
        timescan_predetect = timescan_d + timescan_p

        # Detection
        if params.detector == "ideal":
            timescan = timescan_predetect
        elif params.detector == "non-ideal":
            timescan = np.zeros(params.n_reads)
            for i in range(params.n_reads):
                if timescan_predetect[i] > 0:
                    # SIR (single-ion response): assumed to be lognormal (as
                    # hypothesed in DOI: 10.1039/D4JA00241E)
                    timescan[i] = np.sum(
                        rngo.lognormal(
                            params.sir_mean,
                            params.sir_std,
                            timescan_predetect[i].astype(int),
                        )
                    )

        # Output yaml file
        if n_particle_events > 0:
            output_dict = {
                "datetime": datetime.datetime.now(),
                "id": params.id,
                "code_name": __codename__,
                "code_version": __version__,
                "n_particles_introduced": n_particles,
                "n_particle_events": n_particle_events,
                "particle_event_cumulative_counts_mean": float(
                    "{:.3e}".format(
                        np.mean(n_atoms_particles_detected_nonzero)
                    )
                ),
                "particle_event_cumulative_counts_std": float(
                    "{:.3e}".format(np.std(n_atoms_particles_detected_nonzero))
                ),
                "particle_event_duration_mean": float(
                    "{:.3e}".format(
                        np.mean(particle_event_duration) * params.dwell_time
                    )
                ),
                "particle_event_duration_std": float(
                    "{:.3e}".format(
                        np.std(particle_event_duration) * params.dwell_time
                    )
                ),
            }
        else:
            output_dict = {
                "id": params.id,
                "n_particles_introduced": n_particles,
                "n_particle_events": 0,
                "particle_event_duration_mean": 0.0,
                "particle_event_duration_std": 0.0,
                "particle_event_cumulative_counts_mean": 0.0,
                "particle_event_cumulative_counts_std": 0.0,
            }
        path_output = os.path.join(folder_name, "output.yaml")
        with open(path_output, "w") as f:
            yaml.dump(output_dict, f, sort_keys=False, default_style=None)

        # Dump time scan as a .csv file
        if params.timescan_csv:
            tra_path = os.path.join(
                folder_name, "timescan+" + params.id + ".csv"
            )
            with open(tra_path, mode="w", newline="\n") as file:
                writer = csv.writer(file)
                if params.timescan_csv_header:
                    writer.writerow(["Synthetic time scan,"])
                    writer.writerow(["Intensity Vs Time", "Counts"])
                    writer.writerow(["Generated: numerically,"])
                    writer.writerow(["Time [s],"])
                for index, item in enumerate(timescan):
                    writer.writerow(
                        [
                            round(
                                (index + 1) * params.dwell_time
                                + params.timescan_csv_start,
                                TIME_PREC,
                            ),
                            item,
                        ]
                    )

        # Dump particle event starts and ends in a .csv file
        event_path = os.path.join(folder_name, "events+" + params.id + ".csv")
        i_particle_event = 0  # counter
        with open(event_path, mode="w", newline="\n") as file:
            writer = csv.writer(file)
            for i in range(n_particles):
                if n_atoms_particles_detected[i] > 0:
                    writer.writerow(
                        [
                            i_particle_event + 1,
                            round(
                                particle_event_start[i_particle_event]
                                * params.dwell_time
                                + params.timescan_csv_start,
                                TIME_PREC,
                            ),
                            round(
                                particle_event_end[i_particle_event]
                                * params.dwell_time
                                + params.timescan_csv_start,
                                TIME_PREC,
                            ),
                        ]
                    )
                    i_particle_event += 1  # counter


if __name__ == "__main__":
    main()

