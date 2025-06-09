# -*- coding: utf-8 -*-
"""
Written by Pierre-Emmanuel Peyneau
This version: 2025-06-05
Licence: CC-BY 4.0
"""

"""
Program generating a realistic synthetic time scan having both a dissolved and a particulate component (nanoparticles are assumed to be spherical)
"""

"""
Inputs (through input.yaml file):
    <id>: identifier of the simulation. TYPE: str.
    <n_reads>: number of readings of the synthetic time scan. TYPE: int.
    <dwell_time>: dwell time. TYPE: float. UNIT: s.
    <f_transmission>: probability of ion transmission plasma → detector. TYPE: float. UNIT: none.
    <flux_dissolved>: number of analyte ions produced each second from dissolved species. TYPE: float. UNIT: s^{-1}.
    <flux_particles>: number of nanoparticles entering the plasma each second. TYPE: float. UNIT: s^{-1}.
    <detector>: type of detector. TYPE: str.
    <sir_mean>: mean of the logarithm of the detector response. TYPE: float. UNIT: none.
    <sir_std>: standard deviation of the logarithm of the detector response. TYPE: float. UNIT: none.
    <size_distribution>: probability distribution used to draw nanoparticle diameters. TYPE: str.
    <size_particle_mean>: mean of the nanoparticle diameter distribution. TYPE: float. UNIT: cm.
    <size_particle_std>: standard deviation of the nanoparticle diameter distribution. TYPE: float. UNIT: cm.
    <shape_distribution>: probability distribution for the time dispersion distribution of analyte ions produced by one nanoparticle relative to the start of the corresponding particle event. TYPE: str.
    <mean>: mean transport time of an ion from the plasma to the detector (relative to particle event start). TYPE: float. UNIT: s.
    <scale>: scale parameter when <shape_distribution> = 'inverse-gaussian'. TYPE: float. UNIT: s.
    <seed>: seed used in the random generators. TYPE: int.
    <molar_mass>: molar mass of the monitored isotope. TYPE: float. UNIT: g/mol.
    <isotopic_abundance>: natural isotope abundance of the monitored isotope. TYPE: float. UNIT: none.
    <molar_fraction>: molar fraction of the element monitored in the nanoparticles. TYPE: float. UNIT: none.
    <mass_density>: mass density of the nanoparticle. TYPE: float. UNIT: g/cm^3.
    <timescan_csv>: save synthetic time scan as a CSV file? TYPE: boolean.
    <timescan_csv_header>: include a header in the CSV file containing the synthetic time scan? TYPE: boolean.
    <timescan_csv_start>: start time of the time scan. TYPE: float. UNIT: s.
    
Outputs:
    Data files:
        output.yaml: main characteristics of the simulation
        timescan+*.csv: synthetic time scan.
        events+*.csv: list of synthetic particle events contained in the time scan.
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
    def __init__(self, input_dict):
        for key in input_dict:
            setattr(self, key, input_dict[key])


def transportToInverseGaussian(length, velocity, dispersion):
    """
    Transport parameters → mean and shape parameter of the inverse Gaussian distribution
    """
    return length / velocity, length * length / (2.0 * dispersion)


def inverseGaussianMoments(mu, lmbd):
    """
    Parameters
    ----------
    mu : mean of the inverse Gaussian distribution
    lmbd : scale of the inverse Gaussian distribution

    Returns
    -------
    Mean and standard deviation of the distribution
    """
    return mu, (mu**3 / lmbd) ** (0.5)


def poissonProcess(rate, maxtime):
    pptimes = []
    tmax = 0
    if rate > 0:  # otherwise, pptimes remains = []
        while tmax <= maxtime:
            u = rngp.uniform(0, 1)
            dt = -(1.0 / rate) * math.log(
                1 - u
            )  # simulation of an exponential distribution with parameter <rate>
            tmax += dt
            pptimes.append(tmax)
        pptimes.pop()  # the last element of this list is necessarily > maxtime; it is thus erased
    return pptimes


parser = argparse.ArgumentParser()
parser.add_argument(
    "--folder",
    nargs="+",  # to handle multiple arguments associated to a single keyword
    type=str,
)
args = parser.parse_args()
folder_names_list = args.folder


for folder_name in folder_names_list:

    # Folder
    path_input = os.path.join(folder_name, "input.yaml")

    # Opening the yaml input file
    with open(path_input, "r") as inputfile:
        yaml_input = yaml.full_load(inputfile)

    # Input parameters: dictionary to class
    params = InputParams(yaml_input)

    # Seeding
    # Three independent number generators so that the position of the particle events, the value of <flux_dissolved> and the type of detector are all independent from each other
    rngd = np.random.default_rng(params.seed)  # for the dissolved component
    rngp = np.random.default_rng(params.seed)  # for the particulate component
    rngo = np.random.default_rng(params.seed)  # for the detector (output)

    # Dissolved component of the time scan
    n_atoms_dissolved_introduced = rngd.poisson(
        params.flux_dissolved * params.dwell_time, params.n_reads
    )
    n_ions_detector = rngd.binomial(
        n_atoms_dissolved_introduced, params.f_transmission, params.n_reads
    )
    timescan_d = n_ions_detector

    # Particulate component of the time scan
    # Numbers of :
    #    | introduced particles
    #    | particle events (nonzero cumulated number of counts)
    timescan_p = np.zeros(params.n_reads)
    ppjumps = poissonProcess(
        params.flux_particles * params.dwell_time, params.n_reads
    )
    n_particles = len(ppjumps)
    if params.size_distribution == "dirac":
        diameters = np.ones(n_particles) * params.size_particle_mean
    elif params.size_distribution == "truncated-normal":
        diameters = rngp.normal(
            params.size_particle_mean, params.size_particle_std, n_particles
        )
        diameters[diameters < 0.0] = 0.0  # get rid of < 0 sizes
    elif params.size_distribution == "lognormal":
        mu_ln = 2 * math.log(params.size_particle_mean) - 0.5 * math.log(
            params.size_particle_mean**2 + params.size_particle_std**2
        )
        sigma_ln = math.sqrt(
            math.log(
                1 + params.size_particle_std**2 / params.size_particle_mean**2
            )
        )
        diameters = rngp.lognormal(mu_ln, sigma_ln, n_particles)
    atoms_per_particle = (
        diameters**3
        * (PI / 6)
        * params.mass_density
        * params.molar_fraction
        * params.isotopic_abundance
        * AVOGADRO
        / params.molar_mass
    )
    atoms_per_particle = (atoms_per_particle).astype(
        int
    )  # array of floats -> array of ints
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
                # SIR (single-ion response): assumed to be lognormal (as hypothesed in DOI: 10.1039/D4JA00241E)
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
                "{:.3e}".format(np.mean(n_atoms_particles_detected_nonzero))
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
        tra_path = os.path.join(folder_name, "timescan+" + params.id + ".csv")
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
