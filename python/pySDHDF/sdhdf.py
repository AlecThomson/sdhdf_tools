#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core SDHDF module
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Table
from xarray import DataArray, Variable


def _decode_df(df: pd.DataFrame) -> pd.DataFrame:
    """Decode a pandas dataframe to a string"""
    str_df = df.select_dtypes([np.object])
    str_df = str_df.stack().str.decode("utf-8").unstack()
    for col in str_df:
        df[col] = str_df[col]
    return df


@dataclass
class MetaData:
    """An SDHDF metadata object

    Args:
        filename (Path): Path to the SDHDF file

    Attributes:
        beam_params (DataFrame): The beam parameters as a pandas DataFrame
        history (DataFrame): The history as a pandas DataFrame
        primary_header (DataFrame): The primary header as a pandas DataFrame
        backend_config (DataFrame): The backend configuration as a pandas DataFrame
        cal_backend_config (DataFrame): The calibration backend configuration as a pandas DataFrame

    Methods:
        print_metadata: Quickly list the metadata

    """

    filename: Path

    def __post_init__(self):
        with h5py.File(self.filename, "r") as f:
            meta = f["metadata"]
            self.beam_params = pd.DataFrame(np.array(meta["beam_params"]))
            self.history = pd.DataFrame(np.array(meta["history"]))
            self.primary_header = pd.DataFrame(np.array(meta["primary_header"]))

            config = f["config"]
            self.backend_config = pd.DataFrame(np.array(config["backend_config"]))
            try:
                self.cal_backend_config = pd.DataFrame(
                    np.array(config["cal_backend_config"])
                )
            except KeyError:
                self.cal_backend_config = pd.DataFrame(np.array([]))

            for df in (
                self.beam_params,
                self.history,
                self.primary_header,
                self.backend_config,
            ):
                df = _decode_df(df)

    def print_metadata(self, format: str = "grid") -> None:
        """Quickly list the metadata"""
        for label, df in zip(
            (
                "Beam Parameters",
                "History",
                "Primary Header",
                "Backend Configuration",
                "Calibration Backend Configuration",
            ),
            (
                self.beam_params,
                self.history,
                self.primary_header,
                self.backend_config,
                self.cal_backend_config,
            ),
        ):
            print(f"{label}:")
            print(df.T.to_markdown(tablefmt=format, headers=[]))


@dataclass
class SubBand:
    """An SDHDF sub-band data object

    Args:
        label (str): Sub-band label
        filename (Path): Path to the SDHDF file
        beam_label (str): Beam label
        in_memory (bool, optional): Load the data into memory. Defaults to False.

    Attributes:
        data (DataArray): The sub-band data as an xarray DataArray
        flag (DataArray): The sub-band flag as an xarray DataArray
        meta (DataFrame): The sub-band metadata as a pandas DataFrame

    Methods:
        plot_waterfall: Plot the sub-band data as a waterfall plot
        plot_spectrum: Plot a single spectrum from the sub-band data

    """

    label: str
    filename: Path
    beam_label: str
    in_memory: bool = False

    def __post_init__(self):
        with h5py.File(self.filename, "r") as h5:
            sb_data = f"{self.beam_label}/{self.label}/astronomy_data/data"
            sb_freq = f"{self.beam_label}/{self.label}/astronomy_data/frequency"
            sb_para = f"{self.beam_label}/{self.label}/metadata/obs_params"
            has_flags = (
                "flag" in h5[f"{self.beam_label}/{self.label}/astronomy_data"].keys()
            )
            data = h5[sb_data]
            if has_flags:
                flag = h5[f"{self.beam_label}/{self.label}/astronomy_data/flag"]
                # Ensure flag has same shape as data
                flag_reshape = flag[:].copy()
                for i, s in enumerate(data.shape):
                    if i > len(flag_reshape.shape) - 1:
                        flag_reshape = np.expand_dims(flag_reshape, axis=-1)
                    else:
                        if flag_reshape.shape[i] == s:
                            continue
                        else:
                            flag_reshape = np.expand_dims(flag_reshape, axis=i)
                flag = flag_reshape
            else:
                flag = np.zeros_like(data)
            freq = h5[sb_freq]
            meta = _decode_df(pd.DataFrame(h5[sb_para][:]))
            if self.in_memory:
                data = np.array(data)
                freq = np.array(freq)
                flag = np.array(flag)
            names = meta.columns
            coords = {name: ("time", meta[name]) for name in names}
            coords["frequency"] = Variable(
                dims="frequency", data=freq, attrs={"units": freq.attrs["UNIT"]}
            )
            dims = h5[sb_data].attrs["DIMENSION_LABELS"]

            # Need to isel beam 0 here - it will always be dimension 0
            data_xr = DataArray(
                data,
                dims=dims,
                coords=coords,
                name=f"{self.label}_data",
            ).isel(beam=0)
            data_xr.attrs["units"] = data_xr.UNIT

            flag_xr = DataArray(
                flag,
                dims=dims,
                coords=coords,
                name=f"{self.label}_flag",
            ).isel(beam=0)

            self.data = data_xr
            self.flag = flag_xr
            self.metadata = meta

    def plot_waterfall(
        self,
        polarization: int = 0,
        bin: int = 0,
        flag: bool = False,
        **plot_kwargs,
    ):
        """Waterfall plot of the data

        Args:
            polarization (int, optional): Polarization to select. Defaults to 0.
            bin (int, optional): Bin to select. Defaults to 0.
            flag (bool, optional): Blank flagged data. Defaults to False.
        """
        sub_data = self.data.isel(polarization=polarization, bin=bin)
        if flag:
            sub_flag = self.flag.isel(polarization=polarization, bin=bin)
            sub_data = sub_data.fillna(0).where(sub_flag == 0)
        sub_data.plot(**plot_kwargs)
        ax = plt.gca()
        return ax

    def plot_spectrum(
        self,
        time,
        polarization: int = 0,
        bin: int = 0,
        flag: bool = False,
        **plot_kwargs,
    ):
        sub_data = self.data.isel(time=time, polarization=polarization, bin=bin)
        if flag:
            sub_flag = self.flag.isel(time=time, polarization=polarization, bin=bin)
            sub_data = sub_data.fillna(0).where(sub_flag == 0)
        sub_data.plot(**plot_kwargs)
        ax = plt.gca()
        return ax


@dataclass
class Beam:
    """An SDHDF beam data object

    Args:
        label (str): The beam label
        filename (Path): The SDHDF file
        in_memory (bool, optional): Load data into memory. Defaults to False.

    Attributes:
        subbands (List[SubBand]): A list of subbands

    Methods:
        plot_waterfall: Plot a waterfall plot of the data
        plot_spectrum: Plot a spectrum of the data
        plot_wide: Plot spectra from all subbands

    """

    label: str
    filename: Path
    in_memory: bool = False

    def __post_init__(self):
        with h5py.File(self.filename, "r") as f:
            sb_avail = Table.read(f, path=self.label + "/metadata/band_params")
            self.subbands = [
                SubBand(
                    label=sb,
                    filename=self.filename,
                    beam_label=self.label,
                    in_memory=self.in_memory,
                )
                for sb in sb_avail["LABEL"]
            ]
            for sb in self.subbands:
                self.__dict__[sb.label] = sb

    def plot_waterfall(
        self,
        subband: Union[int, str],
        polarization: int = 0,
        bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(subband, int):
            subband = self.subbands[subband]
        elif isinstance(subband, str):
            subband = self.__dict__[subband]

        ax = subband.plot_waterfall(
            polarization=polarization,
            bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_spectrum(
        self,
        subband: Union[int, str],
        time: int = 0,
        polarization: int = 0,
        bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(subband, int):
            subband = self.subbands[subband]
        elif isinstance(subband, str):
            subband = self.__dict__[subband]

        ax = subband.plot_spectrum(
            time=time,
            polarization=polarization,
            bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_wide(
        self,
        time: int = 0,
        polarization: int = 0,
        bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        fig, ax = plt.subplots()
        for i, sb in enumerate(self.subbands):
            sb.plot_spectrum(
                time=time,
                polarization=polarization,
                bin=bin,
                flag=flag,
                ax=ax,
                label=sb.label,
                **plot_kwargs,
            )
        ax.legend()
        return ax


@dataclass
class SDHDF:
    """An SDHDF data object

    Args:
        filename (Path): Path to the SDHDF file
        in_memory (bool, optional): Load data into memory. Defaults to False.

    Attributes:
        metadata (MetaData): Observation metadata
        beams (List[Beam]): List of beams

    Methods:
        plot_waterfall: Waterfall plot of the data
        plot_spectrum: Spectrum plot of the data
        plot_wide: Plot spectra from all subbands
        print_metadata: List the metadata in the file
        write: Write the data to a new file

    """

    filename: Path
    in_memory: bool = False

    def __post_init__(self):
        self.metadata = MetaData(self.filename)
        with h5py.File(self.filename, "r") as f:
            keys = list(f.keys())
            self.beams = [
                Beam(
                    label=key,
                    filename=self.filename,
                    in_memory=self.in_memory,
                )
                for key in keys
                if "beam_" in key
            ]
            for beam in self.beams:
                self.__dict__[beam.label] = beam

    def plot_waterfall(
        self,
        beam: Union[int, str],
        subband: Union[int, str],
        polarization: int = 0,
        bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        """Waterfall plot of the data

        Args:
            beam (int | str): Beam to select.
            subband (int | str): Subband to select.
            polarization (int, optional): Polarization to select. Defaults to 0.
            bin (int, optional): Bin to select. Defaults to 0.
            flag (bool, optional): Blank flagged data. Defaults to False.
        """
        if isinstance(beam, int):
            beam = self.beams[beam]
        elif isinstance(beam, str):
            beam = self.__dict__[beam]
        ax = beam.plot_waterfall(
            subband=subband,
            polarization=polarization,
            bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_spectrum(
        self,
        beam: Union[int, str],
        subband: Union[int, str],
        time: int = 0,
        polarization: int = 0,
        bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(beam, int):
            beam = self.beams[beam]
        elif isinstance(beam, str):
            beam = self.__dict__[beam]

        ax = beam.plot_spectrum(
            subband=subband,
            time=time,
            polarization=polarization,
            bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def plot_wide(
        self,
        beam: Union[int, str],
        time: int = 0,
        polarization: int = 0,
        bin=0,
        flag: bool = False,
        **plot_kwargs,
    ):
        if isinstance(beam, int):
            beam = self.beams[beam]
        elif isinstance(beam, str):
            beam = self.__dict__[beam]

        ax = beam.plot_wide(
            time=time,
            polarization=polarization,
            bin=bin,
            flag=flag,
            **plot_kwargs,
        )
        return ax

    def print_metadata(self, format: str = "grid"):
        self.metadata.print_metadata(format=format)

    def write(self, filename: Path):
        """Write the SDHDF object to a file.

        Args:
            filename (Path): Filename to write to.
        """
        raise NotImplementedError
