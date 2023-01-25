#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core SDHDF module
"""

from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
import pkg_resources
import inspect
from typing import List, Optional, Tuple, Union
import warnings

from dask.distributed import Client, get_client, get_task_stream, progress
from dask.diagnostics import ProgressBar
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from astropy.table import Table
from tqdm.auto import tqdm
from pySDHDF import flagging, history
from xarray import DataArray, Variable, Dataset
import xarray as xr


def _decode_df(df: pd.DataFrame) -> pd.DataFrame:
    """Decode a pandas dataframe to a string"""
    str_df = df.select_dtypes([np.object])
    str_df = str_df.stack().str.decode("utf-8").unstack()
    for col in str_df:
        df[col] = str_df[col]
    return df

def _get_sdhdf_version(filename: Path) -> str:
    """Get the SDHDF version of a file

    Args:
        filename (Path): Path to the SDHDF file

    Returns:
        str: SDHDF version
    """
    with h5py.File(filename, "r") as f:
        # Need to hardcode this path for now
        # Assuming it always exists
        try:
            version = f["metadata/primary_header"]["HDR_DEFN_VERSION"][0].decode()
        except KeyError:
            raise KeyError(f"SDHDF version not found in file '{filename}'")
    return version

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
        version = _get_sdhdf_version(self.filename)
        defintion_file = pkg_resources.resource_filename(
            "pySDHDF", f"definitions/sdhdf_def_v{version}.json"
        )
        with open(defintion_file, "r") as f:
            self.definition = json.load(f)
        with h5py.File(self.filename, "r") as f:
            # Get the metadata and configs
            for name in ("metadata", "config"):
                for key, val in self.definition[name].items():
                    df = _decode_df(pd.DataFrame(f[name][val][:]))
                    self.__dict__[key] = df

    def print_metadata(self, format: str = "grid") -> None:
        """Print the metadata to the terminal"""
        for label in self.definition["metadata"].keys():
            df = self.__dict__[label]
            print(f"{label}:")
            print(df.T.to_markdown(tablefmt=format, headers=[]))

    def print_config(self, format: str = "grid") -> None:
        """Print the configuration to the terminal"""
        for label in self.definition["config"].keys():
            df = self.__dict__[label]
            print(f"{label}:")
            print(df.T.to_markdown(tablefmt=format, headers=[]))


@dataclass
class SubBand:
    """An SDHDF sub-band data object

    Args:
        label (str): Sub-band label
        filename (Path): Path to the SDHDF file
        definition (dict): SDHDF definition
        beam_label (str): Beam label
        in_memory (bool, optional): Load the data into memory. Defaults to False.
        client (Client, optional): Dask client. Defaults to None.

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
    definition: dict
    beam_label: str
    in_memory: bool = False
    client: Union[Client, None] = None

    def __post_init__(self):
        # Get the astronomy data
        self.astronomy_dataset = self._get_data()
        # Now get the calibrator data
        self.calibrator_dataset = self._get_cal()
        # TODO: Get the calibrator data

    def _get_cal(self):
        return

    def _get_data(self):
        """Get the astronomy sub-band data"""
        astro_def = self.definition["subband"]["astronomy"]
        sb_path = f"{self.beam_label}/{self.label}"

        with h5py.File(self.filename, "r") as h5:
            data_path = f"{sb_path}/{astro_def['data']}"
            freq_path = f"{sb_path}/{astro_def['frequency']}"
            meta_path = f"{sb_path}/{astro_def['metadata']}"

            data = h5[data_path]
            freqs = h5[freq_path]
            meta = pd.DataFrame(h5[meta_path][:])

            # Get the flags (if they exists)
            has_flags = "flags" in astro_def.keys()
            if has_flags:
                flag_path = f"{sb_path}/{astro_def['flags']}"
                flags = h5[flag_path][:]
                # Ensure flag has same shape as data
                flag_reshape = flags[:].copy()
                for i, s in enumerate(data.shape):
                    if i > len(flag_reshape.shape) - 1:
                        flag_reshape = np.expand_dims(flag_reshape, axis=-1)
                    else:
                        if flag_reshape.shape[i] == s:
                            continue
                        else:
                            flag_reshape = np.expand_dims(flag_reshape, axis=i)
                flags = flag_reshape
            else:
                warnings.warn(
                    f"""
                    No flags found for sub-band '{self.label}' in file '{self.filename}'!
                    SDHDF version is {self.definition['version']}.
                    Flags will be set to all zeros.
                    """
                )
                flags = np.zeros_like(data)

            # Load into memory if requested
            if self.in_memory:
                data = np.array(data)
                freqs = np.array(freqs)
                flags = np.array(flags)

            # Process into xarray
            names = meta.columns
            coords = {name: ("time", meta[name]) for name in names}
            coords["frequency"] = Variable(
                dims="frequency", data=freqs, attrs={"units": freqs.attrs["UNIT"]}
            )
            dims = h5[data_path].attrs["DIMENSION_LABELS"]

            # Need to isel beam 0 here - it will always be dimension 0
            data_xr = DataArray(
                data,
                dims=dims,
                coords=coords,
                name=f"{self.label}_data",
            ).isel(beam=0)
            data_xr.attrs["units"] = data_xr.UNIT

            flag_xr = DataArray(
                flags,
                dims=dims,
                coords=coords,
                name=f"{self.label}_flag",
            ).isel(beam=0)

            astronomy_dataset = Dataset(
                {
                    "data": data_xr,
                    "flag": flag_xr,
                    "metadata": (("time", "meta"), meta),
                }
            )
            return astronomy_dataset

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
        sub_data = self.astronomy_dataset.isel(polarization=polarization, bin=bin)
        if flag:
            sub_data = sub_data.where(sub_data.flag == 0)
        sub_data.data.plot(**plot_kwargs)
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
        sub_data = self.astronomy_dataset.isel(time=time, polarization=polarization, bin=bin)
        if flag:
            sub_data = sub_data.where(sub_data.flag == 0)
        sub_data.data.plot(**plot_kwargs)
        ax = plt.gca()
        return ax

    def autoflag(self, sigma=3, n_windows=100):
        """Automatic flagging using rolling sigma clipping"""
        client = self.client
        data_xr_flg = self.astronomy_dataset.data.where(
            ~self.astronomy_dataset.flag.astype(bool)
        )
        # Set chunks for parallel processing
        chunks = {d:1 for d in data_xr_flg.dims}
        chunks["frequency"] = len(self.astronomy_dataset.data.frequency)
        data_xr_flg = data_xr_flg.chunk(chunks)
        mask = xr.apply_ufunc(
            flagging.box_filter,
            data_xr_flg,
            input_core_dims=[["frequency"]],
            output_core_dims=[["frequency"]],
            kwargs={"sigma": sigma, "n_windows": n_windows},
            dask="parallelized",
            vectorize=True,
            output_dtypes=(bool),
        )
        self.astronomy_dataset["flag"] = mask.astype(int).compute()
        hist = history.generate_history_row()
        return hist

@dataclass
class Beam:
    """An SDHDF beam data object

    Args:
        label (str): The beam label
        filename (Path): The SDHDF file
        definition (dict): The SDHDF definition
        in_memory (bool, optional): Load data into memory. Defaults to False.
        client (Client, optional): Dask client. Defaults to None.

    Attributes:
        subbands (List[SubBand]): A list of subbands

    Methods:
        plot_waterfall: Plot a waterfall plot of the data
        plot_spectrum: Plot a spectrum of the data
        plot_wide: Plot spectra from all subbands

    """

    label: str
    filename: Path
    definition: dict
    in_memory: bool = False
    client: Union[Client, None] = None

    def __post_init__(self):
        with h5py.File(self.filename, "r") as f:
            sb_avail = Table.read(f, path=self.label + "/metadata/band_params")
            self.subbands = [
                SubBand(
                    label=sb,
                    filename=self.filename,
                    definition=self.definition,
                    beam_label=self.label,
                    in_memory=self.in_memory,
                    client=self.client,
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

    def autoflag(self, sigma=3, n_windows=100):
        """Automatic flagging using rolling sigma clipping"""
        hists = []
        for sb in tqdm(self.subbands, desc="Flagging subbands"):
            hist = sb.autoflag(sigma=sigma, n_windows=n_windows,)
            hists.append(hist)
        return hists

@dataclass
class SDHDF:
    """An SDHDF data object

    Args:
        filename (Path): Path to the SDHDF file
        in_memory (bool, optional): Load data into memory. Defaults to False.
        parallel (bool, optional): Use dask for parallel processing. Defaults to False.

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
    parallel: bool = False

    def __post_init__(self):
        self.client = Client() if self.parallel else None
        if self.parallel:
            print(f"Dask dashboard at: {self.client.dashboard_link}")
        self.metadata = MetaData(self.filename)
        self.definition = self.metadata.definition
        with h5py.File(self.filename, "r") as f:
            keys = list(f.keys())
            self.beams = [
                Beam(
                    label=key,
                    filename=self.filename,
                    in_memory=self.in_memory,
                    definition=self.definition,
                    client=self.client,
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

    def print_config(self):
        self.metadata.print_config()


    def flag_persistent_rfi(self):
        """Flag persistent RFI in all subbands."""
        telescope = self.metadata.primary_header["TELESCOPE"].values[0]
        rfi = flagging.get_persistent_rfi(telescope=telescope)
        for i, x in tqdm(rfi.iterrows(), desc="Flagging persistent RFI", total=len(rfi)):
            for beam in self.beams:
                for sb in beam.subbands:
                    sb.astronomy_dataset.flag.loc[
                        dict(
                            frequency=slice(
                                x["freq0 MHz"], x["freq1 MHz"]
                            )
                        )
                    ] = 1
        row = history.generate_history_row()
        self.metadata.history = pd.concat([self.metadata.history, row])

    def auto_flag_rfi(self, sigma=3, n_windows=100,):
        """Automatic flagging using rolling sigma clipping"""
        self.flag_persistent_rfi()
        histss = []
        for beam in tqdm(self.beams, desc="Flagging beams"):
            hists = beam.autoflag(sigma=sigma, n_windows=n_windows)
            histss.extend(hists)

        self.metadata.history = pd.concat([self.metadata.history]+ histss)

    def write(self, filename: Path):
        """Write the SDHDF object to a file.

        Args:
            filename (Path): Filename to write to.
        """
        raise NotImplementedError
