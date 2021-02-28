# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c), The AiiDA-CP2K authors.                                      #
# SPDX-License-Identifier: MIT                                                #
# AiiDA-CP2K is hosted on GitHub at https://github.com/aiidateam/aiida-cp2k   #
# For further information on the license, see the LICENSE.txt file.           #
###############################################################################
"""AiiDA-CP2K output parser."""

from aiida.parsers import Parser
from aiida.common import exceptions
from aiida.engine import ExitCode
from aiida.orm import Dict
from aiida.plugins import DataFactory

from aiida_cp2k.utils.parser import parse_cp2k_output, parse_cp2k_trajectory

StructureData = DataFactory('structure')  # pylint: disable=invalid-name
BandsData = DataFactory('array.bands')  # pylint: disable=invalid-name


class Cp2kBaseParser(Parser):
    """Basic AiiDA parser for the output of CP2K."""
    sections = []

    def parse(self, **kwargs):
        """Receives in input a dictionary of retrieved nodes. Does all the logic here."""

        try:
            _ = self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # Check if dictionary is returned and create output-nodes.
        returned = self._parse_stdout()
        if isinstance(returned, dict):
            self._create_output_nodes(returned)
        else:  # in case this is an error code
            return returned

        try:
            returned = self._parse_trajectory()
            if isinstance(returned, StructureData):
                self.out('output_structure', returned)
            else:  # in case this is an error code
                return returned
        except exceptions.NotExistent:
            pass

        return ExitCode(0)

    def _parse_stdout(self):
        """Basic CP2K output file parser."""

        fname = self.node.get_attribute('output_filename')

        if fname not in self.retrieved.list_object_names():
            return self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING

        try:
            output_string = self.retrieved.get_object_content(fname)
        except IOError:
            return self.exit_codes.ERROR_OUTPUT_STDOUT_READ

        result_dict = parse_cp2k_output(output_string, self.sections)

        # All exit_codes are triggered here
        if "geo_not_converged" in result_dict:
            return self.exit_codes.ERROR_GEOMETRY_CONVERGENCE_NOT_REACHED
        if "uks_needed" in result_dict:
            return self.exit_codes.ERROR_UKS_NEEDED
        if "aborted" in result_dict:
            return self.exit_codes.ERROR_OUTPUT_CONTAINS_ABORT

        return result_dict

    def _parse_trajectory(self):
        """CP2K trajectory parser."""

        from ase import Atoms

        fname = self.node.process_class._DEFAULT_RESTART_FILE_NAME  # pylint: disable=protected-access

        # Check if the restart file is present.
        if fname not in self.retrieved.list_object_names():
            raise exceptions.NotExistent("No restart file available, so the output trajectory can't be extracted.")

        # Read the restart file.
        try:
            output_string = self.retrieved.get_object_content(fname)
        except IOError:
            return self.exit_codes.ERROR_OUTPUT_STDOUT_READ

        return StructureData(ase=Atoms(**parse_cp2k_trajectory(output_string)))

    def _create_output_nodes(self, result_dict):
        self.out("output_parameters", Dict(dict=result_dict))


class Cp2kAdvancedParser(Cp2kBaseParser):
    """Advanced AiiDA parser class for the output of CP2K."""
    sections = ['spin_density', 'natoms', 'scf_parameters', 'init_nel', 'kpoint_data', 'motion_info']

    def _create_output_nodes(self, result_dict):
        # Compute the bandgap for Spin1 and Spin2 if eigen was parsed (works also with smearing!)
        if 'eigen_spin1_au' in result_dict:
            if result_dict['dft_type'] == "RKS":
                result_dict['eigen_spin2_au'] = result_dict['eigen_spin1_au']

            lumo_spin1_idx = result_dict['init_nel_spin1']
            lumo_spin2_idx = result_dict['init_nel_spin2']
            if (lumo_spin1_idx > len(result_dict['eigen_spin1_au'])-1) or \
               (lumo_spin2_idx > len(result_dict['eigen_spin2_au'])-1):
                #electrons jumped from spin1 to spin2 (or opposite): assume last eigen is lumo
                lumo_spin1_idx = len(result_dict['eigen_spin1_au']) - 1
                lumo_spin2_idx = len(result_dict['eigen_spin2_au']) - 1
            homo_spin1 = result_dict['eigen_spin1_au'][lumo_spin1_idx - 1]
            homo_spin2 = result_dict['eigen_spin2_au'][lumo_spin2_idx - 1]
            lumo_spin1 = result_dict['eigen_spin1_au'][lumo_spin1_idx]
            lumo_spin2 = result_dict['eigen_spin2_au'][lumo_spin2_idx]
            result_dict['bandgap_spin1_au'] = lumo_spin1 - homo_spin1
            result_dict['bandgap_spin2_au'] = lumo_spin2 - homo_spin2

        if "kpoint_data" in result_dict:
            bnds = BandsData()
            bnds.set_kpoints(result_dict["kpoint_data"]["kpoints"])
            bnds.labels = result_dict["kpoint_data"]["labels"]
            bnds.set_bands(
                result_dict["kpoint_data"]["bands"],
                units=result_dict["kpoint_data"]["bands_unit"],
            )
            self.out("output_bands", bnds)
            del result_dict["kpoint_data"]
        self.out("output_parameters", Dict(dict=result_dict))
